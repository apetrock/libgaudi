#ifndef GAUDI_ENGMAN_SOLVER_HPP
#define GAUDI_ENGMAN_SOLVER_HPP

/// The Engman solver (`scripts/sympy/cr_poisson/MATH.md` section 17): one
/// implicit step over the hepworth blocks (positions, the edge store, ...)
/// with soft energies and hard constraints.
///
/// - Energies are hepworth rows as they are: projection rows contribute
///   `1/2 |A q - p(q)|^2` (`fill_A`, `project`), Gauss-Newton rows
///   `1/2 w^2 |r(q)|^2` (`fill_J`, `add_gradient`, `residual_sq`), and the
///   blocks' inertia `1/2 |q - s|^2_{M / h^2}`, `s` the momentum prediction.
///   Their weights keep their meaning: these are a compromise by intent.
/// - Constraints (`engman::constraint`) are met: each step is a Newton-KKT
///   step of the proximal method of multipliers,
///
///       [ H   J^T     ] [ dq      ]     [ g + J^T lambda ]
///       [ J  -I / rho ] [ dlambda ] = - [ C              ],
///
///   `H` the Gauss-Newton Hessian of the energies. The step meets the
///   linearized constraints and the energies act only in their tangent space,
///   so an energy never trades against a constraint. The matrix is symmetric
///   quasi-definite: one `SimplicialLDLT`, no pivoting. The factor is kept
///   for `refresh_every` iterates with a fresh right-hand side; each step
///   takes a line search on the merit `E + lambda^T C + rho / 2 |C|^2`.

#include "gaudi/engman/constraint.hpp"
#include "gaudi/engman/kkt.hpp"
#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/cache_compute.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {
struct solver_context;
}
} // namespace hepworth

namespace engman {

using matS = Eigen::SparseMatrix<real>;

struct solver_config {
  /// The multiplier penalty. Above 0 it is used as is; otherwise `rho_rel`
  /// times the mean diagonal of `H` over the mean squared constraint row at
  /// the first factor. A conditioning parameter: with the duals on, the
  /// solution does not depend on it. Usable up to `rho_rel` about 1e6; at
  /// 1e9 the quasi-definite factor loses the constraint block and steps
  /// crawl.
  real rho = 0;
  real rho_rel = 1e6;
  /// Off: `lambda` stays as it is and the constraints act as the penalty
  /// `rho / 2 |C|^2` (`rho = w^2` reproduces a `w`-weighted penalty row).
  bool dual_updates = true;
  /// Iterates a factor is kept for (fresh right-hand side each); 1 refactors
  /// every iterate, 0 only when a step cannot descend.
  int refresh_every = 1;
  /// Step halvings tried before a step counts as failed.
  int backtracks = 8;
  /// Refactors after a failed step on a stale factor; a failed step on a
  /// fresh factor ends the step.
  int max_ascent_refreshes = 2;
  /// The step ends when the update is below `tol` relative to `|q|`.
  real tol = 1e-13;
  /// Diagonal floor added to `H` in the factor, relative to its mean
  /// diagonal (variables no energy holds).
  real reg = 1e-12;
  /// Feasibility restoration after the iterates: Newton steps on `C = 0`
  /// alone (the KKT step with no energy gradient: the smallest correction in
  /// the `H` metric), so every step ends on the constraints even when the
  /// iterates ran out descending. Stops when `|C|` stops falling 10x. Only
  /// with `dual_updates` (as penalties the constraints stay soft).
  int restore = 6;
  /// The KKT solve (`kkt_mode`), one of:
  ///   kkt_mode::full       one factor of the whole quasi-definite matrix
  ///                        (default: the same fill as `projection_solver`'s
  ///                        matrix, 3.9 against 4.4 s on 8192 vertices)
  ///   kkt_mode::schur      the coupling's store blocks and multipliers
  ///                        eliminated per edge onto positions (falls back
  ///                        to full when a row couples store edges); 7x
  ///                        smaller but fills more (5.9 s): eliminating the
  ///                        store couples the 2-ring (MATH.md 17.2)
  ///   kkt_mode::automatic  schur when it applies, else full
  kkt_mode kkt = kkt_mode::full;
};

struct solver_report {
  int iterations = 0;
  int refreshes = 0;
  real E_first = 0, E_last = 0;
  real C_first = 0, C_last = 0;  ///< max |C|
  real rho = 0;
  bool stopped = false;  ///< a step on a fresh factor could not descend
  int restorations = 0;  ///< feasibility steps taken after the iterates
  bool schur = false;  ///< the last factor was the Schur reduction
  double factor_seconds = 0;  ///< total over the step's factors
  index_t factored_size = 0, factor_nonzeros = 0;  ///< of the last factor
};

class solver {
public:
  solver_config cfg;
  solver_report last;

  void set_energies(std::vector<hepworth::projection_constraint::ptr> energies) {
    _energies = std::move(energies);
  }
  void set_constraints(std::vector<constraint::ptr> constraints) { _constraints = std::move(constraints); }
  void set_computes(std::vector<hepworth::cache_compute::ptr> computes) {
    _computes = std::move(computes);
  }

  void step(std::vector<hepworth::sim_block::ptr> &blocks, real h, real damping, int iterations,
            hepworth::block::solver_context *rest_ctx = nullptr) {
    vecX q, s, M;
    for (auto &b : blocks)
      b->map_mass(M);
    for (index_t i = 0; i < M.size(); ++i)
      M[i] = std::max(M[i], real(1e-8));
    const vecX Mh = M / (h * h);
    for (auto &b : blocks)
      b->map_to_x(q);
    if (rest_ctx)
      for (auto &e : _energies)
        if (auto bc = std::dynamic_pointer_cast<hepworth::block::block_constraint>(e))
          if (!bc->rest_initialized())
            bc->init_rest(q, *rest_ctx);
    for (auto &b : blocks) {
      b->integrate_inertia(h);
      b->map_to_x(s);
    }
    q = s;
    for (auto &c : _computes)
      if (c)
        c->compute(q);
    const index_t n = q.size();

    // Projection rows are fixed; Gauss-Newton rows are refilled per factor.
    std::vector<int> proj, gn;
    std::vector<trip> at;
    index_t np = 0;
    for (int i = 0; i < static_cast<int>(_energies.size()); ++i) {
      if (_energies[static_cast<size_t>(i)]->gauss_newton()) {
        gn.push_back(i);
      } else {
        proj.push_back(i);
        _energies[static_cast<size_t>(i)]->fill_A(np, at);
      }
    }
    matS Ap(np, n);
    Ap.setFromTriplets(at.begin(), at.end());
    vecX p = vecX::Zero(np);
    auto project_all = [&](const vecX &qk) {
      p.setZero();
      for (int i : proj)
        _energies[static_cast<size_t>(i)]->project(qk, p);
    };

    index_t m = 0;
    std::vector<index_t> c0;
    for (auto &c : _constraints) {
      c0.push_back(m);
      m += c->rows();
    }
    vecX lambda = vecX::Zero(m);
    for (size_t j = 0; j < _constraints.size(); ++j)
      if (_constraints[j]->lambda.size() == _constraints[j]->rows())
        lambda.segment(c0[j], _constraints[j]->rows()) = _constraints[j]->lambda;
    auto eval_C = [&](const vecX &qk, vecX &C, std::vector<trip> *J) {
      C.resize(m);
      for (size_t j = 0; j < _constraints.size(); ++j) {
        vecX Cj;
        std::vector<trip> Jj;
        _constraints[j]->evaluate(qk, Cj, J ? &Jj : nullptr);
        C.segment(c0[j], _constraints[j]->rows()) = Cj;
        if (J)
          for (const trip &t : Jj)
            J->emplace_back(c0[j] + t.row(), t.col(), t.value());
      }
    };
    auto energy = [&](const vecX &qk) {
      project_all(qk);
      const vecX d = qk - s;
      real e = d.dot(Mh.cwiseProduct(d)) + (Ap * qk - p).squaredNorm();
      for (int i : gn)
        e += _energies[static_cast<size_t>(i)]->residual_sq(qk);
      return real(0.5) * e;
    };
    auto gradient = [&](const vecX &qk) {
      project_all(qk);
      vecX g = Mh.cwiseProduct(qk - s) + Ap.transpose() * (Ap * qk - p);
      for (int i : gn)
        _energies[static_cast<size_t>(i)]->add_gradient(qk, g);
      return g;
    };
    real rho = cfg.rho;
    auto merit = [&](const vecX &qk) {
      vecX C;
      eval_C(qk, C, nullptr);
      return energy(qk) + lambda.dot(C) + real(0.5) * rho * C.squaredNorm();
    };

    std::vector<std::pair<index_t, index_t>> pairs;
    for (size_t j = 0; j < _constraints.size(); ++j) {
      std::vector<std::pair<index_t, index_t>> pj;
      _constraints[j]->pairs(pj);
      for (const auto &[r, c] : pj)
        pairs.emplace_back(c0[j] + r, c);
    }
    kkt_system K;
    auto factor = [&](const vecX &qk) {
      std::vector<trip> gt;
      index_t ng = 0;
      for (int i : gn)
        _energies[static_cast<size_t>(i)]->fill_J(qk, ng, gt);
      matS Ag(ng, n);
      Ag.setFromTriplets(gt.begin(), gt.end());
      matS H = matS(Ap.transpose() * Ap) + matS(Ag.transpose() * Ag);
      vecX C;
      std::vector<trip> jt;
      eval_C(qk, C, &jt);
      matS J(m, n);
      J.setFromTriplets(jt.begin(), jt.end());
      real hd = 0;
      for (index_t i = 0; i < n; ++i)
        hd += Mh[i] + H.coeff(i, i);
      hd /= std::max<real>(real(n), 1);
      if (!(rho > real(0))) {
        real jd = 0;
        for (index_t r = 0; r < m; ++r)
          jd += J.row(r).squaredNorm();
        jd /= std::max<real>(real(m), 1);
        // No constraint rows (energies only): rho is unused, and the
        // automatic value would be infinite and turn the merit's
        // `rho / 2 |C|^2` into inf * 0.
        rho = m > 0 && jd > real(0) ? cfg.rho_rel * hd / jd : real(1);
      }
      std::vector<trip> dt;
      for (index_t i = 0; i < n; ++i)
        dt.emplace_back(i, i, Mh[i] + cfg.reg * hd);
      matS Dg(n, n);
      Dg.setFromTriplets(dt.begin(), dt.end());
      K.factor(matS(H + Dg), J, rho, pairs, cfg.kkt);
      ++last.refreshes;
      last.schur = K.schur();
      last.factored_size = K.factored_size();
      last.factor_nonzeros = K.factor_nonzeros();
    };

    last = solver_report{};
    {
      vecX C;
      eval_C(q, C, nullptr);
      last.E_first = energy(q);
      last.C_first = m > 0 ? C.cwiseAbs().maxCoeff() : real(0);
    }
    bool need_factor = true, fresh = false;
    int ascent_refreshes = 0;
    int k = 0;
    for (; k < iterations; ++k) {
      const bool scheduled = cfg.refresh_every > 0 && k % cfg.refresh_every == 0;
      if (need_factor || scheduled) {
        factor(q);
        need_factor = false;
        fresh = true;
      }
      if (!K.ok())
        break;
      vecX C;
      std::vector<trip> jt;
      eval_C(q, C, &jt);
      matS J(m, n);
      J.setFromTriplets(jt.begin(), jt.end());
      vecX dq, dl;
      K.solve(-(gradient(q) + J.transpose() * lambda), -C, dq, dl);
      const real phi0 = merit(q);
      real a = 1;
      bool ok = false;
      vecX q1;
      for (int b = 0; b <= cfg.backtracks; ++b, a *= real(0.5)) {
        q1 = q + a * dq;
        const real phi1 = merit(q1);
        if (std::isfinite(phi1) && phi1 <= phi0 + real(1e-14) * std::abs(phi0)) {
          ok = true;
          break;
        }
      }
      if (!ok) {
        if (fresh || ascent_refreshes >= cfg.max_ascent_refreshes) {
          last.stopped = true;
          break;
        }
        ++ascent_refreshes;
        need_factor = true;
        continue;
      }
      fresh = false;
      q = q1;
      if (cfg.dual_updates)
        lambda += a * dl;
      if (a * dq.cwiseAbs().maxCoeff() <= cfg.tol * (real(1) + q.cwiseAbs().maxCoeff()))
        break;
    }
    last.iterations = std::min(k + 1, iterations);
    for (int r = 0; r < cfg.restore && m > 0 && cfg.dual_updates; ++r) {
      vecX C;
      eval_C(q, C, nullptr);
      const real c0n = C.cwiseAbs().maxCoeff();
      if (!(c0n > real(0)))
        break;
      if (r == 0)
        factor(q);
      if (!K.ok())
        break;
      vecX dq, dl;
      K.solve(vecX::Zero(n), -C, dq, dl);
      const vecX q1 = q + dq;
      vecX C1;
      eval_C(q1, C1, nullptr);
      const real c1n = C1.cwiseAbs().maxCoeff();
      if (!(c1n < c0n))
        break;
      q = q1;
      ++last.restorations;
      if (!(c1n < real(0.1) * c0n))
        break;
    }
    last.rho = rho;
    last.factor_seconds = K.seconds();
    {
      vecX C;
      eval_C(q, C, nullptr);
      last.E_last = energy(q);
      last.C_last = m > 0 ? C.cwiseAbs().maxCoeff() : real(0);
    }
    for (size_t j = 0; j < _constraints.size(); ++j)
      _constraints[j]->lambda = lambda.segment(c0[j], _constraints[j]->rows());
    for (auto &b : blocks)
      b->map_from_x(q, h, damping);
  }

private:
  std::vector<hepworth::projection_constraint::ptr> _energies;
  std::vector<constraint::ptr> _constraints;
  std::vector<hepworth::cache_compute::ptr> _computes;
};

} // namespace engman
} // namespace gaudi

#endif
