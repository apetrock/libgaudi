
#ifndef __HEP_ROD_SOLVER__
#define __HEP_ROD_SOLVER__

#include "gaudi/sparse_solver.h"

#include "projection_diagnostics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "block_constraint.hpp"
#include "cache_compute.hpp"
#include "gaudi/common.h"
#include "gaudi/kusama/laplacian.hpp"
#include "sim_block.hpp"
#include "modal_corrector.hpp"
#include "variational_audit.hpp"

///////////////////////////////////////////////////////
// solver
///////////////////////////////////////////////////////

namespace gaudi {
namespace hepworth {
namespace block {

struct solver_context;

namespace detail {

inline size_t hash_triplets(const std::vector<trip> &triplets) {
  size_t h = triplets.size();
  for (const trip &t : triplets) {
    h ^= std::hash<int>()(static_cast<int>(t.row())) + 0x9e3779b9 + (h << 6) +
         (h >> 2);
    h ^= std::hash<int>()(static_cast<int>(t.col())) + 0x9e3779b9 + (h << 6) +
         (h >> 2);
    h ^= std::hash<double>()(t.value()) + 0x9e3779b9 + (h << 6) + (h >> 2);
  }
  return h;
}

inline real mass_diagonal_sum(const matS &M) {
  real sum = 0.0;
  for (int k = 0; k < M.outerSize(); ++k) {
    for (matS::InnerIterator it(M, k); it; ++it) {
      if (it.row() == it.col())
        sum += it.value();
    }
  }
  return sum;
}

struct linear_system_cache_key {
  int nm = 0;
  index_t a_rows = 0;
  size_t constraint_count = 0;
  size_t triplet_hash = 0;
  real h = 0.0;
  real mass_sum = 0.0;

  bool operator==(const linear_system_cache_key &o) const {
    return nm == o.nm && a_rows == o.a_rows &&
           constraint_count == o.constraint_count &&
           triplet_hash == o.triplet_hash && h == o.h && mass_sum == o.mass_sum;
  }
};

struct p_span {
  std::string name;
  index_t begin = 0;
  index_t end = 0;
};

/// Sum of `p` over each constraint's rows. Printed on the same iterations as
/// the solver trace, and once when `p` or `q` leaves the reals.
inline void log_p_by_name(int k, const vecX &q, const vecX &p,
                          const std::vector<p_span> &spans) {
  struct acc {
    real n2 = 0;
    int n = 0;
  };
  std::map<std::string, acc> by_name;
  for (const p_span &sp : spans) {
    acc &a = by_name[sp.name];
    ++a.n;
    if (sp.end > sp.begin && sp.begin >= 0 && sp.end <= p.size())
      a.n2 += p.segment(sp.begin, sp.end - sp.begin).squaredNorm();
  }
  std::vector<std::pair<std::string, acc>> ranked(by_name.begin(), by_name.end());
  std::sort(ranked.begin(), ranked.end(),
            [](const auto &u, const auto &v) { return u.second.n2 > v.second.n2; });
  gaudi::duchamp::demo_log() << "[m_solver] p iter=" << (k + 1)
                             << " |q|=" << q.norm() << " |p|=" << p.norm();
  const int show = std::min(4, static_cast<int>(ranked.size()));
  for (int i = 0; i < show; ++i) {
    gaudi::duchamp::demo_log() << " " << ranked[static_cast<size_t>(i)].first
                               << "=" << std::sqrt(ranked[static_cast<size_t>(i)].second.n2)
                               << " n=" << ranked[static_cast<size_t>(i)].second.n;
  }
  gaudi::duchamp::demo_log() << std::endl;
  gaudi::duchamp::demo_log().flush();
}

/// Mean diagonal of `A^T A` that each constraint name adds over the columns
/// it touches, beside the mean inertia diagonal `M / h^2` of DOFs with real
/// mass. Printed when the matrix is rebuilt, so row weights can be compared.
inline void log_stiffness_by_name(const std::vector<trip> &triplets,
                                  const std::vector<p_span> &spans,
                                  const matS &Mh) {
  std::vector<int> row_span;
  for (int s = 0; s < static_cast<int>(spans.size()); ++s)
    for (index_t r = spans[s].begin; r < spans[s].end; ++r) {
      if (static_cast<index_t>(row_span.size()) <= r)
        row_span.resize(static_cast<size_t>(r) + 1, -1);
      row_span[static_cast<size_t>(r)] = s;
    }
  struct acc {
    real sum = 0;
    std::map<int, char> cols;
  };
  std::map<std::string, acc> by_name;
  for (const trip &t : triplets) {
    if (t.row() < 0 || static_cast<size_t>(t.row()) >= row_span.size())
      continue;
    const int s = row_span[static_cast<size_t>(t.row())];
    if (s < 0)
      continue;
    acc &a = by_name[spans[static_cast<size_t>(s)].name];
    a.sum += t.value() * t.value();
    a.cols[t.col()] = 1;
  }
  real inertia = 0;
  int n_inertia = 0;
  for (int k = 0; k < Mh.outerSize(); ++k)
    for (matS::InnerIterator it(Mh, k); it; ++it)
      if (it.row() == it.col() && it.value() > real(1e-4)) {
        inertia += it.value();
        ++n_inertia;
      }
  gaudi::duchamp::demo_log() << "[m_solver] stiffness inertia="
                             << (n_inertia > 0 ? inertia / n_inertia : real(0));
  for (const auto &kv : by_name)
    if (!kv.second.cols.empty())
      gaudi::duchamp::demo_log()
          << " " << kv.first << "="
          << kv.second.sum / real(kv.second.cols.size());
  gaudi::duchamp::demo_log() << std::endl;
}

} // namespace detail

class projection_solver {
public:
  projection_solver() {}

  void set_mass(std::vector<sim_block::ptr> &blocks) {

    vecX M;
    for (auto &block : blocks) {
      block->map_mass(M);
    }
    std::vector<trip> triplets;
    for (int i = 0; i < M.rows(); i++) {
      real Mi = M[i] < 1e-8 ? 1e-8 : M[i];
      triplets.push_back(trip(i, i, Mi));
    }
    __M = matS(M.size(), M.size());
    __M.setFromTriplets(triplets.begin(), triplets.end());
  }

  void step(std::vector<sim_block::ptr> &blocks, const real &h = 0.01,
            const real &damping = 0.5, const int &ITS = 50,
            solver_context *rest_ctx = nullptr) {

    vecX q;
    vecX s;
    set_mass(blocks);

    for (auto &block : blocks) {
      block->map_to_x(q);
    }

    if (rest_ctx) {
      for (int i = 0; i < static_cast<int>(_constraints.size()); i++) {
        if (auto bc =
                std::dynamic_pointer_cast<block_constraint>(_constraints[i])) {
          if (!bc->rest_initialized())
            bc->init_rest(q, *rest_ctx);
        }
      }
    }

    for (auto &block : blocks) {
      block->integrate_inertia(h);
      block->map_to_x(s);
    }
    // Soler et al. Alg.1 line 6: start local/global from the momentum prediction.
    // Without this, project() runs on q^t and never sees torque-tipped frames, so
    // bend/twist springs cannot propagate twist to neighbors.
    q = s;
    for (const cache_compute::ptr &compute : _computes) {
      if (compute)
        compute->compute(q);
    }
    int Nm = q.size();

    // Projection rows are fixed for the step; Gauss–Newton rows hold J(q) at
    // a reference iterate and are refilled at the same row offsets.
    std::vector<trip> triplets;
    std::vector<trip> gn_triplets;
    std::vector<int> gn_index;
    index_t id0 = 0;
    std::vector<detail::p_span> spans;
    spans.reserve(_constraints.size());

    for (int i = 0; i < _constraints.size(); i++) {
      auto &constraint = _constraints[i];
      const index_t begin = id0;
      if (constraint->gauss_newton()) {
        gn_index.push_back(i);
        constraint->fill_J(q, id0, gn_triplets);
      } else {
        constraint->fill_A(id0, triplets);
      }
      spans.push_back(detail::p_span{constraint->name(), begin, id0});
    }
    const bool has_gn = !gn_index.empty();
    last_gn_refreshes = 0;

    matS Mh = __M;
    Mh *= 1.0 / h / h;

    const detail::linear_system_cache_key key{
        Nm,
        id0,
        _constraints.size(),
        detail::hash_triplets(triplets),
        h,
        detail::mass_diagonal_sum(__M),
    };

    const bool diag = projection_diag_enabled();
    matS A;
    m_solver *S = nullptr;

    // With Gauss–Newton rows the key covers the projection rows only: a hit
    // reuses the previous factor (and its J_ref) until the first refresh,
    // which is safe because the fixed point does not depend on J_ref.
    auto all_triplets = [&]() {
      std::vector<trip> t = triplets;
      t.insert(t.end(), gn_triplets.begin(), gn_triplets.end());
      return t;
    };
    bool cache_hit = false;
    if (_cache_valid && key == _cache_key) {
      if (diag)
        projection_diag_log("linear_system cache hit");
      A = _cache_A;
      S = _cache_solver.get();
      cache_hit = true;
    } else {
      if (diag)
        projection_diag_log("linear_system cache miss");
      A = matS(id0, Nm);
      if (diag) {
        projection_diag_log("setFromTriplets...");
        std::cout << "Nm: " << Nm << std::endl;
      }
      const std::vector<trip> all = all_triplets();
      A.setFromTriplets(all.begin(), all.end());
      detail::log_stiffness_by_name(all, spans, Mh);
      if (diag)
        log_matrix_report("A", diagnose_sparse(A));
      matS AtA = A.transpose() * A;
      if (diag)
        log_matrix_report("AtA", diagnose_sparse(AtA));
      matS MAtA = Mh + AtA;
      if (diag)
        log_matrix_report("MAtA", diagnose_sparse(MAtA));

      _cache_solver.reset(); // free the old factor before building the new one
      _cache_solver = std::make_unique<m_solver>(MAtA);
      _cache_A = A;
      _cache_AtA = AtA;
      _cache_MAtA = MAtA;
      _cache_key = key;
      _cache_valid = true;
      S = _cache_solver.get();
    }

    if (!S->success()) {
      std::cout << "Solve failed, attempting normalization" << std::endl;
      matS MAtA = _cache_MAtA;
      Eigen::SparseMatrix<double> I(MAtA.rows(), MAtA.cols());
      MAtA += 1e-6 * I;
      _cache_solver.reset(); // free the old factor before building the new one
      _cache_solver = std::make_unique<m_solver>(MAtA);
      _cache_MAtA = MAtA;
      S = _cache_solver.get();
    }
    if (!S->success()) {
      std::cout << "Solve failed" << std::endl;
      return;
    }

    // Refill the Gauss–Newton rows at `qk` and refactor Mh + A^T A.
    auto refresh_gn = [&](const vecX &qk) {
      gn_triplets.clear();
      for (const int i : gn_index) {
        index_t r0 = spans[static_cast<size_t>(i)].begin;
        _constraints[i]->fill_J(qk, r0, gn_triplets);
      }
      const std::vector<trip> all = all_triplets();
      A = matS(id0, Nm);
      A.setFromTriplets(all.begin(), all.end());
      matS MAtA = Mh + matS(A.transpose() * A);
      _cache_solver.reset(); // free the old factor before building the new one
      _cache_solver = std::make_unique<m_solver>(MAtA);
      _cache_A = A;
      _cache_MAtA = MAtA;
      S = _cache_solver.get();
      ++last_gn_refreshes;
    };

    vecX p = vecX::Zero(id0);
    vecX q0 = q;
    bool announced_bad = false;
    // Gauss–Newton bookkeeping: the iterate the last step started from, its
    // energy, and whether that step used a factor fresh at that iterate.
    vecX q_prev = q;
    real E_prev = std::numeric_limits<real>::infinity();
    bool last_step_fresh = !cache_hit;
    bool refreshed_by_backtrack = false;
    auto project_all = [&](const vecX &qk) {
      p.setZero();
      for (int i = 0; i < _constraints.size(); i++) {
        auto &constraint = _constraints[i];
        if (!constraint->gauss_newton())
          constraint->project(qk, p);
      }
    };
    // p rows of GN spans are J_ref q, so A^T p carries J_ref^T J_ref q.
    auto gn_energy = [&](const vecX &qk) {
      const vecX Aq = A * qk;
      for (const int i : gn_index)
        for (index_t r = spans[i].begin; r < spans[i].end; ++r)
          p[r] = Aq[r];
      const vecX dq_s = qk - s;
      real E = dq_s.dot(Mh * dq_s) + (Aq - p).squaredNorm();
      for (const int i : gn_index)
        E += _constraints[i]->residual_sq(qk);
      return E;
    };
    auto energy_at = [&](const vecX &qk) {
      project_all(qk);
      return gn_energy(qk);
    };
    // Energy by constraint name at `qk` ("inertia" for |q - s|^2_Mh), and the
    // single span (constraint instance) with the largest share per name.
    struct part {
      real E = 0;
      real worst = 0;
      int worst_i = -1;
    };
    auto energy_parts = [&](const vecX &qk) {
      project_all(qk);
      const vecX Aq = A * qk;
      std::map<std::string, part> out;
      const vecX dq_s = qk - s;
      out["inertia"].E = dq_s.dot(Mh * dq_s);
      for (int i = 0; i < static_cast<int>(_constraints.size()); ++i) {
        const auto &c = _constraints[static_cast<size_t>(i)];
        const detail::p_span &sp = spans[static_cast<size_t>(i)];
        real e = 0;
        if (c->gauss_newton())
          e = c->residual_sq(qk);
        else if (sp.end > sp.begin)
          e = (Aq.segment(sp.begin, sp.end - sp.begin) - p.segment(sp.begin, sp.end - sp.begin))
                  .squaredNorm();
        part &pt = out[sp.name];
        pt.E += e;
        if (e > pt.worst) {
          pt.worst = e;
          pt.worst_i = i;
        }
      }
      return out;
    };
    // Which constraints the step from `qa` to `qb` raised, largest first.
    auto log_ascent = [&](int k, const vecX &qa, const vecX &qb) {
      const auto a = energy_parts(qa);
      const auto b = energy_parts(qb);
      std::vector<std::pair<std::string, real>> rise;
      for (const auto &[name, pb] : b)
        rise.push_back({name, pb.E - (a.count(name) ? a.at(name).E : real(0))});
      std::sort(rise.begin(), rise.end(), [](const auto &u, const auto &v) { return u.second > v.second; });
      auto &log = gaudi::duchamp::demo_log();
      log << "[m_solver] ascent at iter " << (k + 1) << ":";
      for (size_t r = 0; r < rise.size() && r < 4; ++r) {
        const part &pb = b.at(rise[r].first);
        log << " " << rise[r].first << " " << (rise[r].second >= 0 ? "+" : "") << rise[r].second
            << " (" << (a.count(rise[r].first) ? a.at(rise[r].first).E : real(0)) << " -> " << pb.E;
        if (pb.worst_i >= 0 && !_constraints[static_cast<size_t>(pb.worst_i)]->gauss_newton())
          log << ", worst span " << pb.worst_i << " " << pb.worst;
        log << ")";
      }
      log << std::endl;
    };
    const real E_first = energy_at(q);
    if (run_variational_audit) {
      auto &log = gaudi::duchamp::demo_log();
      log << "[m_solver] variational audit (relative |FD - envelope|, flagged > 1e-3):" << std::endl;
      for (const auto &row : audit_variational(_constraints, q)) {
        log << "  " << row.name << (row.gauss_newton ? " [GN]" : "") << " n=" << row.count
            << " flagged=" << row.flagged << " median=" << row.median << " max=" << row.max;
        if (!row.instances.empty()) {
          log << " instances:";
          for (real m : row.instances)
            log << " " << m;
        }
        log << std::endl;
      }
    }
    // The buckling preconditioner: modes of H v = lambda P v at the step
    // start, H by central differences of the half gradient r = P q - b(q).
    const bool modal_on = preconditioner == step_preconditioner::modal && modal.modes > 0;
    auto residual_at = [&](const vecX &qk) {
      project_all(qk);
      const vecX Aq = A * qk;
      for (const int i : gn_index)
        for (index_t r = spans[i].begin; r < spans[i].end; ++r)
          p[r] = Aq[r];
      vecX g = vecX::Zero(Nm);
      for (const int i : gn_index)
        _constraints[i]->add_gradient(qk, g);
      return vecX(Mh * (qk - s) + A.transpose() * (Aq - p) + g);
    };
    if (modal_on) {
      if (_modes_n != Nm) {
        _modes.reset();
        _modes_n = Nm;
      }
      const vecX q_lin = q;
      const real eps = modal.fd_eps * std::max(real(1), q_lin.lpNorm<Eigen::Infinity>());
      const hessian_vector_fn hv = [&](const vecX &v) {
        const real vn = v.lpNorm<Eigen::Infinity>();
        if (!(vn > real(0)))
          return vecX(vecX::Zero(v.size()));
        const real e = eps / vn;
        return vecX((residual_at(q_lin + e * v) - residual_at(q_lin - e * v)) / (real(2) * e));
      };
      const bool ok = _modes.refresh(*S, hv, modal);
      auto &log = gaudi::duchamp::demo_log();
      log << "[m_solver] modal: " << (ok ? "lambda" : "no modes");
      if (ok)
        for (Eigen::Index i = 0; i < _modes.lambdas().size(); ++i)
          log << " " << _modes.lambdas()[i];
      log << std::endl;
    }
    int its_used = 0;
    real dq_last = 0;
    bool stopped_on_ascent = false;
    // The lowest-energy iterate seen; the step ends there, never worse.
    vecX q_best = q;
    real E_best = E_first;
    for (int k = 0; k < ITS; k++) {
      its_used = k + 1;
      project_all(q);
      vecX g_gn;
      if (has_gn) {
        // E = |q - s|^2_Mh + |A_s q - p_s|^2 + |r|^2 at this iterate.
        real E = gn_energy(q);
        bool fresh = (k == 0 && !cache_hit) || refreshed_by_backtrack;
        refreshed_by_backtrack = false;
        const bool ascent = k > 0 && E > E_prev * (real(1) + gn_ascent_tol);
        if (ascent)
          log_ascent(k, q_prev, q);
        if (gn_on_ascent == ascent_policy::stop && ascent) {
          q = q_prev;
          stopped_on_ascent = true;
          break;
        }
        if (gn_on_ascent == ascent_policy::refresh && ascent && !last_step_fresh) {
          // The step on a stale factor went uphill: take it back and
          // retake it from the previous iterate with a fresh J.
          q = q_prev;
          refresh_gn(q);
          project_all(q);
          E = gn_energy(q);
          fresh = true;
        } else if (k > 0 && gn_refresh_every > 0 &&
                   k % gn_refresh_every == 0) {
          refresh_gn(q);
          E = gn_energy(q);
          fresh = true;
        }
        q_prev = q;
        E_prev = E;
        if (E < E_best) {
          E_best = E;
          q_best = q;
        }
        last_step_fresh = fresh;
        g_gn = vecX::Zero(Nm);
        for (const int i : gn_index)
          _constraints[i]->add_gradient(q, g_gn);
        if (k < 3 || (k + 1) % 10 == 0)
          gaudi::duchamp::demo_log()
              << "[m_solver] gn iter=" << (k + 1) << " E=" << E
              << " |g_gn|=" << g_gn.norm()
              << " refreshes=" << last_gn_refreshes << std::endl;
      }
      const bool bad = !q.allFinite() || !p.allFinite();
      const bool milestone = k < 3 || ((k + 1) % 10 == 0);
      if (milestone || (bad && !announced_bad)) {
        if (bad)
          announced_bad = true;
        detail::log_p_by_name(k, q, p, spans);
      }
      if (!p.allFinite()) {
        detail::log_p_by_name(k, q, p, spans);
        std::cerr << "[m_solver] abort: non-finite p at iter " << (k + 1)
                  << std::endl;
        assert(false && "hepworth projection produced non-finite p");
        return;
      }
      if (!q.allFinite()) {
        std::cerr << "[m_solver] abort: non-finite q at iter " << (k + 1)
                  << std::endl;
        assert(false && "hepworth projection state q is non-finite");
        return;
      }

      if (diag && (k == 0 || k % 10 == 0))
        std::cout << "k: " << k << " -pnorm: " << p.norm() << std::endl;

      vecX b = Mh * s + A.transpose() * p;
      if (has_gn)
        b -= g_gn;
      if (!b.allFinite()) {
        std::cerr << "[m_solver] abort: non-finite rhs at iter " << (k + 1)
                  << std::endl;
        assert(false && "hepworth rhs b is non-finite");
        return;
      }

      const vecX q_start = q;
      q = S->solve(b);
      if (!q.allFinite()) {
        std::cerr << "[m_solver] abort: non-finite solve at iter " << (k + 1)
                  << std::endl;
        assert(false && "hepworth LDLT solve produced non-finite q");
        return;
      }
      if (modal_on && _modes.ready()) {
        // P q_plain = b, so the residual at the iterate's start is
        // P (q_start - q_plain).
        const vecX r_prev = _cache_MAtA * (q_start - q);
        const real scale = std::sqrt(std::max(q.dot(_cache_MAtA * q), real(0)));
        q = _modes.correct(q, r_prev, energy_at, energy_at(q), scale, modal);
      }

      if (has_gn && gn_backtracks > 0) {
        // Halve the step toward the iterate it started from until the
        // energy does not rise. With no descent on a stale factor, refactor
        // there and retake the step next iterate.
        const vecX q_full = q;
        real E1 = energy_at(q);
        real a = 1;
        const real E_ok = E_prev * (real(1) + gn_ascent_tol);
        for (int t = 0; t < gn_backtracks && !(E1 <= E_ok); ++t) {
          a *= real(0.5);
          q = q_prev + a * (q_full - q_prev);
          E1 = energy_at(q);
        }
        if (!(E1 <= E_ok)) {
          log_ascent(k, q_prev, q);
          q = q_prev;
          if (last_step_fresh || gn_on_ascent == ascent_policy::stop) {
            stopped_on_ascent = true;
            break;
          }
          refresh_gn(q);
          refreshed_by_backtrack = true;
          continue;
        }
      }

      real dq = (q - q0).norm();
      dq_last = dq;
      if (dq < 1e-8)
        break;
      q0 = q;

      if (diag && k % 10 == 0)
        std::cout << "k: " << k << " -norms: "
                  << " q-q0: " << dq << std::endl;
    }

    real E_last = energy_at(q);
    bool kept_best = false;
    if (has_gn && E_best < E_last) {
      q = q_best;
      E_last = E_best;
      kept_best = true;
    }
    last_step = step_report{its_used, E_first, E_last, dq_last, q.norm(),
                            last_gn_refreshes, stopped_on_ascent, kept_best};
    if (modal_on) {
      const auto &mr = _modes.last();
      last_step.lambda_min = mr.lambda_min;
      last_step.modal_negative = mr.negative;
      last_step.modal_tried = mr.tried;
      last_step.modal_accepted = mr.accepted;
      gaudi::duchamp::demo_log() << "[m_solver] modal: lambda_min " << mr.lambda_min
                                 << " negative " << mr.negative << " corrections "
                                 << mr.accepted << "/" << mr.tried << std::endl;
    }
    gaudi::duchamp::demo_log() << "[m_solver] step: its " << its_used << " E " << E_first
                               << " -> " << last_step.E_last << " last |dq| " << dq_last
                               << " refreshes " << last_gn_refreshes
                               << (stopped_on_ascent ? " (stopped on ascent)" : "")
                               << (kept_best ? " (kept best iterate)" : "") << std::endl;

    for (auto &block : blocks) {
      block->map_from_x(q, h, damping);
    }
  }

  /// One `step`'s convergence: iterates used, energy at the first and the
  /// returned iterate, the last step size, the final `|q|`, Gauss–Newton
  /// refactors, and whether it ended on an ascent or fell back to the best
  /// iterate.
  struct step_report {
    int its = 0;
    real E_first = 0;
    real E_last = 0;
    real dq_last = 0;
    real q_norm = 0;
    int refreshes = 0;
    bool stopped_on_ascent = false;
    bool kept_best = false;
    /// The buckling preconditioner (NaN / 0 when off): the smallest
    /// generalized eigenvalue found, modes at or below `lambda_neg`, and
    /// corrections accepted of those tried.
    real lambda_min = std::numeric_limits<real>::quiet_NaN();
    int modal_negative = 0;
    int modal_tried = 0;
    int modal_accepted = 0;
  };
  step_report last_step;

  void
  set_constraints(const std::vector<projection_constraint::ptr> &constraints) {
    _cache_valid = false;
    _constraints = constraints;
  }

  void set_computes(const std::vector<cache_compute::ptr> &computes) {
    _computes = computes;
  }

  /// Gauss–Newton rows keep their factor for this many iterates (fresh
  /// gradient every iterate); 0 never refactors on a schedule.
  int gn_refresh_every = 4;
  /// When an iterate raises the energy: `refresh` takes it back and
  /// refactors at the previous iterate; `stop` takes it back and ends the
  /// step there (no refactor; the step keeps an imperfect solve).
  enum class ascent_policy { refresh, stop };
  ascent_policy gn_on_ascent = ascent_policy::refresh;
  /// Relative rise tolerated before an iterate counts as an ascent: a stale
  /// factor's model minimum sits a few percent off the true one, so small
  /// rises are wobble, not divergence. The step still returns its best
  /// iterate.
  real gn_ascent_tol = 0;
  /// Run `audit_variational` on the constraints at the start of each step
  /// and log the table.
  bool run_variational_audit = false;
  /// The plain factored step (`ldlt`) or the same with the low-mode
  /// correction (`modal`, the buckling preconditioner, `modal_corrector`).
  step_preconditioner preconditioner = step_preconditioner::ldlt;
  modal_config modal;
  /// Step halvings allowed when an iterate raises the energy (0: full
  /// steps, the projective-dynamics default).
  int gn_backtracks = 0;
  int last_gn_refreshes = 0;

  matS __M;
  std::vector<projection_constraint::ptr> _constraints;
  std::vector<cache_compute::ptr> _computes;

private:
  bool _cache_valid = false;
  detail::linear_system_cache_key _cache_key;
  matS _cache_A;
  matS _cache_AtA;
  matS _cache_MAtA;
  std::unique_ptr<m_solver> _cache_solver;
  modal_corrector _modes;
  index_t _modes_n = -1;
}; // class projection_solver
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
