#ifndef GAUDI_HEPWORTH_BLOCK_MODAL_CORRECTOR_HPP
#define GAUDI_HEPWORTH_BLOCK_MODAL_CORRECTOR_HPP

/// Low-mode Newton correction for `projection_solver` (the buckling
/// preconditioner). Each iterate is `q_{k+1} = q_k - P^{-1} r(q_k)` with
/// `P = M/h^2 + A^T A` the factored system and `r` half the energy gradient.
/// Along a generalized eigenvector of `H v = lambda P v` (`H` the true half
/// Hessian) that step is too short where `lambda < 1` (slow, low-frequency
/// modes) and has nowhere to go where `lambda <= 0` (negative curvature, a
/// buckling direction a symmetric start cannot leave). The corrector finds
/// the `modes` smallest `lambda` and, per iterate, replaces the plain step
/// along each slow mode by its Newton step `a_i (1 - 1/lambda_i)` and steps
/// out along each negative mode, keeping the result only when its energy is
/// below the plain iterate's (Peng et al. 2018, Algorithm 1).
///
/// The modes come from the existing LDLT factor, `Pi P Pi^T = L D L^T`: with
/// `B = L D^{1/2}`, `H v = lambda P v` is `C y = lambda y` for
/// `C = B^{-1} Pi H Pi^T B^{-T}`, `v = Pi^T B^{-T} y` (P-orthonormal). Spectra
/// Lanczos takes the largest eigenvalues of `x - C x`, the smallest `lambda`.
/// `H v` is a central difference of `r`, valid because every constraint is
/// variational (`audit_variational`).

#include "gaudi/sparse_solver.h"

#include <Spectra/SymEigsSolver.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Which step `projection_solver` takes: the plain factored solve, or the
/// same with the low-mode correction.
enum class step_preconditioner { ldlt, modal };

struct modal_config {
  /// Modes corrected per step (the smallest generalized eigenvalues).
  int modes = 4;
  /// Lanczos subspace size; 0 is `max(2 modes + 1, 20)`.
  int ncv = 0;
  /// Modes with `lambda` below this get the Newton step; at or above it the
  /// plain step is already about right.
  real lambda_max = 0.9;
  /// Modes with `lambda` at or below this are negative (or vanishing)
  /// curvature, buckling directions: stepped out along, both signs tried when
  /// the gradient does not pick one. Above it the Newton step is used.
  real lambda_neg = 1e-3;
  /// Negative-mode trust radius in the P-norm, relative to the iterate's
  /// P-norm (`scale`); halved up to 12 times until the energy drops.
  real tau = 0.05;
  /// Central-difference step for `H v`, relative to `|q|_inf`.
  real fd_eps = 1e-6;
  /// Lanczos restarts.
  int max_iterations = 200;
};

/// `y = H v` (the half Hessian) at the step's linearization point.
using hessian_vector_fn = std::function<vecX(const vecX &)>;

namespace detail {

/// `y = x - C x` on Spectra's raw buffers.
class modal_operator {
public:
  using Scalar = real;
  using LDLT = Eigen::SimplicialLDLT<matS, Eigen::Lower, Eigen::AMDOrdering<int>>;

  modal_operator(const LDLT &f, const hessian_vector_fn &hv, vecX dinv_sqrt)
      : _f(f), _hv(hv), _dinv(std::move(dinv_sqrt)) {}

  Eigen::Index rows() const { return _dinv.size(); }
  Eigen::Index cols() const { return _dinv.size(); }

  /// `v = Pi^T L^{-T} D^{-1/2} y`.
  vecX to_dofs(const vecX &y) const {
    vecX z = y.cwiseProduct(_dinv);
    z = _f.matrixU().solve(z);
    return _f.permutationPinv() * z;
  }

  void perform_op(const Scalar *x_in, Scalar *y_out) const {
    Eigen::Map<const vecX> x(x_in, rows());
    Eigen::Map<vecX> y(y_out, rows());
    vecX u = _f.permutationP() * _hv(to_dofs(x));
    u = _f.matrixL().solve(u);
    y = x - u.cwiseProduct(_dinv);
  }

private:
  const LDLT &_f;
  const hessian_vector_fn &_hv;
  vecX _dinv;
};

} // namespace detail

class modal_corrector {
public:
  /// Per-step counts, for the step report.
  struct report {
    int modes = 0;
    int negative = 0;
    int tried = 0;
    int accepted = 0;
    real lambda_min = std::numeric_limits<real>::quiet_NaN();
  };

  void reset() {
    _V.resize(0, 0);
    _lambda.resize(0);
    _warm.resize(0);
    _rep = {};
  }

  /// Modes of `H v = lambda P v` from the factor of `P`. Returns false (and
  /// holds no modes) when the factor is not a direct LDLT, `D` is not
  /// positive, or no mode converges.
  bool refresh(const m_solver &S, const hessian_vector_fn &hv, const modal_config &cfg) {
    _V.resize(0, 0);
    _lambda.resize(0);
    _rep = {};
    if (!S.direct())
      return false;
    const auto &f = S.__solver;
    const vecX d = f.vectorD();
    const Eigen::Index n = d.size();
    const int k = std::min<int>(cfg.modes, static_cast<int>(n) - 2);
    if (k < 1 || !(d.minCoeff() > real(0)))
      return false;
    const vecX dinv = d.cwiseSqrt().cwiseInverse();
    detail::modal_operator op(f, hv, dinv);
    int ncv = cfg.ncv > 0 ? cfg.ncv : std::max(2 * k + 1, 20);
    ncv = std::min<int>(std::max(ncv, k + 1), static_cast<int>(n));
    Spectra::SymEigsSolver<detail::modal_operator> eigs(op, k, ncv);
    if (_warm.size() == n)
      eigs.init(_warm.data());
    else
      eigs.init();
    // Keep the converged modes even when a degenerate cluster (lambda near 1,
    // where P already matches H) stops the rest from converging.
    const Eigen::Index nconv =
        eigs.compute(Spectra::SortRule::LargestAlge, cfg.max_iterations, real(1e-8));
    if (nconv <= 0)
      return false;
    const vecX mu = eigs.eigenvalues();
    const Eigen::MatrixXd Y = eigs.eigenvectors();
    // Keep a mode only if its Rayleigh quotient v^T H v (v is P-unit) agrees
    // with its eigenvalue: Lanczos can break down on a nearly degenerate
    // operator and return garbage.
    std::vector<Eigen::Index> keep;
    std::vector<vecX> vs;
    std::vector<real> ls;
    for (Eigen::Index i = 0; i < mu.size(); ++i) {
      const real lam = real(1) - mu[i];
      const vecX v = op.to_dofs(Y.col(i));
      if (!std::isfinite(lam) || !v.allFinite())
        continue;
      const real rq = v.dot(hv(v));
      if (!(std::abs(rq - lam) <= real(1e-3) * std::max(real(1), std::abs(lam))))
        continue;
      keep.push_back(i);
      vs.push_back(v);
      ls.push_back(lam);
    }
    if (vs.empty())
      return false;
    _lambda.resize(static_cast<Eigen::Index>(ls.size()));
    _V.resize(n, static_cast<Eigen::Index>(vs.size()));
    _warm = vecX::Zero(n);
    for (size_t j = 0; j < vs.size(); ++j) {
      _lambda[static_cast<Eigen::Index>(j)] = ls[j];
      _V.col(static_cast<Eigen::Index>(j)) = vs[j];
      _warm += Y.col(keep[j]);
    }
    _rep.modes = static_cast<int>(_lambda.size());
    _rep.lambda_min = _lambda.minCoeff();
    for (Eigen::Index i = 0; i < _lambda.size(); ++i)
      _rep.negative += _lambda[i] <= cfg.lambda_neg;
    return true;
  }

  bool ready() const { return _V.cols() > 0; }

  /// The corrected iterate from the plain one `q_plain`. `r_prev` is the
  /// residual at the iterate the plain step started from (`P (q_start -
  /// q_plain)`), `energy` the solver's energy, `E_plain` its value at
  /// `q_plain`, `scale` the iterate's P-norm (the trust radius' unit).
  /// Slow modes take their Newton step; negative modes a step of the trust
  /// radius, halved until the energy drops, both signs when the gradient does
  /// not choose. Returns `q_plain` unless the correction lowers the energy.
  vecX correct(const vecX &q_plain, const vecX &r_prev,
               const std::function<real(const vecX &)> &energy, real E_plain, real scale,
               const modal_config &cfg) {
    if (!ready() || !std::isfinite(E_plain))
      return q_plain;
    vecX delta = vecX::Zero(q_plain.size());
    std::vector<std::pair<Eigen::Index, real>> neg;
    for (Eigen::Index i = 0; i < _V.cols(); ++i) {
      const real lam = _lambda[i];
      const real a = _V.col(i).dot(r_prev);
      if (lam <= cfg.lambda_neg)
        neg.push_back({i, a});
      else if (lam < cfg.lambda_max)
        delta += a * (real(1) - real(1) / lam) * _V.col(i);
    }
    if (delta.squaredNorm() == real(0) && neg.empty())
      return q_plain;
    ++_rep.tried;
    vecX best = q_plain;
    real E_best = E_plain;
    if (delta.squaredNorm() > real(0)) {
      const vecX qc = q_plain + delta;
      const real E = energy(qc);
      if (std::isfinite(E) && E < E_best) {
        E_best = E;
        best = qc;
      }
    }
    // Negative modes: from the best so far, a trust-radius step along each.
    const real t0 = cfg.tau * std::max(scale, real(1));
    for (const auto &[i, a] : neg) {
      const vecX &v = _V.col(i);
      // The descent side when the gradient has one; else both.
      const real lean = std::abs(a) > real(1e-12) * t0 ? (a > real(0) ? real(-1) : real(1)) : real(0);
      const std::vector<real> sides = lean != real(0) ? std::vector<real>{lean}
                                                      : std::vector<real>{real(1), real(-1)};
      vecX q_side = best;
      real E_side = E_best;
      for (real side : sides) {
        real t = t0;
        for (int h = 0; h < 12; ++h, t *= real(0.5)) {
          const vecX qc = best + side * t * v;
          const real E = energy(qc);
          if (std::isfinite(E) && E < E_side) {
            E_side = E;
            q_side = qc;
            break;
          }
        }
      }
      best = q_side;
      E_best = E_side;
    }
    if (E_best < E_plain)
      ++_rep.accepted;
    return best;
  }

  const report &last() const { return _rep; }
  const vecX &lambdas() const { return _lambda; }
  const Eigen::MatrixXd &modes() const { return _V; }

private:
  Eigen::MatrixXd _V;
  vecX _lambda;
  vecX _warm;
  report _rep;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
