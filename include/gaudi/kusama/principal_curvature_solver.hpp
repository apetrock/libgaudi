#ifndef GAUDI_KUSAMA_PRINCIPAL_CURVATURE_SOLVER_HPP
#define GAUDI_KUSAMA_PRINCIPAL_CURVATURE_SOLVER_HPP

/// Curvature-domain shape solve, Eigensatz 2010 §3.3 and §2.4
/// (`reference/eigensatz2010-constraint-based-surface-processing.pdf`).
/// Least squares over vertex positions:
///   k_pc    A_v (kappa_i - kappa_hat_i)^2   per interior vertex, eq. (3.12)
///   k_conf  (alpha - alpha_rest)^2           per triangle corner, eq. (3.11)
///   k_areal (A_f - A_f_rest)^2               per triangle
///   k_close |x - x_rest|^2                   per vertex (the paper's 0*)
/// minimized by Levenberg–Marquardt, one iteration per `iterate` call.
/// Derivatives are the closed forms of Appendix A.1; `analytic = false`
/// switches to central differences for checking.

#include "gaudi/kusama/principal_curvature.hpp"

#include <Eigen/Sparse>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

namespace gaudi {
namespace kusama {

struct principal_curvature_weights {
  real k_pc = 100;
  real k_conf = 1;
  real k_areal = 1e-5;
  real k_close = 1e-5;
};

struct principal_curvature_report {
  real energy0 = 0;
  real energy1 = 0;
  real lambda = 0;
  real step = 0; ///< max |delta x| of the accepted step, 0 if rejected
  bool accepted = false;
  bool reused = false; ///< accepted on the stashed factor, no refactor
};

/// Seconds spent in the last `iterate`.
struct principal_curvature_timing {
  double jacobian = 0;
  double factor = 0;
  double solve = 0;
  double energy = 0;
  double total = 0;
};

inline real principal_curvature_corner_angle(const vec3 &a, const vec3 &b,
                                   const vec3 &c) {
  const vec3 u = b - a;
  const vec3 w = c - a;
  return std::atan2(u.cross(w).norm(), u.dot(w));
}

/// Net force and torque (about the mass centroid) of a per-vertex field `g`
/// (packed `3 v + c`). A rigid-invariant energy has both zero.
inline void principal_curvature_rigid_moments(const std::vector<vec3> &x,
                                    const std::vector<real> &m,
                                    const vecX &g, vec3 &force,
                                    vec3 &torque) {
  real M = 0;
  vec3 c = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    M += m[i];
    c += m[i] * x[i];
  }
  c /= M;
  force = vec3::Zero();
  torque = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    const vec3 gi = g.segment<3>(3 * i);
    force += gi;
    torque += (x[i] - c).cross(gi);
  }
}

/// Removes the rigid-motion part of `g` in the mass metric: the net force
/// (spread by mass) and the torque (`m_i omega x r_i`, `omega = I^{-1} tau`).
/// Leaves the gradient of any rigid-invariant energy unchanged.
inline void principal_curvature_remove_rigid(const std::vector<vec3> &x,
                                   const std::vector<real> &m, vecX &g) {
  real M = 0;
  vec3 c = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    M += m[i];
    c += m[i] * x[i];
  }
  if (!(M > real(0)))
    return;
  c /= M;
  vec3 F = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i)
    F += g.segment<3>(3 * i);
  for (size_t i = 0; i < x.size(); ++i)
    g.segment<3>(3 * i) -= (m[i] / M) * F;
  vec3 tau = vec3::Zero();
  mat3 I = mat3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    const vec3 r = x[i] - c;
    tau += r.cross(vec3(g.segment<3>(3 * i)));
    I += m[i] * (r.squaredNorm() * mat3::Identity() - r * r.transpose());
  }
  const vec3 omega = I.ldlt().solve(tau);
  if (!omega.allFinite())
    return;
  for (size_t i = 0; i < x.size(); ++i)
    g.segment<3>(3 * i) -= m[i] * omega.cross(x[i] - c);
}

class principal_curvature_solver {
public:
  principal_curvature_weights weights;
  real lambda = real(1e-3);
  real lambda_min = real(1e-12);
  real lambda_max = real(1e12);
  /// Closed-form Jacobian (Appendix A.1). Off: central differences.
  bool analytic = true;
  /// Keep the factored `J^T J + lambda D` between calls and step on it with
  /// a fresh gradient; refactor when that step fails to lower energy, or
  /// after `refresh_every` reused steps (a stale factor still descends, but
  /// slowly once the shape has moved).
  bool prefactor = true;
  int refresh_every = 4;
  principal_curvature_timing last_timing;

  void reset(int n, const std::vector<tri3> &tris,
             const std::vector<vec3> &x_rest) {
    _mesh = principal_curvature_mesh::build(n, tris);
    _x_rest = x_rest;
    _alpha_rest.resize(tris.size());
    _area_rest.resize(tris.size());
    for (size_t f = 0; f < tris.size(); ++f) {
      const tri3 &t = tris[f];
      for (int k = 0; k < 3; ++k)
        _alpha_rest[f][k] = principal_curvature_corner_angle(
            x_rest[t[k]], x_rest[t[(k + 1) % 3]], x_rest[t[(k + 2) % 3]]);
      _area_rest[f] = real(0.5) * principal_curvature_face_cross(x_rest, t).norm();
    }
    _vertex_area_rest.resize(static_cast<size_t>(n));
    for (int v = 0; v < n; ++v)
      _vertex_area_rest[v] = principal_curvature_vertex_area(x_rest, _mesh, v);
    real L = 0;
    int nl = 0;
    for (const tri3 &t : tris)
      for (int k = 0; k < 3; ++k, ++nl)
        L += (x_rest[t[k]] - x_rest[t[(k + 1) % 3]]).norm();
    _eps = real(1e-6) * (nl > 0 ? L / nl : real(1));
    _factored = false;
    _reused = 0;
    _pattern_nnz = -1;
  }

  /// Keeps `lambda` and the weights; only topology and rest change.
  void reset_rest(const std::vector<vec3> &x_rest) {
    reset(_mesh.n, _mesh.tris, x_rest);
  }

  const principal_curvature_mesh &mesh() const { return _mesh; }

  /// Principal curvatures; boundary vertices get NaN.
  std::vector<vec2> principal(const std::vector<vec3> &x) const {
    std::vector<vec2> k(static_cast<size_t>(_mesh.n),
                        vec2::Constant(std::numeric_limits<real>::quiet_NaN()));
    for (int v = 0; v < _mesh.n; ++v) {
      const principal_curvatures pr = principal_curvature_at(x, _mesh, v);
      if (pr.ok)
        k[v] = vec2(pr.k1, pr.k2);
    }
    return k;
  }

  /// Sum of squared residuals. Non-finite `kappa_hat` entries are skipped.
  real energy(const std::vector<vec3> &x,
              const std::vector<vec2> &kappa_hat) const {
    vecX r;
    residuals(x, kappa_hat, r, nullptr);
    return r.squaredNorm();
  }

  /// Residuals and Jacobian at `x` (exposed for tests).
  void linearize(const std::vector<vec3> &x,
                 const std::vector<vec2> &kappa_hat, vecX &r, matS &J) const {
    std::vector<trip> jt;
    residuals(x, kappa_hat, r, &jt);
    J.resize(r.size(), 3 * _mesh.n);
    J.setFromTriplets(jt.begin(), jt.end());
  }

  /// One iteration. With a stashed factor (and `prefactor`), a step on it
  /// with the fresh gradient is tried first. Otherwise, or if that step does
  /// not lower the energy, a Levenberg–Marquardt step: refactor at the
  /// current `J`, raise `lambda` on rejection, lower it on acceptance.
  /// On rejection `x` is unchanged.
  principal_curvature_report iterate(std::vector<vec3> &x,
                           const std::vector<vec2> &kappa_hat) {
    using clock = std::chrono::steady_clock;
    auto secs = [](clock::time_point a, clock::time_point b) {
      return std::chrono::duration<double>(b - a).count();
    };
    last_timing = principal_curvature_timing{};
    const auto t_start = clock::now();
    principal_curvature_report rep;
    vecX r;
    matS J;
    linearize(x, kappa_hat, r, J);
    rep.energy0 = r.squaredNorm();
    rep.energy1 = rep.energy0;
    const vecX g = J.transpose() * r;
    const auto t_jac = clock::now();
    last_timing.jacobian = secs(t_start, t_jac);
    auto finish = [&]() {
      rep.lambda = lambda;
      last_timing.total = secs(t_start, clock::now());
      return rep;
    };
    if (!(rep.energy0 > real(0)) || !(g.norm() > real(0)))
      return finish();

    auto try_step = [&](const vecX &d) {
      std::vector<vec3> xt = x;
      real step = 0;
      for (int v = 0; v < _mesh.n; ++v) {
        const vec3 dv(d[3 * v], d[3 * v + 1], d[3 * v + 2]);
        xt[v] += dv;
        step = std::max(step, dv.norm());
      }
      const auto t0 = clock::now();
      const real e1 = energy(xt, kappa_hat);
      last_timing.energy += secs(t0, clock::now());
      if (!(std::isfinite(e1) && e1 < rep.energy0))
        return false;
      x = std::move(xt);
      rep.energy1 = e1;
      rep.step = step;
      rep.accepted = true;
      return true;
    };

    if (prefactor && _factored && _reused < refresh_every) {
      const auto t0 = clock::now();
      const vecX d = _ldlt.solve(-g);
      last_timing.solve += secs(t0, clock::now());
      if (d.allFinite() && try_step(d)) {
        rep.reused = true;
        ++_reused;
        return finish();
      }
    }
    _reused = 0;

    const matS JtJ = J.transpose() * J;
    const int nx = 3 * _mesh.n;
    real dmax = 0;
    for (int i = 0; i < nx; ++i)
      dmax = std::max(dmax, JtJ.coeff(i, i));
    const real floor = std::max(real(1e-12), real(1e-9) * dmax);
    _factored = false;
    for (int attempt = 0; attempt < 8; ++attempt) {
      matS A = JtJ;
      for (int i = 0; i < nx; ++i)
        A.coeffRef(i, i) += lambda * std::max(JtJ.coeff(i, i), floor);
      const auto t0 = clock::now();
      if (A.nonZeros() != _pattern_nnz) {
        _ldlt.analyzePattern(A);
        _pattern_nnz = A.nonZeros();
      }
      _ldlt.factorize(A);
      last_timing.factor += secs(t0, clock::now());
      if (_ldlt.info() != Eigen::Success) {
        lambda = std::min(lambda * real(4), lambda_max);
        continue;
      }
      const auto t1 = clock::now();
      const vecX d = _ldlt.solve(-g);
      last_timing.solve += secs(t1, clock::now());
      if (!d.allFinite()) {
        lambda = std::min(lambda * real(4), lambda_max);
        continue;
      }
      if (try_step(d)) {
        _factored = true;
        lambda = std::max(lambda / real(3), lambda_min);
        break;
      }
      lambda = std::min(lambda * real(4), lambda_max);
    }
    return finish();
  }

private:
  /// Stacks every residual; with `jt`, also the Jacobian triplets.
  void residuals(const std::vector<vec3> &x,
                 const std::vector<vec2> &kappa_hat, vecX &r,
                 std::vector<trip> *jt) const {
    std::vector<real> rv;
    rv.reserve(static_cast<size_t>(2 * _mesh.n + 4 * _mesh.tris.size() +
                                   3 * _mesh.n));
    std::vector<vec3> xs;
    if (jt && !analytic)
      xs = x;

    const real s_pc = std::sqrt(std::max(weights.k_pc, real(0)));
    if (s_pc > real(0)) {
      Eigen::MatrixXd Jv;
      std::vector<mat3> dK;
      for (int v = 0; v < _mesh.n; ++v) {
        if (static_cast<size_t>(v) >= kappa_hat.size() ||
            !kappa_hat[v].allFinite())
          continue;
        const principal_curvatures pr = principal_curvature_at(x, _mesh, v);
        if (!pr.ok) {
          // Keep the row count fixed for callers that refill J in place.
          rv.push_back(real(0));
          rv.push_back(real(0));
          continue;
        }
        // Rest area: a live A_v would put d(A_v)/dx (kappa - kappa_hat) in J.
        const real w = s_pc * std::sqrt(_vertex_area_rest[v]);
        const int row = static_cast<int>(rv.size());
        rv.push_back(w * (pr.k1 - kappa_hat[v][0]));
        rv.push_back(w * (pr.k2 - kappa_hat[v][1]));
        if (!jt)
          continue;
        const std::vector<int> st = principal_curvature_stencil(_mesh, v);
        const bool ok =
            analytic ? principal_curvature_tensor_jacobian(x, _mesh, v, st, dK)
                     : principal_curvature_tensor_jacobian_fd(xs, _mesh, v, st, _eps, dK);
        if (!ok || !principal_curvature_jacobian(pr, dK, Jv))
          continue;
        for (int s = 0; s < static_cast<int>(st.size()); ++s)
          for (int c = 0; c < 3; ++c)
            for (int i = 0; i < 2; ++i)
              jt->emplace_back(row + i, 3 * st[s] + c, w * Jv(i, 3 * s + c));
      }
    }

    const real s_conf = std::sqrt(std::max(weights.k_conf, real(0)));
    const real s_areal = std::sqrt(std::max(weights.k_areal, real(0)));
    if (s_conf > real(0) || s_areal > real(0)) {
      for (size_t f = 0; f < _mesh.tris.size(); ++f) {
        const tri3 &t = _mesh.tris[f];
        auto face_terms = [&](const std::vector<vec3> &p) {
          Eigen::Matrix<real, 4, 1> out;
          for (int k = 0; k < 3; ++k)
            out[k] = s_conf * (principal_curvature_corner_angle(p[t[k]], p[t[(k + 1) % 3]],
                                                      p[t[(k + 2) % 3]]) -
                               _alpha_rest[f][k]);
          out[3] = s_areal *
                   (real(0.5) * principal_curvature_face_cross(p, t).norm() - _area_rest[f]);
          return out;
        };
        const int row = static_cast<int>(rv.size());
        const Eigen::Matrix<real, 4, 1> r0 = face_terms(x);
        for (int i = 0; i < 4; ++i)
          rv.push_back(r0[i]);
        if (!jt)
          continue;
        if (analytic) {
          face_jacobian(x, t, s_conf, s_areal, row, *jt);
          continue;
        }
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 3; ++c) {
            const real x0 = xs[t[k]][c];
            xs[t[k]][c] = x0 + _eps;
            const Eigen::Matrix<real, 4, 1> rp = face_terms(xs);
            xs[t[k]][c] = x0 - _eps;
            const Eigen::Matrix<real, 4, 1> rm = face_terms(xs);
            xs[t[k]][c] = x0;
            const Eigen::Matrix<real, 4, 1> d = (rp - rm) / (real(2) * _eps);
            for (int i = 0; i < 4; ++i)
              jt->emplace_back(row + i, 3 * t[k] + c, d[i]);
          }
      }
    }

    const real s_close = std::sqrt(std::max(weights.k_close, real(0)));
    if (s_close > real(0) && _x_rest.size() == x.size()) {
      for (int v = 0; v < _mesh.n; ++v)
        for (int c = 0; c < 3; ++c) {
          if (jt)
            jt->emplace_back(static_cast<int>(rv.size()), 3 * v + c, s_close);
          rv.push_back(s_close * (x[v][c] - _x_rest[v][c]));
        }
    }

    r = Eigen::Map<const vecX>(rv.data(), static_cast<index_t>(rv.size()));
  }

  /// Corner angles (A.3) and face area (A.1) rows of face `t`, starting at
  /// `row`: three angle rows (corner `k` of `t`), then the area row.
  static void face_jacobian(const std::vector<vec3> &x, const tri3 &t,
                            real s_conf, real s_areal, int row,
                            std::vector<trip> &jt) {
    const vec3 n = principal_curvature_face_cross(x, t).normalized();
    for (int k = 0; k < 3; ++k) {
      const int ia = t[k], ib = t[(k + 1) % 3], ic = t[(k + 2) % 3];
      const vec3 u = x[ib] - x[ia];
      const vec3 w = x[ic] - x[ia];
      const vec3 gb = u.cross(n) / u.squaredNorm();
      const vec3 gc = n.cross(w) / w.squaredNorm();
      const vec3 ga = -(gb + gc);
      for (int c = 0; c < 3; ++c) {
        jt.emplace_back(row + k, 3 * ia + c, s_conf * ga[c]);
        jt.emplace_back(row + k, 3 * ib + c, s_conf * gb[c]);
        jt.emplace_back(row + k, 3 * ic + c, s_conf * gc[c]);
      }
    }
    for (int k = 0; k < 3; ++k) {
      const vec3 gA = principal_curvature_area_gradient(x, t, k);
      for (int c = 0; c < 3; ++c)
        jt.emplace_back(row + 3, 3 * t[k] + c, s_areal * gA[c]);
    }
  }

  principal_curvature_mesh _mesh;
  std::vector<vec3> _x_rest;
  std::vector<std::array<real, 3>> _alpha_rest;
  std::vector<real> _area_rest;
  std::vector<real> _vertex_area_rest;
  real _eps = real(1e-6);
  Eigen::SimplicialLDLT<matS> _ldlt;
  bool _factored = false;
  int _reused = 0;
  index_t _pattern_nnz = -1;
};

} // namespace kusama
} // namespace gaudi

#endif
