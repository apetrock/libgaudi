#ifndef GAUDI_DUCHAMP_SOFT_TP_NEWTON_HPP
#define GAUDI_DUCHAMP_SOFT_TP_NEWTON_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/common.h"

namespace gaudi {
namespace duchamp {

/// Soft TP integrator mode (displacement proposer before CFL).
enum class soft_tp_mode {
  gradient,      ///< δ = h² · w·G  (classical explicit)
  hessian_force, ///< δ = h²·w·G + ½·H·δ₀,  δ₀ = h²·w·G
  newton         ///< CG on (H+λI)δ=−∇E, line search → δ = x'−x
};

/// Matrix-free damped Newton on soft TP; propose δ then f' = δ/h².
enum class soft_tp_newton_mode { single_iter, full_step };

struct soft_tp_newton_config {
  soft_tp_newton_mode mode = soft_tp_newton_mode::single_iter;
  int newton_iters = 4; ///< used when mode == full_step
  int cg_iters = 12;
  real lambda = 1.0e-2; ///< LM damping on (H + λI)
  real tol = 1.0e-8;    ///< stop when |∇E|_rms < tol
  real cg_tol = 1.0e-3;
  int line_search_max = 6;
  real line_search_shrink = 0.5;
  /// Cap |δ|_∞ per Newton iter (≤0 → 0.25 · lavg).
  real delta_max = 0.0;
  bool verbose = true;
};

inline real soft_tp_field_rms(const std::vector<vec3> &f) {
  real s = 0.0;
  for (const vec3 &a : f)
    s += a.squaredNorm();
  return std::sqrt(s / std::max<real>(real(1), real(f.size())));
}

inline real soft_tp_field_inf(const std::vector<vec3> &f) {
  real m = 0.0;
  for (const vec3 &a : f)
    m = std::max(m, a.norm());
  return m;
}

inline real soft_tp_field_dot(const std::vector<vec3> &a,
                              const std::vector<vec3> &b) {
  real s = 0.0;
  const size_t n = std::min(a.size(), b.size());
  for (size_t i = 0; i < n; ++i)
    s += a[i].dot(b[i]);
  return s;
}

inline void soft_tp_clamp_inf(std::vector<vec3> &v, real vmax) {
  if (!(vmax > 0.0))
    return;
  for (vec3 &a : v) {
    const real n = a.norm();
    if (n > vmax)
      a *= (vmax / n);
  }
}

/// Energy-Hessian · v (calder assembly; same scale as ∇E before force negate).
inline std::vector<vec3> soft_tp_hess_vec(asawa::rod::rod &rod,
                                          const std::vector<vec3> &v,
                                          real R_min, real tau, real p,
                                          real w) {
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  auto Hv =
      calder::tangent_point_hessian_vec_soft(rod, x, l, T, v, R_min, tau, p);
  if (w != 1.0) {
    for (auto &a : Hv)
      a *= w;
  }
  return Hv;
}

/// ∇E from soft TP (before force negate), weighted.
inline std::vector<vec3> soft_tp_grad_energy(asawa::rod::rod &rod, real R_min,
                                             real tau, real p, real w) {
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  auto g = calder::tangent_point_gradient_soft(rod, x, l, T, R_min, tau, p,
                                               /*hess_alpha=*/0.0);
  if (w != 1.0) {
    for (auto &a : g)
      a *= w;
  }
  return g;
}

/// Hessian-augmented explicit step: δ = δ₀ + ½·H·δ₀,  δ₀ = h²·(−∇E).
inline std::vector<vec3>
compute_soft_tp_hessian_displacement(asawa::rod::rod &rod, real R_min, real tau,
                                     real p, real w, real h) {
  const real h2 = std::max(h * h, real(1e-20));
  const std::vector<vec3> g = soft_tp_grad_energy(rod, R_min, tau, p, w);
  const size_t n = g.size();
  std::vector<vec3> delta0(n);
  for (size_t i = 0; i < n; ++i)
    delta0[i] = -h2 * g[i];
  const std::vector<vec3> Hdelta0 =
      soft_tp_hess_vec(rod, delta0, R_min, tau, p, w);
  std::vector<vec3> disp(n);
  for (size_t i = 0; i < n; ++i)
    disp[i] = delta0[i] + real(0.5) * Hdelta0[i];
  return disp;
}

/// Explicit descent force f = −∇E (legacy soft TP force).
inline std::vector<vec3> soft_tp_explicit_force(asawa::rod::rod &rod, real R_min,
                                                real tau, real p, real w) {
  auto g = soft_tp_grad_energy(rod, R_min, tau, p, w);
  for (auto &a : g)
    a *= real(-1.0);
  return g;
}

/// CG solve (H+λI) δ = b.
inline bool soft_tp_cg_solve(asawa::rod::rod &rod, const std::vector<vec3> &b,
                             std::vector<vec3> &delta, real R_min, real tau,
                             real p, real w, real lambda, int cg_iters,
                             real cg_tol) {
  const size_t n = b.size();
  delta.assign(n, vec3::Zero());
  std::vector<vec3> r = b;
  std::vector<vec3> pvec = r;
  real rr = soft_tp_field_dot(r, r);
  const real b2 = std::max(soft_tp_field_dot(b, b), real(1e-30));
  if (rr < cg_tol * cg_tol * b2)
    return true;

  for (int it = 0; it < cg_iters; ++it) {
    std::vector<vec3> Ap = soft_tp_hess_vec(rod, pvec, R_min, tau, p, w);
    for (size_t i = 0; i < n; ++i)
      Ap[i] += lambda * pvec[i];
    const real denom = soft_tp_field_dot(pvec, Ap);
    if (!(denom > real(1e-30)))
      return false;
    const real alpha = rr / denom;
    for (size_t i = 0; i < n; ++i) {
      delta[i] += alpha * pvec[i];
      r[i] -= alpha * Ap[i];
    }
    const real rr_new = soft_tp_field_dot(r, r);
    if (rr_new < cg_tol * cg_tol * b2)
      return true;
    const real beta = rr_new / std::max(rr, real(1e-30));
    for (size_t i = 0; i < n; ++i)
      pvec[i] = r[i] + beta * pvec[i];
    rr = rr_new;
  }
  return true;
}

/// Newton on soft TP at current x; restore x; return δ = x'−x₀.
/// On CG fail or stall, falls back to GD proposal δ = h²(−∇E).
inline std::vector<vec3>
compute_soft_tp_newton_displacement(asawa::rod::rod &rod, real R_min, real tau,
                                     real p, real w, real h,
                                     const soft_tp_newton_config &ncfg = {}) {
  const real h2 = std::max(h * h, real(1e-20));
  std::vector<vec3> &x = rod.x();
  const size_t n = x.size();
  std::vector<vec3> x0 = x;

  const real lavg = std::max(rod.lavg(), real(1e-6));
  const real dmax =
      (ncfg.delta_max > 0.0) ? ncfg.delta_max : real(0.25) * lavg;

  const int niters =
      (ncfg.mode == soft_tp_newton_mode::single_iter) ? 1 : ncfg.newton_iters;

  real lambda = std::max(ncfg.lambda, real(0.0));
  real last_g_rms = 0.0;
  real last_d_rms = 0.0;
  real last_alpha = 0.0;
  int last_cg_ok = 1;
  int n_moved = 0;

  for (int k = 0; k < niters; ++k) {
    std::vector<vec3> g = soft_tp_grad_energy(rod, R_min, tau, p, w);
    last_g_rms = soft_tp_field_rms(g);
    std::vector<vec3> b(n);
    for (size_t i = 0; i < n; ++i)
      b[i] = -g[i];
    if (last_g_rms < ncfg.tol)
      break;

    std::vector<vec3> delta;
    bool ok = soft_tp_cg_solve(rod, b, delta, R_min, tau, p, w, lambda,
                               ncfg.cg_iters, ncfg.cg_tol);
    if (!ok) {
      lambda = std::max(lambda * real(10.0), real(1e-6));
      ok = soft_tp_cg_solve(rod, b, delta, R_min, tau, p, w, lambda,
                            ncfg.cg_iters, ncfg.cg_tol);
    }
    last_cg_ok = ok ? 1 : 0;
    if (!ok) {
      // GD proposal at PD scale: δ = h²(−∇E), not −∇E.
      delta.resize(n);
      for (size_t i = 0; i < n; ++i)
        delta[i] = h2 * b[i];
    }

    soft_tp_clamp_inf(delta, dmax);
    last_d_rms = soft_tp_field_rms(delta);

    const real g0 = soft_tp_field_dot(g, g);
    real alpha = 1.0;
    const std::vector<vec3> x_base = x;
    bool moved = false;
    for (int ls = 0; ls < ncfg.line_search_max; ++ls) {
      for (size_t i = 0; i < n; ++i)
        x[i] = x_base[i] + alpha * delta[i];
      std::vector<vec3> g1 = soft_tp_grad_energy(rod, R_min, tau, p, w);
      if (soft_tp_field_dot(g1, g1) <= g0) {
        moved = true;
        last_alpha = alpha;
        break;
      }
      alpha *= ncfg.line_search_shrink;
    }
    if (!moved) {
      last_alpha = std::pow(ncfg.line_search_shrink, ncfg.line_search_max);
      for (size_t i = 0; i < n; ++i)
        x[i] = x_base[i] + last_alpha * delta[i];
      moved = true;
    }
    if (moved)
      ++n_moved;
  }

  std::vector<vec3> disp(n);
  real dx_rms = 0.0;
  for (size_t i = 0; i < n; ++i) {
    disp[i] = x[i] - x0[i];
    dx_rms += disp[i].squaredNorm();
  }
  dx_rms = std::sqrt(dx_rms / std::max<real>(real(1), real(n)));

  const bool newton_moved = dx_rms > real(1e-12) * lavg;
  if (!newton_moved) {
    x = x0;
    auto g = soft_tp_grad_energy(rod, R_min, tau, p, w);
    for (size_t i = 0; i < n; ++i)
      disp[i] = -h2 * g[i];
  }
  x = x0;

  const real disp_rms = soft_tp_field_rms(disp);
  if (ncfg.verbose) {
    std::cout << "[soft TP Newton] g_rms=" << last_g_rms
              << " d_rms=" << last_d_rms << " alpha=" << last_alpha
              << " dx_rms=" << dx_rms << " disp_rms=" << disp_rms
              << " cg_ok=" << last_cg_ok << " moved=" << n_moved << "/"
              << niters << " fallback=" << (newton_moved ? 0 : 1)
              << " mode="
              << (ncfg.mode == soft_tp_newton_mode::single_iter ? "single"
                                                                : "full")
              << " h=" << h << std::endl;
  }
  return disp;
}

/// Hessian-augmented velocity step: v = v₀ + ½·H·v₀,  v₀ = −∇E (weighted).
inline std::vector<vec3>
compute_soft_tp_hessian_velocity(asawa::rod::rod &rod, real R_min, real tau,
                                 real p, real w, real h) {
  const real h_safe = std::max(h, real(1e-20));
  const std::vector<vec3> g = soft_tp_grad_energy(rod, R_min, tau, p, w);
  const size_t n = g.size();
  std::vector<vec3> v0(n);
  for (size_t i = 0; i < n; ++i)
    v0[i] = -g[i];
  std::vector<vec3> delta0(n);
  for (size_t i = 0; i < n; ++i)
    delta0[i] = h_safe * v0[i];
  const std::vector<vec3> Hdelta0 =
      soft_tp_hess_vec(rod, delta0, R_min, tau, p, w);
  std::vector<vec3> v(n);
  for (size_t i = 0; i < n; ++i)
    v[i] = v0[i] + real(0.5) * Hdelta0[i] / h_safe;
  return v;
}

/// Legacy force slot: f' = δ/h².
inline std::vector<vec3>
compute_soft_tp_newton_force(asawa::rod::rod &rod, real R_min, real tau, real p,
                             real w, real h,
                             const soft_tp_newton_config &ncfg = {}) {
  const real h2 = std::max(h * h, real(1e-20));
  std::vector<vec3> disp =
      compute_soft_tp_newton_displacement(rod, R_min, tau, p, w, h, ncfg);
  for (auto &a : disp)
    a /= h2;
  return disp;
}

/// Velocity slot: v = δ/h.
inline std::vector<vec3>
compute_soft_tp_newton_velocity(asawa::rod::rod &rod, real R_min, real tau,
                                real p, real w, real h,
                                const soft_tp_newton_config &ncfg = {}) {
  const real h_safe = std::max(h, real(1e-20));
  std::vector<vec3> disp =
      compute_soft_tp_newton_displacement(rod, R_min, tau, p, w, h, ncfg);
  for (auto &a : disp)
    a /= h_safe;
  return disp;
}

} // namespace duchamp
} // namespace gaudi

#endif // GAUDI_DUCHAMP_SOFT_TP_NEWTON_HPP
