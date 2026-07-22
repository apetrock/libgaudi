#ifndef GAUDI_ALBERS_MEDIAL_SHAPE_SEARCH_HPP
#define GAUDI_ALBERS_MEDIAL_SHAPE_SEARCH_HPP

#include <algorithm>
#include <cmath>
#include <limits>

#include "gaudi/albers/darboux_medial_geometry.hpp"

// Trust-region boundary handling: when a Newton step exits the D<0 interior,
// how do we return to the trust region?
//
//   Default (gradient descent on D): walk back along -dD/dt = -(g·dir) until
//   D < 0. The Darboux field D is mostly convex in the exterior, so dD/dt > 0
//   out there and this monotonically re-enters the trust region using the
//   actual field gradient — adaptive to steepness, guaranteed to converge
//   under convexity.
//
//   GAUDI_MEDIAL_USE_REFLECTOR: use the older 1D Nelder-Mead-style mirror flip
//   across t_cap (t <- 2*t_cap - t) instead. Kept for comparison/A-B.
//
// #define GAUDI_MEDIAL_USE_REFLECTOR

namespace gaudi {
namespace albers {

struct medial_shape_search_params {
  int max_iters = 24;
  real tol = 1e-6;
  real min_denom = 1e-14;
  real lm_lambda = 1e-8;
  real grad_scale = 0.1;
  int stall_iters = 3;
  int bracket_samples = 200;
  real eps = 1e-12;
  // Brake on the Newton step (|dt| capped). Default uncapped; the trust-region
  // return handles overshoot via D<0 rather than relying on tiny steps.
  real max_step = std::numeric_limits<real>::infinity();
  // Step scale for the gradient-descent-on-D trust-region return (active when
  // GAUDI_MEDIAL_USE_REFLECTOR is NOT defined). Scales the dD/dt step; the
  // return step is additionally bounded to 0.25*t_cap for robustness.
  real trust_return_scale = 1.0;
};

struct medial_shape_search_result {
  vec3 center_local = vec3::Zero();
  vec3 center_world = vec3::Zero();
  real t = 0.0;
  real energy = std::numeric_limits<real>::infinity();
  real energy_prime = std::numeric_limits<real>::infinity();
  real travel = 0.0;
  bool converged = false;
};

inline real clamp_t(real t, real t_max) {
  // Inward-only march: the ray dir is aligned inward, so t >= 0 is the
  // interior side and t < 0 is the outward quartic artifact.
  return std::max(real(0.0), std::min(t_max, t));
}

inline bool eval_energy_safe(const darboux_geometry_bundle &geom,
                             const ray_line_bundle &line, real t, real eps,
                             real &E, real &E_p, real &E_pp) {
  const vec3 g = geom.G(line.f(t));
  if (!g.allFinite() || g.norm() < 1e-14) {
    return false;
  }
  eval_medial_energy_at_t(geom, line, t, E, E_p, E_pp, eps);
  return std::isfinite(E) && std::isfinite(E_p) && std::isfinite(E_pp);
}

inline real bracket_min_energy(const darboux_geometry_bundle &geom,
                               const ray_line_bundle &line, real t_min,
                               real t_max, int n_samples, real eps,
                               real *best_e_out = nullptr) {
  real best_t = t_min;
  real best_e = std::numeric_limits<real>::infinity();
  for (int i = 0; i <= n_samples; ++i) {
    const real t = t_min + (t_max - t_min) * real(i) / real(n_samples);
    real E, E_p, E_pp;
    if (!eval_energy_safe(geom, line, t, eps, E, E_p, E_pp)) {
      continue;
    }
    if (E < best_e) {
      best_e = E;
      best_t = t;
    }
  }
  if (best_e_out != nullptr) {
    *best_e_out = best_e;
  }
  return best_t;
}

inline medial_shape_search_result
search_medial_along_ray(const vec14 &Q, const vec3 &foot_world,
                        const vec3 &dir_world, real t_max,
                        const medial_shape_search_params &params = {}) {
  medial_shape_search_result out;
  if (dir_world.norm() < 1e-12 || t_max <= 0.0) {
    return out;
  }

  const darboux_geometry_bundle geom(Q);
  const ray_line_bundle line(vec3::Zero(), dir_world.normalized());

  // D(x) along the ray — the implicit-surface field. The trust region is the
  // FIRST D<0 segment from the surface: [0, t_surf], where t_surf is the first
  // D=0 crossing (the opposite surface). Beyond t_surf the quartic extrapolates
  // and D oscillates — it goes positive (exterior artifact) and then NEGATIVE
  // again in the deep far field (quartic re-flip). So "D<0" alone is NOT a valid
  // trust region; we must cap at the first crossing.
  auto D_at = [&](real tt) { return eval_darboux(Q, line.f(tt)); };

  // Find the FIRST D=0 crossing in (0, t_max] by scanning outward from t=0.
  // We cannot just check D(t_max): the quartic extrapolation re-flips D
  // negative in the deep far field, so D(t_max) may be <0 even though D crossed
  // 0 (went exterior) and came back. The trust region is the FIRST D<0 segment,
  // [0, t_surf], bounded by the opposite surface.
  real t_surf = t_max;
  {
    const int N = 256;
    real t_prev = real(0.0);
    real D_prev = D_at(real(0.0));
    for (int i = 1; i <= N; ++i) {
      const real t_i = t_max * real(i) / real(N);
      const real D_i = D_at(t_i);
      if (D_prev < real(0.0) && D_i >= real(0.0)) {
        real a = t_prev, b = t_i;
        for (int k = 0; k < 48; ++k) {
          const real m = real(0.5) * (a + b);
          if (D_at(m) < real(0.0)) {
            a = m;
          } else {
            b = m;
          }
        }
        t_surf = real(0.5) * (a + b);
        break;
      }
      t_prev = t_i;
      D_prev = D_i;
    }
  }
  const real t_cap = std::min(t_max, t_surf); // hard trust-region bound

  // Brake tied to the trust region: no single step may exceed half the interior
  // segment. Prevents a giant first Newton leap (from a small E'') from clearing
  // the ridge and bouncing around the boundary.
  const real max_step_eff =
      std::isfinite(params.max_step) && params.max_step > real(0.0)
          ? std::min(params.max_step, real(0.5) * t_cap)
          : real(0.5) * t_cap;

  auto damped_step = [&](real t, real e_p, real e_pp) {
    real dt;
    if (std::abs(e_pp) >= params.min_denom) {
      dt = -e_p / e_pp;
    } else if (std::abs(e_p) < real(1e-14)) {
      return t;
    } else {
      dt = -params.grad_scale * e_p;
    }
    if (std::abs(dt) > max_step_eff) {
      dt = std::copysign(max_step_eff, dt);
    }
    return t + dt;
  };

  // Trust-region return: when a Newton step exits the D<0 interior, get back
  // inside. By default we descend D (walk along -dD/dt) — monotone under
  // convexity and adaptive to field steepness. The Nelder-Mead mirror flip is
  // available behind GAUDI_MEDIAL_USE_REFLECTOR for comparison.
  auto return_to_trust = [&](real t_out) -> real {
    real t = std::min(t_out, t_max);
    const real step_cap = real(0.25) * t_cap;
    for (int k = 0; k < 32 && D_at(t) >= real(0.0); ++k) {
      const vec3 g = geom.G(line.f(t));
      if (!g.allFinite() || g.norm() < real(1e-14)) {
        break;
      }
      const real dDdt = g.dot(line.dir); // = grad D . dir
      if (std::abs(dDdt) < real(1e-14)) {
        break; // at a D-critical point; can't descend further
      }
      real dt = -params.trust_return_scale * dDdt;
      if (std::abs(dt) > step_cap) {
        dt = std::copysign(step_cap, dt);
      }
      t = std::max(real(0.0), std::min(t_cap, t + dt));
    }
    return t;
  };

  auto trust_step = [&](real t_prop) -> real {
    t_prop = std::max(real(0.0), std::min(t_max, t_prop));
#ifdef GAUDI_MEDIAL_USE_REFLECTOR
    // Nelder-Mead-style mirror flip across the opposite surface (t_cap).
    if (t_prop > t_cap) {
      t_prop = real(2.0) * t_cap - t_prop;
    }
    if (t_prop < real(0.0)) {
      t_prop = real(0.0); // inward-only
    }
    if (t_prop > t_cap) {
      t_prop = t_cap; // reflection left it past the boundary — clamp
    }
    return t_prop;
#else
    // Accept if inside the trust region; otherwise descend D back inside.
    if (D_at(t_prop) < real(0.0)) {
      return t_prop;
    }
    return return_to_trust(t_prop);
#endif
  };

  // Deepest-interior seed: trapezoidal centroid of -D over [0, t_cap]. The
  // medial ridge is where D is most negative (center of the maximal inscribed
  // ball); the centroid is a robust, derivative-free estimate of that point —
  // well-defined even when E = 1/||W||^2 is too shallow to drive Newton from
  // the surface. Used both as Newton's seed and as a fallback estimate.
  auto centroid_seed = [&]() -> real {
    const int M = 96;
    real W = real(0.0), Mt = real(0.0);
    real t_prev = real(0.0);
    real w_prev = std::max(real(0.0), -D_at(real(0.0)));
    for (int i = 1; i <= M; ++i) {
      const real t_i = t_cap * real(i) / real(M);
      const real w_i = std::max(real(0.0), -D_at(t_i));
      const real dt = t_i - t_prev;
      W += real(0.5) * (w_prev + w_i) * dt;
      Mt += real(0.5) * (t_prev * w_prev + t_i * w_i) * dt;
      t_prev = t_i;
      w_prev = w_i;
    }
    if (W <= real(1e-18)) {
      return real(0.5) * t_cap;
    }
    return Mt / W;
  };
  const real t_seed = centroid_seed();

  // Seed Newton at the deepest-interior point (not t=0). Starting in the ridge
  // basin avoids the shallow/noisy E' regime near the surface that stalls
  // Newton, and the cap keeps it from escaping to the exterior quartic artifact.
  real t = t_seed;
  real E = std::numeric_limits<real>::infinity();
  real E_p = std::numeric_limits<real>::infinity();
  real E_pp = 0.0;
  int stalls = 0;

  for (int iter = 0; iter < params.max_iters; ++iter) {
    if (!eval_energy_safe(geom, line, t, params.eps, E, E_p, E_pp)) {
      break;
    }
    if (std::abs(E_p) < params.tol) {
      out.converged = true;
      break;
    }

    const real t_new = trust_step(damped_step(t, E_p, E_pp));
    if (std::abs(t_new - t) < params.tol * std::max<real>(1.0, std::abs(t))) {
      ++stalls;
    } else {
      stalls = 0;
    }
    if (stalls >= params.stall_iters) {
      break;
    }
    t = t_new;
  }

  // Accept if Newton converged inside the trust region. If it didn't, fall back
  // to the deepest-interior seed (a valid medial estimate whenever D<0 there).
  if (eval_energy_safe(geom, line, t, params.eps, E, E_p, E_pp) &&
      std::abs(E_p) < params.tol && t >= real(0.0) && t <= t_cap) {
    out.converged = true;
  } else if (t_seed > real(0.0) && D_at(t_seed) < real(0.0)) {
    t = t_seed;
    if (eval_energy_safe(geom, line, t, params.eps, E, E_p, E_pp)) {
      out.converged = true;
    }
  } else {
    out.converged = false;
  }

  out.t = t;
  out.energy = E;
  out.energy_prime = E_p;
  out.center_local = line.f(t);
  out.center_world = foot_world + out.center_local;
  out.travel = out.center_local.norm();
  return out;
}

} // namespace albers
} // namespace gaudi

#endif // GAUDI_ALBERS_MEDIAL_SHAPE_SEARCH_HPP
