#ifndef __CALDER_SHAPE_OPERATOR_WEIGHTS__
#define __CALDER_SHAPE_OPERATOR_WEIGHTS__

#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/calder/weight_functions.hpp"
#include "gaudi/common.h"

#include <cmath>
#include <vector>

namespace gaudi {
namespace calder {

/// Anisotropic weight from a principal frame: inv_dist × exp(-M_W).
/// Kept for A/B; not the hybrid stage-1 default.
inline real aniso_inv_dist_weight(const vec3 &dp, real l0, real p,
                                  const albers::principal_curvature_frame &pc,
                                  real eps = 1e-3) {
  if (!pc.valid) {
    return calc_inv_dist(dp, l0, p);
  }
  const real M = calc_mahalanobis_curvature(dp, pc.k_min, pc.k_max, pc.e_min,
                                            pc.e_max, pc.n, eps);
  return calc_inv_dist(dp, l0, p) * std::exp(-M);
}

/// One-sided jet Green double-layer: max(0, dp·n_ref) / (|dp|³ + l0³).
/// n_ref = Gi(x_ij) → Gi(0) → Ni (always a half-space).
inline real jet_green_double_layer_weight(const vec3 &dp, real l0,
                                          const asawa::shell::face_height_jet &jet,
                                          const vec3 &Ni) {
  const real dist = dp.norm();
  if (dist < 1e-12) {
    return 0.0;
  }

  vec3 n_ref = vec3::Zero();
  if (jet.valid) {
    const vec3 g = jet.grad_at_rel(dp);
    if (g.allFinite() && g.norm() > 1e-10) {
      n_ref = g.normalized();
    } else {
      const vec3 g0 = jet.grad_at_rel(vec3::Zero());
      if (g0.allFinite() && g0.norm() > 1e-10) {
        n_ref = g0.normalized();
      }
    }
  }
  if (n_ref.squaredNorm() < 1e-20) {
    if (Ni.norm() > 1e-12) {
      n_ref = Ni.normalized();
    } else {
      return 0.0;
    }
  }

  const real flux = dp.dot(n_ref);
  if (flux <= 0.0) {
    return 0.0;
  }
  const real lp = std::pow(std::max(l0, real(1e-12)), 3.0);
  const real dist3 = dist * dist * dist;
  return flux / (dist3 + lp);
}

/// Soft-convex × Gaussian(σ) × torus-distance falloff (deferred / A/B).
/// Note: uses world query+dp; prefer foot-jet torus_inv_normal_align_weight.
inline real torus_soft_gaussian_weight(const vec3 &dp, const vec3 &Nj,
                                       real sigma, const vec3 &query,
                                       const albers::osculating_torus &T,
                                       real torus_p = 2.0) {
  const real dist = dp.norm();
  if (dist < 1e-12) {
    return 0.0;
  }
  const real convexity = std::max(real(0.0), dp.normalized().dot(Nj));
  const real g = calc_gaussian(dp, std::max(sigma, real(1e-12)));
  const vec3 pj = query + dp;
  const real dT = albers::torus_abs_distance(pj, T);
  const real wT =
      1.0 / (std::pow(dT, torus_p) +
             std::pow(std::max(sigma, real(1e-12)), torus_p));
  return convexity * g * wT;
}

/// Soft Green distance weight: max(0, dp·n̂) / (|dp|³ + l0³).
inline real green_soft_dist_weight(const vec3 &dp, real l0, const vec3 &n_ref) {
  const real dist = dp.norm();
  if (dist < 1e-12 || !(n_ref.norm() > 1e-12)) {
    return 0.0;
  }
  const real flux = dp.dot(n_ref.normalized());
  if (flux <= 0.0) {
    return 0.0;
  }
  const real lp = std::pow(std::max(l0, real(1e-12)), 3.0);
  return flux / (dist * dist * dist + lp);
}

/// True when PAT torus is unusable as a stage-2 prior (invalid or κ_lo≈0 → R~1/eps).
inline bool osculating_torus_degenerate(const albers::osculating_torus &T,
                                       real l0) {
  if (!T.valid || !std::isfinite(T.R) || !std::isfinite(T.r) || T.r < 1e-12) {
    return true;
  }
  const real scale = std::max(T.r, std::max(l0, real(1e-12)));
  return T.R > real(1e6) * scale;
}

/// Conservative soft-convex inv_dist (same gate as soft_inv_convex_weight).
inline real soft_inv_convex_kernel(const vec3 &dp, real l0, real p,
                                   const vec3 &Nj) {
  const real dist = dp.norm();
  if (dist < 1e-12 || !(Nj.norm() > 1e-12)) {
    return 0.0;
  }
  const real convexity = std::max(real(0.0), dp.normalized().dot(Nj.normalized()));
  return convexity * calc_inv_dist(dp, l0, p);
}

/// Soft-convex × foot-jet torus normal-align × radius Gaussian:
///   w = max(0, dp̂·Nⱼ) * max(0, Nⱼ·∇T̂(dp)) * G(dp; σ)
/// Degenerate / invalid torus drops the align factor → soft-convex × G(σ)
/// (same backbone as adaptive gaussian; torus is a prior multiplier only).
inline real soft_torus_align_gaussian_weight(const vec3 &dp, real l0,
                                             const vec3 &Nj,
                                             const albers::osculating_torus &T,
                                             real sigma) {
  const real dist = dp.norm();
  if (dist < 1e-12 || !(Nj.norm() > 1e-12)) {
    return 0.0;
  }
  const real convexity =
      std::max(real(0.0), dp.normalized().dot(Nj.normalized()));
  if (convexity <= 0.0) {
    return 0.0;
  }
  const real sig = std::max(sigma, real(1e-12));
  const real g = calc_gaussian(dp, sig);

  if (osculating_torus_degenerate(T, l0)) {
    return convexity * g;
  }
  const vec3 Gt = albers::torus_sdf_grad(dp, T);
  if (!Gt.allFinite() || Gt.norm() < 1e-10) {
    return convexity * g;
  }
  const real align =
      std::max(real(0.0), Nj.normalized().dot(Gt.normalized()));
  if (align <= 0.0) {
    return 0.0;
  }
  return convexity * align * g;
}

/// Legacy: torus-align × G without soft-convex. Prefer soft_torus_align_gaussian_weight.
inline real torus_normal_align_gaussian_weight(const vec3 &dp, real l0,
                                               const vec3 &Nj,
                                               const albers::osculating_torus &T,
                                               real sigma) {
  const real dist = dp.norm();
  if (dist < 1e-12) {
    return 0.0;
  }
  const real sig = std::max(sigma, real(1e-12));
  const real g = calc_gaussian(dp, sig);

  if (osculating_torus_degenerate(T, l0) || !(Nj.norm() > 1e-12)) {
    return g;
  }

  const vec3 G = albers::torus_sdf_grad(dp, T);
  if (!G.allFinite() || G.norm() < 1e-10) {
    return g;
  }
  const real align =
      std::max(real(0.0), Nj.normalized().dot(G.normalized()));
  if (align <= 0.0) {
    return 0.0;
  }
  return align * g;
}

/// Legacy name: Green soft-dist form of Nj·∇T̂ (no Gaussian). Prefer
/// torus_normal_align_gaussian_weight for stage-2.
inline real torus_inv_normal_align_weight(const vec3 &dp, real l0, real p,
                                          const vec3 &Nj,
                                          const albers::osculating_torus &T) {
  (void)p;
  const real dist = dp.norm();
  if (dist < 1e-12 || !(Nj.norm() > 1e-12)) {
    return 0.0;
  }

  if (osculating_torus_degenerate(T, l0)) {
    if (T.valid) {
      const vec3 G = albers::torus_sdf_grad(dp, T);
      if (G.allFinite() && G.norm() > 1e-10) {
        return green_soft_dist_weight(dp, l0, G);
      }
    }
    return green_soft_dist_weight(dp, l0, Nj);
  }

  const vec3 G = albers::torus_sdf_grad(dp, T);
  if (!G.allFinite() || G.norm() < 1e-10) {
    return green_soft_dist_weight(dp, l0, Nj);
  }
  const real align =
      std::max(real(0.0), Nj.normalized().dot(G.normalized()));
  if (align <= 0.0) {
    return 0.0;
  }
  const real lp = std::pow(std::max(l0, real(1e-12)), 3.0);
  return align / (dist * dist * dist + lp);
}

/// Foot-jet torus Green double-layer:
///   w = max(0, dp · ∇T(dp).normalized()) / (|dp|³ + l0³)
/// Invalid / degenerate torus → Green soft dist with Nj.
inline real torus_green_dp_weight(const vec3 &dp, real l0, real p,
                                  const vec3 &Nj,
                                  const albers::osculating_torus &T) {
  (void)p;
  if (osculating_torus_degenerate(T, l0)) {
    return green_soft_dist_weight(dp, l0, Nj);
  }
  const real dist = dp.norm();
  if (dist < 1e-12) {
    return 0.0;
  }
  const vec3 G = albers::torus_sdf_grad(dp, T);
  if (!G.allFinite() || G.norm() < 1e-10) {
    return green_soft_dist_weight(dp, l0, Nj);
  }
  return green_soft_dist_weight(dp, l0, G);
}

} // namespace calder
} // namespace gaudi

#endif // __CALDER_SHAPE_OPERATOR_WEIGHTS__
