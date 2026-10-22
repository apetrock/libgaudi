#ifndef GAUDI_ALBERS_OSCULATING_TORUS_HPP
#define GAUDI_ALBERS_OSCULATING_TORUS_HPP

#include "gaudi/common.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>

namespace gaudi {
namespace albers {

struct principal_curvature_frame {
  vec3 n = vec3::UnitZ();
  vec3 e_min = vec3::UnitX();
  vec3 e_max = vec3::UnitY();
  real k_min = 0.0;
  real k_max = 0.0;
  bool valid = false;
};

/// Extract the two tangential principal curvatures from \(W\).
/// Uses \p n_hint to identify / flip the near-normal eigenpair.
inline principal_curvature_frame
principal_frame_from_shape_operator(const mat3 &W, const vec3 &n_hint) {
  principal_curvature_frame out;
  if (!W.allFinite()) {
    return out;
  }
  Eigen::SelfAdjointEigenSolver<mat3> es(W);
  if (es.info() != Eigen::Success) {
    return out;
  }
  const vec3 n_ref =
      (n_hint.norm() > 1e-12) ? n_hint.normalized() : vec3::UnitZ();

  int in = 0;
  real best_n = -1.0;
  for (int i = 0; i < 3; ++i) {
    const real a = std::abs(es.eigenvectors().col(i).dot(n_ref));
    if (a > best_n) {
      best_n = a;
      in = i;
    }
  }
  const int i0 = (in + 1) % 3;
  const int i1 = (in + 2) % 3;
  real k0 = es.eigenvalues()[i0];
  real k1 = es.eigenvalues()[i1];
  vec3 e0 = es.eigenvectors().col(i0);
  vec3 e1 = es.eigenvectors().col(i1);
  if (k0 > k1) {
    std::swap(k0, k1);
    std::swap(e0, e1);
  }
  out.k_min = k0;
  out.k_max = k1;
  out.e_min = e0.normalized();
  out.e_max = e1.normalized();
  out.n = es.eigenvectors().col(in);
  if (out.n.dot(n_ref) < 0.0) {
    out.n *= -1.0;
  }
  out.n = out.n.normalized();
  // Right-handed (e_min, e_max, n).
  if (out.e_min.cross(out.e_max).dot(out.n) < 0.0) {
    out.e_max *= -1.0;
  }
  out.valid = true;
  return out;
}

enum class torus_sign_mode {
  pat_eq22,       // follow κ_max < 0 → interior (+tube toward -n if κ_max>0 ...)
  interior_biased,
  unsigned_side,
};

struct osculating_torus {
  vec3 center = vec3::Zero();
  vec3 axis = vec3::UnitZ();
  real R = 1.0; // major
  real r = 0.1; // minor
  int sign = 1; // side of surface (+1 or -1)
  bool valid = false;
  // Tube frame at the foot. `spine` is the tube center nearest the foot, a
  // distance `r` along the normal. `radial` points from the circle's center
  // toward that spine point. `ruling` is the circle tangent there, along the
  // smaller principal direction. `tube` stays set when the circle is refused.
  vec3 spine = vec3::Zero();
  vec3 radial = vec3::UnitZ();
  vec3 ruling = vec3::UnitX();
  bool tube = false;
};

/// Tube radius \(r=1/|\kappa|_{\mathrm{hi}}\), the same center as
/// `cylinder_minor`: the minor spine is `foot + n/κ_hi` on both sheets.
/// `κ_lo` does not move that spine. Same-sign curvature puts the major center
/// further along that normal, at distance \(1/|\kappa|_{\mathrm{lo}}\).
/// Opposite signs put the major center on the other side of the foot.
/// Vanishing \(\kappa_{\mathrm{lo}}\) keeps a large-\(R\) torus (Feng, Gkioulekas
/// & Crane, *Points as Tori*, `feng-points-as-tori.pdf`): cylinder / plane are
/// limits of that torus, not a refused circle.
inline osculating_torus
osculating_torus_from_shape_operator(const vec3 &foot, const vec3 &n_hint,
                                     const mat3 &W,
                                     torus_sign_mode mode = torus_sign_mode::pat_eq22,
                                     real fallback_r = 0.1) {
  osculating_torus T;
  const principal_curvature_frame pc =
      principal_frame_from_shape_operator(W, n_hint);
  if (!pc.valid) {
    return T;
  }

  // Tube from larger-|κ|; major from the other.
  const real a0 = std::abs(pc.k_min);
  const real a1 = std::abs(pc.k_max);
  const bool max_is_kmax = a1 >= a0;
  const real k_hi = max_is_kmax ? pc.k_max : pc.k_min;
  const real k_lo = max_is_kmax ? pc.k_min : pc.k_max;
  const vec3 e_lo = max_is_kmax ? pc.e_min : pc.e_max;

  // Both principal curvatures vanishing is a plane. No tube.
  if (!(std::abs(k_hi) > real(1e-8)) || !pc.n.allFinite() || !e_lo.allFinite())
    return T;

  T.r = 1.0 / std::abs(k_hi);
  if (!std::isfinite(T.r) || T.r < 1e-12)
    T.r = fallback_r;
  if (!(T.r > real(1e-12)))
    return T;

  // sign stays +1 when κ_hi>0 so the Monge height graph does not flip.
  // The spine itself follows 1/κ_hi, which is the other direction.
  int side = 1;
  if (mode == torus_sign_mode::pat_eq22)
    side = (k_hi < 0.0) ? -1 : 1;
  const real s_prod =
      (pc.k_min * pc.k_max >= 0.0) ? real(1.0) : real(-1.0);
  const real sgn = (k_hi >= 0.0) ? real(1.0) : real(-1.0);
  const vec3 n = pc.n.normalized();
  vec3 ruling = e_lo - n * e_lo.dot(n);
  if (!(ruling.squaredNorm() > real(1e-16)))
    return T;
  ruling.normalize();
  T.sign = side;
  T.ruling = ruling;
  // cylinder_minor center. κ_lo does not enter.
  T.spine = foot + sgn * T.r * n;
  // From the spine toward the foot, so spine + r * radial recovers the foot.
  T.radial = -sgn * n;
  T.tube = T.spine.allFinite() && T.radial.allFinite() && T.ruling.allFinite();
  if (!T.tube)
    return T;

  // Floor a vanishing major κ so R is huge but finite (Points as Tori).
  constexpr real k_lo_floor = real(1e-8);
  const real abs_lo = std::max(std::abs(k_lo), k_lo_floor);
  const real inv_lo = 1.0 / abs_lo;
  T.R = inv_lo - s_prod * T.r;
  if (!std::isfinite(T.R))
    return T;
  // Umbilic / sphere limit is R→0 (valid degenerate torus). Do not inflate.
  T.R = std::max(T.R, real(0.0));

  const vec3 revolution = n.cross(e_lo);
  if (!(revolution.squaredNorm() > real(1e-16)))
    return T;
  T.axis = revolution.normalized();
  // Same sign: further along n/κ_hi. Opposite signs: the other side of the foot.
  T.center = foot + sgn * s_prod * inv_lo * n;
  vec3 radial = T.spine - T.center;
  if (radial.squaredNorm() > real(1e-16))
    T.radial = radial.normalized();
  T.valid = T.center.allFinite() && T.axis.allFinite() && std::isfinite(T.R) &&
            std::isfinite(T.r) && std::isfinite(inv_lo);
  return T;
}

/// Nearest point on the spine circle. Measured from `spine`, so a large R
/// does not subtract two far radii. R→∞ is the line through `spine` along
/// `ruling`.
inline bool nearest_circle_spine(const osculating_torus &T, const vec3 &p,
                                 vec3 &ps) {
  if (!T.valid || !T.tube || !p.allFinite())
    return false;
  const vec3 v = p - T.spine;
  const real tang = v.dot(T.ruling);
  const real radial_off = v.dot(T.radial);
  const real along = T.R + radial_off;
  const real rho = std::hypot(along, tang);
  if (!std::isfinite(rho))
    return false;
  if (!(rho > real(1e-18))) {
    ps = T.spine;
    return ps.allFinite();
  }
  const real scale = T.R / rho;
  const real denom = along + rho;
  const real sag =
      (std::abs(denom) > real(1e-18)) ? (scale * tang * tang) / denom : real(0);
  ps = T.spine - sag * T.radial + (scale * tang) * T.ruling;
  return ps.allFinite();
}

/// Nearest point on the minor-radius cylinder's spine: the line through
/// `spine` along `ruling`.
inline bool nearest_line_spine(const osculating_torus &T, const vec3 &p,
                               vec3 &ps) {
  if (!T.tube || !p.allFinite() || !(T.ruling.squaredNorm() > real(1e-16)))
    return false;
  const vec3 v = p - T.spine;
  ps = T.spine + v.dot(T.ruling) * T.ruling;
  return ps.allFinite();
}

/// \(p_s + r\,\mathrm{normalize}(p - p_s)\). Same ray on either side of the tube.
inline bool land_from_spine(const vec3 &ps, real radius, const vec3 &p,
                            vec3 &out) {
  if (!ps.allFinite() || !p.allFinite() || !(radius > real(1e-12)))
    return false;
  const vec3 ray = p - ps;
  const real q = ray.norm();
  if (!(q > real(1e-12)) || !std::isfinite(q))
    return false;
  out = ps + radius * (ray / q);
  return out.allFinite();
}

/// Approximate SDF of an infinite torus of revolution (positive outside).
inline real torus_sdf(const vec3 &p, const osculating_torus &T) {
  if (!T.valid) {
    return std::numeric_limits<real>::infinity();
  }
  const vec3 d = p - T.center;
  const real axial = d.dot(T.axis);
  const vec3 radial = d - axial * T.axis;
  const real rho = radial.norm();
  const real q = std::sqrt((rho - T.R) * (rho - T.R) + axial * axial);
  return q - T.r;
}

/// ∇ of torus_sdf. For foot-jet use: build T with foot=0, evaluate at dp.
inline vec3 torus_sdf_grad(const vec3 &p, const osculating_torus &T) {
  if (!T.valid) {
    return vec3::Zero();
  }
  const vec3 axis = T.axis.normalized();
  const vec3 d = p - T.center;
  const real axial = d.dot(axis);
  const vec3 radial = d - axial * axis;
  const real rho = radial.norm();
  const real u = rho - T.R;
  const real q = std::sqrt(u * u + axial * axial);
  if (q < 1e-12) {
    // On the major circle: any radial direction in the tube plane is fine.
    if (rho > 1e-12) {
      return radial.normalized();
    }
    // Degenerate: fall back to axis-orthogonal from d, else a fixed basis.
    vec3 t = d - axial * axis;
    if (t.squaredNorm() < 1e-20) {
      t = axis.unitOrthogonal();
    }
    return t.normalized();
  }
  vec3 rho_hat = vec3::Zero();
  if (rho > 1e-12) {
    rho_hat = radial / rho;
  } else {
    rho_hat = axis.unitOrthogonal();
  }
  return (u / q) * rho_hat + (axial / q) * axis;
}

inline real torus_abs_distance(const vec3 &p, const osculating_torus &T) {
  return std::abs(torus_sdf(p, T));
}

} // namespace albers
} // namespace gaudi

#endif // GAUDI_ALBERS_OSCULATING_TORUS_HPP
