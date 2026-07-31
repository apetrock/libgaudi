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
};

/// PAT algebra (Feng–Gkioulekas–Crane): \(r=1/|\kappa|_{\mathrm{hi}}\),
/// \(R=1/|\kappa|_{\mathrm{lo}}-\mathrm{sign}(\kappa_+\kappa_-)r\),
/// \(\mathrm{sign}(\mathbb{T})\) from Eq. 22 (~ follow \(\kappa_{\max}<0\)).
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

  const real eps = 1e-12;
  T.r = 1.0 / std::max(std::abs(k_hi), eps);
  if (!std::isfinite(T.r) || T.r < 1e-12) {
    T.r = fallback_r;
  }
  const real s_prod =
      (pc.k_min * pc.k_max >= 0.0) ? real(1.0) : real(-1.0);
  const real inv_lo = 1.0 / std::max(std::abs(k_lo), eps);
  T.R = inv_lo - s_prod * T.r;
  if (!std::isfinite(T.R)) {
    T.R = std::max(T.r, fallback_r);
  } else {
    // Umbilic / sphere limit is R→0 (valid degenerate torus). Do not inflate.
    T.R = std::max(T.R, real(0.0));
  }

  // Offset of geometric center is foot - side*(R+r)*n.
  // With outward normals and κ>0 (convex), side=+1 puts the torus on the
  // interior / medial side (-n). Concave κ_hi<0 → side=-1 → exterior of the
  // local bowl (still the center of curvature).
  int side = 1;
  if (mode == torus_sign_mode::pat_eq22) {
    side = (k_hi < 0.0) ? -1 : 1;
  } else if (mode == torus_sign_mode::interior_biased) {
    side = 1; // always toward -n
  } else {
    side = 1;
  }
  T.sign = side;
  T.axis = e_lo.normalized();
  // Geometric center of revolution (outer contact): foot is at distance R+r
  // from C along ±n. Tube/major-circle point would be foot - side*r*n.
  T.center =
      foot - static_cast<real>(side) * (T.R + T.r) * pc.n;
  T.valid = T.center.allFinite() && T.axis.allFinite() && std::isfinite(T.R) &&
            std::isfinite(T.r);
  return T;
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
