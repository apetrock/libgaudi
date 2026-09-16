#ifndef GAUDI_HEPWORTH_SLAB_PRISM_OPERATORS_HPP
#define GAUDI_HEPWORTH_SLAB_PRISM_OPERATORS_HPP

#include "slab_pierce_stencil.hpp"
#include "slab_prism_config.hpp"
#include "bivector3.hpp"
#include "trivector3.hpp"

#include <Eigen/SVD>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

struct slab_prism_rest {
  real V_rest = 0.0;
  /// Raw SVD σ-product of top-ring covariance at rest.
  real V_sigma_top = 1.0;
  /// Raw SVD σ-product of bottom-ring covariance at rest.
  real V_sigma_bot = 1.0;
  /// Top-ring SVD singular values at rest (bulk length/volume spring target).
  vec3 sigma_rest_top = vec3::Ones();
  /// Bottom-ring SVD singular values at rest.
  vec3 sigma_rest_bot = vec3::Ones();
  trivector3 tau_rest;
  real L_rest = 1.0;
  mat3 R0_rest = mat3::Identity();
  mat3 R1_rest = mat3::Identity();
  quat O_rest = quat::Identity();
  vec3 c_top_rest = vec3::Zero();
  vec3 c_bot_rest = vec3::Zero();
  vec3 d_rest = vec3::Zero();
  std::vector<vec3> d_top_rest;
  std::vector<vec3> d_bot_rest;
  std::vector<real> top_w;
  std::vector<real> bot_w;
};

inline real effective_slab_weight(real w_base, real L_rest, real L_min,
                                  real w_eff_max = 0.0) {
  const real w =
      w_base / std::max(L_rest, std::max(L_min, real(1e-12)));
  if (w_eff_max > real(0.0))
    return std::min(w, w_eff_max);
  return w;
}

inline mat3 fit_procrustes_rotation(const mat3 &M) {
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  mat3 R = U * V.transpose();
  if (R.determinant() < real(0.0)) {
    U.col(2) *= real(-1.0);
    R = U * V.transpose();
  }
  return R;
}

/// Nearest det=+1 affine map from weighted covariance M = Σ w_k x_k x_k^restᵀ.
inline mat3 fit_volume_preserving_F(const mat3 &M) {
  if (M.squaredNorm() < 1e-24)
    return mat3::Identity();
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  vec3 sigma = svd.singularValues();

  mat3 R = U * V.transpose();
  if (R.determinant() < real(0.0)) {
    U.col(2) *= real(-1.0);
    sigma[2] = -sigma[2];
    R = U * V.transpose();
  }

  const real det_sigma = sigma[0] * sigma[1] * sigma[2];
  if (std::abs(det_sigma) < 1e-24)
    return R;

  const mat3 S = V * sigma.asDiagonal() * V.transpose();
  const real lam = std::cbrt(real(1.0) / det_sigma);
  return R * (lam * S);
}

/// Clamp det=+1 stretch singular values; keeps rotation, limits volume-preserving shear response.
inline mat3 clamp_volume_preserving_stretch(mat3 F, real max_stretch) {
  if (max_stretch <= real(1.0) + real(1e-6))
    return fit_procrustes_rotation(F);
  Eigen::JacobiSVD<mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  vec3 sigma = svd.singularValues();
  const real lo = real(1.0) / max_stretch;
  for (int i = 0; i < 3; ++i)
    sigma[i] = std::clamp(sigma[i], lo, max_stretch);
  const real det = sigma[0] * sigma[1] * sigma[2];
  if (std::abs(det) > real(1e-24)) {
    const real lam = std::cbrt(real(1.0) / det);
    sigma *= lam;
  }
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  if ((U * V.transpose()).determinant() < real(0.0))
    U.col(2) *= real(-1.0);
  return U * sigma.asDiagonal() * V.transpose();
}

inline mat3 fit_volume_preserving_F_clamped(const mat3 &M, real max_stretch) {
  return clamp_volume_preserving_stretch(fit_volume_preserving_F(M), max_stretch);
}

inline mat3 ring_covariance(const vec3 &hub, const std::vector<index_t> &ring_ids,
                            const std::vector<real> &weights,
                            const std::vector<vec3> &d_rest,
                            const std::function<vec3(index_t)> &pos_at) {
  mat3 M = mat3::Zero();
  for (size_t k = 0; k < weights.size(); ++k) {
    const real wk = weights[k];
    const vec3 dq = pos_at(ring_ids[k + 1]) - hub;
    M += wk * dq * d_rest[k].transpose();
  }
  return M;
}

inline real svd_sigma_product(const mat3 &M) {
  if (M.squaredNorm() < 1e-24)
    return real(1.0);
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  vec3 sigma = svd.singularValues();
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  if ((U * V.transpose()).determinant() < real(0.0))
    sigma[2] = -sigma[2];
  return sigma[0] * sigma[1] * sigma[2];
}

inline vec3 svd_sigma_vector(const mat3 &M) {
  if (M.squaredNorm() < 1e-24)
    return vec3::Ones();
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  vec3 sigma = svd.singularValues();
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  if ((U * V.transpose()).determinant() < real(0.0))
    sigma[2] = -sigma[2];
  return sigma.cwiseAbs();
}

struct bulk_svd_frame {
  mat3 U = mat3::Identity();
  mat3 V = mat3::Identity();
  vec3 sigma = vec3::Ones();
};

inline bulk_svd_frame bulk_svd_frame_from_M(const mat3 &M) {
  bulk_svd_frame out;
  if (M.squaredNorm() < 1e-24)
    return out;
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  out.U = svd.matrixU();
  out.V = svd.matrixV();
  out.sigma = svd.singularValues();
  if ((out.U * out.V.transpose()).determinant() < real(0.0)) {
    out.U.col(2) *= real(-1.0);
    out.sigma[2] = -out.sigma[2];
  }
  out.sigma = out.sigma.cwiseAbs();
  return out;
}

/// Pole normal from a hub-and-spoke ring (live geometry).
inline vec3 top_ring_pole_normal(
    const vec3 &hub, const std::vector<index_t> &ring_ids,
    const std::function<vec3(index_t)> &pos_at, real L_min = 1e-8) {
  if (ring_ids.size() < 3)
    return vec3::UnitZ();
  vec3 n = vec3::Zero();
  const size_t nsp = ring_ids.size() - 1;
  for (size_t k = 0; k < nsp; ++k) {
    const vec3 a = pos_at(ring_ids[k + 1]) - hub;
    const vec3 b = pos_at(ring_ids[(k + 1) % nsp + 1]) - hub;
    n += a.cross(b);
  }
  const real nn = n.norm();
  if (nn <= L_min)
    return vec3::UnitZ();
  return n / nn;
}

/// Orient ring normal to point outward (opposite the foot→pierce spine).
inline vec3 outward_ring_normal(const vec3 &n_ring, const vec3 &d_spine,
                                real L_min = 1e-8) {
  vec3 n = n_ring;
  if (n.squaredNorm() <= L_min * L_min)
    return vec3::UnitZ();
  n.normalize();
  const vec3 d = d_spine.norm() > L_min ? d_spine.normalized() : vec3::UnitZ();
  if (n.dot(d) > real(0.0))
    n = -n;
  return n;
}

/// Live bulk SVD axis (column of U) with largest alignment to the ring normal.
inline vec3 bulk_live_axis_nearest_normal(const mat3 &M, const vec3 &n_ring,
                                          real L_min = 1e-8) {
  vec3 n = n_ring;
  if (n.squaredNorm() <= L_min * L_min)
    return vec3::UnitZ();
  n.normalize();

  const bulk_svd_frame frame = bulk_svd_frame_from_M(M);
  int best = 0;
  real best_dot = std::abs(frame.U.col(0).dot(n));
  for (int i = 1; i < 3; ++i) {
    const real d = std::abs(frame.U.col(i).dot(n));
    if (d > best_dot) {
      best_dot = d;
      best = i;
    }
  }
  vec3 axis = frame.U.col(best);
  if (axis.dot(n) < real(0.0))
    axis = -axis;
  const real an = axis.norm();
  return an > L_min ? axis / an : n;
}

/// Thickness pump along outward normal: positive δ expands, negative δ compresses inward.
inline real bulk_sigma_volume_scale(const vec3 &sigma_rest, real V_rest,
                                    real L_spine, real L_min) {
  const real L = std::max(L_spine, L_min);
  const real geom = L * L * L;
  const real prod =
      std::abs(sigma_rest[0] * sigma_rest[1] * sigma_rest[2]);
  const real planar =
      std::abs(sigma_rest[0] * sigma_rest[1] *
               std::max(sigma_rest[2], L_min / L));
  return std::max({geom, prod, planar, std::abs(V_rest), real(1e-12)});
}

inline real bulk_volume_pump_delta(real V_live, real V_target,
                                   const vec3 &sigma_rest, real V_rest,
                                   real L_spine, real L_min, real max_step) {
  if (max_step <= real(0.0))
    return real(0.0);
  const real scale =
      bulk_sigma_volume_scale(sigma_rest, V_rest, L_spine, L_min);
  const real mid =
      std::max(real(0.5) * (std::abs(V_live) + std::abs(V_target)),
               scale * real(1e-3));
  const real denom = std::max(scale, mid);
  const real frac = (V_target - V_live) / denom;
  return std::clamp(frac * max_step, -max_step, max_step);
}

inline bool bulk_sigma_targets_valid(const vec3 &sigma_rest, real V_target) {
  return sigma_rest.minCoeff() > real(1e-12) &&
         std::abs(V_target) > real(1e-12);
}

inline vec3 project_bulk_sigma_springs(const vec3 &l0_in, const vec3 &sigma_rest,
                                       real V_target, real max_stretch,
                                       real w_volume, real w_length) {
  vec3 s = l0_in.cwiseAbs();
  if (!bulk_sigma_targets_valid(sigma_rest, V_target))
    return s;

  const vec3 lt = sigma_rest;
  const real vol_tol = 1e-12 * std::max(std::abs(V_target), real(1e-12));
  const real p0 = s[0] * s[1] * s[2];
  if ((s - lt).squaredNorm() < 1e-24 && std::abs(p0 - V_target) < vol_tol)
    return lt;

  if (w_volume <= real(0.0) && w_length <= real(0.0))
    return s;

  if (w_volume <= real(0.0)) {
    s = lt;
  } else if (w_length <= real(0.0)) {
    if (std::abs(p0) > real(1e-24))
      s *= std::cbrt(V_target / p0);
  } else {
    const real ratio = w_volume / w_length;
    for (int iter = 0; iter < 32; ++iter) {
      const real p = s[0] * s[1] * s[2];
      const real dvol = p - V_target;
      vec3 next = lt;
      for (int i = 0; i < 3; ++i) {
        const real denom = std::max(s[i], real(1e-24));
        next[i] = lt[i] - ratio * dvol * p / denom;
      }
      next = next.cwiseMax(real(1e-24));
      if ((next - s).norm() < 1e-10 * std::max(s.norm(), real(1e-12)))
        break;
      s = next;
    }
  }

  if (max_stretch > real(1.0) + real(1e-6)) {
    const real s_max = max_stretch;
    for (int i = 0; i < 3; ++i)
      s[i] = std::clamp(s[i], lt[i] / s_max, lt[i] * s_max);
    const real p = s[0] * s[1] * s[2];
    if (w_volume > real(0.0) && std::abs(p) > real(1e-24) &&
        std::abs(V_target) > real(1e-24))
      s *= std::cbrt(V_target / p);
  }

  return s;
}

struct bulk_ring_map {
  mat3 Z = mat3::Identity();
  mat3 R = mat3::Identity();
  mat3 V = mat3::Identity();
  vec3 l0 = vec3::Ones();
  vec3 l = vec3::Ones();
};

inline vec3 clamp_slab_spoke_target_delta(const vec3 &target, const vec3 &current,
                                          real max_norm) {
  if (max_norm <= real(0.0))
    return target;
  const vec3 delta = target - current;
  const real n = delta.norm();
  if (n <= max_norm || n <= real(1e-24))
    return target;
  return current + (max_norm / n) * delta;
}

/// One vertex sample in the fixed top-foot rest frame.
struct slab_foot_local_sample {
  index_t vi = -1;
  vec3 local_rest = vec3::Zero();
  real weight = 1.0;
};

inline std::vector<slab_foot_local_sample>
build_top_ring_fit_samples(const slab_pierce_stencil &stencil,
                           const std::function<vec3(index_t)> &pos_at,
                           vec3 &foot_out) {
  std::vector<slab_foot_local_sample> samples;
  foot_out = pos_at(stencil.foot);
  samples.reserve(stencil.top_weights.size());
  for (size_t k = 0; k < stencil.top_weights.size(); ++k) {
    slab_foot_local_sample s;
    s.vi = stencil.top_ids[k + 1];
    s.weight = stencil.top_weights[k];
    s.local_rest = pos_at(s.vi) - foot_out;
    samples.push_back(s);
  }
  return samples;
}

inline mat3 foot_frame_cross_covariance(
    index_t foot, const std::vector<slab_foot_local_sample> &samples,
    const std::function<vec3(index_t)> &pos_at) {
  const vec3 c = pos_at(foot);
  mat3 M = mat3::Zero();
  for (const slab_foot_local_sample &s : samples) {
    const vec3 dq = pos_at(s.vi) - c;
    M += s.weight * dq * s.local_rest.transpose();
  }
  return M;
}

inline mat3 foot_frame_cross_covariance_rest(
    const std::vector<slab_foot_local_sample> &samples) {
  mat3 M = mat3::Zero();
  for (const slab_foot_local_sample &s : samples)
    M += s.weight * s.local_rest * s.local_rest.transpose();
  return M;
}

struct foot_frame_bulk_rest {
  vec3 sigma_rest = vec3::Ones();
  /// Rest principal axes (columns); stretch in Z is expressed in this frame.
  mat3 E_rest = mat3::Identity();
  real V_rest = 1.0;
  real L_rest = 1.0;
};

inline foot_frame_bulk_rest
build_foot_frame_bulk_rest(const std::vector<slab_foot_local_sample> &samples,
                           real L_min) {
  foot_frame_bulk_rest out;
  const mat3 Cov = foot_frame_cross_covariance_rest(samples);
  const bulk_svd_frame frame = bulk_svd_frame_from_M(Cov);
  out.sigma_rest = frame.sigma;
  // Covariance is SPD: U and V span the same axes.
  out.E_rest = frame.V;
  out.V_rest = out.sigma_rest[0] * out.sigma_rest[1] * out.sigma_rest[2];
  for (const slab_foot_local_sample &s : samples)
    out.L_rest = std::max(out.L_rest, s.local_rest.norm());
  out.L_rest = std::max(out.L_rest, L_min);
  return out;
}

inline mat3 foot_frame_live_covariance(
    index_t foot, const std::vector<slab_foot_local_sample> &samples,
    const std::function<vec3(index_t)> &pos_at) {
  const vec3 c = pos_at(foot);
  mat3 Cov = mat3::Zero();
  for (const slab_foot_local_sample &s : samples) {
    const vec3 dq = pos_at(s.vi) - c;
    Cov += s.weight * dq * dq.transpose();
  }
  return Cov;
}

inline mat3 bulk_Z_from_ring_map(const bulk_ring_map &m, const mat3 &E_rest,
                                 const vec3 &sigma_rest) {
  const vec3 denom = sigma_rest.cwiseMax(real(1e-24));
  const vec3 ratio = m.l.cwiseQuotient(denom);
  // R follows live deformation; stretch corrects extents in the rest frame.
  return m.R * E_rest * ratio.asDiagonal() * E_rest.transpose();
}

/// @param M_cross  Σ w (x_live − hub)(x_rest − hub)ᵀ  — Procrustes / R
/// @param Cov_live Σ w (x_live − hub)(x_live − hub)ᵀ — geometric live extents
inline bulk_ring_map bulk_procrustes_map(const mat3 &M_cross, const mat3 &Cov_live,
                                         const foot_frame_bulk_rest &rest,
                                         real max_stretch, real w_volume = 1.0,
                                         real w_length = 1.0) {
  bulk_ring_map out;
  if (M_cross.squaredNorm() < 1e-24)
    return out;

  Eigen::JacobiSVD<mat3> svd(M_cross, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  out.V = svd.matrixV();
  out.R = U * out.V.transpose();
  if (out.R.determinant() < real(0.0)) {
    U.col(2) *= real(-1.0);
    out.R = U * out.V.transpose();
  }

  // Live ring sizes from Cov_live — NOT singular values of M_cross.
  out.l0 = svd_sigma_vector(Cov_live);
  const real V_target =
      rest.sigma_rest[0] * rest.sigma_rest[1] * rest.sigma_rest[2];
  out.l = bulk_sigma_targets_valid(rest.sigma_rest, V_target)
              ? project_bulk_sigma_springs(out.l0, rest.sigma_rest, V_target,
                                           max_stretch, w_volume, w_length)
              : out.l0;
  out.Z = bulk_Z_from_ring_map(out, rest.E_rest, rest.sigma_rest);
  return out;
}

/// Legacy overload: treats M as both cross-cov and Cov_live (ok at rest only).
inline bulk_ring_map bulk_procrustes_map(const mat3 &M, const vec3 &sigma_rest,
                                         real V_target, real max_stretch,
                                         real w_volume = 1.0,
                                         real w_length = 1.0) {
  foot_frame_bulk_rest rest;
  rest.sigma_rest = sigma_rest;
  rest.V_rest = V_target;
  return bulk_procrustes_map(M, M, rest, max_stretch, w_volume, w_length);
}

inline mat3 fit_ring_rotation(index_t hub, const smallvec<index_t, 16> &ring_ids,
                                const smallvec<real, 16> &weights,
                                const vec3 &c_curr,
                                const std::vector<vec3> &d_rest,
                                const std::function<vec3(index_t)> &pos_at) {
  mat3 M = mat3::Zero();
  for (size_t k = 0; k < weights.size(); ++k) {
    const real wk = weights[k];
    const vec3 dq = pos_at(ring_ids[k + 1]) - c_curr;
    M += wk * dq * d_rest[k].transpose();
  }
  if (M.squaredNorm() < 1e-24)
    return mat3::Identity();
  (void)hub;
  return fit_procrustes_rotation(M);
}

inline quat mat3_to_quat(const mat3 &R) {
  return quat(R);
}

inline mat3 quat_to_mat3(const quat &q) {
  return q.normalized().toRotationMatrix();
}

inline void split_bend_twist_omega(const quat &dO, const vec3 &d_hat,
                                   vec3 &omega_b, vec3 &omega_t) {
  Eigen::AngleAxis<real> aa(dO);
  const vec3 omega = aa.angle() * aa.axis();
  const vec3 e2 = d_hat.squaredNorm() > 1e-24 ? d_hat.normalized() : vec3::UnitZ();
  omega_t = omega.dot(e2) * e2;
  omega_b = omega - omega_t;
}

inline quat omega_to_quat(const vec3 &omega) {
  const real ang = omega.norm();
  if (ang < 1.0e-16)
    return quat::Identity();
  return quat(Eigen::AngleAxis<real>(ang, omega / ang));
}

inline trivector3 eval_cap_trivector(index_t hub,
                                     const smallvec<index_t, 16> &ring_ids,
                                     const vec3 &d,
                                     const std::function<vec3(index_t)> &pos_at) {
  const vec3 c = pos_at(hub);
  const size_t n = ring_ids.size() - 1;
  trivector3 tau;
  if (n < 2)
    return tau;
  for (size_t k = 0; k < n; ++k) {
    const vec3 e0 = pos_at(ring_ids[k + 1]) - c;
    const vec3 e1 = pos_at(ring_ids[(k + 1) % n + 1]) - c;
    tau += trivector3::wedge(bivector3::wedge(e0, e1), d);
  }
  return tau;
}

inline real eval_cap_flux(index_t hub, const smallvec<index_t, 16> &ring_ids,
                          const vec3 &d, const std::function<vec3(index_t)> &pos_at) {
  return eval_cap_trivector(hub, ring_ids, d, pos_at).coeff;
}

inline trivector3 eval_slab_trivector(const slab_pierce_stencil &stencil,
                                      const std::function<vec3(index_t)> &pos_at) {
  const vec3 c_top = pos_at(stencil.foot);
  const vec3 c_bot = pos_at(stencil.pierce_foot);
  const vec3 d = c_bot - c_top;
  const trivector3 tau_top =
      eval_cap_trivector(stencil.foot, stencil.top_ids, d, pos_at);
  const trivector3 tau_bot =
      eval_cap_trivector(stencil.pierce_foot, stencil.bot_ids, -d, pos_at);
  return real(0.5) * (tau_top + tau_bot);
}

inline real eval_prism_volume(const slab_pierce_stencil &stencil,
                              const std::function<vec3(index_t)> &pos_at) {
  return eval_slab_trivector(stencil, pos_at).coeff;
}

inline real prism_volume_error_frac_from_values(real V_cur, real V_rest) {
  return std::abs(V_cur - V_rest) / std::max(real(1.0), std::abs(V_rest));
}

inline real prism_volume_error_frac(const slab_pierce_stencil &stencil,
                                    const std::function<vec3(index_t)> &pos_at,
                                    real V_rest) {
  return prism_volume_error_frac_from_values(
      eval_prism_volume(stencil, pos_at), V_rest);
}

inline real prism_volume_abs_tol(real V_rest) {
  return real(1e-12) * std::max(real(1.0), std::abs(V_rest));
}

/// When true, SL3 project() copies live hub-spoke edges (p = Aq on top ring).
inline bool sl3_volume_pass_through(real V_cur, real V_rest,
                                    real min_response_frac) {
  const real dv = std::abs(V_cur - V_rest);
  if (dv <= prism_volume_abs_tol(V_rest))
    return true;
  if (min_response_frac > real(0.0) &&
      prism_volume_error_frac_from_values(V_cur, V_rest) < min_response_frac)
    return true;
  return false;
}

inline quat mat3_relative_quat(const mat3 &Ra, const mat3 &Rb) {
  quat O(Ra.transpose() * Rb);
  if (O.w() < real(0.0))
    O.coeffs() *= real(-1.0);
  return O.normalized();
}

inline quat mat3_to_quat_hemisphere(const mat3 &R) {
  quat q(R);
  if (q.w() < real(0.0))
    q.coeffs() *= real(-1.0);
  return q.normalized();
}

inline void apply_twist_to_top_frame(mat3 &R0, mat3 &R1, const mat3 &R0_rest,
                                     const mat3 &R1_rest, const vec3 &d_hat,
                                     bool penalize_bend) {
  const quat O = mat3_relative_quat(R0, R1);
  const quat O_rest = mat3_relative_quat(R0_rest, R1_rest);
  quat dO = O.conjugate() * O_rest;
  if (dO.w() < real(0.0))
    dO.coeffs() *= real(-1.0);

  vec3 omega_b, omega_t;
  split_bend_twist_omega(dO, d_hat, omega_b, omega_t);
  dO = penalize_bend ? omega_to_quat(omega_b + omega_t) : omega_to_quat(omega_t);

  quat q0 = mat3_to_quat_hemisphere(R0);
  quat q1 = mat3_to_quat_hemisphere(R1);
  const quat q0p = (q0 * dO.conjugate()).normalized();
  const quat q1p = (q1 * dO).normalized();
  R0 = quat_to_mat3(q0p);
  R1 = quat_to_mat3(q1p);
}

inline slab_prism_rest
build_slab_prism_rest(const slab_pierce_stencil &stencil,
                      const std::function<vec3(index_t)> &pos_at) {
  slab_prism_rest out;
  out.c_top_rest = pos_at(stencil.foot);
  out.c_bot_rest = pos_at(stencil.pierce_foot);
  out.d_rest = out.c_bot_rest - out.c_top_rest;
  out.L_rest = std::max(out.d_rest.norm(), real(1e-8));

  out.top_w.resize(stencil.top_weights.size());
  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    out.top_w[k] = stencil.top_weights[k];
  out.bot_w.resize(stencil.bot_weights.size());
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    out.bot_w[k] = stencil.bot_weights[k];

  out.d_top_rest.resize(stencil.top_weights.size());
  for (size_t k = 0; k < stencil.top_weights.size(); ++k)
    out.d_top_rest[k] = pos_at(stencil.top_ids[k + 1]) - out.c_top_rest;

  out.d_bot_rest.resize(stencil.bot_weights.size());
  for (size_t k = 0; k < stencil.bot_weights.size(); ++k)
    out.d_bot_rest[k] = pos_at(stencil.bot_ids[k + 1]) - out.c_bot_rest;

  out.tau_rest = eval_slab_trivector(stencil, pos_at);
  out.V_rest = out.tau_rest.coeff;
  {
    const std::vector<index_t> top_ids(stencil.top_ids.begin(),
                                       stencil.top_ids.end());
    const std::vector<real> top_w(stencil.top_weights.begin(),
                                  stencil.top_weights.end());
    const std::vector<index_t> bot_ids(stencil.bot_ids.begin(),
                                       stencil.bot_ids.end());
    const std::vector<real> bot_w(stencil.bot_weights.begin(),
                                  stencil.bot_weights.end());
    const mat3 M_top = ring_covariance(out.c_top_rest, top_ids, top_w,
                                       out.d_top_rest, pos_at);
    const mat3 M_bot = ring_covariance(out.c_bot_rest, bot_ids, bot_w,
                                       out.d_bot_rest, pos_at);
    out.sigma_rest_top = svd_sigma_vector(M_top);
    out.sigma_rest_bot = svd_sigma_vector(M_bot);
    out.V_sigma_top = out.sigma_rest_top[0] * out.sigma_rest_top[1] *
                      out.sigma_rest_top[2];
    out.V_sigma_bot = out.sigma_rest_bot[0] * out.sigma_rest_bot[1] *
                      out.sigma_rest_bot[2];
  }
  out.R0_rest =
      fit_ring_rotation(stencil.foot, stencil.top_ids, stencil.top_weights,
                        out.c_top_rest, out.d_top_rest, pos_at);
  out.R1_rest =
      fit_ring_rotation(stencil.pierce_foot, stencil.bot_ids, stencil.bot_weights,
                        out.c_bot_rest, out.d_bot_rest, pos_at);
  out.O_rest = mat3_relative_quat(out.R0_rest, out.R1_rest);
  return out;
}

inline slab_prism_rest
build_slab_prism_rest(const slab_pierce_stencil &stencil,
                      const std::vector<vec3> &x_rest) {
  const auto pos_at = [&](index_t vi) { return x_rest[static_cast<size_t>(vi)]; };
  return build_slab_prism_rest(stencil, pos_at);
}

struct top_ring_positions {
  std::vector<index_t> ids;
  std::vector<vec3> x;
  std::unordered_map<index_t, size_t> idx;

  vec3 at(index_t vi) const { return x.at(idx.at(vi)); }

  void set(index_t vi, const vec3 &p) { x.at(idx.at(vi)) = p; }

  static top_ring_positions from_stencil(const slab_pierce_stencil &stencil,
                                         const std::function<vec3(index_t)> &pos_at) {
    top_ring_positions out;
    out.ids.reserve(stencil.top_ids.size());
    for (size_t i = 0; i < stencil.top_ids.size(); ++i) {
      out.idx[stencil.top_ids[i]] = out.ids.size();
      out.ids.push_back(stencil.top_ids[i]);
      out.x.push_back(pos_at(stencil.top_ids[i]));
    }
    return out;
  }
};

/// Union of top- and bottom-ring verts (both hubs + one-rings).
struct slab_ring_positions {
  std::vector<index_t> ids;
  std::vector<vec3> x;
  std::unordered_map<index_t, size_t> idx;

  vec3 at(index_t vi) const { return x.at(idx.at(vi)); }

  void set(index_t vi, const vec3 &p) { x.at(idx.at(vi)) = p; }

  static slab_ring_positions from_stencil(const slab_pierce_stencil &stencil,
                                          const std::function<vec3(index_t)> &pos_at) {
    slab_ring_positions out;
    const auto add = [&](index_t vi) {
      if (out.idx.count(vi))
        return;
      out.idx[vi] = out.ids.size();
      out.ids.push_back(vi);
      out.x.push_back(pos_at(vi));
    };
    for (index_t vi : stencil.top_ids)
      add(vi);
    for (index_t vi : stencil.bot_ids)
      add(vi);
    return out;
  }
};

inline trivector3 eval_prism_trivector_top_only(
    const slab_pierce_stencil &stencil, const top_ring_positions &top,
    const std::function<vec3(index_t)> &anchor_at) {
  const vec3 c_bot = anchor_at(stencil.pierce_foot);
  const vec3 c_top = top.at(stencil.foot);
  const vec3 d = c_bot - c_top;
  const auto pos_at = [&](index_t vi) -> vec3 {
    if (vi == stencil.pierce_foot)
      return c_bot;
    if (vi == stencil.foot)
      return c_top;
    const auto it = top.idx.find(vi);
    if (it != top.idx.end())
      return top.x[it->second];
    return anchor_at(vi);
  };
  const trivector3 tau_top =
      eval_cap_trivector(stencil.foot, stencil.top_ids, d, pos_at);
  const trivector3 tau_bot =
      eval_cap_trivector(stencil.pierce_foot, stencil.bot_ids, -d, pos_at);
  return real(0.5) * (tau_top + tau_bot);
}

inline real eval_prism_volume_top_only(const slab_pierce_stencil &stencil,
                                       const top_ring_positions &top,
                                       const std::function<vec3(index_t)> &anchor_at) {
  return eval_prism_trivector_top_only(stencil, top, anchor_at).coeff;
}

inline void prism_trivector_grad_top_analytic(
    const slab_pierce_stencil &stencil, const top_ring_positions &top,
    const std::function<vec3(index_t)> &anchor_at, std::vector<vec3> &grad) {
  grad.assign(top.x.size(), vec3::Zero());
  const vec3 c_bot = anchor_at(stencil.pierce_foot);
  const vec3 c_top = top.at(stencil.foot);
  const vec3 d = c_bot - c_top;

  const auto pos_top = [&](index_t vi) -> vec3 {
    if (vi == stencil.foot)
      return c_top;
    return top.at(vi);
  };

  const auto add = [&](index_t vi, const vec3 &g) {
    const auto it = top.idx.find(vi);
    if (it != top.idx.end())
      grad[it->second] += g;
  };

  const size_t n_top = stencil.top_ids.size() - 1;
  for (size_t k = 0; k < n_top; ++k) {
    const index_t ik = stencil.top_ids[k + 1];
    const index_t ik1 = stencil.top_ids[(k + 1) % n_top + 1];
    const vec3 e0 = pos_top(ik) - c_top;
    const vec3 e1 = pos_top(ik1) - c_top;
    const vec3 g0 = real(0.5) * e1.cross(d);
    const vec3 g1 = real(0.5) * d.cross(e0);
    const vec3 gc =
        real(0.5) * (-e1.cross(d) - d.cross(e0) - e0.cross(e1));
    add(ik, g0);
    add(ik1, g1);
    add(stencil.foot, gc);
  }

  const vec3 c_bot_hub = c_bot;
  const size_t n_bot = stencil.bot_ids.size() - 1;
  for (size_t k = 0; k < n_bot; ++k) {
    const vec3 e0 = anchor_at(stencil.bot_ids[k + 1]) - c_bot_hub;
    const vec3 e1 =
        anchor_at(stencil.bot_ids[(k + 1) % n_bot + 1]) - c_bot_hub;
    add(stencil.foot, real(0.5) * e0.cross(e1));
  }
}

inline void prism_volume_grad_fd(const slab_pierce_stencil &stencil,
                                 const top_ring_positions &top,
                                 const std::function<vec3(index_t)> &anchor_at,
                                 std::vector<vec3> &grad, real eps = 1e-6) {
  grad.assign(top.x.size(), vec3::Zero());
  const auto energy_at = [&](const top_ring_positions &pos) {
    return eval_prism_volume_top_only(stencil, pos, anchor_at);
  };
  top_ring_positions pert = top;
  for (size_t i = 0; i < top.x.size(); ++i) {
    for (int ax = 0; ax < 3; ++ax) {
      pert.x[i] = top.x[i];
      pert.x[i][ax] += eps;
      const real ep = energy_at(pert);
      pert.x[i][ax] -= real(2.0) * eps;
      const real em = energy_at(pert);
      grad[i][ax] = (ep - em) / (real(2.0) * eps);
      pert.x[i] = top.x[i];
    }
  }
}

/// Minimum-norm Newton on the scalar constraint V(x) = V_rest (top ring DOFs only).
/// Matches ARAP-style closed-form local resolve: one linearized step per iter.
inline top_ring_positions
resolve_prism_volume_newton(const slab_pierce_stencil &stencil,
                            const top_ring_positions &x0,
                            const std::function<vec3(index_t)> &anchor_at,
                            real V_rest, int iters) {
  top_ring_positions x = x0;
  const real dv_tol =
      real(1e-12) * std::max(real(1.0), std::abs(V_rest));

  for (int iter = 0; iter < iters; ++iter) {
    const real dv =
        eval_prism_volume_top_only(stencil, x, anchor_at) - V_rest;
    if (std::abs(dv) <= dv_tol)
      break;

    std::vector<vec3> grad;
    prism_trivector_grad_top_analytic(stencil, x, anchor_at, grad);
    real g2 = 0.0;
    for (const vec3 &g : grad)
      g2 += g.squaredNorm();
    if (g2 < 1e-32)
      break;

    real alpha = -dv / g2;
    top_ring_positions trial = x;
    for (int ls = 0; ls < 8; ++ls) {
      for (size_t i = 0; i < x.x.size(); ++i)
        trial.x[i] = x.x[i] + alpha * grad[i];
      const real dv_trial =
          eval_prism_volume_top_only(stencil, trial, anchor_at) - V_rest;
      if (std::abs(dv_trial) < std::abs(dv)) {
        x = trial;
        break;
      }
      alpha *= real(0.5);
    }
  }
  return x;
}

inline top_ring_positions
resolve_prism_volume_gd(const slab_pierce_stencil &stencil,
                        const top_ring_positions &x0,
                        const std::function<vec3(index_t)> &anchor_at,
                        real V_rest, int iters, real step_scale) {
  (void)step_scale;
  return resolve_prism_volume_newton(stencil, x0, anchor_at, V_rest, iters);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PRISM_OPERATORS_HPP
