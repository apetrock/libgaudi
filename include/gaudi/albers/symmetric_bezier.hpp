#ifndef GAUDI_ALBERS_SYMMETRIC_BEZIER_HPP
#define GAUDI_ALBERS_SYMMETRIC_BEZIER_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace albers {

// Local frame for discrete screw / π-flip symmetry about the rod axis.
// Pairing is central inversion through origin: p → 2o - p, which sends
// (s, v_⊥) → (-s, -v_⊥) — ahead↔behind with lateral flip. That is the
// stencil-level model of a half-turn screw about T (coil-like), not a
// planar mirror (arch-like).
struct screw_frame {
  vec3 origin = vec3::Zero();
  vec3 e_T = vec3::UnitX(); // screw / sort axis (local tangent)
  real sigma[3] = {0, 0, 0};
  bool ok = false;
};

inline vec3 screw_flip(const vec3 &p, const vec3 &origin) {
  return 2.0 * origin - p;
}

// π rotation about axis e_T through origin (orientation-preserving half-turn).
// Alone this keeps the T-coordinate; compose with T-reversal (= screw_flip).
inline vec3 rotate_pi_about_axis(const vec3 &p, const vec3 &origin,
                                 const vec3 &e_T) {
  const vec3 T = e_T.normalized();
  const vec3 d = p - origin;
  return origin + 2.0 * d.dot(T) * T - d;
}

inline std::array<real, 4> bernstein3(real t) {
  const real u = 1.0 - t;
  const real u2 = u * u;
  const real u3 = u2 * u;
  const real t2 = t * t;
  const real t3 = t2 * t;
  return {u3, 3.0 * u2 * t, 3.0 * u * t2, t3};
}

inline vec3 eval_cubic_bezier(const std::array<vec3, 4> &P, real t) {
  const auto b = bernstein3(t);
  return b[0] * P[0] + b[1] * P[1] + b[2] * P[2] + b[3] * P[3];
}

// Axis from cached tangent, refined by PCA (pick eigenvector closest to T).
// Rejects near-isotropic clouds where the screw axis is ill-defined.
inline screw_frame screw_frame_from_centroid_pca(const std::vector<vec3> &points,
                                                 const vec3 &T_cached,
                                                 real sep_ratio_min = 1.15) {
  screw_frame f;
  if (points.size() < 3)
    return f;

  vec3 cen = vec3::Zero();
  for (const vec3 &p : points)
    cen += p;
  cen /= real(points.size());
  f.origin = cen;

  mat3 C = mat3::Zero();
  for (const vec3 &p : points) {
    const vec3 d = p - cen;
    C += d * d.transpose();
  }
  C /= real(points.size());

  Eigen::SelfAdjointEigenSolver<mat3> es(C);
  if (es.info() != Eigen::Success)
    return f;

  const vec3 s = es.eigenvalues();
  if (s.hasNaN() || s[2] < 1e-18)
    return f;

  std::array<vec3, 3> e = {es.eigenvectors().col(0), es.eigenvectors().col(1),
                           es.eigenvectors().col(2)};
  f.sigma[0] = s[2];
  f.sigma[1] = s[1];
  f.sigma[2] = s[0];

  // Need a dominant axis to serve as T; reject when the top two σ's are close
  // (e.g. isotropic in-plane pentagon).
  const real r01 = f.sigma[0] / std::max(f.sigma[1], 1e-18);
  if (r01 < sep_ratio_min)
    return f;

  vec3 T = T_cached;
  if (T.norm() < 1e-12)
    T = points.back() - points.front();
  if (T.norm() < 1e-12)
    T = e[2];
  T.normalize();

  int iT = 0;
  real best = -1.0;
  for (int i = 0; i < 3; ++i) {
    const real a = std::abs(e[i].dot(T));
    if (a > best) {
      best = a;
      iT = i;
    }
  }
  f.e_T = (e[iT].dot(T) < 0.0 ? -e[iT] : e[iT]).normalized();
  f.ok = true;
  return f;
}

// Double the cloud by discrete screw pairing (central inversion through origin).
inline std::vector<vec3> screw_symmetrize_points(const std::vector<vec3> &points,
                                                 const screw_frame &f) {
  std::vector<vec3> cloud;
  cloud.reserve(points.size() * 2);
  for (const vec3 &p : points) {
    cloud.push_back(p);
    cloud.push_back(screw_flip(p, f.origin));
  }
  return cloud;
}

inline std::vector<real> assign_t_by_tangent_projection(
    const std::vector<vec3> &cloud, const vec3 &origin, const vec3 &e_T) {
  std::vector<real> t(cloud.size(), 0.0);
  if (cloud.empty())
    return t;

  const vec3 T = e_T.normalized();
  std::vector<std::pair<real, size_t>> keyed;
  keyed.reserve(cloud.size());
  for (size_t i = 0; i < cloud.size(); ++i)
    keyed.push_back({(cloud[i] - origin).dot(T), i});

  std::sort(keyed.begin(), keyed.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  const real s0 = keyed.front().first;
  const real s1 = keyed.back().first;
  if (std::abs(s1 - s0) < 1e-14)
    return {};

  for (size_t rank = 0; rank < keyed.size(); ++rank) {
    const real s = keyed[rank].first;
    t[keyed[rank].second] = (s - s0) / (s1 - s0);
  }
  return t;
}

inline real tangent_t(const vec3 &p, const vec3 &origin, const vec3 &e_T,
                      real s_min, real s_max) {
  if (std::abs(s_max - s_min) < 1e-14)
    return 0.5;
  const real s = (p - origin).dot(e_T.normalized());
  return std::min(1.0, std::max(0.0, (s - s_min) / (s_max - s_min)));
}

inline bool tangent_s_range(const std::vector<vec3> &cloud, const vec3 &origin,
                            const vec3 &e_T, real &s_min, real &s_max) {
  if (cloud.empty())
    return false;
  const vec3 T = e_T.normalized();
  s_min = std::numeric_limits<real>::infinity();
  s_max = -std::numeric_limits<real>::infinity();
  for (const vec3 &q : cloud) {
    const real s = (q - origin).dot(T);
    s_min = std::min(s_min, s);
    s_max = std::max(s_max, s);
  }
  return std::abs(s_max - s_min) > 1e-14;
}

inline bool fit_cubic_bezier(const std::vector<vec3> &cloud,
                             const std::vector<real> &t,
                             std::array<vec3, 4> &P_out) {
  P_out = {vec3::Zero(), vec3::Zero(), vec3::Zero(), vec3::Zero()};
  if (cloud.size() < 4 || cloud.size() != t.size())
    return false;

  Eigen::Matrix<real, 4, 4> AtA = Eigen::Matrix<real, 4, 4>::Zero();
  Eigen::Matrix<real, 4, 3> Atb = Eigen::Matrix<real, 4, 3>::Zero();

  for (size_t j = 0; j < cloud.size(); ++j) {
    const auto b = bernstein3(t[j]);
    for (int a = 0; a < 4; ++a) {
      for (int c = 0; c < 4; ++c)
        AtA(a, c) += b[a] * b[c];
      Atb.row(a) += b[a] * cloud[j].transpose();
    }
  }

  Eigen::LDLT<Eigen::Matrix<real, 4, 4>> ldlt(AtA);
  if (ldlt.info() != Eigen::Success)
    return false;
  Eigen::Matrix<real, 4, 3> X = ldlt.solve(Atb);
  if (!X.allFinite())
    return false;

  for (int i = 0; i < 4; ++i)
    P_out[i] = X.row(i).transpose();
  return true;
}

struct symmetric_bezier_fit {
  screw_frame frame;
  std::array<vec3, 4> P = {vec3::Zero(), vec3::Zero(), vec3::Zero(),
                           vec3::Zero()};
  real s_min = 0.0;
  real s_max = 1.0;
  bool ok = false;

  vec3 eval_at_point(const vec3 &p) const {
    const real tt = tangent_t(p, frame.origin, frame.e_T, s_min, s_max);
    return eval_cubic_bezier(P, tt);
  }
};

inline symmetric_bezier_fit fit_symmetric_cubic_bezier(
    const std::vector<vec3> &points, const vec3 &T_cached,
    real sep_ratio_min = 1.15) {
  symmetric_bezier_fit out;
  out.frame = screw_frame_from_centroid_pca(points, T_cached, sep_ratio_min);
  if (!out.frame.ok)
    return out;

  const std::vector<vec3> cloud = screw_symmetrize_points(points, out.frame);
  if (!tangent_s_range(cloud, out.frame.origin, out.frame.e_T, out.s_min,
                       out.s_max))
    return out;

  const std::vector<real> t =
      assign_t_by_tangent_projection(cloud, out.frame.origin, out.frame.e_T);
  if (t.empty())
    return out;

  if (!fit_cubic_bezier(cloud, t, out.P))
    return out;

  out.ok = true;
  return out;
}

// Alias kept for existing call sites.
inline symmetric_bezier_fit fit_symmetric_quartic_bezier(
    const std::vector<vec3> &points, const vec3 &T_cached,
    real sep_ratio_min = 1.15) {
  return fit_symmetric_cubic_bezier(points, T_cached, sep_ratio_min);
}

} // namespace albers
} // namespace gaudi

#endif
