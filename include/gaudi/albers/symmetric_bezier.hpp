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

struct mirror_frame {
  vec3 origin = vec3::Zero();
  vec3 e_T = vec3::UnitX();      // along-curve
  vec3 n_mirror = vec3::UnitY(); // reflection normal (left/right)
  vec3 n_plane = vec3::UnitZ();  // out-of-plane
  real sigma[3] = {0, 0, 0};
  bool ok = false;
};

inline vec3 reflect_plane(const vec3 &p, const vec3 &origin, const vec3 &n) {
  const vec3 nn = n.normalized();
  const vec3 d = p - origin;
  return p - 2.0 * d.dot(nn) * nn;
}

inline std::array<real, 5> bernstein4(real t) {
  const real u = 1.0 - t;
  const real u2 = u * u;
  const real u3 = u2 * u;
  const real u4 = u3 * u;
  const real t2 = t * t;
  const real t3 = t2 * t;
  const real t4 = t3 * t;
  return {u4, 4.0 * u3 * t, 6.0 * u2 * t2, 4.0 * u * t3, t4};
}

inline vec3 eval_quartic_bezier(const std::array<vec3, 5> &P, real t) {
  const auto b = bernstein4(t);
  return b[0] * P[0] + b[1] * P[1] + b[2] * P[2] + b[3] * P[3] + b[4] * P[4];
}

// Centroid-PCA frame.
// For a mild C-buckle the largest principal axis is left/right = mirror normal,
// and is also the natural axis to sort along for t (open the arc left→right).
// Cached T only disambiguates sign / continuity (and a secondary out-of-plane axis).
inline mirror_frame mirror_plane_from_centroid_pca(
    const std::vector<vec3> &points, const vec3 &T_cached,
    real sep_ratio_min = 1.15) {
  mirror_frame f;
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

  // Ascending eigenvalues: col(0) smallest.
  const vec3 s = es.eigenvalues();
  if (s.hasNaN() || s[2] < 1e-18)
    return f;

  std::array<vec3, 3> e = {es.eigenvectors().col(0), es.eigenvectors().col(1),
                           es.eigenvectors().col(2)};
  f.sigma[0] = s[2];
  f.sigma[1] = s[1];
  f.sigma[2] = s[0];

  const real r01 = f.sigma[0] / std::max(f.sigma[1], 1e-18);
  const real r12 = f.sigma[1] / std::max(f.sigma[2], 1e-18);
  if (r01 < sep_ratio_min || r12 < sep_ratio_min)
    return f; // near-degenerate / isotropic-in-plane regime

  // Largest variance → mirror normal (and sort axis for t).
  f.n_mirror = e[2].normalized();
  // Sign: prefer agreement with chord implied by T, else keep SVD sign.
  vec3 T = T_cached;
  if (T.norm() < 1e-12)
    T = f.n_mirror;
  T.normalize();
  if (f.n_mirror.dot(T) < 0.0)
    f.n_mirror = -f.n_mirror;

  // Sort / "tangent" axis for parameterization = mirror normal for open C.
  f.e_T = f.n_mirror;

  // Out-of-plane = smallest; bisector-ish = middle.
  f.n_plane = e[0].normalized();
  vec3 bitan = e[1].normalized();
  if (f.e_T.cross(bitan).dot(f.n_plane) < 0.0)
    f.n_plane = -f.n_plane;
  (void)bitan;

  f.ok = true;
  return f;
}

inline std::vector<vec3> symmetrize_points(const std::vector<vec3> &points,
                                           const mirror_frame &f) {
  std::vector<vec3> cloud;
  cloud.reserve(points.size() * 2);
  for (const vec3 &p : points) {
    cloud.push_back(p);
    cloud.push_back(reflect_plane(p, f.origin, f.n_mirror));
  }
  return cloud;
}

// Returns t in [0,1] parallel to cloud (same order). Fails → empty.
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
    // Use continuous normalized s (not rank) so duplicate projections share t.
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

// Unconstrained quartic LS: min Σ ||B(t_j) - q_j||^2 over P0..P4 in R^3.
inline bool fit_quartic_bezier(const std::vector<vec3> &cloud,
                               const std::vector<real> &t,
                               std::array<vec3, 5> &P_out) {
  P_out = {vec3::Zero(), vec3::Zero(), vec3::Zero(), vec3::Zero(),
           vec3::Zero()};
  if (cloud.size() < 5 || cloud.size() != t.size())
    return false;

  // Normal equations: (B^T B) X = B^T Y, with 5 controls, 3 independent axes.
  Eigen::Matrix<real, 5, 5> AtA = Eigen::Matrix<real, 5, 5>::Zero();
  Eigen::Matrix<real, 5, 3> Atb = Eigen::Matrix<real, 5, 3>::Zero();

  for (size_t j = 0; j < cloud.size(); ++j) {
    const auto b = bernstein4(t[j]);
    for (int a = 0; a < 5; ++a) {
      for (int c = 0; c < 5; ++c)
        AtA(a, c) += b[a] * b[c];
      Atb.row(a) += b[a] * cloud[j].transpose();
    }
  }

  Eigen::LDLT<Eigen::Matrix<real, 5, 5>> ldlt(AtA);
  if (ldlt.info() != Eigen::Success)
    return false;
  Eigen::Matrix<real, 5, 3> X = ldlt.solve(Atb);
  if (!X.allFinite())
    return false;

  for (int i = 0; i < 5; ++i)
    P_out[i] = X.row(i).transpose();
  return true;
}

struct symmetric_bezier_fit {
  mirror_frame frame;
  std::array<vec3, 5> P = {vec3::Zero(), vec3::Zero(), vec3::Zero(),
                           vec3::Zero(), vec3::Zero()};
  real s_min = 0.0;
  real s_max = 1.0;
  bool ok = false;

  vec3 eval_at_point(const vec3 &p) const {
    const real tt = tangent_t(p, frame.origin, frame.e_T, s_min, s_max);
    return eval_quartic_bezier(P, tt);
  }
};

inline symmetric_bezier_fit fit_symmetric_quartic_bezier(
    const std::vector<vec3> &points, const vec3 &T_cached,
    real sep_ratio_min = 1.15) {
  symmetric_bezier_fit out;
  out.frame = mirror_plane_from_centroid_pca(points, T_cached, sep_ratio_min);
  if (!out.frame.ok)
    return out;

  const std::vector<vec3> cloud = symmetrize_points(points, out.frame);
  if (!tangent_s_range(cloud, out.frame.origin, out.frame.e_T, out.s_min,
                       out.s_max))
    return out;

  const std::vector<real> t =
      assign_t_by_tangent_projection(cloud, out.frame.origin, out.frame.e_T);
  if (t.empty())
    return out;

  if (!fit_quartic_bezier(cloud, t, out.P))
    return out;

  out.ok = true;
  return out;
}

} // namespace albers
} // namespace gaudi

#endif
