#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

#include <Eigen/Eigenvalues>

#include "gaudi/common.h"
#include "gaudi/duchamp/demo_trait.hpp"

namespace gaudi {
namespace vermeer {

struct oriented_bbox {
  vec3 com = vec3::Zero();    // sample mean (look-at)
  vec3 center = vec3::Zero(); // PCA AABB center (extents / sphere)
  // Ascending eigenvalues: axes.col(0) smallest extent direction.
  mat3 axes = mat3::Identity();
  vec3 extents = vec3::Zero(); // half-sizes along axes
  real sphere_radius = 0.0;
  bool ok = false;
};

inline void append_snapshot_points(const std::optional<duchamp::mesh_snapshot> &snap,
                                   std::vector<vec3> &pts) {
  if (!snap)
    return;
  pts.insert(pts.end(), snap->positions.begin(), snap->positions.end());
}

inline oriented_bbox oriented_bbox_from_points(const std::vector<vec3> &pts) {
  oriented_bbox out;
  if (pts.size() < 3)
    return out;

  vec3 mean = vec3::Zero();
  for (const vec3 &p : pts)
    mean += p;
  mean /= real(pts.size());

  mat3 C = mat3::Zero();
  for (const vec3 &p : pts) {
    const vec3 d = p - mean;
    C += d * d.transpose();
  }
  C /= real(pts.size());

  Eigen::SelfAdjointEigenSolver<mat3> es(C);
  if (es.info() != Eigen::Success)
    return out;

  out.com = mean;
  out.axes = es.eigenvectors();
  if (out.axes.determinant() < 0.0)
    out.axes.col(0) = -out.axes.col(0);

  vec3 bmin = vec3::Constant(std::numeric_limits<real>::infinity());
  vec3 bmax = vec3::Constant(-std::numeric_limits<real>::infinity());
  for (const vec3 &p : pts) {
    const vec3 q = out.axes.transpose() * (p - mean);
    bmin = bmin.cwiseMin(q);
    bmax = bmax.cwiseMax(q);
  }
  out.extents = real(0.5) * (bmax - bmin);
  const vec3 local_c = real(0.5) * (bmin + bmax);
  out.center = mean + out.axes * local_c;

  real r2 = 0.0;
  for (const vec3 &p : pts)
    r2 = std::max(r2, (p - out.com).squaredNorm());
  out.sphere_radius = std::sqrt(r2);
  out.ok = out.sphere_radius > 1.0e-8;
  return out;
}

inline oriented_bbox
oriented_bbox_from_snapshots(const std::optional<duchamp::mesh_snapshot> &shell,
                             const std::optional<duchamp::mesh_snapshot> &rod) {
  std::vector<vec3> pts;
  pts.reserve((shell ? shell->positions.size() : 0) +
              (rod ? rod->positions.size() : 0));
  append_snapshot_points(shell, pts);
  append_snapshot_points(rod, pts);
  return oriented_bbox_from_points(pts);
}

} // namespace vermeer
} // namespace gaudi
