#ifndef GAUDI_CALDER_POINT_TANGENT_POINT_HPP
#define GAUDI_CALDER_POINT_TANGENT_POINT_HPP

/// Tangent-point gradient on a point+tangent cloud (no rod topology).
/// Used by hierarchical TP MG at extracted BVH cuts.

#include "gaudi/arp/arp.h"
#include "gaudi/arp/datums.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/calder/weight_functions.hpp"
#include "gaudi/common.h"

#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <vector>

namespace gaudi {
namespace calder {

/// Cauchy TP ∇ on samples (c, T, w). Builds a fresh T1 BVH over \p c.
/// \p bh_eps ≤ 0 ⇒ leaves-only. Branch far field uses AABB center + pyramid data.
inline std::vector<vec3> tangent_point_gradient_points_cauchy(
    const std::vector<vec3> &c, const std::vector<vec3> &T,
    const std::vector<real> &w, real l0, real p = 3.0, real bh_eps = 0.25) {
  const size_t n = c.size();
  std::vector<vec3> out(n, vec3::Zero());
  if (n < 2)
    return out;

  std::vector<index_t> adj(n);
  std::vector<index_t> ids(n);
  for (size_t i = 0; i < n; ++i) {
    adj[i] = static_cast<index_t>(i);
    ids[i] = static_cast<index_t>(i);
  }

  arp::T1::ptr tree = arp::T1::create(adj, c, 12);
  if (!tree || tree->internal_nodes_.empty()) {
    // Degenerate tree: dense fallback for tiny n.
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        const vec3 dp = c[j] - c[i];
        if (dp.squaredNorm() < real(1e-24))
          continue;
        vec3 Ti = T[i];
        vec3 Tj = T[j];
        if (Ti.squaredNorm() > real(1e-24))
          Ti.normalize();
        if (Tj.squaredNorm() > real(1e-24))
          Tj.normalize();
        const vec3 ci = Ti.cross(dp);
        const vec3 cj = Tj.cross(dp);
        if (ci.squaredNorm() < real(1e-32) || cj.squaredNorm() < real(1e-32))
          continue;
        const vec3 Ni = ci.normalized().cross(Ti).normalized();
        const vec3 Nj = cj.normalized().cross(Tj).normalized();
        const vec3 gj = calc_tangent_point_radius_grad_cauchy(dp, Nj, l0, p);
        const vec3 gi = calc_tangent_point_radius_grad_cauchy(-dp, Ni, l0, p);
        const real wi = w[i];
        const real wj = w[j];
        const vec3 f = real(0.5) * (wi * gi - wj * gj);
        out[i] += f;
        out[j] -= f;
      }
    }
    return out;
  }

  fast_summation<arp::T1> sum(*tree);
  sum.bind(scalar_datum::create(ids, w));
  sum.bind(vec3_datum::create(ids, T));

  const real eps_open = bh_eps;

  auto pair_contrib = [l0, p](const vec3 &pi, const vec3 &pj, real wi, real wj,
                              vec3 Ti, vec3 Tj) -> vec3 {
    if (Ti.squaredNorm() > real(1e-24))
      Ti.normalize();
    else
      return vec3::Zero();
    if (Tj.squaredNorm() > real(1e-24))
      Tj.normalize();
    else
      return vec3::Zero();
    const vec3 dp = pj - pi;
    if (dp.squaredNorm() < real(1e-24))
      return vec3::Zero();
    const vec3 ci = Ti.cross(dp);
    const vec3 cj = Tj.cross(dp);
    if (ci.squaredNorm() < real(1e-32) || cj.squaredNorm() < real(1e-32))
      return vec3::Zero();
    const vec3 Ni = ci.normalized().cross(Ti).normalized();
    const vec3 Nj = cj.normalized().cross(Tj).normalized();
    const vec3 gj = calc_tangent_point_radius_grad_cauchy(dp, Nj, l0, p);
    const vec3 gi = calc_tangent_point_radius_grad_cauchy(-dp, Ni, l0, p);
    return real(0.5) * (wi * gi - wj * gj);
  };

  out = sum.calc<vec3>(
      c,
      [&](const index_t i, const index_t j, const vec3 &pi,
          const std::vector<datum::ptr> &data, Node_Type node_type,
          const arp::T1 &tr) -> vec3 {
        (void)node_type;
        const auto s = tr.leaf_simplex(j);
        const vec3 pj = s[0];
        // Skip self: same original index.
        const index_t orig = tr.get_index(j);
        if (orig == i)
          return vec3::Zero();
        const real wi = w[static_cast<size_t>(i)];
        const real wj = get_data<real>(LEAF, j, 0, data);
        const vec3 Tj = get_data<vec3>(LEAF, j, 1, data);
        return pair_contrib(pi, pj, wi, wj, T[static_cast<size_t>(i)], Tj);
      },
      [&](const index_t i, const index_t j, const vec3 &pi,
          const std::vector<datum::ptr> &data, Node_Type node_type,
          const arp::T1 &tr) -> vec3 {
        (void)node_type;
        const ext::extents_t &ext = tr.bvh_.internal[static_cast<size_t>(j)];
        const vec3 pj = real(0.5) * (ext[0] + ext[1]);
        const real wi = w[static_cast<size_t>(i)];
        const real wj = get_data<real>(BRANCH, j, 0, data);
        vec3 Tj = get_data<vec3>(BRANCH, j, 1, data);
        // Pyramid sums tangents; renormalize for direction.
        if (Tj.squaredNorm() > real(1e-24))
          Tj.normalize();
        return pair_contrib(pi, pj, wi, wj, T[static_cast<size_t>(i)], Tj);
      },
      eps_open);

  return out;
}

/// Kabsch SO(3) for point cloud motion x0 → x1 (sim-time tangent transport).
inline mat3 kabsch_sim_rotation(const std::vector<vec3> &x0,
                                const std::vector<vec3> &x1) {
  const size_t n = std::min(x0.size(), x1.size());
  if (n < 2)
    return mat3::Identity();
  vec3 c0 = vec3::Zero();
  vec3 c1 = vec3::Zero();
  for (size_t i = 0; i < n; ++i) {
    c0 += x0[i];
    c1 += x1[i];
  }
  c0 /= real(n);
  c1 /= real(n);
  mat3 M = mat3::Zero();
  for (size_t i = 0; i < n; ++i)
    M += (x1[i] - c1) * (x0[i] - c0).transpose();
  if (M.squaredNorm() < real(1e-24))
    return mat3::Identity();
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

} // namespace calder
} // namespace gaudi

#endif
