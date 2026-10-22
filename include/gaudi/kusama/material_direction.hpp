#ifndef GAUDI_KUSAMA_MATERIAL_DIRECTION_HPP
#define GAUDI_KUSAMA_MATERIAL_DIRECTION_HPP

/// A unit tangent direction per vertex that moves with the surface instead of
/// being re-read from the geometry each frame (a fiber direction). Kept as a
/// vec3 vertex datum so remeshing interpolates it. Per step:
///   `transport` turns each direction with its one ring's rotation over the
///   step; `smooth` is a short vector heat flow (Sharp, Soliman, Crane 2019)
///   with the Stein et al. 2020 vector Dirichlet energy on the CR edge
///   unknowns, which removes the remesh and transport noise.
/// The field is a vector field: interpolation and smoothing assume neighbors
/// agree in sign (`orient` after seeding from a sign-free axis). A line field
/// that cannot be oriented around some loop leaves a seam where the smoothed
/// field cancels.

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/common.h"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"

#include <Eigen/SVD>
#include <deque>
#include <vector>

namespace gaudi {
namespace kusama {
namespace material_direction {

/// `v` tangent to the unit normal `n`, unit length; false (and `v` unusable)
/// when the tangent part is not longer than `min_len`.
inline bool tangent_unit(vec3 &v, const vec3 &n, real min_len = real(1e-12)) {
  if (!v.allFinite() || !n.allFinite())
    return false;
  v -= n * v.dot(n);
  const real len = v.norm();
  if (!(len > min_len))
    return false;
  v /= len;
  return true;
}

inline vec3 unit_vert_normal(const asawa::shell::shell &M, int i, const std::vector<vec3> &x) {
  vec3 n = asawa::shell::vert_normal(M, asawa::shell::vert_id(i), x);
  const real len = n.norm();
  return len > real(1e-20) && n.allFinite() ? vec3(n / len) : vec3::Zero();
}

/// Flips entries so each agrees in sign with the vertex it was reached from,
/// breadth-first over the one rings from the lowest vertex of each component.
/// Zero entries are skipped and do not propagate.
inline void orient(const asawa::shell::shell &M, std::vector<vec3> &u) {
  const int n = static_cast<int>(std::min<size_t>(u.size(), M.vert_count()));
  std::vector<char> seen(static_cast<size_t>(n), 0);
  std::deque<int> queue;
  for (int root = 0; root < n; ++root) {
    if (seen[static_cast<size_t>(root)] || !(u[static_cast<size_t>(root)].squaredNorm() > 0))
      continue;
    seen[static_cast<size_t>(root)] = 1;
    queue.push_back(root);
    while (!queue.empty()) {
      const int i = queue.front();
      queue.pop_front();
      const vec3 ui = u[static_cast<size_t>(i)];
      for (asawa::shell::VertId vj : M.get_one_ring(asawa::shell::vert_id(i))) {
        const int j = static_cast<int>(vj);
        if (j < 0 || j >= n || seen[static_cast<size_t>(j)])
          continue;
        vec3 &uj = u[static_cast<size_t>(j)];
        if (!(uj.squaredNorm() > 0))
          continue;
        if (uj.dot(ui) < real(0))
          uj = -uj;
        seen[static_cast<size_t>(j)] = 1;
        queue.push_back(j);
      }
    }
  }
}

/// Rotation of the polar decomposition of `A` (no reflection); identity when
/// `A` is rank deficient.
inline mat3 polar_rotation(const mat3 &A) {
  Eigen::JacobiSVD<mat3> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
  if (!(svd.singularValues()[1] > real(1e-12) * std::max(svd.singularValues()[0], real(1e-300))))
    return mat3::Identity();
  mat3 D = mat3::Identity();
  D(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant() < real(0) ? real(-1)
                                                                               : real(1);
  return svd.matrixU() * D * svd.matrixV().transpose();
}

/// `u_i <- R_i u_i` with `R_i` the rotation fitting vertex `i`'s spokes from
/// `x0` onto `x1` (`sum_j (x1_j - x1_i)(x0_j - x0_i)^T`, uniform weights),
/// then tangent to the `x1` normal and unit. `x0` and `x1` share the shell's
/// topology. An entry that is zero, or whose ring or normal degenerates,
/// keeps its value.
inline void transport(const asawa::shell::shell &M, const std::vector<vec3> &x0,
                      const std::vector<vec3> &x1, std::vector<vec3> &u) {
  const int n = static_cast<int>(std::min({u.size(), x0.size(), x1.size(),
                                           static_cast<size_t>(M.vert_count())}));
  for (int i = 0; i < n; ++i) {
    const vec3 ui = u[static_cast<size_t>(i)];
    if (!(ui.squaredNorm() > 0) || !ui.allFinite())
      continue;
    mat3 A = mat3::Zero();
    for (asawa::shell::VertId vj : M.get_one_ring(asawa::shell::vert_id(i))) {
      const int j = static_cast<int>(vj);
      if (j < 0 || j >= n)
        continue;
      A += (x1[static_cast<size_t>(j)] - x1[static_cast<size_t>(i)]) *
           (x0[static_cast<size_t>(j)] - x0[static_cast<size_t>(i)]).transpose();
    }
    if (!A.allFinite())
      continue;
    vec3 v = polar_rotation(A) * ui;
    if (tangent_unit(v, unit_vert_normal(M, i, x1)))
      u[static_cast<size_t>(i)] = v;
  }
}

/// `u_i <- normalize((1 - alpha) u_i + alpha g_i)`, `g_i` first flipped to
/// agree with `u_i` (sign free), so `alpha` 0 keeps `u` and 1 takes `g`.
/// Entries where either is zero keep their value.
inline void lean_toward(std::vector<vec3> &u, const std::vector<vec3> &g, real alpha) {
  if (!(alpha > real(0)))
    return;
  alpha = std::min(alpha, real(1));
  for (size_t i = 0; i < u.size() && i < g.size(); ++i) {
    if (!(u[i].squaredNorm() > 0) || !(g[i].squaredNorm() > 0) || !g[i].allFinite())
      continue;
    const vec3 gi = u[i].dot(g[i]) < real(0) ? vec3(-g[i]) : g[i];
    const vec3 v = (real(1) - alpha) * u[i] + alpha * gi;
    if (v.squaredNorm() > real(1e-24))
      u[i] = v.normalized();
  }
}

/// One implicit vector heat step of time `t` (squared length):
/// `(L + M_e / t) w = (M_e / t) g`, `L` the vector Dirichlet stiffness and
/// `M_e` the edge mass. Guidance `g` on an edge is the sum of its two
/// endpoint directions; the solved edge field returns to each vertex as the
/// sum over its incident edges, tangent to the vertex normal, unit. `t <= 0`
/// only makes every entry tangent and unit. A vertex whose sum is under a
/// tenth of the summed edge magnitudes (its edges disagree, near a defect)
/// keeps its value. Needs dense-packed vertices.
inline void smooth(asawa::shell::shell &M, const std::vector<vec3> &x, std::vector<vec3> &u,
                   real t) {
  const int nV = M.vert_count();
  if (static_cast<int>(x.size()) != nV || static_cast<int>(u.size()) < nV)
    return;
  for (int i = 0; i < nV; ++i) {
    vec3 v = u[static_cast<size_t>(i)];
    if (tangent_unit(v, unit_vert_normal(M, i, x)))
      u[static_cast<size_t>(i)] = v;
  }
  if (!(t > real(0)) || !std::isfinite(t) || !M.verts_are_dense_packed())
    return;
  std::vector<int> slot_map;
  const int nE = build_compact_edge_dof_map(M, slot_map);
  if (nE <= 0)
    return;
  const Eigen::SparseMatrix<real> L = laplace::vector_dirichlet::energy(M, x);
  if (L.rows() != 2 * nE)
    return;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  build_dof_to_corner(M, slot_map, nE, dof_to_corner);
  std::vector<std::vector<asawa::shell::FaceId>> incident;
  build_edge_incident_faces(M, slot_map, nE, incident);
  std::vector<vec3> edge_dir(static_cast<size_t>(nE), vec3::Zero());
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    edge_dir[static_cast<size_t>(e)] = u[M.vert(c)] + u[M.vert(M.other(c))];
  }
  const Eigen::VectorXd g = build_edge_vector_guidance(M, x, dof_to_corner, incident, edge_dir);
  const Eigen::VectorXd mass = build_edge_barycentric_mass_diagonal(M, x, nE, dof_to_corner);
  Eigen::VectorXd w;
  try {
    w = solve_guided_vector_dirichlet(L, mass, g, real(1) / t);
  } catch (const std::exception &) {
    return;
  }
  if (w.size() != 2 * nE || !w.allFinite())
    return;
  std::vector<vec3> sum(static_cast<size_t>(nV), vec3::Zero());
  std::vector<real> mag(static_cast<size_t>(nV), real(0));
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    const vec3 n_e = edge_average_normal(M, x, incident[static_cast<size_t>(e)]);
    const vec3 W = edge_guidance_vector_3d(M, c, x, n_e, w(e), w(e + nE));
    for (const int i : {static_cast<int>(M.vert(c)), static_cast<int>(M.vert(M.other(c)))}) {
      sum[static_cast<size_t>(i)] += W;
      mag[static_cast<size_t>(i)] += W.norm();
    }
  }
  for (int i = 0; i < nV; ++i) {
    vec3 v = sum[static_cast<size_t>(i)];
    if (tangent_unit(v, unit_vert_normal(M, i, x), real(0.1) * mag[static_cast<size_t>(i)]))
      u[static_cast<size_t>(i)] = v;
  }
}

} // namespace material_direction
} // namespace kusama
} // namespace gaudi

#endif
