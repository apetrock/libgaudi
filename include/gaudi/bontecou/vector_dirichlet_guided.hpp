#ifndef GAUDI_BONTECOU_VECTOR_DIRICHLET_GUIDED_HPP
#define GAUDI_BONTECOU_VECTOR_DIRICHLET_GUIDED_HPP

/// Guided vector Dirichlet: minimize \f$\tfrac12 u^\top L u + \tfrac{\lambda}{2}(u-g)^\top M (u-g)\f$
/// with diagonal \f$M\f$ from \ref asawa::shell::edge_barycentric_dual_mass (Stein \f$L\f$ from
/// \ref build_vector_dirichlet_energy). Solve \f$(L+\lambda M)u=\lambda M g\f$ via \ref m_solver.

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/bontecou/vector_dirichlet.hpp"
#include "gaudi/common.h"
#include "gaudi/sparse_solver.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gaudi {
namespace bontecou {

struct edge_triangle_link {
  int other_edge = -1;
  asawa::shell::FaceId face = asawa::shell::face_id(-1);
};

/// +parallel direction in 3D for the global Stein DOF on this half-edge (unit).
inline vec3 edge_cr_parallel_3d(const asawa::shell::shell &M,
                                asawa::shell::CornerId c,
                                const std::vector<vec3> &x) {
  vec3 t = x[M.vert(M.next(c))] - x[M.vert(c)];
  real len = t.norm();
  if (len < 1e-20)
    return vec3::UnitX();
  return static_cast<real>(edge_orientation_sign(M, c)) * t / len;
}

/// In-plane perpendicular \f$n\times e_\parallel\f$ (unit); \p n_unit is a unit averaged normal.
inline vec3 edge_cr_perp_3d(const asawa::shell::shell &M,
                          asawa::shell::CornerId c,
                          const std::vector<vec3> &x, const vec3 &n_unit) {
  vec3 e_par = edge_cr_parallel_3d(M, c, x);
  vec3 e_perp = n_unit.cross(e_par);
  real lp = e_perp.norm();
  if (lp < 1e-20)
    return vec3::UnitY();
  return e_perp / lp;
}

inline vec3 edge_guidance_vector_3d(const asawa::shell::shell &M,
                                    asawa::shell::CornerId c,
                                    const std::vector<vec3> &x,
                                    const vec3 &n_unit, real gp, real gq) {
  return gp * edge_cr_parallel_3d(M, c, x) +
         gq * edge_cr_perp_3d(M, c, x, n_unit);
}

inline vec3 project_to_tangent_plane(vec3 v, const vec3 &n_unit) {
  return v - n_unit * v.dot(n_unit);
}

/// Area-weighted average of face normals for faces incident on this edge (unique).
inline vec3 edge_average_normal(
    const asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<asawa::shell::FaceId> &incident_faces) {
  vec3 acc = vec3::Zero();
  for (asawa::shell::FaceId f : incident_faces) {
    if (static_cast<int>(f) < 0)
      continue;
    real a = asawa::shell::face_area(M, f, x);
    acc += a * asawa::shell::face_normal(M, f, x);
  }
  real s = acc.norm();
  if (s < 1e-20)
    return vec3::UnitZ();
  return acc / s;
}

/// For each triangle, record the other two edges (by compact DOF) and the face.
inline void build_edge_triangle_adjacency(
    const asawa::shell::shell &M, const std::vector<int> &slot_map, int nE,
    std::vector<std::vector<edge_triangle_link>> &adj) {
  adj.assign(static_cast<size_t>(nE), {});
  for (asawa::shell::FaceId fi : M.get_face_range(true)) {
    if (M.fsize(fi) != 3)
      continue;
    M.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          const int d0 =
              slot_map[static_cast<size_t>(static_cast<int>(c0) / 2)];
          const int d1 =
              slot_map[static_cast<size_t>(static_cast<int>(c1) / 2)];
          const int d2 =
              slot_map[static_cast<size_t>(static_cast<int>(c2) / 2)];
          adj[static_cast<size_t>(d0)].push_back({d1, fi});
          adj[static_cast<size_t>(d0)].push_back({d2, fi});
          adj[static_cast<size_t>(d1)].push_back({d0, fi});
          adj[static_cast<size_t>(d1)].push_back({d2, fi});
          adj[static_cast<size_t>(d2)].push_back({d0, fi});
          adj[static_cast<size_t>(d2)].push_back({d1, fi});
        });
  }
}

inline void unique_sort_faces(std::vector<asawa::shell::FaceId> &faces) {
  std::sort(faces.begin(), faces.end());
  faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
}

/// Unique incident triangles per physical edge.
inline void build_edge_incident_faces(
    const asawa::shell::shell &M, const std::vector<int> &slot_map, int nE,
    std::vector<std::vector<asawa::shell::FaceId>> &out) {
  out.assign(static_cast<size_t>(nE), {});
  for (asawa::shell::FaceId fi : M.get_face_range(true)) {
    if (M.fsize(fi) != 3)
      continue;
    M.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          out[static_cast<size_t>(
                  slot_map[static_cast<size_t>(static_cast<int>(c0) / 2)])]
              .push_back(fi);
          out[static_cast<size_t>(
                  slot_map[static_cast<size_t>(static_cast<int>(c1) / 2)])]
              .push_back(fi);
          out[static_cast<size_t>(
                  slot_map[static_cast<size_t>(static_cast<int>(c2) / 2)])]
              .push_back(fi);
        });
  }
  for (auto &v : out)
    unique_sort_faces(v);
}

inline void build_dof_to_corner(const asawa::shell::shell &M,
                                const std::vector<int> &slot_map, int nE,
                                std::vector<asawa::shell::CornerId> &dof_to_corner) {
  dof_to_corner.assign(static_cast<size_t>(nE), asawa::shell::corner_id(-1));
  for (asawa::shell::CornerId c : M.get_edge_range()) {
    const int slot = static_cast<int>(c) / 2;
    const int dof = slot_map[static_cast<size_t>(slot)];
    if (dof >= 0 && dof < nE)
      dof_to_corner[static_cast<size_t>(dof)] = c;
  }
}

/// Stacked guidance \f$g\f$: \c g(e) and \c g(e+nE) are parallel / perp Stein components.
inline Eigen::VectorXd build_edge_curvature_guidance(
    const asawa::shell::shell &M, const std::vector<vec3> &x,
    const std::vector<int> &slot_map, int nE,
    const std::vector<asawa::shell::CornerId> &dof_to_corner,
    const std::vector<std::vector<asawa::shell::FaceId>> &edge_incident_faces,
    asawa::shell::face_curvature_stencil stencil) {
  Eigen::VectorXd g = Eigen::VectorXd::Zero(2 * nE);
  for (int e = 0; e < nE; ++e) {
    asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    const vec3 n_avg =
        edge_average_normal(M, x, edge_incident_faces[static_cast<size_t>(e)]);
    real wsum = 0;
    real acc_p = 0, acc_q = 0;
    for (asawa::shell::FaceId f :
         edge_incident_faces[static_cast<size_t>(e)]) {
      if (static_cast<int>(f) < 0)
        continue;
      asawa::shell::face_curvature_frame fc =
          asawa::shell::face_curvature_frame_fit(M, x, f, stencil);
      vec3 d = project_to_tangent_plane(fc.t_min, n_avg);
      real dn = d.norm();
      if (dn < 1e-14)
        continue;
      d /= dn;
      vec3 e_par = edge_cr_parallel_3d(M, c, x);
      vec3 e_perp = edge_cr_perp_3d(M, c, x, n_avg);
      real w = asawa::shell::face_area(M, f, x);
      wsum += w;
      acc_p += w * d.dot(e_par);
      acc_q += w * d.dot(e_perp);
    }
    if (wsum > 1e-20) {
      g(e) = acc_p / wsum;
      g(e + nE) = acc_q / wsum;
    }
  }
  return g;
}

/// BFS from lexicographically smallest vertex pair; flip \f$g\f$ on an edge when its projected
/// 3D direction disagrees with the parent's on their shared triangle.
inline void orient_edge_guidance_sign_coherence(
    const asawa::shell::shell &M, const std::vector<vec3> &x, int nE,
    const std::vector<asawa::shell::CornerId> &dof_to_corner,
    const std::vector<vec3> &edge_n_unit,
    const std::vector<std::vector<edge_triangle_link>> &adj,
    Eigen::VectorXd &g) {
  if (nE <= 0)
    return;
  int root = 0;
  std::pair<int, int> best = {std::numeric_limits<int>::max(),
                              std::numeric_limits<int>::max()};
  for (int e = 0; e < nE; ++e) {
    asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    int a = std::min(static_cast<int>(M.vert(c)),
                     static_cast<int>(M.vert(M.next(c))));
    int b = std::max(static_cast<int>(M.vert(c)),
                     static_cast<int>(M.vert(M.next(c))));
    std::pair<int, int> key = {a, b};
    if (key < best) {
      best = key;
      root = e;
    }
  }

  std::vector<char> vis(static_cast<size_t>(nE), 0);
  std::deque<int> q;
  vis[static_cast<size_t>(root)] = 1;
  q.push_back(root);

  while (!q.empty()) {
    const int e = q.front();
    q.pop_front();
    asawa::shell::CornerId c_e = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c_e) < 0)
      continue;
    vec3 n_e = edge_n_unit[static_cast<size_t>(e)];
    vec3 W_e = edge_guidance_vector_3d(M, c_e, x, n_e, g(e), g(e + nE));

    std::vector<edge_triangle_link> nb = adj[static_cast<size_t>(e)];
    std::sort(nb.begin(), nb.end(), [](const edge_triangle_link &A,
                                       const edge_triangle_link &B) {
      if (A.other_edge != B.other_edge)
        return A.other_edge < B.other_edge;
      return static_cast<int>(A.face) < static_cast<int>(B.face);
    });

    for (const edge_triangle_link &lk : nb) {
      const int f = lk.other_edge;
      if (f < 0 || f >= nE)
        continue;
      if (vis[static_cast<size_t>(f)])
        continue;

      asawa::shell::CornerId c_f = dof_to_corner[static_cast<size_t>(f)];
      if (static_cast<int>(c_f) < 0)
        continue;
      vec3 n_f = edge_n_unit[static_cast<size_t>(f)];
      vec3 W_f =
          edge_guidance_vector_3d(M, c_f, x, n_f, g(f), g(f + nE));

      vec3 nT = asawa::shell::face_normal(M, lk.face, x);
      vec3 We_t = project_to_tangent_plane(W_e, nT);
      vec3 Wf_t = project_to_tangent_plane(W_f, nT);
      real ne = We_t.norm();
      real nf = Wf_t.norm();
      if (ne > 1e-14 && nf > 1e-14 &&
          We_t.dot(Wf_t) < 0) {
        g(f) *= -1;
        g(f + nE) *= -1;
      }

      vis[static_cast<size_t>(f)] = 1;
      q.push_back(f);
    }
  }
}

/// Diagonal entries for \f$2|E|\f$ DOFs (parallel then perp blocks), same mass on both.
inline Eigen::VectorXd build_edge_barycentric_mass_diagonal(
    const asawa::shell::shell &M, const std::vector<vec3> &x, int nE,
    const std::vector<asawa::shell::CornerId> &dof_to_corner) {
  Eigen::VectorXd m = Eigen::VectorXd::Zero(2 * nE);
  for (int e = 0; e < nE; ++e) {
    asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    real me = asawa::shell::edge_barycentric_dual_mass(M, c, x);
    m(e) = me;
    m(e + nE) = me;
  }
  return m;
}

/// \p g layout matches \ref build_vector_dirichlet_energy (parallel block then perp).
inline Eigen::VectorXd solve_guided_vector_dirichlet(
    const Eigen::SparseMatrix<real> &L, const Eigen::VectorXd &mass_diag,
    const Eigen::VectorXd &g, real lambda) {
  const int n = static_cast<int>(L.rows());
  if (L.cols() != n || mass_diag.size() != n || g.size() != n)
    throw std::runtime_error("solve_guided_vector_dirichlet: size mismatch");
  if (lambda <= 0)
    throw std::runtime_error("solve_guided_vector_dirichlet: lambda must be > 0");

  Eigen::VectorXd b = lambda * mass_diag.cwiseProduct(g);
  Eigen::SparseMatrix<real> A = L;
  for (int i = 0; i < n; ++i)
    A.coeffRef(i, i) += lambda * mass_diag(i);

  m_solver solver(A);
  if (!solver.success())
    throw std::runtime_error("solve_guided_vector_dirichlet: factorization failed");
  Eigen::VectorXd b_mut = b;
  return solver.solve(b_mut);
}

/// Face principal direction \c t_min → edge guidance, optional sign pass, then solve.
inline Eigen::VectorXd solve_curvature_guided_vector_dirichlet(
    asawa::shell::shell &M, const std::vector<vec3> &x, real lambda,
    asawa::shell::face_curvature_stencil stencil,
    bool apply_sign_coherence = true) {
  if (!M.verts_are_dense_packed()) {
    throw std::runtime_error(
        "solve_curvature_guided_vector_dirichlet: vertices must be dense-packed");
  }
  std::vector<int> slot_map;
  int nE = build_compact_edge_dof_map(M, slot_map);
  if (nE <= 0)
    return Eigen::VectorXd();

  Eigen::SparseMatrix<real> L = build_vector_dirichlet_energy(M, x, nullptr);
  if (L.rows() != 2 * nE)
    throw std::runtime_error("solve_curvature_guided_vector_dirichlet: L size mismatch");

  std::vector<asawa::shell::CornerId> dof_to_corner;
  build_dof_to_corner(M, slot_map, nE, dof_to_corner);

  std::vector<std::vector<asawa::shell::FaceId>> edge_incident;
  build_edge_incident_faces(M, slot_map, nE, edge_incident);

  Eigen::VectorXd g = build_edge_curvature_guidance(
      M, x, slot_map, nE, dof_to_corner, edge_incident, stencil);

  std::vector<vec3> edge_n(static_cast<size_t>(nE), vec3::UnitZ());
  for (int e = 0; e < nE; ++e)
    edge_n[static_cast<size_t>(e)] =
        edge_average_normal(M, x, edge_incident[static_cast<size_t>(e)]);

  if (apply_sign_coherence) {
    std::vector<std::vector<edge_triangle_link>> adj;
    build_edge_triangle_adjacency(M, slot_map, nE, adj);
    orient_edge_guidance_sign_coherence(M, x, nE, dof_to_corner, edge_n, adj,
                                        g);
  }

  Eigen::VectorXd mass_diag =
      build_edge_barycentric_mass_diagonal(M, x, nE, dof_to_corner);
  return solve_guided_vector_dirichlet(L, mass_diag, g, lambda);
}

} // namespace bontecou
} // namespace gaudi

#endif
