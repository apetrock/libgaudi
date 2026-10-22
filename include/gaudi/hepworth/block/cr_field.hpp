#ifndef GAUDI_HEPWORTH_CR_FIELD_HPP
#define GAUDI_HEPWORTH_CR_FIELD_HPP

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

/// Compact Crouzeix–Raviart edge topology for the Stein vector Dirichlet operator.
struct cr_topology {
  int nE = 0;
  std::vector<int> slot_map;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  std::vector<std::vector<asawa::shell::FaceId>> edge_incident;
  std::vector<std::vector<kusama::edge_triangle_link>> adj;
};

inline bool assemble_cr_topology(asawa::shell::shell &M, cr_topology &topo) {
  topo = cr_topology{};
  if (!M.verts_are_dense_packed())
    return false;
  topo.nE = kusama::build_compact_edge_dof_map(M, topo.slot_map);
  if (topo.nE <= 0)
    return false;
  kusama::build_dof_to_corner(M, topo.slot_map, topo.nE, topo.dof_to_corner);
  kusama::build_edge_incident_faces(M, topo.slot_map, topo.nE,
                                   topo.edge_incident);
  kusama::build_edge_triangle_adjacency(M, topo.slot_map, topo.nE, topo.adj);
  return true;
}

/// Curvature guidance in the Stein layout: parallel block, then perpendicular.
inline Eigen::VectorXd cr_curvature_guidance(const asawa::shell::shell &M,
                                            const std::vector<vec3> &x,
                                            const cr_topology &topo) {
  Eigen::VectorXd g = kusama::build_edge_curvature_guidance(
      M, x, topo.slot_map, topo.nE, topo.dof_to_corner, topo.edge_incident,
      asawa::shell::face_curvature_stencil::one_ring);
  std::vector<vec3> edge_n(static_cast<size_t>(topo.nE), vec3::UnitZ());
  for (int e = 0; e < topo.nE; ++e)
    edge_n[static_cast<size_t>(e)] = kusama::edge_average_normal(
        M, x, topo.edge_incident[static_cast<size_t>(e)]);
  kusama::orient_edge_guidance_sign_coherence(
      M, x, topo.nE, topo.dof_to_corner, edge_n, topo.adj, g);
  return g;
}

/// Parallel and perpendicular units for one edge, the same split the
/// coefficients are stored in.
inline bool cr_edge_frame(const asawa::shell::shell &M, const std::vector<vec3> &x,
                          asawa::shell::CornerId c0,
                          const std::vector<asawa::shell::FaceId> &faces,
                          vec3 &e_par, vec3 &e_perp) {
  if (static_cast<int>(c0) < 0)
    return false;
  e_par = kusama::edge_cr_parallel_3d(M, c0, x);
  const vec3 N = kusama::edge_average_normal(M, x, faces);
  e_perp = kusama::edge_cr_perp_3d(M, c0, x, N);
  return e_par.allFinite() && e_perp.allFinite();
}

/// Polar factor of `Σ cur restᵀ`. Same determinant fix as the diamond ARAP
/// rotation: `R = UVᵀ`, and `U.col(2)` flips when `det R < 0`.
inline bool diamond_polar(const vec3 rest[4], const vec3 cur[4], mat3 &R) {
  mat3 M = mat3::Zero();
  int used = 0;
  for (int k = 0; k < 4; ++k) {
    if (!rest[k].allFinite() || !cur[k].allFinite())
      continue;
    if (!(rest[k].squaredNorm() > real(1e-20)) ||
        !(cur[k].squaredNorm() > real(1e-20)))
      continue;
    M += cur[k] * rest[k].transpose();
    ++used;
  }
  if (used < 3 || !M.allFinite() || !(M.squaredNorm() > real(1e-16)))
    return false;
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  mat3 V = svd.matrixV();
  R = U * V.transpose();
  if (R.determinant() < real(0)) {
    U.col(2) *= real(-1);
    R = U * V.transpose();
  }
  return R.allFinite() && R.determinant() > real(0);
}

inline void diamond_spokes_at(const vec3 &p0, const vec3 &p1, const vec3 &a,
                             const vec3 &b, vec3 spoke[4]) {
  const vec3 mid = real(0.5) * (p0 + p1);
  spoke[0] = p0 - mid;
  spoke[1] = p1 - mid;
  spoke[2] = a - mid;
  spoke[3] = b - mid;
}

/// Child tangents of one split edge. `out[0]` and `out[1]` are the halves
/// toward `p0` and `p1`. `out[2]` and `out[3]` are the spokes toward `a` and
/// `b`. Each is `parent` rotated by the polar factor of that child's diamond
/// against the parent diamond.
inline void split_diamond_tangents(const vec3 &p0, const vec3 &p1, const vec3 &a,
                                  const vec3 &b, real s, const vec3 &parent,
                                  vec3 out[4]) {
  const vec3 n = (real(1) - s) * p0 + s * p1;
  vec3 rest[4];
  diamond_spokes_at(p0, p1, a, b, rest);
  const vec3 u[4] = {p0, n, n, n};
  const vec3 v[4] = {n, p1, a, b};
  const vec3 c[4] = {a, a, p0, p0};
  const vec3 d[4] = {b, b, p1, p1};
  for (int i = 0; i < 4; ++i) {
    vec3 cur[4];
    diamond_spokes_at(u[i], v[i], c[i], d[i], cur);
    mat3 R = mat3::Identity();
    out[i] = diamond_polar(rest, cur, R) ? R * parent : parent;
    if (!out[i].allFinite())
      out[i] = parent;
  }
}

/// Solved Dirichlet tangent, one vector per edge slot. A split writes the
/// polar of the parent diamond onto the two halves and the two spokes.
/// Pack renumbers the slots with the corners.
struct edge_tangent_datum : public asawa::datum_t<vec3> {
  DEFINE_CREATE_FUNC(edge_tangent_datum)

  edge_tangent_datum(asawa::prim_type type, const std::vector<vec3> &data)
      : asawa::datum_t<vec3>(type, data) {}

  void subdivide(const asawa::shell::shell &M, const index_t &i0,
                 const index_t &i1, const index_t &i2, const index_t &i3,
                 const real &s, const index_t &source_corner) override {
    if (__type != asawa::EDGE)
      return;
    auto &d = data();
    auto at = [&](index_t i) -> vec3 {
      if (i < 0 || static_cast<size_t>(i) >= d.size())
        return vec3::Zero();
      return d[static_cast<size_t>(i)];
    };
    const vec3 parent = at(i0);
    auto put = [&](index_t i, const vec3 &v) {
      if (i < 0)
        return;
      if (static_cast<size_t>(i) >= d.size())
        d.resize(static_cast<size_t>(i) + 1, vec3::Zero());
      d[static_cast<size_t>(i)] = v;
    };
    if (!(parent.squaredNorm() > real(1e-16)) || !parent.allFinite()) {
      put(i0, vec3::Zero());
      put(i1, vec3::Zero());
      put(i2, vec3::Zero());
      put(i3, vec3::Zero());
      return;
    }
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    int v0 = -1;
    int v1 = -1;
    int opp0 = -1;
    int opp1 = -1;
    if (source_corner >= 0) {
      const asawa::shell::CornerId c0 = asawa::shell::corner_id(source_corner);
      const asawa::shell::CornerId c1 = M.other(c0);
      v0 = static_cast<int>(M.vert(c0));
      v1 = static_cast<int>(M.vert(c1));
      auto opposite = [&](asawa::shell::CornerId c, int &opp) {
        opp = -1;
        if (static_cast<int>(c) < 0)
          return;
        const asawa::shell::FaceId f = M.face(c);
        if (static_cast<int>(f) < 0 || M.fsize(f) != 3)
          return;
        M.const_for_each_face(
            f, [&](asawa::shell::CornerId ci, const asawa::shell::shell &Ms) {
              const int v = static_cast<int>(Ms.vert(ci));
              if (v != v0 && v != v1)
                opp = v;
            });
      };
      opposite(c0, opp0);
      opposite(c1, opp1);
    }
    const bool verts =
        v0 >= 0 && v1 >= 0 && opp0 >= 0 && opp1 >= 0 &&
        v0 < static_cast<int>(x.size()) && v1 < static_cast<int>(x.size()) &&
        opp0 < static_cast<int>(x.size()) && opp1 < static_cast<int>(x.size());
    vec3 child[4] = {parent, parent, parent, parent};
    if (verts && std::isfinite(s)) {
      split_diamond_tangents(x[static_cast<size_t>(v0)],
                            x[static_cast<size_t>(v1)],
                            x[static_cast<size_t>(opp0)],
                            x[static_cast<size_t>(opp1)], s, parent, child);
    }
    put(i0, child[0]);
    put(i1, child[1]);
    put(i2, child[2]);
    put(i3, child[3]);
  }
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
