#ifndef GAUDI_KUSAMA_SHAPE_OPERATOR_HPP
#define GAUDI_KUSAMA_SHAPE_OPERATOR_HPP

/// Mesh-wide Stein shape operator.
/// Stein, Wardetzky, Jacobson, Grinspun 2020, eq. (7) and Appendix A.2
/// (`reference/stein2020-vector-dirichlet-energy.pdf`).
/// One G per edge, one S per face, then the area average at each vertex
/// in the Dirichlet frame of `coeff`.

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/kusama/vector_dirichlet.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"

#include <vector>

namespace gaudi {
namespace kusama {

/// Writes `(k11, k12, k22)` into `H`. `coeff` is parallel on `[0, nE)` and
/// perpendicular on `[nE, 2 nE)`. `wrote`, when set, is 1 at each vertex
/// that received a frame. Vertices without a frame are left unchanged.
inline bool solve_shape_operator(asawa::shell::shell &shell,
                                 const std::vector<vec3> &x,
                                 const std::vector<real> &coeff,
                                 std::vector<vec3> &H,
                                 std::vector<char> *wrote = nullptr) {
  if (!shell.verts_are_dense_packed())
    return false;
  const int n = shell.vert_count();
  if (static_cast<int>(x.size()) != n)
    return false;
  if (static_cast<int>(H.size()) < n)
    H.resize(static_cast<size_t>(n), vec3::Zero());
  if (wrote)
    wrote->assign(static_cast<size_t>(n), 0);

  std::vector<int> slot_map;
  const int nE = build_compact_edge_dof_map(shell, slot_map);
  if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
    return false;
  std::vector<asawa::shell::CornerId> dof_to_corner;
  build_dof_to_corner(shell, slot_map, nE, dof_to_corner);

  std::vector<hepworth::block::cr_tri> tris;
  std::vector<int> face_local(static_cast<size_t>(shell.face_count()), -1);
  for (asawa::shell::FaceId fi : shell.get_face_range(true)) {
    if (shell.fsize(fi) != 3)
      continue;
    hepworth::block::cr_tri tri;
    bool filled = false;
    shell.const_for_each_face_tri(
        fi, [&](asawa::shell::CornerId c0, asawa::shell::CornerId c1,
                asawa::shell::CornerId c2, const asawa::shell::shell &Ms) {
          if (filled)
            return;
          filled = true;
          const asawa::shell::CornerId c[3] = {c0, c1, c2};
          for (int i = 0; i < 3; ++i) {
            tri.v[static_cast<size_t>(i)] = static_cast<int>(Ms.vert(c[i]));
            const int slot = static_cast<int>(c[i]) / 2;
            tri.e[static_cast<size_t>(i)] =
                slot >= 0 && slot < static_cast<int>(slot_map.size())
                    ? slot_map[static_cast<size_t>(slot)]
                    : -1;
          }
        });
    face_local[static_cast<size_t>(static_cast<int>(fi))] =
        static_cast<int>(tris.size());
    tris.push_back(tri);
  }
  if (tris.empty())
    return false;

  std::vector<hepworth::block::cr_edge> edges(static_cast<size_t>(nE));
  for (int e = 0; e < nE; ++e) {
    const asawa::shell::CornerId c = dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    edges[static_cast<size_t>(e)].v0 = static_cast<int>(shell.vert(c));
    edges[static_cast<size_t>(e)].v1 =
        static_cast<int>(shell.vert(shell.next(c)));
  }

  const std::vector<hepworth::block::cr_face_operator> faceS =
      hepworth::block::cr_face_operators(x, tris, edges);

  for (int i = 0; i < n; ++i) {
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    if (static_cast<int>(shell.vbegin(vi)) < 0)
      continue;
    std::vector<int> faces;
    std::vector<int> edge_ids;
    shell.const_for_each_vertex(
        vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &Ms) {
          const int slot = static_cast<int>(ci) / 2;
          if (slot >= 0 && slot < static_cast<int>(slot_map.size())) {
            const int e = slot_map[static_cast<size_t>(slot)];
            if (e >= 0)
              edge_ids.push_back(e);
          }
          const int lj = static_cast<int>(Ms.face(ci));
          if (lj >= 0 && lj < static_cast<int>(face_local.size())) {
            const int local = face_local[static_cast<size_t>(lj)];
            if (local >= 0)
              faces.push_back(local);
          }
        });
    vec3 N;
    if (!hepworth::block::cr_fan_normal(x, tris, faces, i, N))
      continue;
    vec3 T;
    if (!hepworth::block::cr_dirichlet_tangent(x, edges, coeff, edge_ids, N, T,
                                              nullptr)) {
      const vec3 a =
          std::abs(N.x()) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
      T = N.cross(a);
      if (!(T.squaredNorm() > real(1e-24)))
        continue;
      T.normalize();
    }
    vec3 B = N.cross(T);
    if (!(B.squaredNorm() > real(1e-24)))
      continue;
    B.normalize();
    T = B.cross(N);
    const mat3 S = hepworth::block::cr_vertex_operator(faceS, faces);
    const real k11 = T.dot(S * T);
    const real k12 = T.dot(S * B);
    const real k22 = B.dot(S * B);
    if (!std::isfinite(k11) || !std::isfinite(k12) || !std::isfinite(k22))
      continue;
    H[static_cast<size_t>(i)] = vec3(k11, k12, k22);
    if (wrote)
      (*wrote)[static_cast<size_t>(i)] = 1;
  }
  return true;
}

} // namespace kusama
} // namespace gaudi

#endif
