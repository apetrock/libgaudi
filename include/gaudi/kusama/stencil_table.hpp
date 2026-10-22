#ifndef GAUDI_KUSAMA_STENCIL_TABLE_HPP
#define GAUDI_KUSAMA_STENCIL_TABLE_HPP

/// Flat per-face stencils for the CR curvature kernels, built once per frame
/// from any `asawa::shell::stencil_walker`. Face `f`'s stencil is its three
/// vertices (winding order) and, in slot `3 + i`, the vertex across local edge
/// `i = (P_i, P_{i+1})` in the neighboring triangle. On a boundary edge the
/// slot is a reflected dummy `P_i + P_{i+1} - P_{i+2}`, which keeps the
/// masked neighbor terms of the generated read-out finite (`h_i = 0`).

#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/common.h"

#include <array>
#include <numeric>
#include <vector>

namespace gaudi {
namespace kusama {

struct stencil_table {
  int n_verts = 0;
  int n_faces = 0;
  int n_edges = 0;
  /// Per face: vertices 0..2, then the opposite vertex of local edge i (-1 on a boundary).
  std::vector<std::array<int, 6>> ids;
  /// Per face: compact edge id of local edge i.
  std::vector<std::array<int, 3>> face_edges;
  /// Per edge: its faces and the edge's local index in each (-1 when absent).
  std::vector<std::array<int, 2>> edge_faces;
  std::vector<std::array<int, 2>> edge_local;
  /// Flat `6 n_faces` gather indices for `permuted_simplex_view<6>` (a
  /// boundary slot points at the face's vertex 0 and is replaced after the
  /// gather) and the identity order over faces.
  std::vector<index_t> adjacency;
  std::vector<index_t> order;

  bool boundary(int f, int i) const {
    return ids[static_cast<size_t>(f)][static_cast<size_t>(3 + i)] < 0;
  }
};

template <asawa::shell::stencil_walker W>
inline stencil_table build_stencil_table(const W &w) {
  stencil_table t;
  t.n_verts = w.vertex_count();
  t.n_faces = w.face_count();
  t.n_edges = w.edge_count();
  t.ids.resize(static_cast<size_t>(t.n_faces));
  t.face_edges.resize(static_cast<size_t>(t.n_faces));
  t.edge_faces.assign(static_cast<size_t>(t.n_edges), {-1, -1});
  t.edge_local.assign(static_cast<size_t>(t.n_edges), {-1, -1});
  t.adjacency.resize(static_cast<size_t>(6 * t.n_faces));
  for (int f = 0; f < t.n_faces; ++f) {
    auto &id = t.ids[static_cast<size_t>(f)];
    for (int i = 0; i < 3; ++i) {
      id[static_cast<size_t>(i)] = w.vertex(f, i);
      id[static_cast<size_t>(3 + i)] = w.opposite(f, i);
      const int e = w.edge(f, i);
      t.face_edges[static_cast<size_t>(f)][static_cast<size_t>(i)] = e;
      if (e < 0 || e >= t.n_edges)
        continue;
      const int s = t.edge_faces[static_cast<size_t>(e)][0] < 0 ? 0 : 1;
      t.edge_faces[static_cast<size_t>(e)][static_cast<size_t>(s)] = f;
      t.edge_local[static_cast<size_t>(e)][static_cast<size_t>(s)] = i;
    }
    for (int k = 0; k < 6; ++k)
      t.adjacency[static_cast<size_t>(6 * f + k)] =
          id[static_cast<size_t>(k)] >= 0 ? id[static_cast<size_t>(k)] : id[0];
  }
  t.order.resize(static_cast<size_t>(t.n_faces));
  std::iota(t.order.begin(), t.order.end(), 0);
  return t;
}

/// Face `f`'s six stencil positions (reflected dummies on boundary edges)
/// and the neighbor flags `h` (1 interior, 0 boundary), the inputs of the
/// generated `cr_curvature_stencil` / `cr_curvature_edge_readings`.
template <Vec3View X>
inline std::array<vec3, 6> stencil_positions(const stencil_table &t, const X &x, int f,
                                             vec3 &h) {
  const permuted_simplex_view<6, X, std::vector<index_t>> view(x, t.adjacency, t.order);
  std::array<vec3, 6> P = view[static_cast<size_t>(f)];
  for (int i = 0; i < 3; ++i) {
    const bool open = t.boundary(f, i);
    h[i] = open ? real(0) : real(1);
    if (open)
      P[static_cast<size_t>(3 + i)] = P[static_cast<size_t>(i)] +
                                      P[static_cast<size_t>((i + 1) % 3)] -
                                      P[static_cast<size_t>((i + 2) % 3)];
  }
  return P;
}

} // namespace kusama
} // namespace gaudi

#endif
