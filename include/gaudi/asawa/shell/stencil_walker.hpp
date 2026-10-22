#ifndef GAUDI_ASAWA_SHELL_STENCIL_WALKER_HPP
#define GAUDI_ASAWA_SHELL_STENCIL_WALKER_HPP

/// Triangle topology for per-face stencils: each face's three vertices in
/// winding order, the vertex across each of its edges in the neighboring
/// face, and a compact id per edge. Local edge `i` of a face runs from vertex
/// `i` to vertex `i + 1`. `kusama::build_stencil_table` consumes any walker.

#include "gaudi/asawa/shell/shell.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <tuple>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

template <typename W>
concept stencil_walker = requires(const W &w, int f, int i) {
  { w.vertex_count() } -> std::convertible_to<int>;
  { w.face_count() } -> std::convertible_to<int>;
  { w.edge_count() } -> std::convertible_to<int>;
  /// Vertex `i` of face `f`, winding order.
  { w.vertex(f, i) } -> std::convertible_to<int>;
  /// Vertex across local edge `i` in the neighboring triangle; -1 on a boundary.
  { w.opposite(f, i) } -> std::convertible_to<int>;
  /// Compact id of local edge `i`.
  { w.edge(f, i) } -> std::convertible_to<int>;
};

/// Walks the corner shell directly (`next`, `other`): one pass over the live
/// triangles, no hashing. Edge ids follow the live-pair slot order, the same
/// numbering as `kusama::build_compact_edge_dof_map`, so they match the CR
/// edge dofs. Faces are the live triangles in id order; `source_face` maps
/// back.
class shell_walker {
public:
  explicit shell_walker(const shell &M) {
    _nv = static_cast<int>(M.vert_count());
    std::vector<int> slot_edge(M.corner_count() / 2, -1);
    for (int i = 0; i < static_cast<int>(M.corner_count()); i += 2)
      if (M.next(corner_id(i)) >= 0)
        slot_edge[static_cast<size_t>(i / 2)] = _ne++;
    for (FaceId fi : M.get_face_range(true)) {
      std::array<CornerId, 3> c{corner_id(-1), corner_id(-1), corner_id(-1)};
      bool got = false;
      M.const_for_each_face_tri(
          fi, [&](CornerId c0, CornerId c1, CornerId c2, const shell &) {
            if (got)
              return;
            c = {c0, c1, c2};
            got = true;
          });
      if (!got)
        continue;
      std::array<int, 3> v{}, o{}, e{};
      for (int i = 0; i < 3; ++i) {
        const CornerId ci = c[static_cast<size_t>(i)];
        v[static_cast<size_t>(i)] = static_cast<int>(M.vert(ci));
        e[static_cast<size_t>(i)] = slot_edge[static_cast<size_t>(static_cast<int>(ci) / 2)];
        const CornerId twin = M.other(ci);
        const FaceId g = M.face(twin);
        o[static_cast<size_t>(i)] =
            static_cast<int>(g) >= 0 && static_cast<int>(M.fbegin(g)) >= 0 &&
                    M.fsize(g) == 3
                ? static_cast<int>(M.vert(M.prev(twin)))
                : -1;
      }
      _v.push_back(v);
      _o.push_back(o);
      _e.push_back(e);
      _source.push_back(fi);
    }
  }

  int vertex_count() const { return _nv; }
  int face_count() const { return static_cast<int>(_v.size()); }
  int edge_count() const { return _ne; }
  int vertex(int f, int i) const { return _v[static_cast<size_t>(f)][static_cast<size_t>(i)]; }
  int opposite(int f, int i) const { return _o[static_cast<size_t>(f)][static_cast<size_t>(i)]; }
  int edge(int f, int i) const { return _e[static_cast<size_t>(f)][static_cast<size_t>(i)]; }
  FaceId source_face(int f) const { return _source[static_cast<size_t>(f)]; }

private:
  int _nv = 0;
  int _ne = 0;
  std::vector<std::array<int, 3>> _v, _o, _e;
  std::vector<FaceId> _source;
};

/// A triangle list (consistent winding): edges paired by sorting their
/// endpoint keys. For tests and uncoupled solves without a shell.
class triangle_walker {
public:
  triangle_walker(int n, std::vector<std::array<int, 3>> tris)
      : _nv(n), _v(std::move(tris)) {
    const int nf = static_cast<int>(_v.size());
    _o.assign(static_cast<size_t>(nf), {-1, -1, -1});
    _e.assign(static_cast<size_t>(nf), {-1, -1, -1});
    std::vector<std::tuple<int, int, int, int>> half; // (lo, hi, face, local)
    half.reserve(static_cast<size_t>(3 * nf));
    for (int f = 0; f < nf; ++f)
      for (int i = 0; i < 3; ++i) {
        const int a = _v[static_cast<size_t>(f)][static_cast<size_t>(i)];
        const int b = _v[static_cast<size_t>(f)][static_cast<size_t>((i + 1) % 3)];
        half.emplace_back(std::min(a, b), std::max(a, b), f, i);
      }
    std::sort(half.begin(), half.end());
    for (size_t k = 0; k < half.size();) {
      size_t k1 = k + 1;
      while (k1 < half.size() && std::get<0>(half[k1]) == std::get<0>(half[k]) &&
             std::get<1>(half[k1]) == std::get<1>(half[k]))
        ++k1;
      const int id = _ne++;
      for (size_t a = k; a < k1; ++a)
        _e[static_cast<size_t>(std::get<2>(half[a]))]
          [static_cast<size_t>(std::get<3>(half[a]))] = id;
      if (k1 - k == 2) {
        for (int s = 0; s < 2; ++s) {
          const auto &me = half[k + static_cast<size_t>(s)];
          const auto &nb = half[k + static_cast<size_t>(1 - s)];
          const int g = std::get<2>(nb);
          const int j = std::get<3>(nb);
          _o[static_cast<size_t>(std::get<2>(me))][static_cast<size_t>(std::get<3>(me))] =
              _v[static_cast<size_t>(g)][static_cast<size_t>((j + 2) % 3)];
        }
      }
      k = k1;
    }
  }

  int vertex_count() const { return _nv; }
  int face_count() const { return static_cast<int>(_v.size()); }
  int edge_count() const { return _ne; }
  int vertex(int f, int i) const { return _v[static_cast<size_t>(f)][static_cast<size_t>(i)]; }
  int opposite(int f, int i) const { return _o[static_cast<size_t>(f)][static_cast<size_t>(i)]; }
  int edge(int f, int i) const { return _e[static_cast<size_t>(f)][static_cast<size_t>(i)]; }

private:
  int _nv = 0;
  int _ne = 0;
  std::vector<std::array<int, 3>> _v, _o, _e;
};

static_assert(stencil_walker<shell_walker>);
static_assert(stencil_walker<triangle_walker>);

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
