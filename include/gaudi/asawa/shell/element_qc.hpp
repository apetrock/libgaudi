#ifndef GAUDI_ASAWA_SHELL_ELEMENT_QC_HPP
#define GAUDI_ASAWA_SHELL_ELEMENT_QC_HPP

/// Element quality for the remesher: a triangle metric, a census of bad
/// faces / n-gons / valences, and the two repair operations the census needs
/// beyond collapse and flip.
///
/// - `dissolve_vertex` removes a vertex and triangulates its ring polygon.
/// - `retriangulate_face` triangulates an n-gon in place.
///
/// Both pick the triangulation that maximizes the minimum metric (an
/// interval dynamic program over the polygon), and both follow the surgery
/// rule: gather the stencil while the shell is valid, rewrite with raw link /
/// label calls until it is a valid polygon mesh again, and only then call
/// primitives that traverse (`split_face`). A dissolve reuses the corner
/// pairs and faces it frees, so it allocates nothing and datum sizes hold.

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/face_quality.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace gaudi {
namespace asawa {
namespace shell {

/// A triangle score in [0, 1], higher is better, from three vertex values.
using tri_metric = std::function<real(const vec3 &, const vec3 &, const vec3 &)>;

struct qc_config {
  /// VERTEX vec3 datum the metric reads (positions).
  index_t datum = 0;
  /// Repair faces below this.
  real q_warn = real(0.1);
  /// Report faces below this as hard failures.
  real q_min = real(1e-4);
  /// Interior valence band; outside it a vertex is repaired.
  int val_lo = 4;
  int val_hi = 8;
  /// Repair rounds per remesh step.
  int rounds = 3;
  bool enabled = true;
};

struct qc_census {
  int faces = 0;  ///< triangles below `q_warn`
  int hard = 0;   ///< triangles below `q_min`
  int ngons = 0;  ///< faces with more than three corners
  int val_lo = 0; ///< vertices below `val_lo`
  int val_hi = 0; ///< vertices above `val_hi`
  real worst_q = real(1);
  int worst_face = -1;
};

/// One remesh step's QC: a census after each sub-step, the repairs, and the
/// elements left bad. Ids are pre-pack (the handler runs before `pack`).
struct qc_report {
  qc_census subdivide, collapse, merge, flip, final;
  int collapsed = 0, flipped = 0, dissolved = 0, retriangulated = 0;
  /// Bad elements the repairs tried and could not fix.
  int refused = 0;
  std::vector<FaceId> bad_faces;
  std::vector<VertId> bad_verts;
};

inline qc_census census(const shell &M, const std::vector<vec3> &x,
                        const qc_config &cfg, const tri_metric &q) {
  qc_census r;
  for (int i = 0; i < static_cast<int>(M.face_count()); ++i) {
    const CornerId c = M.fbegin(face_id(i));
    if (c < 0)
      continue;
    if (face_len(M, c) != 3) {
      ++r.ngons;
      continue;
    }
    const real qi = q(x[M.vert(c)], x[M.vert(M.next(c))], x[M.vert(M.prev(c))]);
    r.faces += qi < cfg.q_warn;
    r.hard += qi < cfg.q_min;
    if (qi < r.worst_q) {
      r.worst_q = qi;
      r.worst_face = i;
    }
  }
  for (VertId v : M.get_vert_range()) {
    const int n = static_cast<int>(M.vsize(v));
    r.val_lo += n < cfg.val_lo;
    r.val_hi += n > cfg.val_hi;
  }
  return r;
}

/// A polygon triangulation as triangles of polygon indices `(i, k, j)`,
/// `i < k < j`, in split order (parents before children).
using polygon_triangulation = std::vector<std::array<int, 3>>;

/// Best triangulation of the polygon `p` (ring order, counterclockwise about
/// `n`). A triangle is admissible when it is positively oriented about `n`
/// and no other polygon vertex lies inside it in the plane of `n`; a
/// diagonal `(i, j)` only when `diagonal_ok(i, j)` (e.g. the edge does not
/// already exist). Maximizes the minimum `q`, then the sum. Empty when no
/// admissible triangulation exists; `q_out` gets the minimum.
inline polygon_triangulation
best_triangulation(const std::vector<vec3> &p, const vec3 &n, const tri_metric &q,
                   const std::function<bool(int, int)> &diagonal_ok,
                   real *q_out = nullptr) {
  const int m = static_cast<int>(p.size());
  if (q_out)
    *q_out = real(0);
  if (m < 3)
    return {};
  const vec3 nn = n.normalized();
  const vec3 t0 = (std::abs(nn[0]) < real(0.9) ? vec3(1, 0, 0) : vec3(0, 1, 0)).cross(nn).normalized();
  const vec3 t1 = nn.cross(t0);
  std::vector<std::array<real, 2>> uv(static_cast<size_t>(m));
  for (int i = 0; i < m; ++i)
    uv[static_cast<size_t>(i)] = {p[static_cast<size_t>(i)].dot(t0), p[static_cast<size_t>(i)].dot(t1)};
  auto orient = [&](int a, int b, int c) {
    const auto &A = uv[static_cast<size_t>(a)], &B = uv[static_cast<size_t>(b)],
               &C = uv[static_cast<size_t>(c)];
    return (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
  };
  auto admissible = [&](int i, int k, int j) {
    if (!(orient(i, k, j) > real(0)))
      return false;
    for (int s = 0; s < m; ++s) {
      if (s == i || s == k || s == j)
        continue;
      if (orient(i, k, s) >= 0 && orient(k, j, s) >= 0 && orient(j, i, s) >= 0)
        return false;
    }
    return true;
  };
  auto chord_ok = [&](int a, int b) { return b - a == 1 || (a == 0 && b == m - 1) || diagonal_ok(a, b); };

  struct cell {
    real lo = -1, sum = 0;
    int k = -1;
  };
  std::vector<cell> T(static_cast<size_t>(m * m));
  auto at = [&](int i, int j) -> cell & { return T[static_cast<size_t>(i * m + j)]; };
  for (int i = 0; i + 1 < m; ++i)
    at(i, i + 1) = {std::numeric_limits<real>::infinity(), 0, -1};
  for (int len = 2; len < m; ++len) {
    for (int i = 0; i + len < m; ++i) {
      const int j = i + len;
      if (!chord_ok(i, j))
        continue;
      cell best;
      for (int k = i + 1; k < j; ++k) {
        const cell &a = at(i, k), &b = at(k, j);
        if (a.lo < 0 || b.lo < 0 || !admissible(i, k, j))
          continue;
        const real qk = q(p[static_cast<size_t>(i)], p[static_cast<size_t>(k)], p[static_cast<size_t>(j)]);
        const real lo = std::min({a.lo, b.lo, qk});
        const real sum = a.sum + b.sum + qk;
        if (lo > best.lo || (lo == best.lo && sum > best.sum))
          best = {lo, sum, k};
      }
      at(i, j) = best;
    }
  }
  if (at(0, m - 1).lo < 0)
    return {};
  if (q_out)
    *q_out = at(0, m - 1).lo;
  polygon_triangulation out;
  std::vector<std::array<int, 2>> stack{{0, m - 1}};
  while (!stack.empty()) {
    const auto [i, j] = stack.back();
    stack.pop_back();
    const int k = at(i, j).k;
    out.push_back({i, k, j});
    if (j - k > 1)
      stack.push_back({k, j});
    if (k - i > 1)
      stack.push_back({i, k});
  }
  return out;
}

namespace detail {

/// Split the polygon face whose corner starting at polygon vertex `t` is
/// `at[t]` along `tri`, reusing dead pairs / faces from the pools (appending
/// when they run out). Recursive over sub-polygons; every split leaves a
/// valid polygon mesh.
inline void apply_triangulation(shell &M, const polygon_triangulation &tri, size_t &next_tri,
                                std::vector<CornerId> &at, std::vector<CornerId> &pairs,
                                std::vector<FaceId> &faces) {
  const auto [i, k, j] = tri[next_tri++];
  auto take_pair = [&]() {
    if (pairs.empty())
      return corner_id(-1);
    const CornerId c = pairs.back();
    pairs.pop_back();
    return c;
  };
  auto take_face = [&]() {
    if (faces.empty())
      return face_id(-1);
    const FaceId f = faces.back();
    faces.pop_back();
    return f;
  };
  // split_face(c0, c1): the new face holds c0 .. prev(c1) plus other(c0i);
  // the old face keeps c1 .. prev(c0) plus c0i.
  if (k - i > 1) {
    const CornerId c0i = split_face(M, at[static_cast<size_t>(i)], at[static_cast<size_t>(k)],
                                    take_pair(), take_face());
    const CornerId keep_k = at[static_cast<size_t>(k)];
    at[static_cast<size_t>(k)] = M.other(c0i);
    apply_triangulation(M, tri, next_tri, at, pairs, faces);
    at[static_cast<size_t>(k)] = keep_k;
    at[static_cast<size_t>(i)] = c0i;
  }
  if (j - k > 1) {
    const CornerId c0i = split_face(M, at[static_cast<size_t>(k)], at[static_cast<size_t>(j)],
                                    take_pair(), take_face());
    const CornerId keep_j = at[static_cast<size_t>(j)];
    at[static_cast<size_t>(j)] = M.other(c0i);
    apply_triangulation(M, tri, next_tri, at, pairs, faces);
    at[static_cast<size_t>(j)] = keep_j;
    at[static_cast<size_t>(k)] = c0i;
  }
}

} // namespace detail

/// Remove interior vertex `v` and triangulate its ring with the best
/// admissible triangulation. Refuses (returns false, shell untouched) when a
/// fan face is not a triangle, the ring repeats a vertex or wraps onto
/// itself, a ring vertex would drop below valence 3, a needed diagonal
/// already exists, or the best minimum `q` is below `q_floor`.
/// dV = -1, dE = -3, dF = -2.
inline bool dissolve_vertex(shell &M, VertId v, const std::vector<vec3> &x,
                            const tri_metric &q, real q_floor = real(0)) {
  // Gather.
  std::vector<CornerId> spoke;
  M.const_for_each_vertex(v, [&](CornerId c, const shell &) { spoke.push_back(c); });
  const int n = static_cast<int>(spoke.size());
  if (n < 3)
    return false;
  std::vector<CornerId> link(static_cast<size_t>(n), corner_id(-1));
  std::vector<VertId> ring(static_cast<size_t>(n), vert_id(-1));
  std::vector<FaceId> fan(static_cast<size_t>(n), face_id(-1));
  for (int k = 0; k < n; ++k) {
    const CornerId s = spoke[static_cast<size_t>(k)];
    if (face_len(M, s) != 3)
      return false;
    link[static_cast<size_t>(k)] = M.next(s);
    ring[static_cast<size_t>(k)] = M.vert(M.other(s));
    fan[static_cast<size_t>(k)] = M.face(s);
  }
  for (int k = 0; k < n; ++k) {
    const VertId w = ring[static_cast<size_t>(k)];
    if (w == v || M.vsize(w) < 4)
      return false;
    if (std::count(ring.begin(), ring.end(), w) > 1)
      return false;
    const FaceId across = M.face(M.other(link[static_cast<size_t>(k)]));
    if (std::find(fan.begin(), fan.end(), across) != fan.end())
      return false;
  }
  std::vector<vec3> p(static_cast<size_t>(n));
  vec3 normal = vec3::Zero();
  for (int k = 0; k < n; ++k) {
    p[static_cast<size_t>(k)] = x[ring[static_cast<size_t>(k)]];
    normal += (x[ring[static_cast<size_t>(k)]] - x[v]).cross(x[ring[static_cast<size_t>((k + 1) % n)]] - x[v]);
  }
  if (!(normal.squaredNorm() > real(0)))
    return false;
  real q_best = 0;
  const polygon_triangulation tri = best_triangulation(
      p, normal, q,
      [&](int a, int b) {
        return M.find_edge_from_verts(ring[static_cast<size_t>(a)], ring[static_cast<size_t>(b)]) < 0;
      },
      &q_best);
  if (tri.empty() || q_best < q_floor)
    return false;

  for (auto &d : M.get_data())
    d->dissolve(M, static_cast<index_t>(v));

  // Rewrite to one n-gon: link the ring corners, drop spokes, fan faces but
  // one, and v. No traversals until the polygon closes.
  const FaceId f0 = fan[0];
  for (int k = 0; k < n; ++k)
    M.link(link[static_cast<size_t>(k)], link[static_cast<size_t>((k + 1) % n)]);
  for (int k = 0; k < n; ++k)
    M.set_face(link[static_cast<size_t>(k)], f0);
  for (int k = 0; k < n; ++k) {
    M.remove_edge_pair(spoke[static_cast<size_t>(k)]);
    M.set_vbegin(ring[static_cast<size_t>(k)], link[static_cast<size_t>(k)]);
  }
  for (int k = 1; k < n; ++k)
    M.remove_face(fan[static_cast<size_t>(k)]);
  M.remove_vertex(v);

  // Triangulate, reusing the freed spoke pairs and fan faces.
  std::vector<CornerId> pairs(spoke.begin(), spoke.end());
  std::vector<FaceId> faces(fan.begin() + 1, fan.end());
  size_t next_tri = 0;
  detail::apply_triangulation(M, tri, next_tri, link, pairs, faces);
  return true;
}

/// Triangulate face `f` in place with the best admissible triangulation
/// (normal from the polygon's vector area). New pairs and faces are
/// appended, which does not grow EDGE / FACE datums. Returns false (shell
/// untouched) when `f` is already a triangle, its cycle is broken, or no
/// admissible triangulation exists.
inline bool retriangulate_face(shell &M, FaceId f, const std::vector<vec3> &x,
                               const tri_metric &q) {
  const CornerId c0 = M.fbegin(f);
  const int n = c0 < 0 ? -1 : face_len(M, c0);
  if (n <= 3)
    return false;
  std::vector<CornerId> at(static_cast<size_t>(n), corner_id(-1));
  std::vector<vec3> p(static_cast<size_t>(n));
  std::vector<VertId> vs(static_cast<size_t>(n), vert_id(-1));
  CornerId c = c0;
  vec3 normal = vec3::Zero();
  for (int k = 0; k < n; ++k, c = M.next(c)) {
    at[static_cast<size_t>(k)] = c;
    vs[static_cast<size_t>(k)] = M.vert(c);
    p[static_cast<size_t>(k)] = x[M.vert(c)];
  }
  for (int k = 0; k < n; ++k)
    normal += p[static_cast<size_t>(k)].cross(p[static_cast<size_t>((k + 1) % n)]);
  if (!(normal.squaredNorm() > real(0)))
    return false;
  const polygon_triangulation tri = best_triangulation(p, normal, q, [&](int a, int b) {
    return vs[static_cast<size_t>(a)] != vs[static_cast<size_t>(b)] &&
           M.find_edge_from_verts(vs[static_cast<size_t>(a)], vs[static_cast<size_t>(b)]) < 0;
  });
  if (tri.empty())
    return false;
  std::vector<CornerId> pairs;
  std::vector<FaceId> faces;
  size_t next_tri = 0;
  detail::apply_triangulation(M, tri, next_tri, at, pairs, faces);
  return true;
}

/// What collapsing edge `c` would leave, before touching the shell: the
/// minimum `q` over the surviving one-ring with `vert(other(c))` at the edge
/// midpoint (the vertex datum's collapse), the minimum before, whether any
/// of those faces turns over, and the merged valence.
struct collapse_prediction {
  real q_after = 1, q_before = 1;
  bool inverts = false;
  int valence = 0;
};

inline collapse_prediction predict_collapse(const shell &M, CornerId c,
                                            const std::vector<vec3> &x,
                                            const tri_metric &q) {
  collapse_prediction r;
  const VertId v0 = M.vert(c), v1 = M.vert(M.other(c));
  const FaceId fa = M.face(c), fb = M.face(M.other(c));
  const vec3 mid = real(0.5) * (x[v0] + x[v1]);
  r.valence = static_cast<int>(M.vsize(v0) + M.vsize(v1)) - 4;
  auto ring = [&](VertId v) {
    M.const_for_each_vertex(v, [&](CornerId s, const shell &m) {
      const FaceId f = m.face(s);
      if (f == fa || f == fb)
        return;
      const vec3 &a = x[m.vert(m.next(s))], &b = x[m.vert(m.prev(s))];
      const vec3 n0 = (a - x[v]).cross(b - x[v]);
      const vec3 n1 = (a - mid).cross(b - mid);
      r.q_before = std::min(r.q_before, q(x[v], a, b));
      r.q_after = std::min(r.q_after, q(mid, a, b));
      r.inverts |= !(n0.dot(n1) > real(0));
    });
  };
  ring(v0);
  ring(v1);
  return r;
}

/// What flipping edge `c` would leave: the two new faces' minimum `q`, the
/// two old faces' minimum, and whether a new face turns against the quad.
struct flip_prediction {
  real q_after = 0, q_before = 0;
  bool inverts = false;
};

inline flip_prediction predict_flip(const shell &M, CornerId c, const std::vector<vec3> &x,
                                    const tri_metric &q) {
  const vec3 &p0 = x[M.vert(c)], &p1 = x[M.vert(M.other(c))];
  const vec3 &a = x[M.vert(M.prev(c))], &b = x[M.vert(M.prev(M.other(c)))];
  // Old faces (v0, v1, a) and (v1, v0, b); new faces (v1, a, b), (a, v0, b).
  const vec3 n_old = (p1 - p0).cross(a - p0) + (p0 - p1).cross(b - p1);
  flip_prediction r;
  r.q_before = std::min(q(p0, p1, a), q(p1, p0, b));
  r.q_after = std::min(q(p1, a, b), q(a, p0, b));
  r.inverts = !((a - p1).cross(b - p1).dot(n_old) > real(0)) ||
              !((p0 - a).cross(b - a).dot(n_old) > real(0));
  return r;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
