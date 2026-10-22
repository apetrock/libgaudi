#ifndef GAUDI_HEPWORTH_FACE_ELEMENT_HPP
#define GAUDI_HEPWORTH_FACE_ELEMENT_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "gaudi/hepworth/block/cache_compute.hpp"
#include "gaudi/hepworth/block/cr_dirichlet.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// One face plus the triangle across each of its edges. Six vertices, nine
/// edges, when every neighbor is present. `vid[0..2]` is the face winding.
struct face_patch {
  int face = -1;
  std::array<int, 3> neighbor{-1, -1, -1};
  /// Edge dofs of the face, in winding order. `-1` when the patch has no chart.
  std::array<int, 3> edge{-1, -1, -1};
  std::array<int, 6> vid{};
  int nvert = 0;
  std::array<std::array<int, 2>, 9> loc{};
  int nedge = 0;
};

inline bool bake_face_patch(const asawa::shell::shell &M,
                            asawa::shell::FaceId face, face_patch &out) {
  out = face_patch{};
  if (M.fsize(face) != 3)
    return false;
  const std::array<asawa::shell::VertId, 3> tri = M.get_tri(face);
  out.face = static_cast<int>(face);
  auto add_vert = [&](int v) -> int {
    for (int i = 0; i < out.nvert; ++i) {
      if (out.vid[static_cast<size_t>(i)] == v)
        return i;
    }
    if (out.nvert >= 6)
      return -1;
    out.vid[static_cast<size_t>(out.nvert)] = v;
    return out.nvert++;
  };
  auto add_edge = [&](int a, int b) {
    const int ia = add_vert(a);
    const int ib = add_vert(b);
    if (ia < 0 || ib < 0)
      return;
    for (int i = 0; i < out.nedge; ++i) {
      const int u = out.loc[static_cast<size_t>(i)][0];
      const int v = out.loc[static_cast<size_t>(i)][1];
      if ((u == ia && v == ib) || (u == ib && v == ia))
        return;
    }
    if (out.nedge >= 9)
      return;
    out.loc[static_cast<size_t>(out.nedge)] = {ia, ib};
    ++out.nedge;
  };
  int gv[3];
  for (int i = 0; i < 3; ++i)
    gv[i] = static_cast<int>(tri[static_cast<size_t>(i)]);
  for (int i = 0; i < 3; ++i)
    add_edge(gv[i], gv[(i + 1) % 3]);
  for (int i = 0; i < 3; ++i) {
    const asawa::shell::CornerId c = M.find_edge_from_verts(
        asawa::shell::vert_id(gv[i]), asawa::shell::vert_id(gv[(i + 1) % 3]));
    if (static_cast<int>(c) < 0)
      continue;
    for (const asawa::shell::FaceId nb : M.dihedral_face_ids(c)) {
      if (nb == face || M.fsize(nb) != 3)
        continue;
      out.neighbor[static_cast<size_t>(i)] = static_cast<int>(nb);
      const std::array<asawa::shell::VertId, 3> nt = M.get_tri(nb);
      for (int k = 0; k < 3; ++k)
        add_edge(static_cast<int>(nt[static_cast<size_t>(k)]),
                 static_cast<int>(nt[static_cast<size_t>((k + 1) % 3)]));
      break;
    }
  }
  return out.nvert >= 3 && out.nedge >= 3;
}

inline int patch_edge_dof(const asawa::shell::shell &M, const cr_topology &topo,
                          int a, int b) {
  const asawa::shell::CornerId c = M.find_edge_from_verts(
      asawa::shell::vert_id(a), asawa::shell::vert_id(b));
  if (static_cast<int>(c) < 0)
    return -1;
  const int slot = static_cast<int>(c) / 2;
  if (slot < 0 || slot >= static_cast<int>(topo.slot_map.size()))
    return -1;
  return topo.slot_map[static_cast<size_t>(slot)];
}

inline bool face_tri_frame(const vec3 &x0, const vec3 &x1, const vec3 &x2,
                           vec3 &N, vec3 &T, vec3 &B) {
  const vec3 cr = (x1 - x0).cross(x2 - x0);
  if (!(cr.squaredNorm() > real(1e-24)) || !cr.allFinite())
    return false;
  N = cr.normalized();
  vec3 t = x1 - x0;
  t -= N * t.dot(N);
  if (!(t.squaredNorm() > real(1e-16)) || !t.allFinite())
    return false;
  t.normalize();
  vec3 b = N.cross(t);
  if (!(b.squaredNorm() > real(1e-16)) || !b.allFinite())
    return false;
  b.normalize();
  T = t;
  B = b;
  return true;
}

/// Four spokes from the edge midpoint. The ARAP polar of those spokes is the
/// rotation that carries the tangent.
inline void diamond_mid_spokes(const vec3 &p0, const vec3 &p1, const vec3 &a,
                               const vec3 &b, vec3 spoke[4]) {
  const vec3 mid = real(0.5) * (p0 + p1);
  spoke[0] = p0 - mid;
  spoke[1] = p1 - mid;
  spoke[2] = a - mid;
  spoke[3] = b - mid;
}

inline mat3 spoke_covariance(const vec3 cur[4], const vec3 rest[4]) {
  mat3 M = mat3::Zero();
  for (int k = 0; k < 4; ++k)
    M += cur[k] * rest[k].transpose();
  return M;
}

/// Rotates each carried edge tangent by the ARAP spoke rotation of its
/// diamond, about the edge midpoint. Writes the tangent cache. Empty
/// tangents are left empty.
class tangent_cache_compute : public cache_compute {
public:
  struct edge {
    int id = -1;
    int v0 = -1;
    int v1 = -1;
    int o0 = -1;
    int o1 = -1;
    std::array<int, 3> t0{};
    std::array<int, 3> t1{};
    real par_sign = 1;
    real perp_sign = 1;
    vec3 rest[4] = {vec3::Zero(), vec3::Zero(), vec3::Zero(), vec3::Zero()};
    mat3 R0 = mat3::Identity();
    vec3 T0 = vec3::Zero();
  };

  tangent_cache_compute(sim_block::ptr pos, std::vector<real> *u, int nE,
                        std::vector<edge> edges)
      : _pos(std::move(pos)), _u(u), _nE(nE), _edges(std::move(edges)) {}

  void compute(vecX &q) override {
    if (!_pos || !_u || _nE <= 0 || static_cast<int>(_u->size()) < 2 * _nE)
      return;
    for (const edge &e : _edges) {
      if (e.id < 0 || e.id >= _nE)
        continue;
      real &gp = (*_u)[static_cast<size_t>(e.id)];
      real &gq = (*_u)[static_cast<size_t>(e.id + _nE)];
      if (!(gp * gp + gq * gq > real(1e-16)))
        continue;
      const vec3 q0 = _pos->get_vec3(e.v0, q);
      const vec3 p1 = _pos->get_vec3(e.v1, q);
      const vec3 a = _pos->get_vec3(e.o0, q);
      const vec3 b = _pos->get_vec3(e.o1, q);
      vec3 spoke[4];
      diamond_mid_spokes(q0, p1, a, b, spoke);
      mat3 R = mat3::Identity();
      if (!diamond_rotation(spoke_covariance(spoke, e.rest), R))
        continue;
      const vec3 T = R * e.R0.transpose() * e.T0;
      vec3 edge_v = p1 - q0;
      const real el = edge_v.norm();
      vec3 N;
      if (!T.allFinite() || !(el > real(1e-20)) ||
          !diamond_normal(q0, p1, a, b, e.t0, e.t1, e.v0, e.v1, e.o0, e.o1, N))
        continue;
      const vec3 e_par = e.par_sign * edge_v / el;
      vec3 e_perp = N.cross(e_par);
      if (!(e_perp.squaredNorm() > real(1e-24)) || !e_perp.allFinite())
        continue;
      e_perp.normalize();
      e_perp *= e.perp_sign;
      const real ngp = T.dot(e_par);
      const real ngq = T.dot(e_perp);
      if (!std::isfinite(ngp) || !std::isfinite(ngq))
        continue;
      gp = ngp;
      gq = ngq;
    }
  }

  sim_block::ptr _pos;
  std::vector<real> *_u = nullptr;
  int _nE = 0;
  std::vector<edge> _edges;
};

/// `G = M^{-1} D(u)` from the positions in `q`, into the edge cache.
class gradient_cache_compute : public cache_compute {
public:
  gradient_cache_compute(sim_block::ptr pos, std::vector<vec3> *par,
                         std::vector<vec3> *perp, std::vector<cr_tri> tris,
                         std::vector<cr_edge> edges, int nV)
      : _pos(std::move(pos)), _par(par), _perp(perp), _tris(std::move(tris)),
        _edges(std::move(edges)), _nV(nV) {}

  void compute(vecX &q) override {
    if (!_pos || !_par || !_perp || _nV <= 0)
      return;
    std::vector<vec3> x(static_cast<size_t>(_nV), vec3::Zero());
    for (int i = 0; i < _nV; ++i)
      x[static_cast<size_t>(i)] = _pos->get_vec3(i, q);
    const cr_edge_gradient G = cr_coordinate_gradients(x, _tris, _edges);
    *_par = G.par;
    *_perp = G.perp;
  }

  sim_block::ptr _pos;
  std::vector<vec3> *_par = nullptr;
  std::vector<vec3> *_perp = nullptr;
  std::vector<cr_tri> _tris;
  std::vector<cr_edge> _edges;
  int _nV = 0;
};

struct face_pin_slot {
  int face = -1;
  /// Index into the complex triangle list.
  int tri = -1;
  std::array<int, 3> v{-1, -1, -1};
  std::array<int, 3> e{-1, -1, -1};
};

/// One face, as the pin last wrote it. `e1` is the first principal direction.
struct projected_face_frame {
  vec3 origin = vec3::Zero();
  vec3 e1 = vec3::Zero();
  vec3 e2 = vec3::Zero();
  vec3 N = vec3::Zero();
  real k1 = 0;
  real k2 = 0;
  bool ok = false;
};

/// Identity pin. Rebuilds `G` from the positions in `q`, then writes the
/// face operator of that same mesh.
class face_shape_pin : public block_constraint {
public:
  DEFINE_CREATE_FUNC(face_shape_pin)

  face_shape_pin(std::vector<face_pin_slot> faces, std::vector<cr_tri> tris,
                 std::vector<cr_edge> edges, std::vector<vec3> *gpar,
                 std::vector<vec3> *gperp, int nV, real w,
                 std::vector<sim_block::ptr> blocks,
                 const std::vector<real> *coeff = nullptr)
      : block_constraint({}, w, std::move(blocks)), _faces(std::move(faces)),
        _tris(std::move(tris)), _edges(std::move(edges)), _gpar(gpar),
        _gperp(gperp), _nV(nV), _u(coeff) {}

  std::string name() override { return "face_shape_pin"; }

  /// Principal frames of the positions in `q`. Does not write the pin target.
  void record(const vecX &q) { rebuild(q, nullptr); }

  void project(const vecX &q, vecX &p) override { rebuild(q, &p); }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.size() < 2 || !_blocks[1])
      return;
    for (const face_pin_slot &face : _faces) {
      const index_t col = _blocks[1]->get_offset_idx(face.face);
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, col + ax, _w));
      id0 += 3;
    }
  }

  std::vector<face_pin_slot> _faces;
  /// Overwritten by `project` and `record`.
  std::vector<projected_face_frame> _projected;
  std::vector<cr_tri> _tris;
  std::vector<cr_edge> _edges;
  std::vector<vec3> *_gpar = nullptr;
  std::vector<vec3> *_gperp = nullptr;
  int _nV = 0;
  /// Parallel on `[0, nE)`, perpendicular on `[nE, 2 nE)`. Empty skips the chart.
  const std::vector<real> *_u = nullptr;

private:
  void rebuild(const vecX &q, vecX *p) {
    if (_blocks.empty() || !_blocks[0] || !_gpar || !_gperp || _nV <= 0)
      return;
    if (p && (_blocks.size() < 2 || !_blocks[1]))
      return;
    std::vector<vec3> x(static_cast<size_t>(_nV), vec3::Zero());
    for (int i = 0; i < _nV; ++i)
      x[static_cast<size_t>(i)] = _blocks[0]->get_vec3(i, q);
    const int nE = static_cast<int>(_edges.size());
    const bool cached =
        _gpar && _gperp && nE > 0 &&
        static_cast<int>(_gpar->size()) >= nE &&
        static_cast<int>(_gperp->size()) >= nE;
    if (!cached) {
      const cr_edge_gradient G = cr_coordinate_gradients(x, _tris, _edges);
      *_gpar = G.par;
      *_gperp = G.perp;
    }
    _projected.assign(_faces.size(), projected_face_frame{});
    for (int k = 0; k < static_cast<int>(_faces.size()); ++k) {
      const face_pin_slot &face = _faces[static_cast<size_t>(k)];
      vec3 coeff = vec3::Zero();
      if (p && _blocks.size() > 1 && _blocks[1])
        coeff = _blocks[1]->get_vec3(face.face, q);
      vec3 N, T, B;
      vec3 next;
      const bool got = face_coeff_from_gradients(
          x, face.v, face.e, _edges, *_gpar, *_gperp, next, N, T, B);
      if (got && next.allFinite())
        coeff = next;
      if (!coeff.allFinite())
        coeff = vec3::Zero();
      if (p)
        p->block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = _w * coeff;
      if (!got || !next.allFinite())
        continue;
      const shape_w W = shape_with_axes(shape_w_from_coeffs(coeff));
      if (!W.ok)
        continue;
      vec3 e1 = W.e1x * T + W.e1y * B;
      vec3 e2 = W.e2x * T + W.e2y * B;
      if (!(e1.squaredNorm() > real(1e-16)) || !(e2.squaredNorm() > real(1e-16)))
        continue;
      e1.normalize();
      e2.normalize();
      const vec3 &x0 = x[static_cast<size_t>(face.v[0])];
      const vec3 &x1 = x[static_cast<size_t>(face.v[1])];
      const vec3 &x2 = x[static_cast<size_t>(face.v[2])];
      projected_face_frame &drawn = _projected[static_cast<size_t>(k)];
      drawn.origin = (x0 + x1 + x2) / real(3);
      drawn.e1 = e1;
      drawn.e2 = e2;
      drawn.N = N;
      drawn.k1 = W.k1;
      drawn.k2 = W.k2;
      drawn.ok = true;
    }
  }
};

/// `p[0..2]` is the face. The foot is their centroid, on the triangle, the
/// same role as the vertex in `frame_fair`. The opposite vertices stay in
/// `p` and are the spokes.
inline bool patch_positions(const face_patch &patch, const sim_block &pos,
                            const vecX &q, std::vector<vec3> &p, vec3 &center,
                            vec3 &N, vec3 &T, vec3 &B) {
  if (patch.nvert < 3)
    return false;
  p.resize(static_cast<size_t>(patch.nvert));
  for (int i = 0; i < patch.nvert; ++i)
    p[static_cast<size_t>(i)] = pos.get_vec3(patch.vid[static_cast<size_t>(i)], q);
  center = (p[0] + p[1] + p[2]) / real(3);
  return face_tri_frame(p[0], p[1], p[2], N, T, B);
}

inline void write_spoke_targets(const face_patch &patch,
                                const std::vector<vec3> &mapped, index_t id0,
                                real w, vecX &p) {
  for (int k = 0; k < patch.nedge; ++k) {
    const int ia = patch.loc[static_cast<size_t>(k)][0];
    const int ib = patch.loc[static_cast<size_t>(k)][1];
    const vec3 d =
        mapped[static_cast<size_t>(ib)] - mapped[static_cast<size_t>(ia)];
    p.block(id0 + static_cast<index_t>(3 * k), 0, 3, 1) = w * d;
  }
}

inline void fill_spoke_rows(const face_patch &patch, const sim_block &pos,
                            index_t &id0, real w, std::vector<trip> &triplets) {
  for (int k = 0; k < patch.nedge; ++k) {
    const int ia = patch.loc[static_cast<size_t>(k)][0];
    const int ib = patch.loc[static_cast<size_t>(k)][1];
    const index_t ca = pos.get_offset_idx(patch.vid[static_cast<size_t>(ia)]);
    const index_t cb = pos.get_offset_idx(patch.vid[static_cast<size_t>(ib)]);
    for (int ax = 0; ax < 3; ++ax) {
      triplets.push_back(trip(id0 + ax, ca + ax, -w));
      triplets.push_back(trip(id0 + ax, cb + ax, w));
    }
    id0 += 3;
  }
}

/// Canal or Monge landing of the opposite vertices, then one ARAP spoke
/// per patch edge. The face triangle is the foot and is not projected.
class face_fair : public block_constraint {
public:
  DEFINE_CREATE_FUNC(face_fair)

  face_fair(face_patch patch, real w, std::vector<objective> objectives,
            std::vector<sim_block::ptr> blocks,
            const std::vector<real> *coeff = nullptr,
            std::shared_ptr<const std::vector<cr_edge>> edges = nullptr)
      : block_constraint({}, w, std::move(blocks)), _patch(patch),
        _objectives(std::move(objectives)), _u(coeff),
        _edges(std::move(edges)) {
    _spines = spine_set_of(_objectives);
    for (const objective &o : _objectives) {
      if (std::holds_alternative<monge_fairing>(o))
        _monge = true;
    }
  }

  std::string name() override { return "face_fair"; }

  void project(const vecX &q, vecX &p) override {
    if (_blocks.size() < 2 || !_blocks[0] || !_blocks[1])
      return;
    std::vector<vec3> pos;
    vec3 center, N, T, B;
    if (!patch_positions(_patch, *_blocks[0], q, pos, center, N, T, B))
      return;
    const vec3 coeff = _blocks[1]->get_vec3(_patch.face, q);
    const shape_w goal = shape_with_axes(shape_w_from_coeffs(coeff));
    std::vector<vec3> mapped = pos;
    if (goal.ok) {
      std::vector<vec3> outer;
      std::vector<int> outer_at;
      outer.reserve(static_cast<size_t>(std::max(_patch.nvert - 3, 0)));
      for (int i = 3; i < _patch.nvert; ++i) {
        outer.push_back(pos[static_cast<size_t>(i)] - center);
        outer_at.push_back(i);
      }
      bool spined = false;
      if (_spines.any() && !outer.empty()) {
        const monge_axes ax = monge_axes_of(goal);
        spine_curvature sc;
        sc.tube = ax.hx * T + ax.hy * B;
        sc.ruling = ax.lx * T + ax.ly * B;
        sc.normal = N;
        sc.k_hi = ax.k_hi;
        sc.k_lo = ax.k_lo;
        const spine_landing land = choose_spine_shapes(sc, outer, _spines);
        if (land.fit != spine_kind::none) {
          for (int k = 0; k < static_cast<int>(outer_at.size()) &&
                          k < static_cast<int>(land.hit.size());
               ++k) {
            if (k < static_cast<int>(land.ok.size()) &&
                land.ok[static_cast<size_t>(k)] &&
                land.hit[static_cast<size_t>(k)].allFinite())
              mapped[static_cast<size_t>(outer_at[static_cast<size_t>(k)])] =
                  center + land.hit[static_cast<size_t>(k)];
          }
          spined = true;
        }
      }
      if (!spined && _monge) {
        for (int k = 0; k < static_cast<int>(outer.size()); ++k) {
          vec3 hit = outer[static_cast<size_t>(k)];
          if (select_monge_spoke(_objectives, outer[static_cast<size_t>(k)], goal,
                                 nullptr, T, B, N, hit) &&
              hit.allFinite())
            mapped[static_cast<size_t>(outer_at[static_cast<size_t>(k)])] =
                center + hit;
        }
      }
    }
    write_spoke_targets(_patch, mapped, _id0, _w, p);
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0])
      return;
    fill_spoke_rows(_patch, *_blocks[0], id0, _w, triplets);
  }

  face_patch _patch;
  std::vector<objective> _objectives;
  spine_set _spines;
  bool _monge = false;
  const std::vector<real> *_u = nullptr;
  std::shared_ptr<const std::vector<cr_edge>> _edges;
};

/// Scales the patch positions on the face axes, then one ARAP spoke per edge.
class face_stretch : public block_constraint {
public:
  DEFINE_CREATE_FUNC(face_stretch)

  face_stretch(face_patch patch, real angle, real sx, real sy, real power,
               real w, std::vector<sim_block::ptr> blocks,
               const std::vector<real> *coeff = nullptr,
               std::shared_ptr<const std::vector<cr_edge>> edges = nullptr)
      : block_constraint({}, w, std::move(blocks)), _patch(patch), _angle(angle),
        _sx(sx), _sy(sy), _power(power), _u(coeff), _edges(std::move(edges)) {}

  std::string name() override { return "face_stretch"; }

  void project(const vecX &q, vecX &p) override {
    if (_blocks.size() < 2 || !_blocks[0] || !_blocks[1])
      return;
    std::vector<vec3> pos;
    vec3 center, N, T, B;
    if (!patch_positions(_patch, *_blocks[0], q, pos, center, N, T, B))
      return;
    const vec3 coeff = _blocks[1]->get_vec3(_patch.face, q);
    const shape_w goal = shape_with_axes(shape_w_from_coeffs(coeff));
    std::vector<vec3> mapped = pos;
    if (goal.ok) {
      const monge_axes ax = monge_axes_of(goal);
      const vec3 minor = ax.hx * T + ax.hy * B;
      const vec3 major = ax.lx * T + ax.ly * B;
      vec3 u, v;
      stretch_frame(minor, major, _angle, u, v);
      const real c = std::cos(_angle);
      const real s = std::sin(_angle);
      const real ux = c * ax.hx + s * ax.lx;
      const real uy = c * ax.hy + s * ax.ly;
      const real vx = -s * ax.hx + c * ax.lx;
      const real vy = -s * ax.hy + c * ax.ly;
      const real sx = stretch_axis_scale(goal, ux, uy, _sx, _power);
      const real sy = stretch_axis_scale(goal, vx, vy, _sy, _power);
      if (u.allFinite() && v.allFinite()) {
        for (int i = 0; i < _patch.nvert; ++i) {
          vec3 d = pos[static_cast<size_t>(i)] - center;
          d = scale_spoke(scale_spoke(d, u, sx), v, sy);
          mapped[static_cast<size_t>(i)] = center + d;
        }
      }
    }
    write_spoke_targets(_patch, mapped, _id0, _w, p);
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0])
      return;
    fill_spoke_rows(_patch, *_blocks[0], id0, _w, triplets);
  }

  face_patch _patch;
  real _angle = 0;
  real _sx = 0;
  real _sy = 0;
  real _power = 0;
  const std::vector<real> *_u = nullptr;
  std::shared_ptr<const std::vector<cr_edge>> _edges;
};

/// Vertex polar frame from the carried edge tangents and the cached gradients.
/// With `stars` and the edge store, `h` is the edge store's local polar
/// (`kusama::weingarten::vertex_tensor`) instead of the face operators.
class vertex_frame_cache_compute : public cache_compute {
public:
  vertex_frame_cache_compute(asawa::shell::shell *shell, sim_block::ptr pos,
                             std::vector<real> *u, std::vector<vec3> *gpar,
                             std::vector<vec3> *gperp,
                             std::vector<vertex_chart> *charts,
                             const kusama::weingarten::vertex_stars *stars = nullptr,
                             sim_block::ptr store = nullptr)
      : _shell(shell), _pos(std::move(pos)), _u(u), _gpar(gpar), _gperp(gperp),
        _charts(charts), _stars(stars), _store(std::move(store)) {}

  void compute(vecX &q) override {
    if (!_shell || !_pos || !_u || !_charts)
      return;
    const int n = _shell->vert_count();
    if (n <= 0)
      return;
    std::vector<vec3> x(static_cast<size_t>(n), vec3::Zero());
    for (int i = 0; i < n; ++i)
      x[static_cast<size_t>(i)] = _pos->get_vec3(i, q);
    std::vector<vec3> H(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> wrote(static_cast<size_t>(n), 0);
    vertex_weingarten(*_shell, x, *_u, H, tang, axes, wrote, _charts, _gpar,
                      _gperp);
    if (_stars && _store && static_cast<int>(_stars->at.size()) == n) {
      auto pos = [&](int i) { return x[static_cast<size_t>(i)]; };
      auto store = [&](int e) { return _store->get_vec3(e, q); };
      for (int i = 0; i < n && i < static_cast<int>(_charts->size()); ++i) {
        vertex_chart &c = (*_charts)[static_cast<size_t>(i)];
        if (!c.ok)
          continue;
        const vec3 h = kusama::weingarten::vertex_tensor(*_stars, i, pos, store, c.T, c.B, c.N);
        if (h.allFinite())
          c.h = h;
        else
          c.ok = false;
      }
    }
  }

  asawa::shell::shell *_shell = nullptr;
  sim_block::ptr _pos;
  std::vector<real> *_u = nullptr;
  std::vector<vec3> *_gpar = nullptr;
  std::vector<vec3> *_gperp = nullptr;
  std::vector<vertex_chart> *_charts = nullptr;
  const kusama::weingarten::vertex_stars *_stars = nullptr;
  sim_block::ptr _store;
};

inline void init_shell_caches(asawa::shell::shell &shell, shell_cr_block &cr,
                              sim_block::ptr positions,
                              std::vector<cache_compute::ptr> &computes) {
  if (!positions || cr.topology.nE <= 0)
    return;
  const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);
  std::vector<tangent_cache_compute::edge> edges;
  for (int e = 0; e < cr.topology.nE; ++e) {
    const auto &faces = cr.topology.edge_incident[static_cast<size_t>(e)];
    if (faces.size() != 2)
      continue;
    const asawa::shell::FaceId f0 = faces[0];
    const asawa::shell::FaceId f1 = faces[1];
    if (shell.fsize(f0) != 3 || shell.fsize(f1) != 3)
      continue;
    const asawa::shell::CornerId c = cr.topology.dof_to_corner[static_cast<size_t>(e)];
    if (static_cast<int>(c) < 0)
      continue;
    const int v0 = static_cast<int>(shell.vert(c));
    const int v1 = static_cast<int>(shell.vert(shell.next(c)));
    const std::array<asawa::shell::VertId, 3> t0 = shell.get_tri(f0);
    const std::array<asawa::shell::VertId, 3> t1 = shell.get_tri(f1);
    int o0 = -1;
    int o1 = -1;
    for (int i = 0; i < 3; ++i) {
      const int a = static_cast<int>(t0[static_cast<size_t>(i)]);
      const int b = static_cast<int>(t1[static_cast<size_t>(i)]);
      if (a != v0 && a != v1)
        o0 = a;
      if (b != v0 && b != v1)
        o1 = b;
    }
    if (o0 < 0 || o1 < 0 || v0 < 0 || v1 < 0 ||
        v0 >= static_cast<int>(x.size()) || v1 >= static_cast<int>(x.size()) ||
        o0 >= static_cast<int>(x.size()) || o1 >= static_cast<int>(x.size()))
      continue;
    if (static_cast<int>(cr.u.size()) < 2 * cr.topology.nE)
      continue;
    const real gp = cr.u[static_cast<size_t>(e)];
    const real gq = cr.u[static_cast<size_t>(e + cr.topology.nE)];
    if (!(gp * gp + gq * gq > real(1e-16)))
      continue;
    vec3 e_par, e_perp;
    if (!cr_edge_frame(shell, x, c, faces, e_par, e_perp))
      continue;
    tangent_cache_compute::edge row;
    row.id = e;
    row.v0 = v0;
    row.v1 = v1;
    row.o0 = o0;
    row.o1 = o1;
    for (int i = 0; i < 3; ++i) {
      row.t0[static_cast<size_t>(i)] = static_cast<int>(t0[static_cast<size_t>(i)]);
      row.t1[static_cast<size_t>(i)] = static_cast<int>(t1[static_cast<size_t>(i)]);
    }
    const vec3 raw = x[static_cast<size_t>(v1)] - x[static_cast<size_t>(v0)];
    const real plen = raw.norm();
    if (!(plen > real(1e-20)))
      continue;
    row.par_sign = e_par.dot(raw) >= real(0) ? real(1) : real(-1);
    vec3 N;
    if (!diamond_normal(x[static_cast<size_t>(v0)], x[static_cast<size_t>(v1)],
                        x[static_cast<size_t>(o0)], x[static_cast<size_t>(o1)],
                        row.t0, row.t1, v0, v1, o0, o1, N))
      continue;
    vec3 trial = N.cross(e_par);
    if (!(trial.squaredNorm() > real(1e-24)) || !trial.allFinite())
      continue;
    trial.normalize();
    row.perp_sign = e_perp.dot(trial) >= real(0) ? real(1) : real(-1);
    diamond_mid_spokes(x[static_cast<size_t>(v0)], x[static_cast<size_t>(v1)],
                       x[static_cast<size_t>(o0)], x[static_cast<size_t>(o1)],
                       row.rest);
    if (!diamond_rotation(spoke_covariance(row.rest, row.rest), row.R0))
      row.R0 = mat3::Identity();
    row.T0 = gp * e_par + gq * e_perp;
    edges.push_back(row);
  }
  computes.push_back(std::make_shared<tangent_cache_compute>(
      positions, &cr.u, cr.topology.nE, std::move(edges)));

  cr_complex cx;
  if (!build_cr_complex(shell, cx))
    return;
  computes.push_back(std::make_shared<gradient_cache_compute>(
      positions, &cr.g_par, &cr.g_perp, cx.tris, cx.edges, shell.vert_count()));
  computes.push_back(std::make_shared<vertex_frame_cache_compute>(
      &shell, positions, &cr.u, &cr.g_par, &cr.g_perp, &cr.vertex_frames, &cr.stars,
      cr.shape_dof()));
}

inline void init_face_shape_pin(asawa::shell::shell &shell, shell_cr_block &cr,
                                sim_block::ptr positions,
                                std::vector<projection_constraint::ptr> &constraints,
                                real w) {
  if (!(w > real(0)) || !positions || !cr.face_dof())
    return;
  cr_complex cx;
  if (!build_cr_complex(shell, cx))
    return;
  std::vector<face_pin_slot> faces;
  for (asawa::shell::FaceId f : shell.get_face_range(true)) {
    if (shell.fsize(f) != 3)
      continue;
    const std::array<asawa::shell::VertId, 3> tri = shell.get_tri(f);
    face_pin_slot slot;
    slot.face = static_cast<int>(f);
    if (slot.face < 0 || slot.face >= static_cast<int>(cr.face_shape.size()) ||
        slot.face >= static_cast<int>(cx.face_local.size()))
      continue;
    slot.tri = cx.face_local[static_cast<size_t>(slot.face)];
    bool ok = true;
    for (int i = 0; i < 3; ++i) {
      slot.v[static_cast<size_t>(i)] = static_cast<int>(tri[static_cast<size_t>(i)]);
      slot.e[static_cast<size_t>(i)] = patch_edge_dof(
          shell, cr.topology, slot.v[static_cast<size_t>(i)],
          static_cast<int>(tri[static_cast<size_t>((i + 1) % 3)]));
      if (slot.e[static_cast<size_t>(i)] < 0)
        ok = false;
    }
    if (ok)
      faces.push_back(slot);
  }
  if (faces.empty())
    return;
  constraints.push_back(face_shape_pin::create(
      std::move(faces), cx.tris, cx.edges, &cr.g_par, &cr.g_perp,
      shell.vert_count(), w,
      std::vector<sim_block::ptr>{positions, cr.face_dof()}, &cr.u));
}

inline std::shared_ptr<const std::vector<cr_edge>>
face_chart_edges(asawa::shell::shell &shell) {
  cr_complex cx;
  if (!build_cr_complex(shell, cx))
    return nullptr;
  return std::make_shared<const std::vector<cr_edge>>(std::move(cx.edges));
}

inline void stamp_patch_edges(const asawa::shell::shell &shell,
                              const cr_topology &topo, face_patch &patch) {
  if (patch.nvert < 3)
    return;
  for (int i = 0; i < 3; ++i) {
    patch.edge[static_cast<size_t>(i)] = patch_edge_dof(
        shell, topo, patch.vid[static_cast<size_t>(i)],
        patch.vid[static_cast<size_t>((i + 1) % 3)]);
  }
}

inline void init_face_stretch(const asawa::shell::shell &shell,
                              const shell_cr_block &cr, sim_block::ptr positions,
                              std::vector<projection_constraint::ptr> &constraints,
                              real w, const std::vector<stretch_target> &targets) {
  if (!(w > real(0)) || !positions || !cr.face_dof() || targets.empty())
    return;
  const auto edges =
      face_chart_edges(const_cast<asawa::shell::shell &>(shell));
  for (asawa::shell::FaceId f : shell.get_face_range(true)) {
    face_patch patch;
    if (!bake_face_patch(shell, f, patch))
      continue;
    if (patch.face < 0 || patch.face >= static_cast<int>(cr.face_shape.size()))
      continue;
    stamp_patch_edges(shell, cr.topology, patch);
    for (const stretch_target &spec : targets) {
      if (!std::isfinite(spec.x) || !std::isfinite(spec.y))
        continue;
      if (!(std::abs(spec.x) > real(1e-16) || std::abs(spec.y) > real(1e-16)))
        continue;
      constraints.push_back(face_stretch::create(
          patch, spec.angle, spec.x, spec.y, spec.power, w,
          std::vector<sim_block::ptr>{positions, cr.face_dof()}, &cr.u,
          edges));
    }
  }
}

inline void init_face_fair(const asawa::shell::shell &shell, const shell_cr_block &cr,
                           sim_block::ptr positions,
                           std::vector<projection_constraint::ptr> &constraints,
                           real w, const std::vector<objective> &objectives) {
  if (!(w > real(0)) || !positions || !cr.face_dof() || objectives.empty())
    return;
  bool want = false;
  for (const objective &o : objectives) {
    if (std::holds_alternative<monge_fairing>(o) ||
        std::holds_alternative<spine_fairing>(o))
      want = true;
  }
  if (!want)
    return;
  const auto edges =
      face_chart_edges(const_cast<asawa::shell::shell &>(shell));
  for (asawa::shell::FaceId f : shell.get_face_range(true)) {
    face_patch patch;
    if (!bake_face_patch(shell, f, patch))
      continue;
    stamp_patch_edges(shell, cr.topology, patch);
    constraints.push_back(face_fair::create(
        patch, w, objectives,
        std::vector<sim_block::ptr>{positions, cr.face_dof()}, &cr.u, edges));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
