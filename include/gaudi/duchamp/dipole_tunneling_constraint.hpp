#pragma once

#include <array>
#include <cassert>
#include <functional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"
#include "gaudi/duchamp/dipole_tunneling_geometry.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace duchamp {

using dipole_Nr_fn = std::function<std::vector<vec3>()>;

namespace dipole_tunneling {

// Dipole weld: mesh point ↔ rod segment onto the dipole cylinder along mesh N.
// ids = {vm, vr0, vr1}; blocks = {shell, rod}.
// Target cylinder is rod-relative (cen = xr + s·r·Nr, axis T); ray dir is Nm.
class dipole_weld : public hepworth::block::block_constraint {
public:
  DEFINE_CREATE_FUNC(dipole_weld)

  dipole_weld(const std::vector<index_t> &ids, const vec3 &Nr0, const vec3 &Nr1,
              const vec3 &Nm, real r, real wm, real wr,
              std::vector<hepworth::sim_block::ptr> blocks)
      : hepworth::block::block_constraint(ids, 0.0, blocks),
        _Nr0(Nr0),
        _Nr1(Nr1),
        _Nm(Nm),
        _r(r),
        _wm(wm),
        _wr(wr) {
    assert(ids.size() == 3);
    assert(blocks.size() == 2);
  }

  virtual std::string name() { return "dipole_weld"; }

  static real segment_s(const vec3 &pm, const vec3 &r0, const vec3 &r1) {
    const vec3 dr = r1 - r0;
    const real dr2 = dr.squaredNorm();
    real s = dr2 > 1e-24 ? (pm - r0).dot(dr) / dr2 : 0.5;
    return va::clamp(s, 0.0, 1.0);
  }

  virtual void project(const vecX &q, vecX &p) {
    const vec3 pm = _blocks[0]->get_vec3(_ids[0], q);
    const vec3 r0 = _blocks[1]->get_vec3(_ids[1], q);
    const vec3 r1 = _blocks[1]->get_vec3(_ids[2], q);

    const real s = segment_s(pm, r0, r1);
    const real a0 = 1.0 - s;
    const real a1 = s;
    const vec3 pr = va::mix(s, r0, r1);

    vec3 dX = vec3::Zero();
    segment_frame f;
    const vec3 Ni = interpolate_rod_normal(_Nr0, _Nr1, s);
    if (try_make_segment_frame(r0, r1, Ni, f) && _r > 0.0 &&
        is_finite_vec(pm) && is_finite_vec(pr)) {
      const vec3 p_proj =
          proj_to_dipole_cylinder_along_normal(pm, pr, f, _r, _Nm);
      if (is_finite_vec(p_proj)) {
        // Displacement = (−side)·(r−dist)·N  (upper −N, lower +N).
        dX = p_proj - pm;
        // Debug: mesh point → SDF step target.
        //if (dX.squaredNorm() > 1e-18)
        //  geometry_logger::line(pm, p_proj, vec4(0.0, 0.4, 0.7, 1.0));
      }
      if (!is_finite_vec(dX))
        dX = vec3::Zero();
    }

    p.block(_id0 + 0, 0, 3, 1) = _wr * (r0 - a0 * dX) - _wm * pm;
    p.block(_id0 + 3, 0, 3, 1) = _wr * (r1 - a1 * dX) - _wm * pm;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t im = _blocks[0]->get_offset_idx(_ids[0]);
    const index_t ir0 = _blocks[1]->get_offset_idx(_ids[1]);
    const index_t ir1 = _blocks[1]->get_offset_idx(_ids[2]);

    for (int ax = 0; ax < 3; ++ax) {
      triplets.push_back(trip(_id0 + 0 + ax, ir0 + ax, _wr));
      triplets.push_back(trip(_id0 + 0 + ax, im + ax, -_wm));
      triplets.push_back(trip(_id0 + 3 + ax, ir1 + ax, _wr));
      triplets.push_back(trip(_id0 + 3 + ax, im + ax, -_wm));
    }
    id0 += 6;
  }

  vec3 _Nr0{0.0, 0.0, 0.0};
  vec3 _Nr1{0.0, 0.0, 0.0};
  vec3 _Nm{0.0, 0.0, 1.0};
  real _r = 0.0;
  real _wm = 1.0;
  real _wr = 1.0;
};

// Face-normal alignment onto the dipole sphere.
// ids = {vs0, vs1, vs2}; blocks = {shell, rod}; rod endpoints _rid0/_rid1.
class triangle_dipole_tunneling : public hepworth::block::block_constraint {
public:
  DEFINE_CREATE_FUNC(triangle_dipole_tunneling)

  triangle_dipole_tunneling(const std::vector<index_t> &ids, const std::vector<vec3> &x,
                            index_t rid0, index_t rid1, const vec3 &Nr0, const vec3 &Nr1,
                            real r, real w, std::vector<hepworth::sim_block::ptr> blocks)
      : hepworth::block::block_constraint(ids, w, blocks),
        _Nr0(Nr0),
        _Nr1(Nr1),
        _r(r),
        _rid0(rid0),
        _rid1(rid1) {
    assert(ids.size() == 3);
    assert(blocks.size() == 2);
    const vec3 &p0 = x[ids[0]];
    const vec3 &p1 = x[ids[1]];
    const vec3 &p2 = x[ids[2]];
    _d0 = p1 - p0;
    _d1 = p2 - p0;
    _N0 = safe_normalized(_d0.cross(_d1), vec3::UnitZ());
  }

  virtual std::string name() { return "triangle_dipole_tunneling"; }

  vec3 calc_proj_N(const vecX &q) const {
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    const vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    const vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);
    const vec3 xf = (q0 + q1 + q2) / 3.0;

    const vec3 xr0 = _blocks[1]->get_vec3(_rid0, q);
    const vec3 xr1 = _blocks[1]->get_vec3(_rid1, q);
    const real s = dipole_weld::segment_s(xf, xr0, xr1);
    const vec3 Ni = interpolate_rod_normal(_Nr0, _Nr1, s);
    segment_frame f;
    if (!try_make_segment_frame(xr0, xr1, Ni, f) || !(_r > 0.0) ||
        !is_finite_vec(xf))
      return _N0;
    const auto &[xr0f, xr1f, N, T, B] = f;
    (void)xr0f;
    (void)xr1f;
    (void)T;
    (void)B;

    const vec3 xr = safe_project_on_line(xr0, xr1, xf);
    const vec3 dp = xf - xr;
    const real sg = va::sgn(N.dot(dp));
    const real side = (sg == 0.0) ? 1.0 : sg;
    const vec3 cen = dipole_sphere_center(xr, dp, N, _r);
    if (!is_finite_vec(cen))
      return _N0;
    const vec3 dpcen = xf - cen;
    const vec3 n_proj = dpcen.squaredNorm() > k_dipole_eps2
                            ? safe_normalized(dpcen, side * N)
                            : side * N;
    // Flip vs raw sphere outward: +N side (upper) tunnels / opens clearance;
    // −N side (lower) expands onto the tube.
    const vec3 out = -side * n_proj;
    return is_finite_vec(out) && out.squaredNorm() > k_dipole_eps2 ? out : _N0;
  }

  virtual void project(const vecX &q, vecX &p) {
    const vec3 Np = calc_proj_N(q);
    if (!(_N0.squaredNorm() > k_dipole_eps2) ||
        !(Np.squaredNorm() > k_dipole_eps2)) {
      p.block(_id0 + 0, 0, 3, 1) = _w * _d0;
      p.block(_id0 + 3, 0, 3, 1) = _w * _d1;
      return;
    }
    const quat qN = quat::FromTwoVectors(_N0, Np);
    p.block(_id0 + 0, 0, 3, 1) = _w * (qN * _d0);
    p.block(_id0 + 3, 0, 3, 1) = _w * (qN * _d1);
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t i0 = _blocks[0]->get_offset_idx(_ids[0]);
    const index_t i1 = _blocks[0]->get_offset_idx(_ids[1]);
    const index_t i2 = _blocks[0]->get_offset_idx(_ids[2]);

    for (int ax = 0; ax < 3; ++ax)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, -_w));
    for (int ax = 0; ax < 3; ++ax)
      triplets.push_back(trip(_id0 + 0 + ax, i1 + ax, _w));

    for (int ax = 0; ax < 3; ++ax)
      triplets.push_back(trip(_id0 + 3 + ax, i0 + ax, -_w));
    for (int ax = 0; ax < 3; ++ax)
      triplets.push_back(trip(_id0 + 3 + ax, i2 + ax, _w));

    id0 += 6;
  }

  vec3 _Nr0{0.0, 0.0, 0.0};
  vec3 _Nr1{0.0, 0.0, 0.0};
  vec3 _d0{0.0, 0.0, 0.0};
  vec3 _d1{0.0, 0.0, 0.0};
  vec3 _N0{0.0, 0.0, 0.0};
  real _r = 0.0;
  index_t _rid0 = 0;
  index_t _rid1 = 0;
};

struct dipole_edge_key {
  index_t a = 0;
  index_t b = 0;
  bool operator==(const dipole_edge_key &o) const { return a == o.a && b == o.b; }
};

struct dipole_edge_key_hash {
  size_t operator()(const dipole_edge_key &k) const {
    return (size_t(k.a) << 32) ^ size_t(k.b);
  }
};

inline dipole_edge_key make_edge_key(index_t a, index_t b) {
  return a < b ? dipole_edge_key{a, b} : dipole_edge_key{b, a};
}

using dipole_face_capture_map =
    std::unordered_map<dipole_edge_key,
                       std::pair<std::vector<asawa::shell::FaceId>,
                                 std::vector<asawa::shell::FaceId>>,
                       dipole_edge_key_hash>;

// Rod AABB prefilter (get_collisions), then dipole classify.
inline dipole_face_capture_map
gather_dipole_faces(asawa::shell::shell &M, const std::vector<vec3> &x,
                    const std::vector<vec3> &x_r, const std::vector<vec3> &Nr,
                    asawa::rod::dynamic &rod_d, real r) {
  dipole_face_capture_map out;

  const real tol = 2.0 * (r + dipole_margin(r));
  std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
  const auto collisions = rod_d.get_collisions(edge_verts_M, x, tol);

  std::unordered_map<dipole_edge_key, std::set<int>, dipole_edge_key_hash>
      candidates;
  for (const auto &c : collisions) {
    if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[3] < 0)
      continue;
    const index_t vr0 = c[2];
    const index_t vr1 = c[3];
    if (vr0 < 0 || vr1 < 0 || static_cast<size_t>(vr0) >= x_r.size() ||
        static_cast<size_t>(vr1) >= x_r.size())
      continue;

    const asawa::shell::CornerId ec = M.find_edge_from_verts(
        asawa::shell::vert_id(c[0]), asawa::shell::vert_id(c[1]));
    if (static_cast<int>(ec) < 0)
      continue;

    auto &face_set = candidates[make_edge_key(vr0, vr1)];
    for (const asawa::shell::FaceId fid : M.dihedral_face_ids(ec))
      face_set.insert(static_cast<int>(fid));
  }

  for (auto &[key, face_set] : candidates) {
    const index_t vr0 = key.a;
    const index_t vr1 = key.b;
    if (static_cast<size_t>(vr0) >= Nr.size() ||
        static_cast<size_t>(vr1) >= Nr.size())
      continue;
    segment_frame f;
    if (!try_make_segment_frame(
            x_r[vr0], x_r[vr1],
            interpolate_rod_normal(Nr[vr0], Nr[vr1], 0.5), f))
      continue;
    std::vector<asawa::shell::FaceId> upper;
    std::vector<asawa::shell::FaceId> lower;
    std::set<int> upper_seen;
    std::set<int> lower_seen;
    for (const int fid_i : face_set) {
      const asawa::shell::FaceId face = asawa::shell::face_id(fid_i);
      const int zone = classify_face_on_segment(M, face, x, f, r);
      if (zone > 0 && upper_seen.insert(fid_i).second)
        upper.push_back(face);
      else if (zone < 0 && lower_seen.insert(fid_i).second)
        lower.push_back(face);
    }
    if (!upper.empty() || !lower.empty())
      out.emplace(key, std::make_pair(std::move(upper), std::move(lower)));
  }
  return out;
}

// Same classify as gather_dipole_faces, but sweeps every mesh face per true rod
// edge (edge_verts from rod::get_edge_vert_ids). Audit / no-dynamic fallback.
inline dipole_face_capture_map
gather_dipole_faces_brute(const asawa::shell::shell &M, const std::vector<vec3> &x,
                          const std::vector<vec3> &x_r, const std::vector<vec3> &Nr,
                          const std::vector<index_t> &edge_verts, real r) {
  dipole_face_capture_map out;
  assert(edge_verts.size() % 2 == 0);
  for (size_t e = 0; e + 1 < edge_verts.size(); e += 2) {
    const index_t v0 = edge_verts[e];
    const index_t v1 = edge_verts[e + 1];
    if (v0 < 0 || v1 < 0 || static_cast<size_t>(v0) >= x_r.size() ||
        static_cast<size_t>(v1) >= x_r.size())
      continue;
    segment_frame f;
    if (!try_make_segment_frame(
            x_r[v0], x_r[v1],
            interpolate_rod_normal(Nr[v0], Nr[v1], 0.5), f))
      continue;
    auto [upper, lower] = gather_faces_in_segment(M, x, f, r);
    if (!upper.empty() || !lower.empty())
      out.emplace(make_edge_key(v0, v1),
                  std::make_pair(std::move(upper), std::move(lower)));
  }
  return out;
}

// Returns {n_a_only_faces, n_b_only_faces, n_shared_faces}.
inline std::array<index_t, 3>
compare_dipole_gather_face_sets(const dipole_face_capture_map &a,
                                const dipole_face_capture_map &b) {
  auto flatten = [](const dipole_face_capture_map &m) {
    std::set<std::tuple<index_t, index_t, int>> s;
    for (const auto &[key, zones] : m) {
      for (const auto fid : zones.first)
        s.insert({key.a, key.b, static_cast<int>(fid)});
      for (const auto fid : zones.second)
        s.insert({key.a, key.b, static_cast<int>(fid)});
    }
    return s;
  };
  const auto A = flatten(a);
  const auto B = flatten(b);
  index_t shared = 0;
  for (const auto &t : A)
    if (B.count(t))
      ++shared;
  return {static_cast<index_t>(A.size()) - shared,
          static_cast<index_t>(B.size()) - shared, shared};
}

// Shared capture → face-orientation + vertex weld constraints.
// Gather once; emit triangle_dipole_tunneling and/or dipole_weld from the
// same face set. Weld targets the rod-relative dipole cylinder along mesh N.
inline void init_dipole_clearance(
    asawa::shell::shell &shell, asawa::rod::dynamic &rod_d,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, real r, real w_tunnel, real wm, real wr,
    std::vector<hepworth::sim_block::ptr> blocks) {
  const bool do_tunnel = w_tunnel > 0.0;
  const bool do_weld = wm > 0.0 || wr > 0.0;
  if (!do_tunnel && !do_weld)
    return;
  assert(Nr.size() >= x_r.size());
  const auto captured = gather_dipole_faces(shell, x, x_r, Nr, rod_d, r);
  std::set<std::tuple<index_t, index_t, index_t>> weld_seen;
  for (const auto &[key, zones] : captured) {
    const index_t vr0 = key.a;
    const index_t vr1 = key.b;
    auto add_face = [&](asawa::shell::FaceId fid) {
      if (shell.fsize(fid) != 3)
        return;
      if (asawa::shell::face_area(shell, fid, x) < 1e-8)
        return;
      const auto tri = shell.get_tri(fid);
      if (do_tunnel) {
        constraints.push_back(triangle_dipole_tunneling::create(
            std::vector<index_t>({tri[0], tri[1], tri[2]}), x, vr0, vr1,
            Nr[vr0], Nr[vr1], r, w_tunnel, blocks));
      }
      if (do_weld) {
        for (index_t vi : tri) {
          if (!weld_seen.insert({vi, vr0, vr1}).second)
            continue;
          constraints.push_back(dipole_weld::create(
              std::vector<index_t>({vi, vr0, vr1}), Nr[vr0], Nr[vr1],
              asawa::shell::vert_normal(shell, asawa::shell::vert_id(vi), x), r,
              wm, wr, blocks));
        }
      }
    };
    for (const auto fid : zones.first)
      add_face(fid);
    for (const auto fid : zones.second)
      add_face(fid);
  }
}

// Thin wrappers: same capture path, one constraint family only.
inline void init_dipole_tunneling(
    asawa::shell::shell &shell, asawa::rod::dynamic &rod_d,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, real r, real w,
    std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_clearance(shell, rod_d, constraints, x, x_r, Nr, r, w, 0.0, 0.0,
                        blocks);
}

inline void init_dipole_weld(
    asawa::shell::shell &shell, asawa::rod::dynamic &rod_d,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, real r, real wm, real wr,
    std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_clearance(shell, rod_d, constraints, x, x_r, Nr, r, 0.0, wm, wr,
                        blocks);
}

// Brute-force fallback: edge_verts must be rod::get_edge_vert_ids() pairs.
inline void init_dipole_clearance_brute(
    const asawa::shell::shell &shell,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, const std::vector<index_t> &edge_verts, real r,
    real w_tunnel, real wm, real wr,
    std::vector<hepworth::sim_block::ptr> blocks) {
  const bool do_tunnel = w_tunnel > 0.0;
  const bool do_weld = wm > 0.0 || wr > 0.0;
  if (!do_tunnel && !do_weld)
    return;
  assert(Nr.size() >= x_r.size());
  assert(edge_verts.size() % 2 == 0);
  std::set<std::tuple<index_t, index_t, index_t>> weld_seen;

  for (size_t e = 0; e + 1 < edge_verts.size(); e += 2) {
    const index_t v0 = edge_verts[e];
    const index_t v1 = edge_verts[e + 1];
    segment_frame f;
    if (!try_make_segment_frame(x_r[v0], x_r[v1],
                                interpolate_rod_normal(Nr[v0], Nr[v1], 0.5), f))
      continue;
    const auto [upper, lower] = gather_faces_in_segment(shell, x, f, r);

    auto add_faces = [&](const std::vector<asawa::shell::FaceId> &faces) {
      for (const asawa::shell::FaceId fid : faces) {
        if (shell.fsize(fid) != 3)
          continue;
        if (asawa::shell::face_area(shell, fid, x) < 1e-8)
          continue;
        const auto tri = shell.get_tri(fid);
        if (do_tunnel) {
          constraints.push_back(triangle_dipole_tunneling::create(
              std::vector<index_t>({tri[0], tri[1], tri[2]}), x, v0, v1, Nr[v0],
              Nr[v1], r, w_tunnel, blocks));
        }
        if (do_weld) {
          for (index_t vi : tri) {
            if (!weld_seen.insert({vi, v0, v1}).second)
              continue;
            constraints.push_back(dipole_weld::create(
                std::vector<index_t>({vi, v0, v1}), Nr[v0], Nr[v1],
                asawa::shell::vert_normal(shell, asawa::shell::vert_id(vi), x),
                r, wm, wr, blocks));
          }
        }
      }
    };

    add_faces(upper);
    add_faces(lower);
  }
}

inline void init_dipole_weld_brute(
    const asawa::shell::shell &shell,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, const std::vector<index_t> &edge_verts, real r,
    real wm, real wr, std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_clearance_brute(shell, constraints, x, x_r, Nr, edge_verts, r, 0.0,
                              wm, wr, blocks);
}

inline void init_dipole_weld_brute(
    const asawa::shell::shell &shell,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r, const vec3 &dipole_N,
    const std::vector<index_t> &edge_verts, real r, real wm, real wr,
    std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_weld_brute(shell, constraints, x, x_r,
                         std::vector<vec3>(x_r.size(), dipole_N), edge_verts, r,
                         wm, wr, blocks);
}

inline void init_dipole_tunneling_brute(
    const asawa::shell::shell &shell,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r,
    const std::vector<vec3> &Nr, const std::vector<index_t> &edge_verts, real r,
    real w, std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_clearance_brute(shell, constraints, x, x_r, Nr, edge_verts, r, w,
                              0.0, 0.0, blocks);
}

inline void init_dipole_tunneling_brute(
    const asawa::shell::shell &shell,
    std::vector<hepworth::projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &x_r, const vec3 &dipole_N,
    const std::vector<index_t> &edge_verts, real r, real w,
    std::vector<hepworth::sim_block::ptr> blocks) {
  init_dipole_tunneling_brute(shell, constraints, x, x_r,
                              std::vector<vec3>(x_r.size(), dipole_N), edge_verts,
                              r, w, blocks);
}

} // namespace dipole_tunneling

template <size_t... Is>
inline hepworth::block::constraint_recompute_fn make_dipole_weld_recompute(
    hepworth::block::shell_position_block::ptr shell,
    hepworth::block::rod_position_block::ptr rod, dipole_Nr_fn Nr_fn, real disk_r,
    real wm, real wr) {
  static_assert(sizeof...(Is) == 2, "dipole weld routes shell then rod");
  return [shell, rod, Nr_fn, disk_r, wm, wr](hepworth::block::solver_context &ctx) {
    if (rod->dynamic) {
      dipole_tunneling::init_dipole_weld(
          *shell->mesh, *rod->dynamic, ctx.constraints, shell->xs->get(),
          rod->rod->x(), Nr_fn(), disk_r, wm, wr,
          hepworth::block::select_blocks<Is...>(ctx));
    } else {
      dipole_tunneling::init_dipole_weld_brute(
          *shell->mesh, ctx.constraints, shell->xs->get(), rod->rod->x(), Nr_fn(),
          rod->rod->get_edge_vert_ids(), disk_r, wm, wr,
          hepworth::block::select_blocks<Is...>(ctx));
    }
  };
}

template <size_t... Is>
inline hepworth::block::constraint_bundle
make_dipole_weld_bundle(hepworth::block::shell_position_block::ptr shell,
                        hepworth::block::rod_position_block::ptr rod,
                        dipole_Nr_fn Nr_fn, real disk_r, real wm, real wr) {
  return {make_dipole_weld_recompute<Is...>(shell, rod, std::move(Nr_fn), disk_r, wm,
                                            wr)};
}

template <size_t... Is>
inline hepworth::block::constraint_bundle
make_dipole_weld_bundle(hepworth::block::shell_position_block::ptr shell,
                        hepworth::block::rod_position_block::ptr rod,
                        const vec3 &dipole_N, real disk_r, real wm, real wr) {
  return make_dipole_weld_bundle<Is...>(
      shell, rod,
      [rod, dipole_N]() {
        return std::vector<vec3>(rod->rod->x().size(), dipole_N);
      },
      disk_r, wm, wr);
}

template <size_t... Is>
inline hepworth::block::constraint_recompute_fn make_dipole_tunneling_recompute(
    hepworth::block::shell_position_block::ptr shell,
    hepworth::block::rod_position_block::ptr rod, dipole_Nr_fn Nr_fn, real disk_r,
    real w) {
  static_assert(sizeof...(Is) == 2, "dipole tunneling routes shell then rod");
  return [shell, rod, Nr_fn, disk_r, w](hepworth::block::solver_context &ctx) {
    if (rod->dynamic) {
      dipole_tunneling::init_dipole_tunneling(
          *shell->mesh, *rod->dynamic, ctx.constraints, shell->xs->get(),
          rod->rod->x(), Nr_fn(), disk_r, w,
          hepworth::block::select_blocks<Is...>(ctx));
    } else {
      dipole_tunneling::init_dipole_tunneling_brute(
          *shell->mesh, ctx.constraints, shell->xs->get(), rod->rod->x(), Nr_fn(),
          rod->rod->get_edge_vert_ids(), disk_r, w,
          hepworth::block::select_blocks<Is...>(ctx));
    }
  };
}

template <size_t... Is>
inline hepworth::block::constraint_bundle
make_dipole_tunneling_bundle(hepworth::block::shell_position_block::ptr shell,
                             hepworth::block::rod_position_block::ptr rod,
                             dipole_Nr_fn Nr_fn, real disk_r, real w) {
  return {make_dipole_tunneling_recompute<Is...>(shell, rod, std::move(Nr_fn), disk_r,
                                                 w)};
}

template <size_t... Is>
inline hepworth::block::constraint_bundle
make_dipole_tunneling_bundle(hepworth::block::shell_position_block::ptr shell,
                             hepworth::block::rod_position_block::ptr rod,
                             const vec3 &dipole_N, real disk_r, real w) {
  return make_dipole_tunneling_bundle<Is...>(
      shell, rod,
      [rod, dipole_N]() {
        return std::vector<vec3>(rod->rod->x().size(), dipole_N);
      },
      disk_r, w);
}

// Combined tunnel+weld: one gather, both constraint families.
template <size_t... Is>
inline hepworth::block::constraint_recompute_fn make_dipole_clearance_recompute(
    hepworth::block::shell_position_block::ptr shell,
    hepworth::block::rod_position_block::ptr rod, dipole_Nr_fn Nr_fn,
    std::function<real()> disk_r_fn, std::function<real()> w_tunnel_fn,
    std::function<real()> wm_fn, std::function<real()> wr_fn) {
  static_assert(sizeof...(Is) == 2, "dipole clearance routes shell then rod");
  return [shell, rod, Nr_fn, disk_r_fn, w_tunnel_fn, wm_fn,
          wr_fn](hepworth::block::solver_context &ctx) {
    const real w_tunnel = w_tunnel_fn();
    const real wm = wm_fn();
    const real wr = wr_fn();
    if (w_tunnel <= 0.0 && wm <= 0.0 && wr <= 0.0)
      return;
    const real disk_r = disk_r_fn();
    const std::vector<vec3> Nr = Nr_fn();
    auto blocks = hepworth::block::select_blocks<Is...>(ctx);
    if (rod->dynamic) {
      dipole_tunneling::init_dipole_clearance(
          *shell->mesh, *rod->dynamic, ctx.constraints, shell->xs->get(),
          rod->rod->x(), Nr, disk_r, w_tunnel, wm, wr, blocks);
    } else {
      dipole_tunneling::init_dipole_clearance_brute(
          *shell->mesh, ctx.constraints, shell->xs->get(), rod->rod->x(), Nr,
          rod->rod->get_edge_vert_ids(), disk_r, w_tunnel, wm, wr, blocks);
    }
  };
}

// Convenience: constant weights (captured once at build time).
template <size_t... Is>
inline hepworth::block::constraint_recompute_fn make_dipole_clearance_recompute(
    hepworth::block::shell_position_block::ptr shell,
    hepworth::block::rod_position_block::ptr rod, dipole_Nr_fn Nr_fn,
    std::function<real()> disk_r_fn, real w_tunnel, real wm, real wr) {
  return make_dipole_clearance_recompute<Is...>(
      shell, rod, std::move(Nr_fn), std::move(disk_r_fn),
      [w_tunnel]() { return w_tunnel; }, [wm]() { return wm; },
      [wr]() { return wr; });
}

} // namespace duchamp
} // namespace gaudi
