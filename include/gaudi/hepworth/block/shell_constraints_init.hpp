#ifndef __HEP_SHELL_BLOCK_CONSTRAINTS_INIT__
#define __HEP_SHELL_BLOCK_CONSTRAINTS_INIT__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/hepworth/block/block_geometry.hpp"
#include "shell_collision_constraint.hpp"
#include "shell_constraints.hpp"
#include "gaudi/hepworth/block/cr_dirichlet.hpp"
#include "gaudi/geometry_logger.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

void init_edge_strain(const asawa::shell::shell &shell,
                      std::vector<projection_constraint::ptr> &constraints,
                      const std::vector<vec3> &x, const real &w,
                      std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();

  for (auto c0 : range) {
    asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
    int i = shell.vert(cid);
    int j = shell.vert(shell.other(cid));
    real l0 = asawa::shell::edge_length(shell, cid, x);

    if (std::isnan(l0))
      continue;
    if (l0 < 1e-8)
      continue;

    real ct = asawa::shell::cotan(shell, cid, x);
    real Ai = asawa::shell::face_area(shell, shell.face(cid), x);
    real Aj =
        asawa::shell::face_area(shell, shell.face(shell.other(cid)), x);
    if (Ai < 1e-8)
      continue;
    if (Aj < 1e-8)
      continue;

    constraints.push_back(
        edge_strain::create(std::vector<index_t>({i, j}), w, blocks));
  }
}

void init_edge_strain(const asawa::shell::shell &shell,
                      std::vector<projection_constraint::ptr> &constraints,
                      const std::vector<vec3> &x, const std::vector<real> &l0,
                      const real &w, std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();

  for (auto c0 : range) {
    asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
    int i = shell.vert(cid);
    int j = shell.vert(shell.other(cid));

    if (std::isnan(l0[c0 / 2]))
      continue;
    if (l0[c0 / 2] < 1e-8)
      continue;

    real ct = asawa::shell::cotan(shell, cid, x);
    real Ai = asawa::shell::face_area(shell, shell.face(cid), x);
    real Aj =
        asawa::shell::face_area(shell, shell.face(shell.other(cid)), x);
    if (Ai < 1e-8)
      continue;
    if (Aj < 1e-8)
      continue;

    edge_strain::ptr edge_strain =
        edge_strain::create(std::vector<index_t>({i, j}), w, blocks);
    edge_strain->set_goal_length(l0[c0 / 2]);
    constraints.push_back(edge_strain);
  }
}

void init_triangle_strain(const asawa::shell::shell &shell,
                          std::vector<projection_constraint::ptr> &constraints,
                          const std::vector<vec3> &x, real w,
                          std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_face_range();
  int i = 0;
  for (auto fi : range) {
    asawa::shell::FaceId fid = asawa::shell::face_id(fi);
    if (shell.fsize(fid) != 3)
      continue;
    if (asawa::shell::face_area(shell, fid, x) < 1e-8)
      continue;
    auto tri = shell.get_tri(fid);
    constraints.push_back( //
        triangle_strain::create(std::vector<index_t>({
                                    tri[0],
                                    tri[1],
                                    tri[2],
                                }),
                                0.5, 1.5, w, blocks));
  }
}

void init_triangle_conformal(const asawa::shell::shell &shell,
                             std::vector<projection_constraint::ptr> &constraints,
                             const std::vector<vec3> &x, real w,
                             std::vector<sim_block::ptr> blocks) {
  if (!(w > real(0)))
    return;
  for (auto fi : shell.get_face_range()) {
    asawa::shell::FaceId fid = asawa::shell::face_id(fi);
    if (shell.fsize(fid) != 3)
      continue;
    if (asawa::shell::face_area(shell, fid, x) < 1e-8)
      continue;
    auto tri = shell.get_tri(fid);
    constraints.push_back(triangle_conformal::create(
        std::vector<index_t>({tri[0], tri[1], tri[2]}), x, w, blocks));
  }
}

void init_anisotropic_triangle_strain(
    const asawa::shell::shell &shell,
    std::vector<projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const std::vector<vec3> &axial_face,
    const vec2 &growth, real w, std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_face_range();
  for (auto fi : range) {
    asawa::shell::FaceId fid = asawa::shell::face_id(fi);
    if (shell.fsize(fid) != 3)
      continue;
    if (asawa::shell::face_area(shell, fid, x) < 1e-8)
      continue;
    vec3 axial = fi < static_cast<int>(axial_face.size()) ? axial_face[fi]
                                                          : vec3::Zero();
    if (!axial.allFinite() || axial.squaredNorm() < 1e-24)
      continue;
    auto tri = shell.get_tri(fid);
    constraints.push_back(anisotropic_triangle_strain::create(
        std::vector<index_t>({tri[0], tri[1], tri[2]}), x, axial, growth, w,
        blocks));
  }
}

/// One `area` row per triangle, holding its area within `band` of its
/// step-start area.
void init_area(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const std::vector<real> &w,
               std::vector<sim_block::ptr> blocks,
               area_mode mode = area_mode::rest, area_band band = area_config{}.band) {

  auto range = shell.get_face_range();
  int i = 0;
  const bool zero = mode == area_mode::zero;
  for (auto fi : range) {
    asawa::shell::FaceId fid = asawa::shell::face_id(fi);
    if (shell.fsize(fid) != 3)
      continue;
    if (asawa::shell::face_area(shell, fid, x) < 1e-8)
      continue;
    auto tri = shell.get_tri(fid);
    constraints.push_back(area::create(std::vector<index_t>({
                                           tri[0],
                                           tri[1],
                                           tri[2],
                                       }),
                                       band.lo, band.hi, w[i++], blocks, zero));
  }
}

void init_area(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const real &w,
               std::vector<sim_block::ptr> blocks,
               area_mode mode = area_mode::rest, area_band band = area_config{}.band) {

  auto range = shell.get_face_range();
  std::vector<real> ws(range.size(), w);
  init_area(shell, constraints, x, ws, blocks, mode, band);
}

void init_area(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const area_config &cfg,
               std::vector<sim_block::ptr> blocks, area_mode mode = area_mode::rest) {
  init_area(shell, constraints, x, cfg.w, std::move(blocks), mode, cfg.band);
}

void init_bending(const asawa::shell::shell &shell,
                  std::vector<projection_constraint::ptr> &constraints,
                  const std::vector<vec3> &x, const real &w,
                  std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    real A =
        asawa::shell::vert_area(shell, asawa::shell::vert_id(iv), x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(asawa::shell::vert_id(iv));
    std::vector<index_t> idx(ring.begin(), ring.end());
    idx.insert(idx.begin(), iv);
    std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, asawa::shell::vert_id(iv), x);
    constraints.push_back(bending::create(idx, weights, w, blocks));
  }
}

void init_arap(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const real &w,
               std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(vid);
    if (ring.empty())
      continue;
    std::vector<index_t> idx(ring.begin(), ring.end());
    idx.insert(idx.begin(), iv);
    std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, vid, x);
    if (weights.size() + 1 != idx.size())
      continue;
    constraints.push_back(arap::create(idx, weights, w, blocks));
  }
}

inline std::vector<vec3>
quantize_normals_axis(const std::vector<vec3> &N) {
  static const vec3 axes[6] = {
      vec3::UnitX(),  -vec3::UnitX(), vec3::UnitY(),
      -vec3::UnitY(), vec3::UnitZ(),  -vec3::UnitZ(),
  };
  std::vector<vec3> out(N.size(), vec3::UnitZ());
  for (size_t i = 0; i < N.size(); ++i) {
    if (N[i].squaredNorm() < real(1e-24)) {
      out[i] = vec3::UnitZ();
      continue;
    }
    const vec3 n = N[i].normalized();
    real best = -real(2.0);
    vec3 g = vec3::UnitZ();
    for (const vec3 &a : axes) {
      const real d = n.dot(a);
      if (d > best) {
        best = d;
        g = a;
      }
    }
    out[i] = g;
  }
  return out;
}

void init_cubification(const asawa::shell::shell &shell,
                       std::vector<projection_constraint::ptr> &constraints,
                       const std::vector<vec3> &x,
                       const std::vector<vec3> &N_rest,
                       const std::vector<vec3> &N_goal, const real &w,
                       real lambda, std::vector<sim_block::ptr> blocks) {
  const size_t nv = static_cast<size_t>(shell.vert_count());
  if (N_rest.size() < nv || N_goal.size() < nv)
    return;

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(vid);
    if (ring.empty())
      continue;
    std::vector<index_t> idx(ring.begin(), ring.end());
    idx.insert(idx.begin(), iv);
    std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, vid, x);
    if (weights.size() + 1 != idx.size())
      continue;
    const real lambda_tip = lambda * A;
    constraints.push_back(arap_normal_align::create(
        idx, weights, N_rest[static_cast<size_t>(iv)],
        N_goal[static_cast<size_t>(iv)], w, lambda_tip, blocks));
  }
}

void init_willmore(const asawa::shell::shell &shell,
                   std::vector<projection_constraint::ptr> &constraints,
                   const std::vector<vec3> &x, const std::vector<real> &w,
                   std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    real A =
        asawa::shell::vert_area(shell, asawa::shell::vert_id(iv), x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(asawa::shell::vert_id(iv));
    std::vector<index_t> idx(ring.begin(), ring.end());
    idx.insert(idx.begin(), iv);
    std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, asawa::shell::vert_id(iv), x);
    real aw = A * w[iv];
    auto constraint = willmore::create(idx, weights, aw, blocks);
    std::vector<char> seen(static_cast<size_t>(shell.vert_count()), 0);
    if (iv >= 0 && static_cast<size_t>(iv) < seen.size())
      seen[static_cast<size_t>(iv)] = 1;
    for (index_t j : idx) {
      if (j >= 0 && static_cast<size_t>(j) < seen.size())
        seen[static_cast<size_t>(j)] = 1;
    }
    std::vector<index_t> outer;
    for (size_t k = 1; k < idx.size(); ++k) {
      const index_t j = idx[k];
      if (j < 0)
        continue;
      for (asawa::shell::VertId vk :
           shell.get_one_ring(asawa::shell::vert_id(j))) {
        const index_t id = static_cast<index_t>(vk);
        if (id < 0 || static_cast<size_t>(id) >= seen.size() ||
            seen[static_cast<size_t>(id)])
          continue;
        seen[static_cast<size_t>(id)] = 1;
        outer.push_back(id);
      }
    }
    constraint->set_fit_ids(outer);
    constraints.push_back(constraint);
  }
}

void init_willmore(const asawa::shell::shell &shell,
                   std::vector<projection_constraint::ptr> &constraints,
                   const std::vector<vec3> &x, const real &w,
                   std::vector<sim_block::ptr> blocks) {

  std::vector<real> ws = std::vector<real>(shell.vert_count(), w);
  init_willmore(shell, constraints, x, ws, blocks);
}

inline std::vector<real>
laplacian_stencil_weights(const asawa::shell::shell &shell,
                          asawa::shell::VertId vi, const std::vector<vec3> &x,
                          laplacian_stencil stencil) {
  switch (stencil) {
  case laplacian_stencil::unitary:
    return asawa::shell::vert_unitary_weights(shell, vi, x);
  case laplacian_stencil::cotan:
    return asawa::shell::vert_cotan_weights(shell, vi, x);
  case laplacian_stencil::angle:
    return asawa::shell::vert_angle_weights(shell, vi, x);
  }
  return asawa::shell::vert_cotan_weights(shell, vi, x);
}

void init_laplacian(const asawa::shell::shell &shell,
                    std::vector<projection_constraint::ptr> &constraints,
                    const std::vector<vec3> &x, laplacian_mode mode,
                    laplacian_stencil stencil, const real &w,
                    std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(vid);
    std::vector<index_t> idx(ring.begin(), ring.end());
    idx.insert(idx.begin(), iv);
    const std::vector<real> weights =
        laplacian_stencil_weights(shell, vid, x, stencil);
    if (weights.size() + 1 != idx.size())
      continue;

    switch (mode) {
    case laplacian_mode::null: {
      // The vertex area weights the row geometry, so `_w` stays the knob.
      std::vector<real> area_weights = weights;
      for (real &k : area_weights)
        k *= A;
      constraints.push_back(laplacian::create(idx, area_weights, w, blocks));
      break;
    }
    case laplacian_mode::smooth:
      constraints.push_back(bending::create(idx, weights, w, blocks));
      break;
    }
  }
}

void init_edge_willmore(const asawa::shell::shell &shell,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &w,
                        std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();
  const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

  int k = 0;
  for (auto c0 : range) {
    asawa::shell::CornerId cid = asawa::shell::corner_id(c0);
    asawa::shell::CornerId c1 = shell.other(cid);
    int i = shell.vert(cid);
    int j = shell.vert(shell.other(cid));

    int i0 = shell.vert(cid);
    int ip = shell.vert(shell.prev(cid));
    int j0 = shell.vert(c1);
    int jp = shell.vert(shell.prev(c1));
    vec3 N = asawa::shell::edge_normal(shell, cid, x);
#if 0
    vec3 xi = asawa::shell::edge_center(shell, c0, x);
    gg::geometry_logger::line(xi, xi + 0.05 * w[k] * N,
                              vec4(0.0, 0.5, 1.0, 1.0));

#endif

    auto constraint = edge_willmore::create(
        std::vector<index_t>({i0, ip, j0, jp}), w[c0 / 2], blocks);
    constraints.push_back(constraint);
    k++;
  }
}

void init_edge_willmore(const asawa::shell::shell &shell,
                        std::vector<projection_constraint::ptr> &constraints,
                        const real &w, std::vector<sim_block::ptr> blocks) {

  std::vector<real> ws = std::vector<real>(shell.edge_count(), w);
  init_edge_willmore(shell, constraints, ws, blocks);
}

void init_edge_edge_collisions(
    asawa::shell::shell &M, asawa::shell::dynamic &dynamic,
    std::vector<projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const real &r, const real &eps, const real &w,
    std::vector<sim_block::ptr> blocks) {

  vector<std::array<index_t, 2>> collisions =
      dynamic.get_internal_edge_edge_collisions(r);
  for (auto &c : collisions) {

    if (c[0] < 0 || c[1] < 0)
      continue;

    asawa::shell::CornerId c00 = asawa::shell::corner_id(c[0]);
    asawa::shell::CornerId c01 = M.other(c00);
    asawa::shell::CornerId c10 = asawa::shell::corner_id(c[1]);
    asawa::shell::CornerId c11 = M.other(c10);

    if (asawa::shell::edges_share_one_ring(M, c00, c10))
      continue;

    index_t v00 = M.vert(c00);
    index_t v00p = M.vert(M.prev(c00));
    index_t v01 = M.vert(c01);
    index_t v01p = M.vert(M.prev(c01));
    index_t v10 = M.vert(c10);
    index_t v10p = M.vert(M.prev(c10));
    index_t v11 = M.vert(c11);
    index_t v11p = M.vert(M.prev(c11));

    std::array<real, 3> d =
        va::distance_Segment_Segment(x[v00], x[v01], x[v10], x[v11]);
    if (d[0] > r)
      continue;

    vec3 NA = asawa::shell::edge_normal(M, c00, x);
    vec3 NB = asawa::shell::edge_normal(M, c10, x);
    real angle = va::dot(NA, NB);
    if (std::abs(angle) < 0.95)
      continue;

    constraints.push_back(edge_edge_normal_collision::create(
        {v00, v00p, v01, v01p, v10, v10p, v11, v11p}, eps, w, blocks));
  }
}

void init_pnt_tri_collisions(
    asawa::shell::shell &M, asawa::shell::dynamic &dynamic,
    std::vector<projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const real &r, const real &eps, const real &w,
    std::vector<sim_block::ptr> blocks, real hardness = 1.0) {

  vector<std::array<index_t, 2>> collisions =
      dynamic.get_internal_pnt_tri_collisions(r);
  for (auto &c : collisions) {
    index_t iv = c[0];
    index_t it = c[1];

    if (it < 0)
      continue;

    std::array<asawa::shell::VertId, 3> tri = M.get_tri(asawa::shell::face_id(it));
    vec3 Nf = asawa::shell::face_normal(M, asawa::shell::face_id(it), x);
    vec3 Nv = asawa::shell::vert_normal(M, asawa::shell::vert_id(iv), x);
    vec3 x0 = x[iv];

    std::array<real, 4> cp =
        va::closest_point({x[tri[0]], x[tri[1]], x[tri[2]]}, x0);
    vec3 xT = cp[1] * x[tri[0]] + cp[2] * x[tri[1]] +
              cp[3] * x[tri[2]];
    vec3 dx = x0 - xT;
    if (dx.squaredNorm() < 1e-24)
      continue;
    real Nfddx = Nf.dot(dx.normalized());
    real NfdNv = Nf.dot(Nv);

    if (abs(Nfddx) < 0.75)
      continue;

    if (abs(NfdNv) < 0.75)
      continue;

    if (cp[0] > r)
      continue;

    if (!asawa::shell::vert_id_valid(tri[0]) ||
        !asawa::shell::vert_id_valid(tri[1]) ||
        !asawa::shell::vert_id_valid(tri[2]))
      continue;

    constraints.push_back(pnt_tri_collision::create(
        {iv, tri[0], tri[1], tri[2]}, eps, w, blocks, hardness));
  }
}

std::vector<fair_ring> build_frame_rings(const asawa::shell::shell &shell,
                                        const std::vector<vec3> &x,
                                        const std::vector<char> &valid) {
  std::vector<fair_ring> rings;
  for (int iv = 0; iv < shell.vert_count(); ++iv) {
    if (iv >= static_cast<int>(valid.size()) || !valid[static_cast<size_t>(iv)])
      continue;
    if (static_cast<int>(shell.vbegin(asawa::shell::vert_id(iv))) < 0)
      continue;
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    auto ring = shell.get_one_ring(vid);
    if (ring.empty())
      continue;
    std::vector<index_t> idx(ring.begin(), ring.end());
    const std::vector<real> weights =
        asawa::shell::vert_cotan_weights(shell, vid, x);
    if (weights.size() != idx.size())
      continue;
    fair_ring fr;
    fr.i = iv;
    fr.spokes.reserve(weights.size());
    for (size_t k = 0; k < weights.size(); ++k) {
      fair_spoke sp;
      sp.j = idx[k];
      sp.K = weights[k];
      const index_t j = sp.j;
      if (j >= 0 && static_cast<size_t>(j) < valid.size() &&
          valid[static_cast<size_t>(j)] && std::isfinite(weights[k]) &&
          weights[k] != 0.0)
        sp.ok = 1;
      fr.spokes.push_back(sp);
    }
    rings.push_back(std::move(fr));
  }
  return rings;
}

/// Vertex ids in the closed 2-ring of `vi` (1-ring ∪ neighbors-of-neighbors),
/// excluding the center. Used as osculating-fit samples for spine projectors.
inline std::vector<index_t>
vertex_two_ring_ids(const asawa::shell::shell &shell, index_t vi) {
  std::vector<index_t> out;
  const int nV = shell.vert_count();
  if (vi < 0 || vi >= nV)
    return out;
  std::vector<char> seen(static_cast<size_t>(nV), 0);
  seen[static_cast<size_t>(vi)] = 1;
  const auto ring1 = shell.get_one_ring(asawa::shell::vert_id(vi));
  out.reserve(ring1.size() * 4);
  for (const asawa::shell::VertId nj : ring1) {
    const int j = static_cast<int>(nj);
    if (j < 0 || j >= nV || seen[static_cast<size_t>(j)])
      continue;
    seen[static_cast<size_t>(j)] = 1;
    out.push_back(j);
  }
  for (const asawa::shell::VertId nj : ring1) {
    const int j = static_cast<int>(nj);
    if (j < 0)
      continue;
    for (const asawa::shell::VertId nk :
         shell.get_one_ring(asawa::shell::vert_id(j))) {
      const int k = static_cast<int>(nk);
      if (k < 0 || k >= nV || seen[static_cast<size_t>(k)])
        continue;
      seen[static_cast<size_t>(k)] = 1;
      out.push_back(k);
    }
  }
  return out;
}

void init_frame_fair(const asawa::shell::shell &shell,
                     std::vector<projection_constraint::ptr> &constraints,
                     const std::vector<vec3> &x, const std::vector<char> &valid,
                     const real &w, std::vector<sim_block::ptr> blocks,
                     const std::vector<objective> &objectives = {},
                     int debug_vert = -1,
                     const std::vector<vertex_chart> *charts = nullptr,
                     sim_block::ptr edge_store = nullptr,
                     const kusama::weingarten::vertex_stars *stars = nullptr,
                     real spine_step = real(1)) {
  if (w == 0.0)
    return;
  std::vector<fair_ring> rings = build_frame_rings(shell, x, valid);
  if (rings.empty())
    return;
  int coeff = -1;
  int nE = 0;
  bool stein = false;
  cr_complex cx;
  cr_coeff_block::ptr cr;
  if (!charts)
    cr = cr_coeff_of(blocks, &coeff);
  if (!charts && cr && cr->nE > 0 && coeff >= 0 && shell.verts_are_dense_packed()) {
    auto &mesh = const_cast<asawa::shell::shell &>(shell);
    if (build_cr_complex(mesh, cx) &&
        static_cast<int>(cx.edges.size()) == cr->nE) {
      for (fair_ring &ring : rings)
        stamp_cr_vertex(shell, cx, ring.i, ring.faces, ring.edges);
      nE = cr->nE;
      stein = true;
    } else {
      coeff = -1;
    }
  } else {
    coeff = -1;
  }
  int shape = -1;
  for (int i = 1; i < static_cast<int>(blocks.size()); ++i) {
    if (std::dynamic_pointer_cast<vec3_block>(blocks[static_cast<size_t>(i)]))
      shape = i;
  }
  std::vector<sim_block::ptr> row_blocks;
  if (charts) {
    row_blocks.push_back(blocks[0]);
    if (edge_store && stars)
      row_blocks.push_back(edge_store);
    else
      stars = nullptr;
  } else if (stein) {
    row_blocks.push_back(blocks[0]);
    row_blocks.push_back(blocks[static_cast<size_t>(coeff)]);
    if (shape >= 0)
      row_blocks.push_back(blocks[static_cast<size_t>(shape)]);
  } else {
    row_blocks = blocks;
  }
  std::vector<int> at(static_cast<size_t>(shell.vert_count()), -1);
  for (int r = 0; r < static_cast<int>(rings.size()); ++r)
    at[static_cast<size_t>(rings[static_cast<size_t>(r)].i)] = r;
  for (const fair_ring &ring : rings) {
    std::vector<index_t> ids{ring.i};
    std::vector<real> weights;
    ids.reserve(ring.spokes.size() + 1);
    weights.reserve(ring.spokes.size());
    std::vector<nbr_fan> nbrs;
    for (const fair_spoke &sp : ring.spokes) {
      ids.push_back(sp.j);
      weights.push_back(sp.K);
      if (sp.j < 0 || sp.j >= static_cast<int>(at.size()))
        continue;
      const int jr = at[static_cast<size_t>(sp.j)];
      if (jr < 0)
        continue;
      const fair_ring &nb = rings[static_cast<size_t>(jr)];
      nbr_fan fan;
      fan.j = nb.i;
      fan.tangent = nb.edges;
      if (stein)
        fan_from_ring(cx, nb.faces, nb.edges, fan.faces, fan.support, fan.edges);
      for (const fair_spoke &nsp : nb.spokes) {
        fan.spokes.push_back(nsp.j);
        fan.spoke_k.push_back(nsp.K);
        fan.spoke_ok.push_back(nsp.ok);
      }
      nbrs.push_back(std::move(fan));
    }
    std::vector<fan_face> faces;
    std::vector<fan_face> support;
    std::vector<fan_edge> edges;
    std::vector<int> tangent = ring.edges;
    if (stein)
      fan_from_ring(cx, ring.faces, ring.edges, faces, support, edges);
    std::vector<index_t> fit_ids = vertex_two_ring_ids(shell, ring.i);
    constraints.push_back(frame_fair::create(
        std::move(ids), std::move(weights), std::move(faces), std::move(support),
        std::move(edges), std::move(tangent), std::move(nbrs), objectives, stein,
        nE, debug_vert, w, row_blocks, charts, std::move(fit_ids), spine_step,
        stars));
  }
}

/// One cotangent row per vertex. Blocks are positions only. The fit reads the
/// ring points and the normals of the fan triangles.
void init_albers_fair(const asawa::shell::shell &shell,
                      std::vector<projection_constraint::ptr> &constraints,
                      const std::vector<vec3> &x, const std::vector<char> &valid,
                      const real &w, std::vector<sim_block::ptr> blocks,
                      const std::vector<objective> &objectives) {
  if (!(w > real(0)) || !std::isfinite(w) || blocks.empty() || !blocks[0] ||
      !albers_mask_of(objectives).any())
    return;
  const std::vector<fair_ring> rings = build_frame_rings(shell, x, valid);
  std::vector<sim_block::ptr> pos{blocks[0]};
  for (const fair_ring &ring : rings) {
    std::vector<index_t> ids{ring.i};
    std::vector<real> weights;
    ids.reserve(ring.spokes.size() + 1);
    weights.reserve(ring.spokes.size());
    for (const fair_spoke &sp : ring.spokes) {
      ids.push_back(sp.j);
      weights.push_back(sp.K);
    }
    constraints.push_back(albers_ring_fair::create(
        std::move(ids), std::move(weights), objectives, w, pos));
  }
}

void init_frame_minimizer(const asawa::shell::shell &shell,
                          std::vector<projection_constraint::ptr> &constraints,
                          const std::vector<vec3> &x,
                          const std::vector<char> &valid, const real &w,
                          std::vector<sim_block::ptr> blocks) {
  if (w == 0.0)
    return;
  std::vector<fair_ring> rings = build_frame_rings(shell, x, valid);
  if (rings.empty())
    return;
  std::vector<int> at(static_cast<size_t>(shell.vert_count()), -1);
  for (int r = 0; r < static_cast<int>(rings.size()); ++r)
    at[static_cast<size_t>(rings[static_cast<size_t>(r)].i)] = r;
  for (const fair_ring &self : rings) {
    std::vector<int> need;
    need.push_back(self.i);
    for (const fair_spoke &sp : self.spokes) {
      if (sp.j < 0 || sp.j >= static_cast<int>(at.size()))
        continue;
      const int jr = at[static_cast<size_t>(sp.j)];
      if (jr < 0)
        continue;
      need.push_back(sp.j);
      for (const fair_spoke &nsp : rings[static_cast<size_t>(jr)].spokes) {
        if (nsp.j < 0 || nsp.j >= static_cast<int>(at.size()))
          continue;
        if (at[static_cast<size_t>(nsp.j)] >= 0)
          need.push_back(nsp.j);
      }
    }
    std::sort(need.begin(), need.end());
    need.erase(std::unique(need.begin(), need.end()), need.end());
    std::vector<fair_ring> local;
    local.reserve(need.size());
    for (const int id : need)
      local.push_back(rings[static_cast<size_t>(at[static_cast<size_t>(id)])]);
    constraints.push_back(
        frame_minimizer::create(std::move(local), self.i, w, blocks));
  }
}
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif