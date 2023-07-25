#ifndef __HEP_SHELL_BLOCK_CONSTRAINTS_INIT__
#define __HEP_SHELL_BLOCK_CONSTRAINTS_INIT__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#include "gaudi/common.h"
#include "shell_collision_constraint.hpp"
#include "shell_constraints.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

void init_edge_growth(const asawa::shell::shell &shell,
                      std::vector<projection_constraint::ptr> &constraints,
                      const std::vector<vec3> &x, const std::vector<real> &l0,
                      const real &g, const real &w,
                      std::vector<sim_block::ptr> blocks) {

  std::vector<real> ct = asawa::shell::edge_cotan_weights(shell, x);

  auto range = shell.get_edge_range();
  int k = 0;
  for (auto c0 : range) {
    int i = shell.vert(c0);
    int j = shell.vert(shell.other(c0));

    if (std::isnan(l0[k]))
      continue;

    real cti = ct[k];
    if (std::isnan(cti))
      continue;
    cti = va::clamp(cti, 0.0, 10.0);

    real Ai = asawa::shell::face_area(shell, shell.face(c0), x);
    real Aj = asawa::shell::face_area(shell, shell.face(shell.other(c0)), x);
    if (Ai < 1e-8)
      continue;
    if (Aj < 1e-8)
      continue;

    constraints.push_back(
        edge_strain::create(std::vector<index_t>({i, j}), 0.85 * l0[c0 / 2], g,
                            w * pow(cti, 2.0) / (Ai + Aj), blocks));
    k++;
  }
}

void init_edge_strain(const asawa::shell::shell &shell,
                      std::vector<projection_constraint::ptr> &constraints,
                      const std::vector<vec3> &x, const std::vector<real> &l0,
                      const real &w, std::vector<sim_block::ptr> blocks) {

  std::vector<real> ct = asawa::shell::edge_cotan_weights(shell, x);

  auto range = shell.get_edge_range();
  int i = 0;
  for (auto c0 : range) {
    int i = shell.vert(c0);
    int j = shell.vert(shell.other(c0));

    if (std::isnan(l0[c0 / 2]))
      continue;
    if (l0[c0 / 2] < 1e-8)
      continue;

    real ct = asawa::shell::cotan(shell, c0, x);
    real Ai = asawa::shell::face_area(shell, shell.face(c0), x);
    real Aj = asawa::shell::face_area(shell, shell.face(shell.other(c0)), x);
    if (Ai < 1e-8)
      continue;
    if (Aj < 1e-8)
      continue;

    constraints.push_back(edge_strain::create(std::vector<index_t>({i, j}),
                                              l0[c0 / 2], 1.0, w, blocks));
  }
}

void init_bending(const asawa::shell::shell &shell,
                  std::vector<projection_constraint::ptr> &constraints,
                  const std::vector<vec3> &x, const real &w,
                  std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    real A = asawa::shell::vert_area(shell, iv, x);
    if (A < 1e-8)
      continue;
    std::vector<index_t> idx = shell.get_one_ring(iv);
    idx.insert(idx.begin(), iv);
    std::vector<real> weights = asawa::shell::vert_cotan_weights(shell, iv, x);
    constraints.push_back(bending::create(idx, weights, x, w, blocks));
  }
}

void init_willmore(const asawa::shell::shell &shell,
                   std::vector<projection_constraint::ptr> &constraints,
                   const std::vector<vec3> &x, const std::vector<real> &w,
                   std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    real A = asawa::shell::vert_area(shell, iv, x);
    if (A < 1e-8)
      continue;
    std::vector<index_t> idx = shell.get_one_ring(iv);
    idx.insert(idx.begin(), iv);
    std::vector<real> weights = asawa::shell::vert_cotan_weights(shell, iv, x);
    real aw = A * w[iv];
    constraints.push_back(willmore::create(idx, weights, x, aw, blocks));
  }
}

void init_willmore(const asawa::shell::shell &shell,
                   std::vector<projection_constraint::ptr> &constraints,
                   const std::vector<vec3> &x, const real &w,
                   std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();
  std::vector<real> ws = std::vector<real>(range.size(), w);
  init_willmore(shell, constraints, x, ws, blocks);
}

void init_area(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const std::vector<real> &w,
               std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_face_range();
  int i = 0;
  for (auto fi : range) {
    if (shell.fsize(fi) != 3)
      continue;
    if (asawa::shell::face_area(shell, fi, x) < 1e-8)
      continue;
    auto tri = shell.get_tri(fi);
    constraints.push_back(area::create(std::vector<index_t>({
                                           tri[0],
                                           tri[1],
                                           tri[2],
                                       }),
                                       x, 0.98, 1.02, w[i++], blocks));
  }
}

void init_area(const asawa::shell::shell &shell,
               std::vector<projection_constraint::ptr> &constraints,
               const std::vector<vec3> &x, const real &w,
               std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_face_range();
  std::vector<real> ws(range.size(), w);
  init_area(shell, constraints, x, ws, blocks);
}

void init_laplacian(const asawa::shell::shell &shell,
                    std::vector<projection_constraint::ptr> &constraints,
                    const std::vector<vec3> &x, int weight_type, const real &w,
                    std::vector<sim_block::ptr> blocks) {

  for (int iv = 0; iv < shell.vert_count(); iv++) {
    real A = asawa::shell::vert_area(shell, iv, x);
    if (A < 1e-8)
      continue;
    std::vector<index_t> idx = shell.get_one_ring(iv);
    idx.insert(idx.begin(), iv);

    std::vector<real> weights;
    switch (weight_type) {
    case 0:
      weights = asawa::shell::vert_unitary_weights(shell, iv, x);
    case 1:
      weights = asawa::shell::vert_cotan_weights(shell, iv, x);
    case 2:
      weights = asawa::shell::vert_angle_weights(shell, iv, x);
    }

    constraints.push_back(laplacian::create(idx, weights, x, w * A, blocks));
  }
}

void init_edge_willmore(const asawa::shell::shell &shell,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &w,
                        std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();
  const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

  int i = 0;
  for (auto c0 : range) {
    index_t c1 = shell.other(c0);
    int i = shell.vert(c0);
    int j = shell.vert(shell.other(c0));

    int i0 = shell.vert(c0);
    int ip = shell.vert(shell.prev(c0));
    int j0 = shell.vert(c1);
    int jp = shell.vert(shell.prev(c1));
    vec3 N = asawa::shell::edge_normal(shell, c0, x);

    auto constraint = edge_willmore::create(
        std::vector<index_t>({i0, ip, j0, jp}), w[i++], blocks);
    constraint->set_align_normal(N);
    constraints.push_back(constraint);
  }
}

void init_edge_willmore(const asawa::shell::shell &shell,
                        std::vector<projection_constraint::ptr> &constraints,
                        const real &w, std::vector<sim_block::ptr> blocks) {

  auto range = shell.get_edge_range();
  std::vector<real> ws = std::vector<real>(range.size(), w);
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

    index_t c00 = c[0];
    index_t c01 = M.other(c[0]);
    index_t c10 = c[1];
    index_t c11 = M.other(c[1]);

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

    vec3 NA = asawa::shell::edge_normal(M, c[0], x);
    vec3 NB = asawa::shell::edge_normal(M, c[1], x);
    real angle = va::dot(NA, NB);

    if (abs(angle) < 0.95)
      continue;

    constraints.push_back(edge_edge_normal_collision::create(
        {v00, v00p, v01, v01p, v10, v10p, v11, v11p}, eps, w, blocks));
  }
}

void init_pnt_tri_collisions(
    asawa::shell::shell &M, asawa::shell::dynamic &dynamic,
    std::vector<projection_constraint::ptr> &constraints,
    const std::vector<vec3> &x, const real &r, const real &eps, const real &w,
    std::vector<sim_block::ptr> blocks) {

  vector<std::array<index_t, 2>> collisions =
      dynamic.get_internal_pnt_tri_collisions(r);
  for (auto &c : collisions) {
    index_t iv = c[0];
    index_t it = c[1];

    if (it < 0)
      continue;

    std::array<index_t, 3> tri = M.get_tri(it);
    vec3 cen = asawa::shell::face_center(M, it, x);
    vec3 Nf = asawa::shell::face_normal(M, it, x);
    vec3 Nv = asawa::shell::vert_normal(M, iv, x);
    vec3 x0 = x[iv];

    std::array<real, 4> cp =
        va::closest_point({x[tri[0]], x[tri[1]], x[tri[2]]}, x0);
    vec3 xT = cp[1] * x[tri[0]] + cp[2] * x[tri[1]] + cp[3] * x[tri[2]];
    vec3 dx = x0 - xT;
    real Nfddx = Nf.dot(dx.normalized());
    real NfdNv = Nf.dot(Nv);

    if (abs(Nfddx) < 0.75)
      continue;

    if (abs(NfdNv) < 0.75)
      continue;

    if (cp[0] > r)
      continue;

    if (tri[0] * tri[1] * tri[2] < 0)
      continue;

    // gg::geometry_logger::line(x0, xT, vec4(1.0, r / dx.norm(), 1.0, 1.0));
    // gg::geometry_logger::line(cen, cen + r * N, vec4(1.0, 0.0, 1.0, 1.0));

    // gg::geometry_logger::line(x0, cen, vec4(0.0, 1.0, 0.0, 1.0));
    constraints.push_back(pnt_tri_collision::create(
        {iv, tri[0], tri[1], tri[2]}, Nv, eps, w, blocks));
  }
}
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif