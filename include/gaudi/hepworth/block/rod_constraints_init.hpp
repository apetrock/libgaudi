#ifndef __HEP_ROD_BLOCK_CONSTRAINTS_INIT__
#define __HEP_ROD_BLOCK_CONSTRAINTS_INIT__

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
#include "rod_collision_constraint.hpp"
#include "rod_constraints.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

void init_smooth(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const real &w, std::vector<sim_block::ptr> blocks) {
  // block is just a single rod position block Xr
  for (int i = 0; i < rod.corner_count(); i++) {
    asawa::rod::consec_t c = rod.consec(i);
    constraints.push_back(smooth::create({c[1], c[0], c[2]}, w, blocks));
  }
}

#if 1
void init_helicity(const asawa::rod::rod &rod,
                   std::vector<projection_constraint::ptr> &constraints,
                   const real &w, std::vector<sim_block::ptr> blocks) {

  for (int i = 0; i < rod.corner_count(); i++) {
    index_t ip0 = rod.prev(i);
    if (ip0 < 0)
      continue;
    index_t ip1 = rod.prev(ip0);
    if (ip1 < 0)
      continue;
    index_t ip2 = rod.prev(ip1);
    if (ip2 < 0)
      continue;
    index_t ip3 = rod.prev(ip2);
    if (ip3 < 0)
      continue;
    index_t in0 = rod.next(i);
    if (in0 < 0)
      continue;
    index_t in1 = rod.next(in0);
    if (in1 < 0)
      continue;
    index_t in2 = rod.next(in1);
    if (in2 < 0)
      continue;
    index_t in3 = rod.next(in2);
    if (in3 < 0)
      continue;
    constraints.push_back(helicitiy::create(
        {i, ip3, ip2, ip1, ip0, i, in0, in1, in2, in3}, 1.25, blocks));
  }
}
#endif

void init_stretch_shear(const asawa::rod::rod &R,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &l0, const real &w,
                        std::vector<sim_block::ptr> blocks) {
  int Ni = R.corner_count();
  std::vector<index_t> verts = R.get_vert_range();
  for (auto i : verts) {
    asawa::rod::consec_t c = R.consec(i);
    if (l0[i] < 1e-10)
      continue;

    constraints.push_back(
        stretch_shear::create({c[1], c[2], c[1], Ni}, w, l0[i], blocks));
  }
}

void init_bend_twist(const asawa::rod::rod &R,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w, std::vector<sim_block::ptr> blocks,
                     bool skip = false) {
  int Ni = R.corner_count();
  std::vector<index_t> verts = R.get_vert_range();
  int N = skip ? verts.size() - 1 : verts.size();
  for (int i = 0; i < N; i++) {
    if (R.length(i) < 1e-8)
      continue;
    asawa::rod::consec_t c = R.consec(verts[i]);
    constraints.push_back(bend_twist::create({c[1], c[2], Ni}, w, blocks));
  }
}

void init_angle(const asawa::rod::rod &R,
                std::vector<projection_constraint::ptr> &constraints,
                const vec3 &z, const real &phi, const real &w,
                std::vector<sim_block::ptr> blocks) {
  int Ni = R.corner_count();
  const std::vector<vec3> &x = R.__x;
  std::vector<index_t> verts = R.get_vert_range();
  for (auto i : verts) {
    if (R.length(i) < 1e-8)
      continue;
    asawa::rod::consec_t c = R.consec(i);
    constraints.push_back(angle::create({c[1], c[2], Ni}, z, phi, w, blocks));
  }
}

void init_collisions(asawa::rod::rod &R, asawa::rod::dynamic &dynamic,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w, std::vector<sim_block::ptr> blocks,
                     real K = 1.0) {
  const std::vector<vec3> &x = R.__x;
  vector<std::array<index_t, 4>> collisions =
      dynamic.get_internal_collisions(K);
  for (auto &c : collisions) {
    if (c[0] > -1) {
      vec3 xA0 = x[c[0]];
      vec3 xA1 = x[c[1]];
      vec3 xB0 = x[c[2]];
      vec3 xB1 = x[c[3]];

      if (R.length(c[0]) < 1e-8)
        continue;
      if (R.length(c[1]) < 1e-8)
        continue;
      if (R.length(c[2]) < 1e-8)
        continue;
      if (R.length(c[3]) < 1e-8)
        continue;

      if (R.prev(c[0]) == c[2])
        continue;
      if (R.next(c[1]) == c[3])
        continue;
      if (R.prev(c[0]) == c[3])
        continue;
      if (R.next(c[1]) == c[2])
        continue;

      constraints.push_back(
          rod_collision::create({c[0], c[1], c[2], c[3]}, w, K * R._r, blocks));
    }
  }
}
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif