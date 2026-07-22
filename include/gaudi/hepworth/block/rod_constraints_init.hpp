#ifndef __HEP_ROD_BLOCK_CONSTRAINTS_INIT__
#define __HEP_ROD_BLOCK_CONSTRAINTS_INIT__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <random>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "rod_collision_constraint.hpp"
#include "rod_constraints.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

void init_smooth(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const real &w, std::vector<sim_block::ptr> blocks) {
  for (int i = 0; i < rod.corner_count(); i++) {
    asawa::rod::consec_t c = rod.consec(asawa::rod::corner_id(i));
    constraints.push_back(smooth::create({c[1], c[0], c[2]}, w, blocks));
  }
}

#if 1
void init_helicity(const asawa::rod::rod &rod,
                   std::vector<projection_constraint::ptr> &constraints,
                   const real &w, std::vector<sim_block::ptr> blocks) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;

  for (int i = 0; i < rod.corner_count(); i++) {
    CornerId ci = corner_id(i);
    CornerId ip0 = rod.prev(ci);
    if (ip0 < 0)
      continue;
    CornerId ip1 = rod.prev(ip0);
    if (ip1 < 0)
      continue;
    CornerId ip2 = rod.prev(ip1);
    if (ip2 < 0)
      continue;
    CornerId ip3 = rod.prev(ip2);
    if (ip3 < 0)
      continue;
    CornerId in0 = rod.next(ci);
    if (in0 < 0)
      continue;
    CornerId in1 = rod.next(in0);
    if (in1 < 0)
      continue;
    CornerId in2 = rod.next(in1);
    if (in2 < 0)
      continue;
    CornerId in3 = rod.next(in2);
    if (in3 < 0)
      continue;
    constraints.push_back(helicitiy::create(
        {ci, ip3, ip2, ip1, ip0, ci, in0, in1, in2, in3}, w, blocks));
  }
}
#endif

void init_stretch_shear(const asawa::rod::rod &R,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &l0, const real &w,
                        std::vector<sim_block::ptr> blocks) {
  int Ni = R.corner_count();
  auto verts = R.get_vert_range();
  for (auto i : verts) {
    asawa::rod::consec_t c = R.consec(i);
    if (l0[i] < 1e-6)
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
  auto verts = R.get_vert_range();
  const std::vector<quat> &q = R.__u;
  int N = skip ? verts.size() - 1 : verts.size();
  for (int i = 0; i < N; i++) {
    if (R.length(verts[i]) < 1e-6)
      continue;
    asawa::rod::consec_t c = R.consec(verts[i]);
    constraints.push_back(bend_twist::create({c[1], c[2], Ni}, q, w, blocks));
  }
}

void init_straight(const asawa::rod::rod &R,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w, std::vector<sim_block::ptr> blocks,
                     bool skip = false) {
  int Ni = R.corner_count();
  auto verts = R.get_vert_range();
  int N = skip ? verts.size() - 1 : verts.size();
  for (int i = 0; i < N; i++) {
    if (R.length(verts[i]) < 1e-6)
      continue;
    asawa::rod::consec_t c = R.consec(verts[i]);
    constraints.push_back(straight::create({c[1], c[2], Ni}, w, blocks));
  }
}

void init_angle(const asawa::rod::rod &R,
                std::vector<projection_constraint::ptr> &constraints,
                const vec3 &z, const real &phi, const real &w,
                std::vector<sim_block::ptr> blocks) {
  int Ni = R.corner_count();
  const std::vector<vec3> &x = R.__x;
  auto verts = R.get_vert_range();
  for (auto i : verts) {
    if (R.length(i) < 1e-8)
      continue;
    asawa::rod::consec_t c = R.consec(i);
    constraints.push_back(angle::create({c[1], c[2], Ni}, z, phi, w, blocks));
  }
}

// Skip rod-rod pairs whose endpoints are within `min_sep` chain hops.
// min_sep=1 ≈ old prev/next neighbor filter; 2 also drops one-hop pairs
// (AB vs CD) that otherwise light up along gently curved rods.
inline bool rod_edges_chain_neighbors(const asawa::rod::rod &R, index_t a0,
                                      index_t a1, index_t b0, index_t b1,
                                      int min_sep = 2) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;
  if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
    return true;
  if (min_sep <= 0)
    return false;

  auto reaches = [&](index_t from, index_t target) {
    CornerId c = corner_id(from);
    for (int i = 0; i < min_sep; ++i) {
      c = R.next(c);
      if (c < 0)
        break;
      if (static_cast<index_t>(c) == target)
        return true;
    }
    c = corner_id(from);
    for (int i = 0; i < min_sep; ++i) {
      c = R.prev(c);
      if (c < 0)
        break;
      if (static_cast<index_t>(c) == target)
        return true;
    }
    return false;
  };

  return reaches(a0, b0) || reaches(a0, b1) || reaches(a1, b0) ||
         reaches(a1, b1);
}

void init_collisions(asawa::rod::rod &R, asawa::rod::dynamic &dynamic,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w, std::vector<sim_block::ptr> blocks,
                     real K = 1.0, int chain_sep = 2, real geom_margin = 0.1,
                     real r_override = -1.0) {
  using asawa::rod::corner_id;
  const std::vector<vec3> &x = R.__x;
  // True seg-seg gate (AABB false positives can have midpoints ≫ R apart).
  // Note: va::distance_Segment_Segment returns *squared* distance in d[0].
  const real r = (r_override < 0.0) ? R._r : r_override;
  const real max_sep = (1.0 + geom_margin) * r;
  const real max_sep2 = max_sep * max_sep;
  vector<std::array<index_t, 2>> collisions =
      dynamic.get_internal_collisions(K);

  for (auto &c : collisions) {
    if (c[0] > -1) {
      std::array<index_t, 4> c4 = dynamic.get_collision_ids(c);
      vec3 xA0 = x[c4[0]];
      vec3 xA1 = x[c4[1]];
      vec3 xB0 = x[c4[2]];
      vec3 xB1 = x[c4[3]];

      if (R.length(corner_id(c[0])) < 1e-8)
        continue;
      if (R.length(corner_id(c[1])) < 1e-8)
        continue;

      if (rod_edges_chain_neighbors(R, c4[0], c4[1], c4[2], c4[3], chain_sep))
        continue;

      const std::array<real, 3> d =
          va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
      if (!(d[0] < max_sep2))
        continue;

      constraints.push_back(
          rod_collision::create({c4[0], c4[1], c4[2], c4[3]}, w, K * R._r, blocks));
    }
  }
}
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif