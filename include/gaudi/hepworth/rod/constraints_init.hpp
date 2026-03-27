#ifndef __HEP_ROD_CONSTRAINTS_INIT__
#define __HEP_ROD_CONSTRAINTS_INIT__

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

#include "../collision_constraint.hpp"
#include "constraints.hpp"
#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace rod {

void init_smooth(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const real &w) {
  using asawa::rod::corner_id;

  for (int i = 0; i < rod.corner_count(); i++) {
    auto ci = corner_id(i);
    auto in = rod.next(ci);
    auto ip = rod.prev(ci);
    constraints.push_back(smooth::create({ci, ip, in}, w));
  }
}
#if 1
void init_helicity(const asawa::rod::rod &rod,
                   std::vector<projection_constraint::ptr> &constraints,
                   const real &w) {
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
    CornerId in0 = rod.next(ci);
    if (in0 < 0)
      continue;
    CornerId in1 = rod.next(in0);
    if (in1 < 0)
      continue;
    CornerId in2 = rod.next(in1);
    if (in2 < 0)
      continue;
    constraints.push_back(
        cylinder::create({ci, ip2, ip1, ip0, ci, in0, in1, in2}, 1.25));
  }
}
#endif

void init_stretch_shear(const asawa::rod::rod &rod,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &l0, const real &w) {
  using asawa::rod::corner_id;
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    auto ci = corner_id(i);
    auto j = rod.next(ci);
    constraints.push_back(stretch_shear::create({ci, j, ci, Ni}, w, l0[i]));
  }
}

void init_straight(const asawa::rod::rod &rod,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w) {
  using asawa::rod::corner_id;
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    auto ci = corner_id(i);
    auto j = rod.next(ci);
    constraints.push_back(straight::create({ci, j, Ni}, w));
  }
}

void init_angle(const asawa::rod::rod &rod,
                std::vector<projection_constraint::ptr> &constraints,
                const vec3 &z, const real &phi, const real &w) {
  using asawa::rod::corner_id;
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    auto ci = corner_id(i);
    auto j = rod.next(ci);
    constraints.push_back(angle::create({ci, j, Ni}, z, phi, w));
  }
}

void init_collisions(asawa::rod::rod &rod, asawa::rod::dynamic &dynamic,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w) {
  const std::vector<vec3> &x = rod.__x;
  vector<std::array<index_t, 2>> collisions = dynamic.get_internal_collisions();
  for (auto &c2 : collisions) {
    if (c2[0] > -1) {
      std::array<index_t, 4> c4 = dynamic.get_collision_ids(c2);
      vec3 xA0 = x[c4[0]];
      vec3 xA1 = x[c4[1]];
      vec3 xB0 = x[c4[2]];
      vec3 xB1 = x[c4[3]];
      if (rod.prev(asawa::rod::corner_id(c4[0])) == c4[2])
        continue;
      if (rod.next(asawa::rod::corner_id(c4[1])) == c4[3])
        continue;
      if (rod.prev(asawa::rod::corner_id(c4[0])) == c4[3])
        continue;
      if (rod.next(asawa::rod::corner_id(c4[1])) == c4[2])
        continue;
      // std::cout << c[0] << " " << c[1] << " - " << c[2] << " " << c[3]
      //           << std::endl;

      /*
    std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, xA0, xA1);
    vec3 xB = va::mix(t, xB0, xB1);
    vec3 xAB = xB - xA;

    real l = xAB.norm();
    real dl = rod._r - l;
    if (dl < 0)
      continue;
*/
      constraints.push_back(hepworth::edge_edge_collision::create(
          {c4[0], c4[1], c4[2], c4[3]}, w, 1.0 * rod._r));
    }
  }
}
} // namespace rod
} // namespace hepworth
} // namespace gaudi

#endif