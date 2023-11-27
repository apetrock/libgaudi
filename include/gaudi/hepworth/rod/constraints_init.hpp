#ifndef __HEP_ROD_CONSTRAINTS_INIT__
#define __HEP_ROD_CONSTRAINTS_INIT__

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

#include "../collision_constraint.hpp"
#include "constraints.hpp"
#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace rod {

void init_smooth(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const real &w) {

  for (int i = 0; i < rod.corner_count(); i++) {
    index_t in = rod.next(i);
    index_t ip = rod.prev(i);
    constraints.push_back(smooth::create({i, ip, in}, w));
  }
}
#if 1
void init_helicity(const asawa::rod::rod &rod,
                   std::vector<projection_constraint::ptr> &constraints,
                   const real &w) {

  for (int i = 0; i < rod.corner_count(); i++) {
    index_t ip0 = rod.prev(i);
    if (!ip0)
      continue;
    index_t ip1 = rod.prev(ip0);
    if (!ip1)
      continue;
    index_t ip2 = rod.prev(ip1);
    if (!ip2)
      continue;
    index_t in0 = rod.next(i);
    if (!in0)
      continue;
    index_t in1 = rod.next(in0);
    if (!in1)
      continue;
    index_t in2 = rod.next(in1);
    if (!in2)
      continue;
    constraints.push_back(
        cylinder::create({i, ip2, ip1, ip0, i, in0, in1, in2}, 1.25));
  }
}
#endif

void init_stretch_shear(const asawa::rod::rod &rod,
                        std::vector<projection_constraint::ptr> &constraints,
                        const std::vector<real> &l0, const real &w) {
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    index_t j = rod.next(i);
    constraints.push_back(stretch_shear::create({i, j, i, Ni}, w, l0[i]));
  }
}

void init_bend_twist(const asawa::rod::rod &rod,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w) {
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    index_t j = rod.next(i);
    constraints.push_back(bend_twist::create({i, j, Ni}, w));
  }
}

void init_angle(const asawa::rod::rod &rod,
                std::vector<projection_constraint::ptr> &constraints,
                const vec3 &z, const real &phi, const real &w) {
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    index_t j = rod.next(i);
    constraints.push_back(angle::create({i, j, Ni}, z, phi, w));
  }
}

void init_collisions(asawa::rod::rod &rod, asawa::rod::dynamic &dynamic,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w) {
  const std::vector<vec3> &x = rod.__x;
  vector<std::array<index_t, 4>> collisions = dynamic.get_internal_collisions();
  for (auto &c : collisions) {
    if (c[0] > -1) {
      vec3 xA0 = x[c[0]];
      vec3 xA1 = x[c[1]];
      vec3 xB0 = x[c[2]];
      vec3 xB1 = x[c[3]];
      if (rod.prev(c[0]) == c[2])
        continue;
      if (rod.next(c[1]) == c[3])
        continue;
      if (rod.prev(c[0]) == c[3])
        continue;
      if (rod.next(c[1]) == c[2])
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
          {c[0], c[1], c[2], c[3]}, w, 1.0 * rod._r));
    }
  }
}
} // namespace rod
} // namespace hepworth
} // namespace gaudi

#endif