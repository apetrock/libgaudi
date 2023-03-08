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

#include "constraints.hpp"
#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace rod {
void init_length(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const std::vector<real> &l0, const real &w) {
  for (int i = 0; i < rod.corner_count(); i++) {
    index_t j = rod.next(i);
    // l0[i] *= 1.01;
    constraints.push_back(length::create({i, j}, w, l0[i]));
  }
} // namespace
  // asawa::rod&rod,std::vector<projection_constraint::ptr>&constraints,conststd::vector<real>&l0,constreal&w)

void init_smooth(const asawa::rod::rod &rod,
                 std::vector<projection_constraint::ptr> &constraints,
                 const real &w) {

  for (int i = 0; i < rod.corner_count(); i++) {
    index_t in = rod.next(i);
    index_t ip = rod.prev(i);
    constraints.push_back(smooth::create({i, ip, in}, w));
  }
}

void init_cylinder(const asawa::rod::rod &rod,
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
        cylinder::create({ip2, ip1, ip0, i, in0, in1, in2}, 1.25));
  }
}
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
} // namespace rod
} // namespace hepworth
} // namespace gaudi

#endif