#ifndef __HEP_ROD_COUPING_BLOCK_CONSTRAINTS_INIT__
#define __HEP_ROD_COUPING_BLOCK_CONSTRAINTS_INIT__

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

void init_coupling_collisions(
    asawa::rod::rod &rod0, asawa::rod::dynamic &dynamic0, asawa::rod::rod &rod1,
    asawa::rod::dynamic &dynamic1,
    std::vector<projection_constraint::ptr> &constraints, const real &w,
    std::vector<sim_block::ptr> blocks) {

  const std::vector<vec3> &x0 = rod0.__x;
  const std::vector<vec3> &x1 = rod1.__x;
  std::vector<index_t> edge_verts_0 = rod0.get_edge_vert_ids();

  vector<std::array<index_t, 4>> collisions =
      dynamic1.get_collisions(edge_verts_0, x0, 0.5 * rod0._r);
  for (auto &c : collisions) {
    if (c[0] > -1) {
      vec3 xA0 = x0[c[0]];
      vec3 xA1 = x0[c[1]];
      vec3 xB0 = x1[c[2]];
      vec3 xB1 = x1[c[3]];
      gg::geometry_logger::line(0.5 * (xA0 + xA1), 0.5 * (xB0 + xB1),
                                vec4(0.0, 1.0, 1.0, 1.0));
      constraints.push_back(rod_collision::create({c[0], c[1], c[2], c[3]}, w,
                                                  1.0 * rod0._r, blocks));
    }
  }
}
} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif