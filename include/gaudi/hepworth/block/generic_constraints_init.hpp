#ifndef __HEP_GENERIC_BLOCK_CONSTRAINTS_INIT__
#define __HEP_GENERIC_BLOCK_CONSTRAINTS_INIT__

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

void init_pinned(const asawa::shell::shell &shell,
                 std::vector<projection_constraint::ptr> &constraints,
                 const std::vector<vec3> &x, const real &w,
                 std::vector<sim_block::ptr> blocks) {
  for (int iv = 0; iv < shell.vert_count(); iv++) {
    constraints.push_back(
        pinned::create(std::vector<index_t>({iv}), x[iv], w, blocks));
  }
}

void init_pinned(const asawa::rod::rod &R,
                 std::vector<projection_constraint::ptr> &constraints,
                 const std::vector<vec3> &x, const real &w,
                 std::vector<sim_block::ptr> blocks) {
  for (int iv = 0; iv < R.corner_count(); iv++) {
    constraints.push_back(
        pinned::create(std::vector<index_t>({iv}), x[iv], w, blocks));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif