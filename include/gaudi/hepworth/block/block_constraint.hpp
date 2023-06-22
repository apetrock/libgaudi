
#ifndef __HEP_BLOCK_CONSTRAINTS__
#define __HEP_BLOCK_CONSTRAINTS__

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

#include "../projection_constraint.hpp"
#include "sim_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

class block_constraint : public projection_constraint {
public:
  typedef std::shared_ptr<block_constraint> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<block_constraint>(ids, w, blocks);
  }

  block_constraint(const std::vector<index_t> &ids, const real &w,
                   std::vector<sim_block::ptr> blocks)
      : _blocks(blocks), projection_constraint(ids, w) {}

  std::vector<sim_block::ptr> _blocks;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif