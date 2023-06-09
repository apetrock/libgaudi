#ifndef __HEP_WELD_BLOCK_CONSTRAINTS_INIT__
#define __HEP_WELD_BLOCK_CONSTRAINTS_INIT__

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

#include "block_constraint.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

class pnt_pnt_weld : public block_constraint {
public:
  typedef std::shared_ptr<pnt_pnt_weld> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<pnt_pnt_weld>(ids, w, blocks);
  }

  pnt_pnt_weld(const std::vector<index_t> &ids, const real &w,
               std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i0 = this->_ids[0];
    index_t i1 = this->_ids[1];

    vec3 q0 = _blocks[0]->get_vec3(i0, q);
    vec3 q1 = _blocks[1]->get_vec3(i1, q);

    p.block(_id0 + 0, 0, 3, 1) = _w * q1;
    p.block(_id0 + 3, 0, 3, 1) = _w * q0;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[1]->get_offset_idx(this->_ids[1]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, i1 + ax, _w));

    id0 += 6;
  }
  vec3 _p;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif