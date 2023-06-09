
#ifndef __HEP_ROD_BLOCK_COLLISION__
#define __HEP_ROD_BLOCK_COLLISION__

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
#include "sim_block.hpp"
namespace gaudi {
namespace hepworth {
namespace block {

class rod_collision : public block_constraint {
public:
  typedef std::shared_ptr<rod_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<rod_collision>(ids, w, l, blocks);
  }

  rod_collision(const std::vector<index_t> &ids, const real &w,
                      const real &l, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _l(l) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[1];

    index_t iB0 = this->_ids[2];
    index_t iB1 = this->_ids[3];

    vec3 xA0 = _blocks[0]->get_vec3(iA0, q);
    vec3 xA1 = _blocks[0]->get_vec3(iA1, q);
    vec3 xB0 = _blocks[1]->get_vec3(iB0, q);
    vec3 xB1 = _blocks[1]->get_vec3(iB1, q);

    std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, xA0, xA1);
    vec3 xB = va::mix(t, xB0, xB1);
    vec3 dA = (xA1 - xA0).normalized();
    vec3 dB = (xB1 - xB0).normalized();

    vec3 xAB = xB - xA;

    real l = xAB.norm();
    real dl = (_l - l);
    xAB /= l;

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;
    if (dl > 0) {
      vec3 dX = 1.0 * dl * xAB;
      p.block(_id0 + 0, 0, 3, 1) = _w * (xA0 - a0 * dX);
      p.block(_id0 + 3, 0, 3, 1) = _w * (xA1 - a1 * dX);
      p.block(_id0 + 6, 0, 3, 1) = _w * (xB0 + b0 * dX);
      p.block(_id0 + 9, 0, 3, 1) = _w * (xB1 + b1 * dX);
    } else {
      p.block(_id0 + 0, 0, 3, 1) = _w * xA0;
      p.block(_id0 + 3, 0, 3, 1) = _w * xA1;
      p.block(_id0 + 6, 0, 3, 1) = _w * xB0;
      p.block(_id0 + 9, 0, 3, 1) = _w * xB1;
    }
    //    p.block(3 * j, 0, 3, 1) = -_w * dq / l;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
 
    index_t iA0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t iA1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t iB0 = _blocks[1]->get_offset_idx(this->_ids[2]);
    index_t iB1 = _blocks[1]->get_offset_idx(this->_ids[3]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, iA0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, iA1 + ax, _w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, iB0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, iB1 + ax, _w));
    id0 += 12;
  }

  real _l = 1.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif