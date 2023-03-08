
#ifndef __HEP_ROD_COLLISION__
#define __HEP_ROD_COLLISION__

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

namespace gaudi {
namespace hepworth {
namespace rod {

class collision : public projection_constraint {
public:
  typedef std::shared_ptr<collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l) {
    return std::make_shared<collision>(ids, w, l);
  }

  collision(const std::vector<index_t> &ids, const real &w, const real &l)
      : projection_constraint(ids, w), _l(l) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    vec3 q0 = q.block(3 * i, 0, 3, 1);
    vec3 q1 = q.block(3 * j, 0, 3, 1);

    vec3 dq = q1 - q0;
    real l = dq.norm();

    // gg::geometry_logger::line(q0, q0 + _w * dq, vec4(1.0, 0.0, 1.0, 1.0));

    // gg::geometry_logger::line(q1, q1 - _w * dq, vec4(1.0, 0.0, 1.0, 1.0));

    p.block(3 * i, 0, 3, 1) = _w * dq / l;
    // p.block(3 * j, 0, 3, 1) = -_w * dq / l;
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * i + ax, -_w / _l));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * j + ax, _w / _l));
    // for (int ax = 0; ax < 3; ax++)
    //   triplets.push_back(trip(3 * j + ax, 3 * j + ax, -_w / _l));
    // for (int ax = 0; ax < 3; ax++)
    //   triplets.push_back(trip(3 * j + ax, 3 * i + ax, _w / _l));
  }
  real _l = 1.0;
};

void init_collisions(const asawa::rod &rod,
                     std::vector<projection_constraint::ptr> &constraints,
                     const std::vector<real> &l0, const real &w) {
  int Ni = rod.corner_count();
  for (int i = 0; i < rod.corner_count(); i++) {
    index_t j = rod.next(i);
    constraints.push_back(stretch_shear::create({i, j, i, Ni}, w, l0[i]));
  }
}
} // namespace rod
} // namespace hepworth
} // namespace gaudi
#endif