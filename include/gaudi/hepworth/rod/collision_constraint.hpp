
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
                    const real &l, const real &s, const real &t) {
    return std::make_shared<collision>(ids, w, l, s, t);
  }

  collision(const std::vector<index_t> &ids, const real &w, const real &l,
            const real &s, const real &t)
      : projection_constraint(ids, w), _l(l), _s(s), _t(t) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[1];

    index_t iB0 = this->_ids[2];
    index_t iB1 = this->_ids[3];

    vec3 xA0 = q.block(3 * iA0, 0, 3, 1);
    vec3 xA1 = q.block(3 * iA1, 0, 3, 1);
    vec3 xB0 = q.block(3 * iB0, 0, 3, 1);
    vec3 xB1 = q.block(3 * iB1, 0, 3, 1);

    std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, xA0, xA1);
    vec3 xB = va::mix(t, xB0, xB1);
    vec3 dA = (xA1 - xA0).normalized();
    vec3 dB = (xB1 - xB0).normalized();

    vec3 xAB = xB - xA;

    real l = xAB.norm();
    real dl = (_l - l) / l;

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;

    // gg::geometry_logger::line(xA, xB, vec4(1.0, 0.0, 1.0, 1.0));
    // if (l < _l)
    /*
    gg::geometry_logger::line(xA0, xA0 - a0 * dl * xAB,
                              vec4(1.0, 0.0, 0.0, 1.0));
    gg::geometry_logger::line(xB0, xB0 + b0 * dl * xAB,
                              vec4(0.0, 1.0, 0.0, 1.0));
    gg::geometry_logger::line(xA1, xA1 - a1 * dl * xAB,
                              vec4(1.0, 0.0, 1.0, 1.0));
    gg::geometry_logger::line(xB1, xB1 + b1 * dl * xAB,
                              vec4(0.0, 1.0, 1.0, 1.0));
*/
    p.block(3 * iA0, 0, 3, 1) += _w * (xA0 - 0.5 * a0 * dl * xAB);
    p.block(3 * iA1, 0, 3, 1) += _w * (xA1 - 0.5 * a1 * dl * xAB);
    p.block(3 * iB0, 0, 3, 1) += _w * (xB0 + 0.5 * b0 * dl * xAB);
    p.block(3 * iB1, 0, 3, 1) += _w * (xB1 + 0.5 * b1 * dl * xAB);

    //    p.block(3 * j, 0, 3, 1) = -_w * dq / l;
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[1];

    index_t iB0 = this->_ids[2];
    index_t iB1 = this->_ids[3];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * iA0 + ax, 3 * iA0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * iA1 + ax, 3 * iA1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * iB0 + ax, 3 * iB0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * iB1 + ax, 3 * iB1 + ax, _w));
  }

  real _t = 0.0;
  real _s = 0.0;
  real _l = 1.0;
};

void init_collisions(asawa::rod::rod &rod, asawa::rod::dynamic &dynamic,
                     std::vector<projection_constraint::ptr> &constraints,
                     const real &w) {
  const std::vector<vec3> &x = rod.__x;
  vector<std::array<index_t, 4>> collisions = dynamic.get_collisions();
  for (auto &c : collisions) {
    if (c[0] > 0) {
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

      gg::geometry_logger::line(xA, xB, vec4(1.0, 0.0, 1.0, 1.0));
      constraints.push_back(hepworth::rod::collision::create(
          {c[0], c[1], c[2], c[3]}, w, rod._r, s, t));
    }
  }
}

} // namespace rod
} // namespace hepworth
} // namespace gaudi
#endif