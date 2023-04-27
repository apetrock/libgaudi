
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

#include "projection_constraint.hpp"

namespace gaudi {
namespace hepworth {

class edge_edge_collision : public projection_constraint {
public:
  typedef std::shared_ptr<edge_edge_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l) {
    return std::make_shared<edge_edge_collision>(ids, w, l);
  }

  edge_edge_collision(const std::vector<index_t> &ids, const real &w,
                      const real &l)
      : projection_constraint(ids, w), _l(l) {}

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
    real dl = (_l - l);
    xAB /= l;

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;
    if (dl > 0) {
      // dl = max(dl, 0.25 * _l);
      gg::geometry_logger::line(xA, xB, vec4(1.0, 1.0, 0.0, 1.0));
      gg::geometry_logger::line(xA0, xA1, vec4(0.0, 0.8, 0.8, 1.0));
      gg::geometry_logger::line(xA0 - 0.5 * a0 * dl * xAB,
                                xA1 - 0.5 * a1 * dl * xAB,
                                vec4(0.0, 1.0, 1.0, 1.0));

      gg::geometry_logger::line(xB0, xB1, vec4(0.8, 0.0, 0.8, 1.0));
      gg::geometry_logger::line(xB0 + 0.5 * b0 * dl * xAB,
                                xB1 + 0.5 * b1 * dl * xAB,
                                vec4(1.0, 0.0, 1.0, 1.0));
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
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[1];

    index_t iB0 = this->_ids[2];
    index_t iB1 = this->_ids[3];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, 3 * iA0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, 3 * iA1 + ax, _w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, 3 * iB0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, 3 * iB1 + ax, _w));
    id0 += 12;
  }

  real _l = 1.0;
};

class edge_edge_normal_collision : public projection_constraint {
public:
  typedef std::shared_ptr<edge_edge_normal_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &eps,
                    const real &w) {
    return std::make_shared<edge_edge_normal_collision>(ids, eps, w);
  }

  edge_edge_normal_collision(const std::vector<index_t> &ids, const real &eps,
                             const real &w)
      : projection_constraint(ids, w), _eps(eps) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i00 = this->_ids[0];
    index_t i00p = this->_ids[1];
    index_t i01 = this->_ids[2];
    index_t i01p = this->_ids[3];

    index_t i10 = this->_ids[4];
    index_t i10p = this->_ids[5];
    index_t i11 = this->_ids[6];
    index_t i11p = this->_ids[7];
    vec3 x00 = from3(i00, q);
    vec3 x00p = from3(i00p, q);
    vec3 x01 = from3(i01, q);
    vec3 x01p = from3(i01p, q);

    vec3 x10 = from3(i10, q);
    vec3 x10p = from3(i10p, q);
    vec3 x11 = from3(i11, q);
    vec3 x11p = from3(i11p, q);

    std::array<real, 3> d = va::distance_Segment_Segment(x00, x01, x10, x11);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, x00, x01);
    vec3 xB = va::mix(t, x10, x11);
    vec3 dA = (x01 - x00).normalized();
    vec3 dB = (x11 - x10).normalized();

    vec3 xAB = xA - xB;

    real l = xAB.norm();

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;
    vec3 N10 = (x10p - x10).cross(x11 - x10).normalized();
    vec3 N11 = (x11p - x11).cross(x10 - x11).normalized();
    vec3 N1 = (N10 + N11).normalized();

    real h = N1.dot(xAB) - _eps;
    if (h < 0) {
      real dp = 1e-7;
      real db = pow(l - dp, 2) * log(l / dp); // barrier potential...
      // vec3 dX = db * xAB / l;
      vec3 dX = xAB - _eps * xAB / l;

      gg::geometry_logger::line(xA, xB, vec4(1.0, 0.0, 0.0, 1.0));
      gg::geometry_logger::line(xA, xA - dX, vec4(0.5, 0.2, 0.0, 1.0));
      gg::geometry_logger::line(xB, xB + dX, vec4(0.5, 0.0, 0.2, 1.0));

      gg::geometry_logger::line(x00, x01, vec4(0.0, 1.0, 0.5, 1.0));
      gg::geometry_logger::line(x10, x11, vec4(1.0, 0.5, 1.0, 1.0));

      p.block(_id0 + 0, 0, 3, 1) = _w * (x00 - a0 * dX);
      p.block(_id0 + 3, 0, 3, 1) = _w * (x01 - a1 * dX);
      p.block(_id0 + 6, 0, 3, 1) = _w * (x10 + b0 * dX);
      p.block(_id0 + 9, 0, 3, 1) = _w * (x11 + b1 * dX);
    } else {

      gg::geometry_logger::line(xA, xB, vec4(0.0, 1.0, 0.0, 1.0));
      p.block(_id0 + 0, 0, 3, 1) = _w * x00;
      p.block(_id0 + 3, 0, 3, 1) = _w * x01;
      p.block(_id0 + 6, 0, 3, 1) = _w * x10;
      p.block(_id0 + 9, 0, 3, 1) = _w * x11;
    }
    //    p.block(3 * j, 0, 3, 1) = -_w * dq / l;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[2];

    index_t iB0 = this->_ids[4];
    index_t iB1 = this->_ids[6];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, 3 * iA0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, 3 * iA1 + ax, _w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, 3 * iB0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, 3 * iB1 + ax, _w));
    id0 += 12;
  }
  real _eps;
};

class pnt_tri_collision : public projection_constraint {
public:
  typedef std::shared_ptr<pnt_tri_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w) {
    return std::make_shared<pnt_tri_collision>(ids, Nv, eps, w);
  }

  pnt_tri_collision(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w)
      : projection_constraint(ids, w), _Nv(Nv), _eps(eps) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t iP = this->_ids[0];

    index_t iT0 = this->_ids[1];
    index_t iT1 = this->_ids[2];
    index_t iT2 = this->_ids[3];
    // std::cout << iP << " " << iT0 << " " << iT1 << " " << iT2 << std::endl;
    vec3 xP = from3(iP, q);
    vec3 xT0 = from3(iT0, q);
    vec3 xT1 = from3(iT1, q);
    vec3 xT2 = from3(iT2, q);
    vec3 c = 0.3333 * (xT0 + xT1 + xT2);
    // vec3 xN;
    // real d0 = va::distance_from_triangle({xT0, xT1, xT2}, xP, xN);
    std::array<real, 4> cp = va::closest_point({xT0, xT1, xT2}, xP);
    real u = cp[1];
    real v = cp[2];
    real w = cp[3];

    vec3 xN = u * xT0 + v * xT1 + w * xT2;

    vec3 N = (xT1 - xT0).cross(xT2 - xT0).normalized();
    vec3 dx = (xP - xN);
    real d = dx.norm();
    real h = -va::sgn(_Nv.dot(dx)) * N.dot(dx);
    real db = pow(h - _eps, 2) * log(h / _eps); // barrier potential...

    if (h < _eps + db) {
      db = d + _eps;
      real dh = (_eps - h) + db;
      // gg::geometry_logger::line(xP, xP + db * N, vec4(1.0, 0.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + dx, vec4(1.0, 0.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + _eps * N, vec4(1.0, 1.0, 0.0, 1.0));
      // gg::geometry_logger::line(xT0, xT1, vec4(0.75, 0.0, 0.25, 1.0));
      // gg::geometry_logger::line(xT1, xT2, vec4(0.75, 0.0, 0.25, 1.0));
      // gg::geometry_logger::line(xT2, xT0, vec4(0.75, 0.0, 0.25, 1.0));
      vec3 Nt = dx.normalized();
      p.block(_id0 + 0, 0, 3, 1) = _w * (xN + dh * Nt);
      p.block(_id0 + 3, 0, 3, 1) = _w * (xT0 - u * dh * Nt);
      p.block(_id0 + 6, 0, 3, 1) = _w * (xT1 - v * dh * Nt);
      p.block(_id0 + 9, 0, 3, 1) = _w * (xT2 - w * dh * Nt);

    } else {
      // gg::geometry_logger::line(xN, xN + dx, vec4(0.0, 1.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + _eps * N, vec4(1.0, 1.0, 0.0, 1.0));

      p.block(_id0 + 0, 0, 3, 1) = _w * xP;
      p.block(_id0 + 3, 0, 3, 1) = _w * xT0;
      p.block(_id0 + 6, 0, 3, 1) = _w * xT1;
      p.block(_id0 + 9, 0, 3, 1) = _w * xT2;
    }
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t iP = this->_ids[0];

    index_t iT0 = this->_ids[1];
    index_t iT1 = this->_ids[2];
    index_t iT2 = this->_ids[3];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, 3 * iP + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, 3 * iT0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, 3 * iT1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, 3 * iT2 + ax, _w));
    id0 += 12;
  }
  vec3 _Nv;
  real _eps;
};

} // namespace hepworth
} // namespace gaudi
#endif