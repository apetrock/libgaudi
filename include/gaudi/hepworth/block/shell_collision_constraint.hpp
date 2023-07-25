
#ifndef __HEP_SHELL_BLOCK_COLLISION__
#define __HEP_SHELL_BLOCK_COLLISION__

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

class edge_edge_normal_collision : public block_constraint {
public:
  typedef std::shared_ptr<edge_edge_normal_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &eps,
                    const real &w, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<edge_edge_normal_collision>(ids, eps, w, blocks);
  }

  edge_edge_normal_collision(const std::vector<index_t> &ids, const real &eps,
                             const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _eps(eps) {}

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

class pnt_tri_collision : public block_constraint {
public:
  typedef std::shared_ptr<pnt_tri_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<pnt_tri_collision>(ids, Nv, eps, w, blocks);
  }

  pnt_tri_collision(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w,
                    std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _Nv(Nv), _eps(eps) {}

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
    real d_s = N.dot(dx);
    real hst = va::sgn(d_s);
    real h = abs(d_s);
    real db = pow(h - _eps, 2.0) * log(h / _eps); // barrier potential...
    db = va::clamp(db, 0.0, 1.0);
    if (h < _eps) {
      // db = d + _eps;
      real dh = 0.5 * hst * (_eps - h - db);
      // real dh = 0.5 * hst * db;

      gg::geometry_logger::line(xP, xP + dh * N, vec4(1.0, 0.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + dx, vec4(1.0, 0.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + _eps * N, vec4(1.0, 1.0, 0.0, 1.0));
      // gg::geometry_logger::line(xT0, xT1, vec4(0.75, 0.0, 0.25, 1.0));
      // gg::geometry_logger::line(xT1, xT2, vec4(0.75, 0.0, 0.25, 1.0));
      // gg::geometry_logger::line(xT2, xT0, vec4(0.75, 0.0, 0.25, 1.0));
      vec3 Nt = dx.normalized();
      /* //original formulation
      p.block(_id0 + 0, 0, 3, 1) = _w * (xP + dh * N);
      p.block(_id0 + 3, 0, 3, 1) = _w * (xT0 - u * dh * N);
      p.block(_id0 + 6, 0, 3, 1) = _w * (xT1 - v * dh * N);
      p.block(_id0 + 9, 0, 3, 1) = _w * (xT2 - w * dh * N);
      */
      vec3 xPp = xP + dh * N;
      p.block(_id0 + 0, 0, 3, 1) = _w * (xT0 - u * dh * N - xPp);
      p.block(_id0 + 3, 0, 3, 1) = _w * (xT1 - v * dh * N - xPp);
      p.block(_id0 + 6, 0, 3, 1) = _w * (xT2 - w * dh * N - xPp);
    } else {
      // gg::geometry_logger::line(xN, xN + dx, vec4(0.0, 1.0, 0.0, 1.0));
      // gg::geometry_logger::line(xN, xN + _eps * N, vec4(1.0, 1.0, 0.0, 1.0));

      p.block(_id0 + 0, 0, 3, 1) = _w * (xT0 - xP);
      p.block(_id0 + 3, 0, 3, 1) = _w * (xT1 - xP);
      p.block(_id0 + 6, 0, 3, 1) = _w * (xT2 - xP);
    }
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t iP = this->_ids[0];

    index_t iT0 = this->_ids[1];
    index_t iT1 = this->_ids[2];
    index_t iT2 = this->_ids[3];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, 3 * iT0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, 3 * iP + ax, -_w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, 3 * iT1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, 3 * iP + ax, -_w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, 3 * iT2 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, 3 * iP + ax, -_w));

    id0 += 9;
  }
  vec3 _Nv;
  real _eps;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif