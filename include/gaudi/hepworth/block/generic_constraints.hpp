
#ifndef __HEP_GENERIC_BLOCK_CONSTRAINTS__
#define __HEP_GENERIC_BLOCK_CONSTRAINTS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <functional>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "block_constraint.hpp"
#include "gaudi/common.h"
#include "shell_constraints.hpp"
#include "sim_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class pinned : public block_constraint {
public:
  DEFINE_CREATE_FUNC(pinned)

  pinned(const std::vector<index_t> &ids, const vec3 &p, const real &w,
         std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _p(p) {}
  virtual std::string name() { return typeid(*this).name(); }
  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    vec3 q0 = _blocks[0]->get_vec3(i, q);
    p.block(_id0, 0, 3, 1) = _w * _p;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, _w));
    id0 += 3;
  }
  vec3 _p;
};

class pnt_pnt_weld : public block_constraint {
public:
  DEFINE_CREATE_FUNC(pnt_pnt_weld)

  pnt_pnt_weld(const std::vector<index_t> &ids, const real &w,
               std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}
  virtual std::string name() { return typeid(*this).name(); }
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

vec3 slide(vec3 xA0, vec3 xA1, vec3 xB, real t0, real t1) {
  vec3 xA_t0 = va::mix(t0, xA0, xA1);
  vec3 xA_t1 = va::mix(t1, xA0, xA1);
  vec3 dXA = xA_t1 - xA_t0;
  // gg::geometry_logger::line(xB, xB + dXA, vec4(0.0, 1.0, 1.0, 1.0));
  return xB + dXA;
};

vec3 rotate_to(vec3 xA, vec3 xB, vec3 N0, vec3 N1, real t) {
  quat q0 = quat::FromTwoVectors(N0, N1);
  double angle = acos(N0.dot(N1) / (N0.norm() * N1.norm()));
  Eigen::Vector3d axis = N0.cross(N1).normalized();
  angle *= t;
  Eigen::AngleAxisd R(angle, axis);
  // gg::geometry_logger::line(xB, xB + dXA, vec4(0.0, 1.0, 1.0, 1.0));
  return xA + R * (xB - xA);
};

class edge_edge_weld : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_edge_weld)

  edge_edge_weld(const std::vector<index_t> &ids, const real &w0,
                 const real &w1, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, 0.0, blocks), _w0(w0), _w1(w1) {}
  virtual std::string name() { return typeid(*this).name(); }
  virtual void project(const vecX &q, vecX &p) {
    index_t iA0 = this->_ids[0];
    index_t iA1 = this->_ids[1];

    index_t iB0 = this->_ids[2];
    index_t iB1 = this->_ids[3];

    vec3 xA0 = _blocks[0]->get_vec3(iA0, q);
    vec3 xA1 = _blocks[0]->get_vec3(iA1, q);
    vec3 xB0 = _blocks[1]->get_vec3(iB0, q);
    vec3 xB1 = _blocks[1]->get_vec3(iB1, q);

    std::array<real, 3> d =
        va::distance_Segment_Segment(xA0, xA1, xB0, xB1, false, true);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, xA0, xA1);
    vec3 xB = va::mix(t, xB0, xB1);
    vec3 dA = (xA1 - xA0).normalized();
    vec3 dB = (xB1 - xB0).normalized();

    vec3 xAB = xB - xA;

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;
    vec3 dX = xAB;
    if (xAB.hasNaN() || std::isnan(a1) || std::isnan(b1)) {
      std::cout << "NAN: " << __PRETTY_FUNCTION__ << std::endl;
      exit(0);
    }
    if (_slide) {
      vec3 xB0t = xB0;
      vec3 xB1t = xB1;

      xB0 = slide(xA0, xA1, xB0, s, _t1);
      xB1 = slide(xA0, xA1, xB1, s, _t1);
    }

    if (_rotate_to) {
      xA0 = rotate_to(xAB, xA0, _N_rot_0, _N_rot_1, _t_rot);
      xB0 = rotate_to(xAB, xB0, _N_rot_0, _N_rot_1, _t_rot);
    }
    // gg::geometry_logger::line(xA, xB, vec4(0.0, 1.0, 1.0, 1.0));
// #define WELD_4
#ifdef WELD_4
    p.block(_id0 + 0, 0, 3, 1) = _w0 * (xA0 + 0.0 * a0 * dX);
    p.block(_id0 + 3, 0, 3, 1) = _w0 * (xA1 + 0.0 * a1 * dX);
    p.block(_id0 + 6, 0, 3, 1) = _w1 * (xB0 - 1.0 * b0 * dX);
    p.block(_id0 + 9, 0, 3, 1) = _w1 * (xB1 - 1.0 * b1 * dX);
#else
    p.block(_id0 + 0, 0, 3, 1) = _w0 * (xA1 - xA0);
    p.block(_id0 + 3, 0, 3, 1) = _w1 * (xB0 - 1.0 * b0 * dX - xA0);
    p.block(_id0 + 6, 0, 3, 1) = _w1 * (xB1 - 1.0 * b1 * dX - xA0);
    //    p.block(3 * j, 0, 3, 1) = -_w * dq / l;
#endif
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;

    index_t iA0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t iA1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t iB0 = _blocks[1]->get_offset_idx(this->_ids[2]);
    index_t iB1 = _blocks[1]->get_offset_idx(this->_ids[3]);
#ifdef WELD_4
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, iA0 + ax, _w0));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, iA1 + ax, _w0));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, iB0 + ax, _w1));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, iB1 + ax, _w1));
    id0 += 12;
#else
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, iA0 + ax, -_w0));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, iA1 + ax, _w0));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, iB0 + ax, _w1));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, iA0 + ax, -_w1));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, iB1 + ax, _w1));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, iA0 + ax, -_w1));
    id0 += 9;

#endif
  }

  // maybe update to slide along either edge
  void set_slide(real t1) {
    _t1 = t1;
    _t1 = va::clamp(_t1, -2.0, 2.0);
    _slide = true;
  }

  void set_rotate_to(real t_rot, vec3 N_rot_0, vec3 N_rot_1) {
    _t_rot = t_rot;
    _N_rot_0 = N_rot_0;
    _N_rot_1 = N_rot_1;
    _rotate_to = true;
  }

  bool _rotate_to = false;
  real _t_rot = 0.0;
  vec3 _N_rot_0;
  vec3 _N_rot_1;

  bool _slide = false;
  real _t1;

  real _w0;
  real _w1;
};

real smin_cubic(real a, real b, real k) {
  real h = max(k - abs(a - b), 0.0) / k;
  return min(a, b) - h * h * h * k * (1.0 / 6.0);
}

real smin_log(real a, real b, real k) {
  real res = exp2(-k * a) + exp2(-k * b);
  return -log2(res) / k;
}

real tangent_point_radius(const vec3 &dp, const vec3 &N) {
  real ndp = dp.squaredNorm();
  real nPdp = (N * N.transpose() * dp).norm();
  return 0.5 * ndp / nPdp;
};

class point_edge_creep : public block_constraint {
public:
  DEFINE_CREATE_FUNC(point_edge_creep)

  point_edge_creep(const std::vector<index_t> &ids, //
                   const real &r,                   //
                   const vec3 &Nr0, const vec3 &Nr1,
                   const vec3 &Ns, //
                   const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _r(r), _Nr0(Nr0), _Nr1(Nr1), _Ns(Ns) {
  }
  virtual std::string name() { return typeid(*this).name(); }
  virtual void project(const vecX &q, vecX &p) {
    index_t iA = this->_ids[0];
    index_t iB0 = this->_ids[1];
    index_t iB1 = this->_ids[2];

    vec3 xA = _blocks[0]->get_vec3(iA, q);
    vec3 xB0 = _blocks[1]->get_vec3(iB0, q);
    vec3 xB1 = _blocks[1]->get_vec3(iB1, q);

    vec3 dB = xB1 - xB0;
    real t0 = (xA - xB0).dot(dB) / dB.dot(dB);
    vec3 xAp = va::mix(t0, xB0, xB1);
    vec3 dXA = xA - xAp;

    vec3 Nr = va::mix(t0, _Nr0, _Nr1).normalized();
    real sgn = va::sgn(Nr.dot(dXA));

    real rA = tangent_point_radius(dXA, Nr);

    real s_max = 4.0 * _r;
    real s = dXA.norm() / s_max;
    s = min(s, 1.0);

    real r = max(_r, va::mix(s, 0.0, rA));
    r *= sgn > 0.0 ? 2.0 : 1.0;

    vec3 x_cen = xAp + sgn * Nr * r;

    vec3 dX = xA - x_cen;
    real d_s = dX.norm() - r;

    vec3 x_proj = x_cen + r * dX.normalized();
    // gg::geometry_logger::line(xA, x_proj, vec4(0.0, 1.0, 0.5, 1.0));

    if (_slide) {
      vec3 xp = x_proj;
      x_proj = slide(xB0, xB1, x_proj, t0, _t1);
    }

    if (_rotate_to) {
      x_proj = rotate_to(xAp, x_proj, _N_rot_0, _N_rot_1, _t_rot);
    }

    if (x_proj.hasNaN()) {
      std::cout << "NAN: " << __PRETTY_FUNCTION__ << std::endl;
      exit(0);
    }
    // p.block(_id0 + 0, 0, 3, 1) = _w * x_proj;
    p.block(_id0 + 0, 0, 3, 1) = _w * x_proj;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;

    index_t iA = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t iB0 = _blocks[1]->get_offset_idx(this->_ids[1]);
    index_t iB1 = _blocks[1]->get_offset_idx(this->_ids[2]);
    // return;
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, iA + ax, _w));
    id0 += 3;
  }

  void set_slide(real t1) {
    _t1 = va::clamp(_t1, -2.0, 2.0);
    _slide = true;
  }

  void set_rotate_to(real t_rot, vec3 N_rot_0, vec3 N_rot_1) {
    _t_rot = t_rot;
    _N_rot_0 = N_rot_0;
    _N_rot_1 = N_rot_1;
    _rotate_to = true;
  }

  bool _rotate_to = false;
  real _t_rot = 0.0;
  vec3 _N_rot_0, _N_rot_1;

  bool _slide = false;
  vec3 _Nr0, _Nr1, _Ns;
  real _t1 = 0.0, _r = 0.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif