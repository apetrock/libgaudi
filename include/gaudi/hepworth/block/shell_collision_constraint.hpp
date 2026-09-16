
#ifndef __HEP_SHELL_BLOCK_COLLISION__
#define __HEP_SHELL_BLOCK_COLLISION__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "block_constraint.hpp"
#include "block_geometry.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/geometry_logger.hpp"
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
    // Position-space targets: Aq = w*[x00,x01,x10,x11]; inactive p matches,
    // active p separates along the closest-approach axis to gap >= _eps.
    index_t i00 = this->_ids[0];
    index_t i00p = this->_ids[1];
    index_t i01 = this->_ids[2];
    index_t i01p = this->_ids[3];
    index_t i10 = this->_ids[4];
    index_t i10p = this->_ids[5];
    index_t i11 = this->_ids[6];
    index_t i11p = this->_ids[7];

    vec3 x00 = _blocks[0]->get_vec3(i00, q);
    vec3 x00p = _blocks[0]->get_vec3(i00p, q);
    vec3 x01 = _blocks[0]->get_vec3(i01, q);
    vec3 x01p = _blocks[0]->get_vec3(i01p, q);
    vec3 x10 = _blocks[0]->get_vec3(i10, q);
    vec3 x10p = _blocks[0]->get_vec3(i10p, q);
    vec3 x11 = _blocks[0]->get_vec3(i11, q);
    vec3 x11p = _blocks[0]->get_vec3(i11p, q);

    std::array<real, 3> d = va::distance_Segment_Segment(x00, x01, x10, x11);
    real s = d[1];
    real t = d[2];

    vec3 xA = va::mix(s, x00, x01);
    vec3 xB = va::mix(t, x10, x11);
    vec3 xAB = xA - xB;
    real l = xAB.norm();

    real a0 = 1.0 - s;
    real a1 = s;
    real b0 = 1.0 - t;
    real b1 = t;

    vec3 N10 = (x10p - x10).cross(x11 - x10);
    vec3 N11 = (x11p - x11).cross(x10 - x11);
    vec3 N1 = N10 + N11;
    if (N1.squaredNorm() > 1e-24)
      N1.normalize();
    else if (l > 1e-12)
      N1 = xAB / l;
    else
      N1 = vec3::UnitZ();

    const real h = N1.dot(xAB) - _eps;
    if (h < 0.0 && l > 1e-12) {
      const vec3 dX = xAB - _eps * (xAB / l);
      geometry_logger::line(xA, xB, vec4(1.0, 0.0, 0.0, 1.0));
      geometry_logger::line(xA, xA - dX, vec4(0.5, 0.2, 0.0, 1.0));
      geometry_logger::line(xB, xB + dX, vec4(0.5, 0.0, 0.2, 1.0));

      p.block(_id0 + 0, 0, 3, 1) = _w * (x00 - a0 * dX);
      p.block(_id0 + 3, 0, 3, 1) = _w * (x01 - a1 * dX);
      p.block(_id0 + 6, 0, 3, 1) = _w * (x10 + b0 * dX);
      p.block(_id0 + 9, 0, 3, 1) = _w * (x11 + b1 * dX);
    } else {
      p.block(_id0 + 0, 0, 3, 1) = _w * x00;
      p.block(_id0 + 3, 0, 3, 1) = _w * x01;
      p.block(_id0 + 6, 0, 3, 1) = _w * x10;
      p.block(_id0 + 9, 0, 3, 1) = _w * x11;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t iA0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t iA1 = _blocks[0]->get_offset_idx(this->_ids[2]);
    index_t iB0 = _blocks[0]->get_offset_idx(this->_ids[4]);
    index_t iB1 = _blocks[0]->get_offset_idx(this->_ids[6]);
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
  real _eps;
};

class pnt_tri_collision : public block_constraint {
public:
  typedef std::shared_ptr<pnt_tri_collision> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &eps,
                    const real &w, std::vector<sim_block::ptr> blocks,
                    real hardness = 1.0) {
    return std::make_shared<pnt_tri_collision>(ids, eps, w, blocks, hardness);
  }

  static ptr create(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w,
                    std::vector<sim_block::ptr> blocks, real hardness = 1.0) {
    return std::make_shared<pnt_tri_collision>(ids, Nv, eps, w, blocks,
                                              hardness);
  }

  pnt_tri_collision(const std::vector<index_t> &ids, const real &eps,
                    const real &w, std::vector<sim_block::ptr> blocks,
                    real hardness = 1.0)
      : block_constraint(ids, w, blocks), _eps(eps),
        _hardness(std::max(hardness, real(0.0))) {}

  pnt_tri_collision(const std::vector<index_t> &ids, const vec3 &Nv,
                    const real &eps, const real &w,
                    std::vector<sim_block::ptr> blocks, real hardness = 1.0)
      : block_constraint(ids, w, blocks), _Nv(Nv), _eps(eps),
        _hardness(std::max(hardness, real(0.0))) {
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    if (_rest_initialized || !ctx.shell_mesh)
      return;
    const std::vector<vec3> x =
        gather_shell_positions(*ctx.shell_mesh, _blocks[0], q);
    _Nv = asawa::shell::vert_normal(*ctx.shell_mesh,
                                    asawa::shell::vert_id(_ids[0]), x);
    _rest_initialized = true;
  }

  virtual void project(const vecX &q, vecX &p) {
    index_t iP = this->_ids[0];
    index_t iT0 = this->_ids[1];
    index_t iT1 = this->_ids[2];
    index_t iT2 = this->_ids[3];

    vec3 xP = _blocks[0]->get_vec3(iP, q);
    vec3 xT0 = _blocks[0]->get_vec3(iT0, q);
    vec3 xT1 = _blocks[0]->get_vec3(iT1, q);
    vec3 xT2 = _blocks[0]->get_vec3(iT2, q);

    std::array<real, 4> cp = va::closest_point({xT0, xT1, xT2}, xP);
    real u = cp[1];
    real v = cp[2];
    real w = cp[3];
    vec3 xN = u * xT0 + v * xT1 + w * xT2;

    vec3 N = (xT1 - xT0).cross(xT2 - xT0);
    if (N.squaredNorm() < 1e-24) {
      p.block(_id0 + 0, 0, 3, 1) = _w * xP;
      p.block(_id0 + 3, 0, 3, 1) = _w * xT0;
      p.block(_id0 + 6, 0, 3, 1) = _w * xT1;
      p.block(_id0 + 9, 0, 3, 1) = _w * xT2;
      return;
    }
    N.normalize();

    // Signed height along the triangle normal. Positive = in front.
    // One-sided: if the point is inside the barrier band *or behind the
    // plane*, always recover along +N. Using abs(s)*sgn(s) (the old path)
    // is two-sided and *increases* penetration once s < 0.
    const real s = N.dot(xP - xN);
    if (s >= _eps) {
      p.block(_id0 + 0, 0, 3, 1) = _w * xP;
      p.block(_id0 + 3, 0, 3, 1) = _w * xT0;
      p.block(_id0 + 6, 0, 3, 1) = _w * xT1;
      p.block(_id0 + 9, 0, 3, 1) = _w * xT2;
      return;
    }

    // IPC-style barrier boost on (0, ε). _hardness amplifies that boost and
    // deep-penetration pull-out so contact can win against soft materials.
    // dh is half of the relative gap target because point and triangle both
    // move (relative Δgap ≈ 2 dh).
    real barrier = 0.0;
    if (s > 1e-12 && _eps > 1e-12) {
      barrier = -std::pow(s - _eps, 2.0) * std::log(s / _eps);
      barrier = va::clamp(barrier, 0.0, _hardness * _eps);
    } else if (_eps > 1e-12) {
      // On/behind the plane: restore to ε plus penetration depth.
      barrier = _eps + std::max(real(0.0), -s);
    }

    const real dh = 0.5 * ((_eps - s) + _hardness * barrier);
    geometry_logger::line(xP, xP + dh * N, vec4(1.0, 0.0, 0.0, 1.0));
    p.block(_id0 + 0, 0, 3, 1) = _w * (xP + dh * N);
    p.block(_id0 + 3, 0, 3, 1) = _w * (xT0 - u * dh * N);
    p.block(_id0 + 6, 0, 3, 1) = _w * (xT1 - v * dh * N);
    p.block(_id0 + 9, 0, 3, 1) = _w * (xT2 - w * dh * N);
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t iP = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t iT0 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t iT1 = _blocks[0]->get_offset_idx(this->_ids[2]);
    index_t iT2 = _blocks[0]->get_offset_idx(this->_ids[3]);
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, iP + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, iT0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, iT1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 9 + ax, iT2 + ax, _w));
    id0 += 12;
  }
  vec3 _Nv;
  real _eps;
  /// Scales barrier / deep-pen separation beyond the geometric gap restore.
  real _hardness = 1.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif