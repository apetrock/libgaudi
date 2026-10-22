
#ifndef __HEP_SHELL_BLOCK_CONSTRAINTS__
#define __HEP_SHELL_BLOCK_CONSTRAINTS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "block_constraint.hpp"
#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/common.h"
#include "albers_fairing.hpp"
#include "hopf_sweep.hpp"
#include "shape_objectives.hpp"
#include "spine_shapes.hpp"
#include "sim_block.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/hepworth/block/cr_shape.hpp"
#include "gaudi/kusama/weingarten.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

// Stencil weights for one-ring Laplacian-type constraints.
enum class laplacian_stencil { unitary, cotan, angle };

// null  — laplacian constraint, projection target p = 0 (Δx = 0, minimal surface).
// smooth — bending constraint, projection toward rest one-ring gradient (membrane).
enum class laplacian_mode { null, smooth };

// rest — per-triangle area constraint toward rest configuration.
// zero — per-triangle area constraint with p = 0 (minimize face area).
enum class area_mode { rest, zero };

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class edge_strain : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_strain)

  edge_strain(const std::vector<index_t> &ids, const real &w,
              std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}

  edge_strain(const std::vector<index_t> &ids, const std::vector<vec3> &x,
              const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    _l = (x[_ids[1]] - x[_ids[0]]).norm();
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    const vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    _l = (q1 - q0).norm();
    _rest_initialized = true;
  }

  virtual std::string name() { return "edge_strain"; }

  virtual void set_goal_length(const real &l) {
    _limit = false;
    _l = l;
  }

  virtual void project(const vecX &q, vecX &p) {

    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    vec3 q0 = _blocks[0]->get_vec3(i, q);
    vec3 q1 = _blocks[0]->get_vec3(j, q);
    vec3 qc = 0.5 * (q0 + q1);
    vec3 dq = q1 - q0;

    real l = dq.norm();
    l = std::max(l, 1e-4);
    dq /= l;
    real r = _limit ? std::clamp(l / _l, 0.1, 2.0) : 1.0;
    // std::cout << "l " << l << " " << _l << " " << l / _l << std::endl;
    // real r = std::clamp(l / _l, 0.1, 2.0);
    // l = std::clamp(_l / l, 0.1, 4.0);
    // geometry_logger::line(qc - 0.5 * l * dq, qc + 0.5 * l * dq,
    //             vec4(0.0, 0.0, 1.0, 1.0));
    p.block(_id0, 0, 3, 1) = _w * r * dq;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {

    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, -_w / _l));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, j + ax, _w / _l));
    id0 += 3;
  }
  bool _limit = true;
  real _l = -1.0;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class triangle_strain : public block_constraint {
public:
  DEFINE_CREATE_FUNC(triangle_strain)

  triangle_strain(const std::vector<index_t> &ids, real rangeMin, real rangeMax,
                  real w, std::vector<sim_block::ptr> blocks, bool zero = false)
      : block_constraint(ids, w, blocks), _rangeMin(rangeMin),
        _rangeMax(rangeMax), _zero(zero) {
    assert(ids.size() == 3);
  }

  triangle_strain(const std::vector<index_t> &ids, const std::vector<vec3> &x,
                  real rangeMin, real rangeMax, real w,
                  std::vector<sim_block::ptr> blocks, bool zero = false)
      : block_constraint(ids, w, blocks), _rangeMin(rangeMin),
        _rangeMax(rangeMax), _zero(zero) {
    assert(ids.size() == 3);
    bake_rest_from_pos(
        [&](index_t vi) { return x[static_cast<size_t>(vi)]; });
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    bake_rest_from_pos([&](index_t vi) { return _blocks[0]->get_vec3(vi, q); });
    _rest_initialized = true;
  }

  virtual std::string name() { return "triangle_strain"; }

  virtual void project(const vecX &q, vecX &p) {
    mat32 edges, P;
    vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);

    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();

    mat2 F = P.transpose() * edges * _rest;
    Eigen::JacobiSVD<mat2> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    vec2 S = svd.singularValues();
    S(0) = clamp(S(0), _rangeMin, _rangeMax);
    S(1) = clamp(S(1), _rangeMin, _rangeMax);
    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();

    mat32 PF = P * F;

    p.block(_id0 + 0, 0, 3, 1) = _w * _wa * PF.col(0);
    p.block(_id0 + 3, 0, 3, 1) = _w * _wa * PF.col(1);
  }

  ///////////////////////////////////////////////////////////////////////////////
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t i2 = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int i = 0; i < n; ++i) {
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i0 + ax,
                                -_w * _wa * (_rest(0, i) + _rest(1, i))));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i1 + ax, _w * _wa * _rest(0, i)));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i2 + ax, _w * _wa * _rest(1, i)));
    }
    id0 += 3 * n;
  }

  real _rangeMax = 1.0;
  real _rangeMin = 0.0;
  bool _zero = false;

  mat2 _rest;
  /// sqrt(rest area): the row geometry, so `_w` stays the knob.
  real _wa = 1;

private:
  template <typename PosAt>
  void bake_rest_from_pos(const PosAt &pos_at) {
    mat32 edges, P;
    edges.col(0) = pos_at(_ids[1]) - pos_at(_ids[0]);
    edges.col(1) = pos_at(_ids[2]) - pos_at(_ids[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    _rest = (P.transpose() * edges).inverse();
    real A = (P.transpose() * edges).determinant() / 2.0f;
    _wa = std::sqrt(std::abs(A));
  }
};

/// Per-triangle conformal energy (Eigensatz 2010 §3.3, `E_conf`): the
/// deformation gradient is projected onto its closest similarity, so uniform
/// scale and rotation are free and only shear pays.
class triangle_conformal : public triangle_strain {
public:
  DEFINE_CREATE_FUNC(triangle_conformal)

  triangle_conformal(const std::vector<index_t> &ids,
                     const std::vector<vec3> &x, real w,
                     std::vector<sim_block::ptr> blocks)
      : triangle_strain(ids, x, real(0), real(0), w, std::move(blocks)) {}

  std::string name() override { return "triangle_conformal"; }

  void project(const vecX &q, vecX &p) override {
    mat32 edges, P;
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    const vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    const vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);
    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    mat2 F = P.transpose() * edges * _rest;
    Eigen::JacobiSVD<mat2> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const vec2 S = svd.singularValues();
    const real s = real(0.5) * (S(0) + S(1));
    F = s * svd.matrixU() * svd.matrixV().transpose();
    const mat32 PF = P * F;
    if (!PF.allFinite())
      return;
    p.block(_id0 + 0, 0, 3, 1) = _w * _wa * PF.col(0);
    p.block(_id0 + 3, 0, 3, 1) = _w * _wa * PF.col(1);
  }
};
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class anisotropic_triangle_strain : public block_constraint {
public:
  DEFINE_CREATE_FUNC(anisotropic_triangle_strain)

  anisotropic_triangle_strain(const std::vector<index_t> &ids,
                              const vec3 &axial_dir, const vec2 &growth, real w,
                              std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _axial_dir(axial_dir),
        _growth(growth) {
    assert(ids.size() == 3);
  }

  anisotropic_triangle_strain(const std::vector<index_t> &ids,
                              const std::vector<vec3> &x, const vec3 &axial_dir,
                              const vec2 &growth, real w,
                              std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _axial_dir(axial_dir),
        _growth(growth) {
    assert(ids.size() == 3);
    bake_rest_from_pos(
        [&](index_t vi) { return x[static_cast<size_t>(vi)]; });
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    bake_rest_from_pos([&](index_t vi) { return _blocks[0]->get_vec3(vi, q); });
    _rest_initialized = true;
  }

  virtual std::string name() { return "anisotropic_triangle_strain"; }

  // Same Aq/p convention as triangle_strain: both live in F-space.
  // Aq ≈ w * P * F_current, p = w * P * F_target with
  // F_target = R_aniso * diag(sx,sy) * R_aniso^T (identity when sx=sy=1).
  virtual void project(const vecX &q, vecX &p) {
    mat32 P;
    vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);
    if (!triangle_uv_frame(q0, q1, q2, P)) {
      p.block(_id0 + 0, 0, 3, 1).setZero();
      p.block(_id0 + 3, 0, 3, 1).setZero();
      return;
    }

    mat2 F_target = material_stretch(P);
    mat32 PF = P * F_target;
    p.block(_id0 + 0, 0, 3, 1) = _w * _wa * PF.col(0);
    p.block(_id0 + 3, 0, 3, 1) = _w * _wa * PF.col(1);
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t i2 = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int i = 0; i < n; ++i) {
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i0 + ax,
                                -_w * _wa * (_rest(0, i) + _rest(1, i))));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i1 + ax, _w * _wa * _rest(0, i)));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i2 + ax, _w * _wa * _rest(1, i)));
    }
    id0 += 3 * n;
  }

  vec3 _axial_dir = vec3::UnitX();
  vec2 _growth = vec2::Ones(); // absolute stretch sx,sy
  mat2 _rest = mat2::Identity(); // (P^T edges)^{-1}, same as triangle_strain
  /// sqrt(rest area): the row geometry, so `_w` stays the knob.
  real _wa = 1;

private:
  static bool triangle_uv_frame(const vec3 &q0, const vec3 &q1, const vec3 &q2,
                                mat32 &P) {
    vec3 e0 = q1 - q0;
    vec3 e1 = q2 - q0;
    if (e0.norm() < 1e-12 || e1.norm() < 1e-12)
      return false;
    P.col(0) = e0.normalized();
    vec3 e1t = e1 - e1.dot(P.col(0)) * P.col(0);
    if (e1t.norm() < 1e-12)
      return false;
    P.col(1) = e1t.normalized();
    return true;
  }

  mat2 material_stretch(const mat32 &P) const {
    const real sx = std::clamp(_growth.x(), real(0.5), real(3.0));
    const real sy = std::clamp(_growth.y(), real(0.5), real(3.0));
    if (std::abs(sx - 1.0) < 1e-12 && std::abs(sy - 1.0) < 1e-12)
      return mat2::Identity();

    vec3 n = P.col(0).cross(P.col(1));
    if (n.squaredNorm() < 1e-24)
      return mat2::Identity();
    n.normalize();
    vec3 axial = _axial_dir - _axial_dir.dot(n) * n;
    if (axial.squaredNorm() < 1e-24)
      axial = P.col(0);
    else
      axial.normalize();
    vec3 hoop = n.cross(axial);
    if (hoop.squaredNorm() < 1e-24)
      hoop = P.col(1);
    else
      hoop.normalize();

    mat2 R_aniso;
    R_aniso.col(0) << axial.dot(P.col(0)), axial.dot(P.col(1));
    R_aniso.col(1) << hoop.dot(P.col(0)), hoop.dot(P.col(1));
    if (std::abs(R_aniso.determinant()) < 1e-12)
      return mat2::Identity();

    mat2 D = mat2::Identity();
    D(0, 0) = sx;
    D(1, 1) = sy;
    return R_aniso * D * R_aniso.transpose();
  }

  template <typename PosAt>
  void bake_rest_from_pos(const PosAt &pos_at) {
    mat32 edges, P;
    const vec3 q0 = pos_at(_ids[0]);
    const vec3 q1 = pos_at(_ids[1]);
    const vec3 q2 = pos_at(_ids[2]);
    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;
    if (!triangle_uv_frame(q0, q1, q2, P)) {
      _rest = mat2::Identity();
      return;
    }
    mat2 G = P.transpose() * edges;
    if (std::abs(G.determinant()) < 1e-12) {
      _rest = mat2::Identity();
      return;
    }
    _rest = G.inverse();
    real A = G.determinant() / 2.0f;
    _wa = std::sqrt(std::abs(A));
  }
};


#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1

class area : public block_constraint {
public:
  DEFINE_CREATE_FUNC(area)

  area(const std::vector<index_t> &ids, real rangeMin, real rangeMax, real w,
       std::vector<sim_block::ptr> blocks, bool zero = false)
      : block_constraint(ids, w, blocks), _rangeMin(rangeMin),
        _rangeMax(rangeMax), _zero(zero) {
    assert(ids.size() == 3);
  }

  area(const std::vector<index_t> &ids, const std::vector<vec3> &x,
       real rangeMin, real rangeMax, real w, std::vector<sim_block::ptr> blocks,
       bool zero = false)
      : block_constraint(ids, w, blocks), _rangeMin(rangeMin),
        _rangeMax(rangeMax), _zero(zero) {
    assert(ids.size() == 3);
    bake_rest_from_pos(
        [&](index_t vi) { return x[static_cast<size_t>(vi)]; });
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    bake_rest_from_pos([&](index_t vi) { return _blocks[0]->get_vec3(vi, q); });
    _rest_initialized = true;
  }

  virtual std::string name() { return "area"; }

  virtual void project(const vecX &q, vecX &p) {
    mat32 edges, P;
    vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);

    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();

    mat2 F = P.transpose() * edges * _rest;
    Eigen::JacobiSVD<mat2> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    vec2 S = svd.singularValues();
    vec2 d(0.0f, 0.0f);
    for (int i = 0; i < 20; ++i) {
      real v = S(0) * S(1);
      real f = v - clamp(v, _rangeMin, _rangeMax);

      vec2 g(S(1), S(0));
      d = -((f - g.dot(d)) / g.dot(g)) * g;
      S = svd.singularValues() + d;
    }

    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    // p.block<3, 2>(0, _id0) = (_w * P * F);
    // p.block(_id0, 0, 3, 2) = (_w * P * F);
    mat32 PF = P * F;

    // zero: minimize triangle area (p = 0). Otherwise project toward rest area.
    if (_zero) {
      p.block(_id0 + 0, 0, 3, 1) = vec3::Zero();
      p.block(_id0 + 3, 0, 3, 1) = vec3::Zero();
    } else {
      p.block(_id0 + 0, 0, 3, 1) = _w * _wa * PF.col(0);
      p.block(_id0 + 3, 0, 3, 1) = _w * _wa * PF.col(1);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t i2 = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int i = 0; i < n; ++i) {
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i0 + ax,
                                -_w * _wa * (_rest(0, i) + _rest(1, i))));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i1 + ax, _w * _wa * _rest(0, i)));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i2 + ax, _w * _wa * _rest(1, i)));
    }
    id0 += 3 * n;
  }
  bool _zero = false;
  real _rangeMax = 1.0;
  real _rangeMin = 0.0;

  mat2 _rest;
  /// sqrt(rest area): the row geometry, so `_w` stays the knob.
  real _wa = 1;

private:
  template <typename PosAt>
  void bake_rest_from_pos(const PosAt &pos_at) {
    mat32 edges, P;
    edges.col(0) = pos_at(_ids[1]) - pos_at(_ids[0]);
    edges.col(1) = pos_at(_ids[2]) - pos_at(_ids[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    _rest = (P.transpose() * edges).inverse();
    real A = (P.transpose() * edges).determinant() / 2.0f;
    _wa = std::sqrt(std::abs(A));
  }
};
#endif

///////////////////////////////////////////////////////////////////////////////
vec3 get_vec3(const std::vector<vec3> x, index_t i) { return x[i]; }
vec3 get_vec3(const vecX q, index_t i) {
  return vec3(q(i + 0), q(i + 1), q(i + 2));
}

typedef std::function<void(const index_t &j, const real &k, const vec3 &dq,
                           const vec3 &N)>
    edge_fcn;

void one_ring_iter(const std::vector<index_t> &ids,
                   const std::vector<real> &weights, //
                   std::vector<vec3> x,              //
                   edge_fcn fcn) {
  index_t i = ids[0];
  vec3 qi0 = x[i];
  for (int k = 1; k < ids.size(); k++) {
    index_t kp = k + 1;
    index_t km = k - 1;
    if (kp == ids.size())
      kp = 1;
    // k == 1 steps onto ids[0], the center. The previous spoke is the last.
    if (km == 0)
      km = static_cast<index_t>(ids.size()) - 1;
    index_t j = ids[k];
    index_t jp = ids[kp];
    index_t jm = ids[km];

    real K = weights[k - 1];
    vec3 qi = x[j];
    vec3 qip = x[jp];
    vec3 qim = x[jm];

    vec3 dqi = qi - qi0;
    vec3 dqip = qip - qi0;
    vec3 Np = dqi.cross(dqip);
    vec3 dqim = qim - qi0;
    vec3 Nm = dqim.cross(dqi);
    vec3 N = Np + Nm;
    N.normalize();
    fcn(j, K, dqi, N);
  }
}

void one_ring_iter(const std::vector<index_t> &ids,
                   const std::vector<real> &weights, //
                   const vecX &q, sim_block &block,  //
                   edge_fcn fcn) {
  index_t i = ids[0];
  vec3 qi0 = block.get_vec3(i, q);

  for (int k = 1; k < ids.size(); k++) {

    index_t kp = k + 1;
    index_t km = k - 1;
    if (kp == ids.size())
      kp = 1;
    if (km == 0)
      km = static_cast<index_t>(ids.size()) - 1;
    index_t j = ids[k];
    index_t jp = ids[kp];
    index_t jm = ids[km];

    real K = weights[k - 1];
    vec3 qi = block.get_vec3(j, q);
    vec3 qip = block.get_vec3(jp, q);
    vec3 qim = block.get_vec3(jm, q);

    vec3 dqi = qi - qi0;
    vec3 dqip = qip - qi0;
    vec3 Np = dqi.cross(dqip);
    vec3 dqim = qim - qi0;
    vec3 Nm = dqim.cross(dqi);
    vec3 N = Np + Nm;
    N.normalize();
    fcn(j, K, dqi, N);
  }
}

void one_ring_fill_A(const index_t &id0, const real &w, //
                     const std::vector<index_t> ids,
                     const std::vector<real> weights,
                     const sim_block &block, //
                     std::vector<trip> &triplets) {
  index_t i = block.get_offset_idx(ids[0]);
  real Km = 0.0;

  for (int k = 1; k < ids.size(); k++) {

    index_t j = block.get_offset_idx(ids[k]);

    real K = weights[k - 1];
    if (std::isnan(K))
      abort();

    // K = 1.0;
    Km += K;
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(id0 + ax, j + ax, w * K));
  }

  for (int ax = 0; ax < 3; ax++)
    triplets.push_back(trip(id0 + ax, i + ax, -w * Km));
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// Bouaziz-style ARAP shortcut (closed form, no SVD): one-ring weighted edge
/// sum vf = Σ k dq. The target is vf rotated onto the rest laplacian, so
/// ||target|| = ||_vg||: p ∝ vf * ||_vg|| / ||vf||. Same signed cotan stencil
/// as the cotangent Laplacian.
class bending : public block_constraint {
public:
  DEFINE_CREATE_FUNC(bending)

  bending(const std::vector<index_t> &ids, const std::vector<real> &edge_weights,
          const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {}

  bending(const std::vector<index_t> &ids,
          const std::vector<real> &edge_weights, const std::vector<vec3> &x,
          const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    _vg = vec3::Zero();
    one_ring_iter(ids, edge_weights, x,
                  [&](const index_t &j, const real &k, const vec3 &dq,
                      const vec3 &N) { _vg += k * dq; });
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    _vg = vec3::Zero();
    one_ring_iter(_ids, _edge_weights, q, *_blocks[0],
                  [&](const index_t &j, const real &k, const vec3 &dq,
                      const vec3 &N) { _vg += k * dq; });
    _rest_initialized = true;
  }

  virtual std::string name() { return "bending"; }

  virtual void project(const vecX &q, vecX &p) {

    vec3 vf = vec3::Zero();
    one_ring_iter(
        _ids, _edge_weights, q, *_blocks[0],
        [&](const index_t &j, const real &k, const vec3 &dq, const vec3 &N) {
          if (dq.hasNaN() || std::isnan(k)) {
            std::cout << __PRETTY_FUNCTION__ << " d2 is nan" << std::endl;
            exit(0);
          }
          vf += k * dq;
        });

    const real vf2 = vf.squaredNorm();
    const real vg2 = _vg.squaredNorm();
    // ||vf|| * ||vg||² / ||vf||² = ||vg||² / ||vf||, which grows as the
    // ring flattens and drives a shallow dent into a pit.
    if (vf2 > 1e-12 && vg2 > 1e-24)
      p.block(_id0, 0, 3, 1) = _w * vf * std::sqrt(vg2 / vf2);
    else
      p.block(_id0, 0, 3, 1) = vec3::Zero();
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    one_ring_fill_A(_id0, _w, _ids, _edge_weights, *_blocks[0], triplets);
    id0 += 3;
  }

  std::vector<real> _edge_weights;
  vec3 _vg;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// Classic per-vertex ARAP (Sorkine–Alexa / Chao spokes): fit R ∈ SO(3) by
/// Procrustes SVD on weighted spoke edges, then project each spoke to R d_rest.
/// Cubic stylization extends the local step with ℓ¹ on R n̂ (ADMM); use this
/// as the base ARAP rows.
class arap : public block_constraint {
public:
  DEFINE_CREATE_FUNC(arap)

  arap(const std::vector<index_t> &ids, const std::vector<real> &edge_weights,
       const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    assert(ids.size() >= 2);
    assert(edge_weights.size() + 1 == ids.size());
    _d_rest.resize(edge_weights.size());
  }

  arap(const std::vector<index_t> &ids, const std::vector<real> &edge_weights,
       const std::vector<vec3> &x, const real &w,
       std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    assert(ids.size() >= 2);
    assert(edge_weights.size() + 1 == ids.size());
    _d_rest.resize(edge_weights.size());
    bake_spokes_from_pos(
        [&](index_t vi) { return x[static_cast<size_t>(vi)]; });
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    bake_spokes_from_pos([&](index_t vi) { return _blocks[0]->get_vec3(vi, q); });
    _rest_initialized = true;
  }

  virtual std::string name() { return "arap"; }

  virtual void project(const vecX &q, vecX &p) {
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const vec3 dq = _blocks[0]->get_vec3(_ids[k + 1], q) - q0;
      M += wk * dq * _d_rest[k].transpose();
    }

    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    mat3 R = U * V.transpose();
    if (R.determinant() < real(0.0)) {
      U.col(2) *= real(-1.0);
      R = U * V.transpose();
    }

    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) =
          _w * wk * (R * _d_rest[k]);
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      const index_t row = id0;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * wk));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * wk));
      id0 += 3;
    }
  }

  std::vector<real> _edge_weights;
  std::vector<vec3> _d_rest;

private:
  template <typename PosAt>
  void bake_spokes_from_pos(const PosAt &pos_at) {
    const vec3 x0 = pos_at(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k)
      _d_rest[k] = pos_at(_ids[k + 1]) - x0;
  }
};

/// Turns the Dirichlet axes by `angle`. Positive moves `u` from T toward B,
/// the same sense as a curve-arap spoke. A non-finite angle is 0.
inline void stretch_frame(const vec3 &T, const vec3 &B, real angle, vec3 &u,
                          vec3 &v) {
  if (!std::isfinite(angle))
    angle = real(0);
  const real c = std::cos(angle);
  const real s = std::sin(angle);
  u = c * T + s * B;
  v = -s * T + c * B;
}

/// The rotation nearest `M` (Procrustes: `argmax_R tr(R^T M)` over SO(3)).
inline mat3 polar_rotation(const mat3 &M) {
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  const mat3 V = svd.matrixV();
  mat3 R = U * V.transpose();
  if (R.determinant() < real(0)) {
    U.col(2) *= real(-1);
    R = U * V.transpose();
  }
  return R;
}

/// Scales the component of `d` along `u` by `1 + alpha`. A zero direction
/// leaves `d` unchanged. Reversing `u` does not change the result.
inline vec3 scale_spoke(const vec3 &d, const vec3 &u, real alpha) {
  if (!std::isfinite(alpha) || !(u.squaredNorm() > real(1e-16)) || !d.allFinite())
    return d;
  const vec3 uh = u.normalized();
  const real along = uh.dot(d);
  if (!std::isfinite(along))
    return d;
  return d + alpha * along * uh;
}

/// Ids of one fan, owned by the row. Local vertex order is global-id order,
/// so the parallel axis sign(v0 < v1) matches the mesh the fan was copied from.
struct fan_local {
  std::vector<int> vids;
  int center = -1;
  std::vector<cr_edge> ledges;
  std::vector<int> dof_of;
  int n_incident = 0;
  std::vector<cr_tri> tris;
  std::vector<int> tangent_e;
};

/// Spoke ARAP. The unscaled spokes are the edges at the start of the step.
/// Each projection scales them on the live minor and major tangents turned
/// by `angle`, or along `dir` (a world tangent, sign free) and across it
/// when `dir` is nonzero. Columns are vertex positions only.
class directed_stretch : public block_constraint {
public:
  DEFINE_CREATE_FUNC(directed_stretch)

  directed_stretch(std::vector<index_t> ids, std::vector<real> edge_weights,
                   const std::vector<vec3> &x, fan_local fan, int nE, real angle,
                   real sx, real sy, real power, real w,
                   std::vector<sim_block::ptr> blocks,
                   const std::vector<vertex_chart> *charts = nullptr,
                   const vec3 &dir = vec3::Zero())
      : block_constraint(std::move(ids), w, std::move(blocks)),
        _edge_weights(std::move(edge_weights)), _fan(std::move(fan)), _nE(nE),
        _angle(std::isfinite(angle) ? angle : real(0)), _sx(sx), _sy(sy),
        _power(std::isfinite(power) ? power : real(0)), _dir(dir.allFinite() ? dir : vec3::Zero()),
        _charts(charts) {
    _d0.resize(_edge_weights.size());
    if (!_ids.empty()) {
      const vec3 x0 = x[static_cast<size_t>(_ids[0])];
      for (size_t k = 0; k < _edge_weights.size(); ++k) {
        const index_t j = _ids[k + 1];
        _d0[k] = x[static_cast<size_t>(j)] - x0;
      }
    }
    _d_rest = _d0;
    _rest_initialized = true;
  }

  std::string name() override { return "directed_stretch"; }

  const std::vector<vec3> &rest() const { return _d_rest; }

  void project(const vecX &q, vecX &p) override;

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      const index_t row = id0;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * wk));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * wk));
      id0 += 3;
    }
  }

  std::vector<real> _edge_weights;
  std::vector<vec3> _d0;
  std::vector<vec3> _d_rest;
  const std::vector<vertex_chart> *_charts = nullptr;
  fan_local _fan;
  int _nE = 0;
  real _angle = 0;
  real _sx = 0;
  real _sy = 0;
  real _power = 0;
  /// World stretch axis held for the step; zero reads `_angle`.
  vec3 _dir = vec3::Zero();
};

/// One sign per vertex for a principal axis field (`tilt_sign::propagate`):
/// breadth-first over `nbr` from the lowest vertex of each component, each
/// vertex flipped to agree with the one it was reached from. Filled on the
/// first `orient` call of the frame, shared by every row of one target.
struct tilt_orientation {
  std::vector<std::vector<int>> nbr;
  std::vector<real> sign;
  principal_axis axis = principal_axis::k1;
  bool ready = false;

  /// World axis of vertex `i` from its chart, before the sign; zero without
  /// a chart.
  static vec3 chart_axis(const std::vector<vertex_chart> *charts, int i, principal_axis axis);
  void orient(const std::vector<vertex_chart> *charts);
};

/// Spoke rows toward the frame-start ring turned about a principal axis or
/// the normal (`tilt_target`). The axis and the turn are read from the
/// vertex charts at the first projection and held for the frame. No live
/// polar, so the ring advances by the turn each frame instead of settling.
/// Columns are vertex positions only.
class directed_tilt : public block_constraint {
public:
  DEFINE_CREATE_FUNC(directed_tilt)

  directed_tilt(std::vector<index_t> ids, std::vector<real> edge_weights,
                const std::vector<vec3> &x, tilt_target spec, real w,
                std::vector<sim_block::ptr> blocks,
                const std::vector<vertex_chart> *charts,
                std::shared_ptr<tilt_orientation> field = nullptr)
      : block_constraint(std::move(ids), w, std::move(blocks)),
        _edge_weights(std::move(edge_weights)), _spec(spec), _charts(charts),
        _field(std::move(field)) {
    _d0.resize(_edge_weights.size());
    if (!_ids.empty())
      for (size_t k = 0; k < _edge_weights.size(); ++k)
        _d0[k] = x[static_cast<size_t>(_ids[k + 1])] - x[static_cast<size_t>(_ids[0])];
    _target = _d0;
  }

  std::string name() override { return "directed_tilt"; }

  void project(const vecX &q, vecX &p) override;

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      for (int ax = 0; ax < 3; ax++) {
        triplets.push_back(trip(id0 + ax, ii + ax, -_w * wk));
        triplets.push_back(trip(id0 + ax, jj + ax, _w * wk));
      }
      id0 += 3;
    }
  }

  /// Turns `_target` from the chart; leaves the frame-start ring without one.
  void hold_turn();

  std::vector<real> _edge_weights;
  std::vector<vec3> _d0;
  std::vector<vec3> _target;
  tilt_target _spec;
  const std::vector<vertex_chart> *_charts = nullptr;
  std::shared_ptr<tilt_orientation> _field;
  bool _held = false;
};

/// One ring, cotan-weighted spoke ARAP whose target fixes one ratio of the
/// ring's 2D stretch and leaves the rest free. In the frame-start chart,
/// with an axis `u` and `u_perp = N x u`, the target spokes are
/// `R (U t_k + n_k)`: `t_k` and `n_k` the tangent and normal parts of the
/// frame-start spoke, `U = [[a, b], [b, c]]` in `(u, u_perp)` with
/// `a = t c`. The skew of `u` (`b`), the rotation `R` and, inside the band
/// `area` of the frame-start ring's area (`det U`), the area
/// are fitted (alternating, each exact; outside the band the fit runs on its
/// edge), so the row resists only a change of that ratio and, past the band,
/// of area. The default band `{1, 1}` keeps the ring's area;
/// `area_band::free()` leaves it free.
/// `absolute`: `t` makes the spoke rms along `u` over that across equal to
/// `ratio` at zero skew (a lock on the ring's covariant anisotropy);
/// otherwise `t = ratio` against the frame-start ring (a drive per frame).
/// The axis is `angle` from the chart's k_hi tangent toward its k_lo
/// tangent. Without a chart the row is ARAP. `label` names the row
/// (`axis_ratio`, `shear`). Columns are vertex positions only.
class axis_ratio : public block_constraint {
public:
  DEFINE_CREATE_FUNC(axis_ratio)

  axis_ratio(std::vector<index_t> ids, std::vector<real> edge_weights,
             const std::vector<vec3> &x, real ratio, bool absolute, real angle, real w,
             std::vector<sim_block::ptr> blocks, const std::vector<vertex_chart> *charts,
             std::string label, area_band area = {})
      : block_constraint(std::move(ids), w, std::move(blocks)),
        _edge_weights(std::move(edge_weights)), _charts(charts), _ratio(ratio),
        _absolute(absolute), _angle(std::isfinite(angle) ? angle : real(0)),
        _area{std::max(area.lo, real(0)), std::max(area.hi, area.lo)},
        _label(std::move(label)) {
    _d0.resize(_edge_weights.size());
    if (!_ids.empty())
      for (size_t k = 0; k < _edge_weights.size(); ++k)
        _d0[k] = x[static_cast<size_t>(_ids[k + 1])] - x[static_cast<size_t>(_ids[0])];
    _d_rest = _d0;
    _rest_initialized = true;
  }

  std::string name() override { return _label; }

  const std::vector<vec3> &rest() const { return _d_rest; }

  void project(const vecX &q, vecX &p) override;

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      for (int ax = 0; ax < 3; ax++) {
        triplets.push_back(trip(id0 + ax, ii + ax, -_w * wk));
        triplets.push_back(trip(id0 + ax, jj + ax, _w * wk));
      }
      id0 += 3;
    }
  }

  std::vector<real> _edge_weights;
  std::vector<vec3> _d0;
  std::vector<vec3> _d_rest;
  const std::vector<vertex_chart> *_charts = nullptr;
  real _ratio = 1;
  bool _absolute = true;
  real _angle = 0;
  area_band _area{};
  std::string _label;
};

/// One triangle held at its frame-start shape mapped by a fixed in-plane
/// stretch `U` (2x2, symmetric, in the face basis `u = (X1 - X0) / |.|`,
/// `v = n x u`), any rotation: the target edges are `R U e0_k`, `R` the
/// polar fit to the live edges, so the projection is exact and the row
/// variational. `U = I` is ARAP on the face. Rows are weighted `sqrt(A /
/// l^2)` (`l^2` the mean squared rest edge), the triangle ARAP energy's
/// scale, so `w` compares with the cotan spoke rows. Columns are vertex positions
/// only. Used by `init_face_shear`.
class face_stretch_target : public block_constraint {
public:
  DEFINE_CREATE_FUNC(face_stretch_target)

  face_stretch_target(const std::array<index_t, 3> &ids, const std::vector<vec3> &x,
                      const Eigen::Matrix<real, 2, 2> &U, real w,
                      std::vector<sim_block::ptr> blocks, std::string label)
      : block_constraint({ids[0], ids[1], ids[2]}, w, std::move(blocks)), _label(std::move(label)) {
    const vec3 X0 = x[static_cast<size_t>(ids[0])];
    const vec3 e1 = x[static_cast<size_t>(ids[1])] - X0, e2 = x[static_cast<size_t>(ids[2])] - X0;
    const vec3 n = e1.cross(e2);
    // Triangle ARAP energy `A |F - R|^2` in edge vectors is `A / l^2 |de|^2`:
    // a dimensionless row scale, as a cotan spoke weight.
    const real l2 = real(0.5) * (e1.squaredNorm() + e2.squaredNorm());
    _s = l2 > real(0) ? std::sqrt(real(0.5) * n.norm() / l2) : real(0);
    const vec3 u = e1.normalized(), v = n.normalized().cross(u);
    for (int k = 0; k < 2; ++k) {
      const vec3 &e = k == 0 ? e1 : e2;
      const Eigen::Matrix<real, 2, 1> c = U * Eigen::Matrix<real, 2, 1>(e.dot(u), e.dot(v));
      _c[static_cast<size_t>(k)] = c[0] * u + c[1] * v;
    }
    _rest_initialized = true;
  }

  std::string name() override { return _label; }

  void project(const vecX &q, vecX &p) override {
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    mat3 M = mat3::Zero();
    for (int k = 0; k < 2; ++k)
      M += (_blocks[0]->get_vec3(_ids[static_cast<size_t>(k + 1)], q) - q0) *
           _c[static_cast<size_t>(k)].transpose();
    const mat3 R = polar_rotation(M);
    for (int k = 0; k < 2; ++k)
      p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = _w * _s * (R * _c[static_cast<size_t>(k)]);
  }

  void fill_A(index_t &id0, std::vector<trip> &triplets) override {
    _id0 = id0;
    const index_t i0 = _blocks[0]->get_offset_idx(_ids[0]);
    for (int k = 0; k < 2; ++k) {
      const index_t ik = _blocks[0]->get_offset_idx(_ids[static_cast<size_t>(k + 1)]);
      for (int ax = 0; ax < 3; ax++) {
        triplets.push_back(trip(id0 + ax, i0 + ax, -_w * _s));
        triplets.push_back(trip(id0 + ax, ik + ax, _w * _s));
      }
      id0 += 3;
    }
  }

  std::array<vec3, 2> _c;
  real _s = 1;
  std::string _label;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// ARAP spokes + normal tip-in for cubic stylization (Liu–Jacobson style).
/// Local Procrustes:
///   M = Σ w_ij d̃ dᵀ + λ_tip N_goal N_restᵀ
/// then p_j = w · w_ij · R · d_rest. Goals are typically axis-quantized once.
class arap_normal_align : public block_constraint {
public:
  DEFINE_CREATE_FUNC(arap_normal_align)

  arap_normal_align(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights, const vec3 &N_rest,
                    const vec3 &N_goal, const real &w, real lambda_tip,
                    std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights),
        _N_rest(N_rest), _N_goal(N_goal), _lambda_tip(lambda_tip) {
    assert(ids.size() >= 2);
    assert(edge_weights.size() + 1 == ids.size());
    _d_rest.resize(edge_weights.size());
    if (_N_rest.squaredNorm() > real(1e-24))
      _N_rest.normalize();
    if (_N_goal.squaredNorm() > real(1e-24))
      _N_goal.normalize();
  }

  arap_normal_align(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights,
                    const std::vector<vec3> &x, const vec3 &N_rest,
                    const vec3 &N_goal, const real &w, real lambda_tip,
                    std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights),
        _N_rest(N_rest), _N_goal(N_goal), _lambda_tip(lambda_tip) {
    assert(ids.size() >= 2);
    assert(edge_weights.size() + 1 == ids.size());
    _d_rest.resize(edge_weights.size());
    const vec3 x0 = x[ids[0]];
    for (size_t k = 0; k < edge_weights.size(); ++k)
      _d_rest[k] = x[ids[k + 1]] - x0;
    if (_N_rest.squaredNorm() > real(1e-24))
      _N_rest.normalize();
    if (_N_goal.squaredNorm() > real(1e-24))
      _N_goal.normalize();
    _rest_initialized = true;
  }

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const vec3 x0 = _blocks[0]->get_vec3(_ids[0], q);
    for (size_t k = 0; k < _edge_weights.size(); ++k)
      _d_rest[k] = _blocks[0]->get_vec3(_ids[k + 1], q) - x0;
    _rest_initialized = true;
  }

  virtual std::string name() { return "arap_normal_align"; }

  virtual void project(const vecX &q, vecX &p) {
    const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const vec3 dq = _blocks[0]->get_vec3(_ids[k + 1], q) - q0;
      M += wk * dq * _d_rest[k].transpose();
    }
    if (_lambda_tip > real(0.0) && _N_rest.squaredNorm() > real(1e-24) &&
        _N_goal.squaredNorm() > real(1e-24)) {
      M += _lambda_tip * _N_goal * _N_rest.transpose();
    }

    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    mat3 R = U * V.transpose();
    if (R.determinant() < real(0.0)) {
      U.col(2) *= real(-1.0);
      R = U * V.transpose();
    }

    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) =
          _w * wk * (R * _d_rest[k]);
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    for (size_t k = 0; k < _edge_weights.size(); ++k) {
      const real wk = _edge_weights[k];
      const index_t jj = _blocks[0]->get_offset_idx(_ids[k + 1]);
      const index_t row = id0;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * wk));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * wk));
      id0 += 3;
    }
  }

  std::vector<real> _edge_weights;
  std::vector<vec3> _d_rest;
  vec3 _N_rest = vec3::UnitZ();
  vec3 _N_goal = vec3::UnitZ();
  real _lambda_tip = 0.0;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class laplacian : public block_constraint {
public:
  DEFINE_CREATE_FUNC(laplacian)

  static ptr create(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<laplacian>(ids, edge_weights, w, blocks);
  }

  static ptr create(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights,
                    const std::vector<vec3> &x, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<laplacian>(ids, edge_weights, x, w, blocks);
  }

  laplacian(const std::vector<index_t> &ids,
            const std::vector<real> &edge_weights, const real &w,
            std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    _rest_initialized = true;
  }

  laplacian(const std::vector<index_t> &ids,
            const std::vector<real> &edge_weights, const std::vector<vec3> &x,
            const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    (void)x;
    _rest_initialized = true;
  }
  virtual std::string name() { return "laplacian"; }

  virtual void project(const vecX &q, vecX &p) {
    p.block(_id0, 0, 3, 1) = vec3::Zero();
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    one_ring_fill_A(_id0, _w, _ids, _edge_weights, *_blocks[0], triplets);
    id0 += 3;
  }
  std::vector<real> _edge_weights;
};

/// `k0/k1` of one shape vector. `k0` is the larger eigenvalue. The gradient is
/// with respect to `(k11, k12, k22)`. A near-zero `k1` leaves the ratio unset.
inline bool principal_curvature_ratio(const vec3 &h, real &ratio, vec3 &grad) {
  const real a = h.x();
  const real b = h.y();
  const real c = h.z();
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c))
    return false;
  const real tr = a + c;
  const real half = real(0.5) * tr;
  const real inside = half * half - (a * c - b * b);
  const real disc = std::sqrt(std::max(inside, real(0)));
  const real k0 = half + disc;
  const real k1 = half - disc;
  if (!(std::abs(k1) > real(1e-4)) || !std::isfinite(k0) || !std::isfinite(k1))
    return false;
  ratio = k0 / k1;
  if (!std::isfinite(ratio) || std::abs(ratio) > real(50))
    return false;
  vec3 ddisc = vec3::Zero();
  if (disc > real(1e-14)) {
    ddisc.x() = (half - c) / (real(2) * disc);
    ddisc.y() = b / disc;
    ddisc.z() = (half - a) / (real(2) * disc);
  }
  const real denom = k1 * k1;
  grad.x() = (-real(2) * disc * real(0.5) + tr * ddisc.x()) / denom;
  grad.y() = tr * ddisc.y() / denom;
  grad.z() = (-real(2) * disc * real(0.5) + tr * ddisc.z()) / denom;
  return grad.allFinite();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// One-ring spokes, projected onto a Weingarten quadric
/// \(h = a s^2 + b s t + c t^2\). The quadric is fit from `_fit_ids`
/// (the two-ring) when those are set; otherwise from the one-ring itself.
/// Mean curvature is the trace, \(H = a + c\). The fit keeps the
/// off-diagonal term as well, so the target is the full shape operator.
class willmore : public block_constraint {
public:
  DEFINE_CREATE_FUNC(willmore)

  willmore(const std::vector<index_t> &ids,
           const std::vector<real> &edge_weights, const real &w,
           std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    _rest_initialized = true;
  }

  willmore(const std::vector<index_t> &ids,
           const std::vector<real> &edge_weights, const std::vector<vec3> &x,
           const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    (void)x;
    _rest_initialized = true;
  }
  virtual std::string name() { return "willmore"; }

  void set_align_normal(vec3 N) {
    _align = true;
    _N_align = N;
  }

  /// Vertex ids used only to fit the quadric. Not constraint rows.
  void set_fit_ids(const std::vector<index_t> &ids) { _fit_ids = ids; }

  virtual void project(const vecX &q, vecX &p) {
    const vec3 qi = _blocks[0]->get_vec3(_ids[0], q);
    const int nring = static_cast<int>(_ids.size()) - 1;
    std::vector<vec3> dq(static_cast<size_t>(std::max(nring, 0)));
    vec3 N = vec3::Zero();
    for (int k = 0; k < nring; ++k) {
      dq[static_cast<size_t>(k)] =
          _blocks[0]->get_vec3(_ids[k + 1], q) - qi;
      const int kp = (k + 1 < nring) ? (k + 1) : 0;
      const vec3 dnext = _blocks[0]->get_vec3(_ids[kp + 1], q) - qi;
      N += dq[static_cast<size_t>(k)].cross(dnext);
    }

    vec3 u = vec3::UnitX();
    vec3 v = vec3::UnitY();
    real a = 0.0;
    real b = 0.0;
    real c = 0.0;
    bool quadric = false;
    if (N.squaredNorm() > 1e-24) {
      N.normalize();
      if (nring >= 3) {
        const vec3 tmp =
            std::abs(N[0]) > 0.5 ? vec3(0, 1, 0) : vec3(1, 0, 0);
        u = N.cross(tmp);
        if (u.squaredNorm() < 1e-24)
          u = vec3(0, 0, 1).cross(N);
        u.normalize();
        v = N.cross(u);

        mat3 ATA = mat3::Zero();
        vec3 Atb = vec3::Zero();
        const bool use_outer = static_cast<int>(_fit_ids.size()) >= 3;
        const int nfit = use_outer ? static_cast<int>(_fit_ids.size()) : nring;
        for (int k = 0; k < nfit; ++k) {
          const vec3 d =
              use_outer ? (_blocks[0]->get_vec3(_fit_ids[static_cast<size_t>(k)],
                                                q) -
                           qi)
                        : dq[static_cast<size_t>(k)];
          const real s = d.dot(u);
          const real t = d.dot(v);
          const real h = d.dot(N);
          const vec3 row(s * s, s * t, t * t);
          ATA += row * row.transpose();
          Atb += h * row;
        }
        Eigen::JacobiSVD<mat3> svd(ATA, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const vec3 coef = svd.solve(Atb);
        a = coef[0];
        b = coef[1];
        c = coef[2];
        quadric = std::isfinite(a) && std::isfinite(b) && std::isfinite(c);
      }
    } else {
      N = vec3::UnitZ();
    }

    quat Ra = quat::Identity();
    if (_align && _N_align.squaredNorm() > 1e-24)
      Ra.setFromTwoVectors(N, _N_align.normalized());

    for (int k = 0; k < nring; ++k) {
      const real wk = _edge_weights[static_cast<size_t>(k)];
      const vec3 &d = dq[static_cast<size_t>(k)];
      const vec3 tng = d - d.dot(N) * N;
      real h = 0.0;
      if (quadric) {
        const real s = tng.dot(u);
        const real tv = tng.dot(v);
        h = a * s * s + b * s * tv + c * tv * tv;
        if (!std::isfinite(h))
          h = 0.0;
      }
      vec3 rdq = tng + h * N;
      if (_align && _N_align.squaredNorm() > 1e-24)
        rdq = Ra.normalized() * rdq;
      p.block(_id0 + 3 * k, 0, 3, 1) = _w * wk * rdq;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t ii = _blocks[0]->get_offset_idx(_ids[0]);

    for (int j = 1; j < _ids.size(); j++) {
      real K = _edge_weights[j - 1];
      index_t jj = _blocks[0]->get_offset_idx(_ids[j]);
      const index_t row = id0;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, ii + ax, -_w * K));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(row + ax, jj + ax, _w * K));
      id0 += 3;
    }
  }

  vec3 _N_align; // Goal N
  bool _align = false;
  std::vector<real> _edge_weights;
  std::vector<index_t> _fit_ids;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class normal : public block_constraint {
public:
  DEFINE_CREATE_FUNC(normal)

  normal(const std::vector<index_t> &ids, const std::vector<vec3> &x, vec3 N0,
         real w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _N0(N0) {}

  virtual void project(const vecX &q, vecX &p) {
    mat32 edges, P;
    vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);
    vec3 qc = 1.0 / 3.0 * (q0 + q1 + q2);
    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;

    // gg::geometry_geometry_logger::line(q0, q0 + edges.col(0), vec4(0.0, 1.0,
    // 0.0, 1.0)); gg::geometry_geometry_logger::line(q0, q0 + edges.col(1),
    // vec4(0.0, 1.0, 0.0, 1.0));

    vec3 N = edges.col(0).cross(edges.col(1)).normalized();
    real d = N.dot(_N0);
    quat qN = quat::FromTwoVectors(N, _N0);
    edges.col(0) = qN * edges.col(0);
    edges.col(1) = qN * edges.col(1);

    // gg::geometry_geometry_logger::line(q0, q0 + edges.col(0), vec4(1.0, 0.0,
    // 0.0, 1.0)); gg::geometry_geometry_logger::line(q0, q0 + edges.col(1), vec4(1.0,
    // 0.0, 0.0, 1.0));

    // gg::geometry_geometry_logger::line(qc, qc + 0.025 * N, vec4(0.0, 1.0, 0.5, 1.0));
    // gg::geometry_geometry_logger::line(qc, qc + 0.025 * _N0, vec4(1.0, 0.0,
    // 0.5, 1.0));

    p.block(_id0 + 0, 0, 3, 1) = _w * edges.col(0);
    p.block(_id0 + 3, 0, 3, 1) = _w * edges.col(1);
  }

  ///////////////////////////////////////////////////////////////////////////////
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t i2 = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i1 + ax, _w));

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, i0 + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, i2 + ax, _w));

    id0 += 6;
  }
  vec3 _N0;
};
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class face_normal_flow : public block_constraint {
//this will be here untested, since I'm ultimately going to do this for edges, instead

public:
  DEFINE_CREATE_FUNC(face_normal_flow)

  face_normal_flow(const std::vector<index_t> &ids, const std::vector<vec3> &x, real d,
         real w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _d(d) {
        _N0 = (x[_ids[1]] - x[_ids[0]]).cross(x[_ids[2]] - x[_ids[0]]).normalized(); 
    _cen = 1.0 / 3.0 * (x[_ids[0]] + x[_ids[1]] + x[_ids[2]]);
      }

  virtual void project(const vecX &q, vecX &p) {
    mat32 edges, P;
    vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 q1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 q2 = _blocks[0]->get_vec3(_ids[2], q);
    vec3 qc = 1.0 / 3.0 * (q0 + q1 + q2);
    edges.col(0) = q1 - q0;
    edges.col(1) = q2 - q0;


    vec3 N = edges.col(0).cross(edges.col(1)).normalized();
    //compute distances from initial plane defined by _N0 and _cen
    real d0 = (q0 - _cen).dot(_N0);
    real d1 = (q1 - _cen).dot(_N0);
    real d2 = (q2 - _cen).dot(_N0);
    
    //compute the flow
    q0 += (_d - d0) * N;
    q1 += (_d - d1) * N;
    q2 += (_d - d2) * N;

    p.block(_id0 + 0, 0, 3, 1) = _w * edges.col(0);
    p.block(_id0 + 3, 0, 3, 1) = _w * edges.col(1);
    p.block(_id0 + 6, 0, 3, 1) = _w * edges.col(1);
  }

  ///////////////////////////////////////////////////////////////////////////////
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t i2 = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, i1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 6 + ax, i2 + ax, _w));

    id0 += 9;
  }
  vec3 _N0;
  vec3 _cen;
  real _d;
};
#endif


class edge_willmore : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_willmore)

  edge_willmore(const std::vector<index_t> &ids, const real &w,
                std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _align(false),
        _N_align(vec3::Zero()) {}

  void init_rest(const vecX &q, solver_context &ctx) override {
    (void)ctx;
    if (_rest_initialized)
      return;
    const vec3 qi0 = _blocks[0]->get_vec3(_ids[0], q);
    const vec3 qi1 = _blocks[0]->get_vec3(_ids[1], q);
    const vec3 qj0 = _blocks[0]->get_vec3(_ids[2], q);
    const vec3 qj1 = _blocks[0]->get_vec3(_ids[3], q);
    const vec3 d00 = qi1 - qi0;
    const vec3 d01 = qi1 - qj0;
    const vec3 d10 = qj1 - qi0;
    const vec3 d11 = qj1 - qj0;
    const vec3 N0 = d00.cross(d01).normalized();
    const vec3 N1 = d11.cross(d10).normalized();
    _N_align = (N0 + N1).normalized();
    _align = true;
    _rest_initialized = true;
  }

  void set_align_normal(vec3 N) {
    _align = true;
    _N_align = N;
  }

  virtual void project(const vecX &q, vecX &p) {
    // d00   i0   j1 d11
    //       /|1/
    // d01  /0|/     d10
    // i1    j0

    vec3 qi0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 qi1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 qj0 = _blocks[0]->get_vec3(_ids[2], q);
    vec3 qj1 = _blocks[0]->get_vec3(_ids[3], q);

    vec3 c0 = 0.333 * (qi0 + qi1 + qj0);
    vec3 c1 = 0.333 * (qj0 + qj1 + qi0);

    vec3 d0 = (qi0 - qj0).normalized();

    vec3 d00 = qi1 - qi0;
    vec3 d01 = qi1 - qj0;

    vec3 d10 = qj1 - qi0;
    vec3 d11 = qj1 - qj0;

    vec3 N0 = d00.cross(d01).normalized();
    vec3 N1 = d11.cross(d10).normalized();
    real cs = va::cos(N0, N1, d0);
    real sn = va::sin(N0, N1, d0);
    real thet = 0.5 * atan2(sn, cs);
    // real thet = atan(va::tan(N0, N1, d0));
    vec3 N = (N0 + N1).normalized();
    // mat3 R = va::rejection_matrix(N);

    Eigen::AngleAxis<real> R(thet, d0);
    vec3 d00p = R * d00;
    vec3 d01p = R * d01;
    vec3 d10p = R.inverse() * d10;
    vec3 d11p = R.inverse() * d11;
    if (_align) {
      vec3 Ni = (N0 + N1).normalized();
      quat Ra;
      Ra.setFromTwoVectors(Ni, _N_align);
      d00p = Ra * d00p;
      d01p = Ra * d01p;
      d10p = Ra * d10p;
      d11p = Ra * d11p;
#if 0
      if (!_debugged) {
        _debugged = true;
        gg::geometry_geometry_logger::line(c0, c0 + 0.01 * _N_align,
                                  vec4(0.0, 1.0, 0.0, 1.0));
        gg::geometry_geometry_logger::line(qi0, qi0 + d00p, vec4(1.0, 0.5, 0.0, 1.0));
        gg::geometry_geometry_logger::line(qj0, qj0 + d01p, vec4(0.0, 1.0, 0.5, 1.0));
        gg::geometry_geometry_logger::line(qi0, qi0 + d10p, vec4(0.0, 0.5, 1.0, 1.0));
        gg::geometry_geometry_logger::line(qj0, qj0 + d11p, vec4(0.5, 1.0, 0.0, 1.0));
      }
#endif
    }

    if (d00p.hasNaN() || d01p.hasNaN() || d10p.hasNaN() || d11p.hasNaN()) {
      std::cout << __PRETTY_FUNCTION__ << " d2 is nan" << std::endl;
      exit(0);
    }
    p.block(_id0 + 0, 0, 3, 1) = _w * d00p;
    // p.block(_id0 + 3, 0, 3, 1) = _w * d01p;
    // p.block(_id0 + 6, 0, 3, 1) = _w * d10p;
    p.block(_id0 + 3, 0, 3, 1) = _w * d11p;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;

    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t j0 = _blocks[0]->get_offset_idx(this->_ids[2]);
    index_t j1 = _blocks[0]->get_offset_idx(this->_ids[3]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, -_w));
    /*
        for (int ax = 0; ax < 3; ax++)
          triplets.push_back(trip(_id0 + 3 + ax, i1 + ax, _w));
        for (int ax = 0; ax < 3; ax++)
          triplets.push_back(trip(_id0 + 3 + ax, j0 + ax, -_w));

        for (int ax = 0; ax < 3; ax++)
          triplets.push_back(trip(_id0 + 6 + ax, j1 + ax, _w));
        for (int ax = 0; ax < 3; ax++)
          triplets.push_back(trip(_id0 + 6 + ax, i0 + ax, -_w));
    */
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, j1 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, j0 + ax, -_w));

    id0 += 6;
  }
  vec3 _N_align; // Goal N
  bool _align = false;
  bool _debugged = false;
};


class edge_normal_flow : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_normal_flow)

  edge_normal_flow(const std::vector<index_t> &ids, const vector<vec3> & x, real d, const real &w,
                std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _d(d) {

    vec3 pi0 = x[_ids[0]];
    vec3 pi1 = x[_ids[1]];
    vec3 pj0 = x[_ids[2]];
    vec3 pj1 = x[_ids[3]];
    vec3 e0 = pi0 - pj0;
    vec3 N0 = (pi1 - pi0).cross(e0).normalized();
    vec3 N1 = -(pj1 - pj0).cross(e0).normalized();
    _N0 = (N0 + N1).normalized();
    _cen = 0.5 * (pi0 + pj0);
  }


  virtual void project(const vecX &q, vecX &p) {
    // d00   i0   j1 d11
    //       /|1/
    // d01  /0|/     d10
    // i1    j0

    vec3 qi0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 qi1 = _blocks[0]->get_vec3(_ids[1], q);
    vec3 qj0 = _blocks[0]->get_vec3(_ids[2], q);
    vec3 qj1 = _blocks[0]->get_vec3(_ids[3], q);

    vec3 c0 = 0.333 * (qi0 + qi1 + qj0);
    vec3 c1 = 0.333 * (qj0 + qj1 + qi0);

    vec3 e0 = (qi0 - qj0).normalized();

    vec3 d00 = qi1 - qi0;
    vec3 d01 = qi1 - qj0;

    vec3 d10 = qj1 - qi0;
    vec3 d11 = qj1 - qj0;

    vec3 N0 = (qi1 - qi0).cross(e0).normalized();
    vec3 N1 = -(qj1 - qj0).cross(e0).normalized();
    vec3 N = (N0 + N1).normalized();

    real di = (qi0 - _cen).dot(_N0);
    real dj = (qj0 - _cen).dot(_N0);
    vec3 qi0p = qi0 + (_d - di) * N;
    vec3 qj0p = qj0 + (_d - dj) * N;
    //geometry_logger::line(qi0, qi0p, vec4(0.0, 1.0, 0.5, 1.0));
    //geometry_logger::line(qj0, qj0p, vec4(0.0, 1.0, 0.5, 1.0));
    
    p.block(_id0 + 0, 0, 3, 1) = _w * qi0p;
    p.block(_id0 + 3, 0, 3, 1) = _w * qj0p;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;

    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t i1 = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t j0 = _blocks[0]->get_offset_idx(this->_ids[2]);
    index_t j1 = _blocks[0]->get_offset_idx(this->_ids[3]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i0 + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, j0 + ax, _w));

    id0 += 6;
  }
  vec3 _N0; // Goal N
  vec3 _cen;
  real _d;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// Body log of a unit quaternion. NaN when the input is not a finite rotation.
inline vec3 frame_quat_log(quat q) {
  const real n2 = q.squaredNorm();
  if (!(n2 > 1e-24) || !std::isfinite(n2))
    return vec3::Constant(std::numeric_limits<real>::quiet_NaN());
  q.normalize();
  if (q.w() < 0.0)
    q.coeffs() *= -1.0;
  const vec3 v = q.vec();
  const real vn = v.norm();
  if (vn < 1e-14)
    return vec3::Zero();
  const real w = std::clamp(q.w(), real(-1.0), real(1.0));
  const real angle = real(2.0) * std::atan2(vn, w);
  return (angle / vn) * v;
}

inline bool frame_quat_ok(const quat &q) {
  const real n2 = q.squaredNorm();
  return std::isfinite(n2) && n2 > 0.25;
}

/// Symmetric shape operator W t = s, unknowns (κ11, κ12, κ22).
struct shape_sample {
  real tx = 0;
  real ty = 0;
  real sx = 0;
  real sy = 0;
};

struct shape_w {
  real k11 = 0;
  real k12 = 0;
  real k22 = 0;
  real k1 = 0;
  real k2 = 0;
  real e1x = 1;
  real e1y = 0;
  real e2x = 0;
  real e2y = 1;
  bool ok = false;
};

/// `fraction * |κ|^power` along the chart direction `(tx, ty)`. `κ` is the
/// normal curvature there, one value for every spoke of that axis. `power` 0
/// returns `fraction`, including where κ is 0. `|κ|` below 1e-3 is treated as
/// 1e-3 so a negative power stays finite. A missing operator leaves `fraction`.
inline real stretch_axis_scale(const shape_w &W, real tx, real ty, real fraction,
                              real power) {
  if (!std::isfinite(fraction))
    return real(0);
  if (!W.ok || !std::isfinite(power) || !(std::abs(power) > real(1e-16)) ||
      !std::isfinite(tx) || !std::isfinite(ty))
    return fraction;
  const real kappa =
      W.k11 * tx * tx + real(2) * W.k12 * tx * ty + W.k22 * ty * ty;
  if (!std::isfinite(kappa))
    return fraction;
  real ak = std::abs(kappa);
  if (!(ak > real(1e-3)))
    ak = real(1e-3);
  const real gain = std::pow(ak, power);
  if (!std::isfinite(gain))
    return fraction;
  return fraction * gain;
}

inline void sym2_eig(real a, real b, real c, real &k1, real &k2, real &x1,
                     real &y1, real &x2, real &y2) {
  const real tr = a + c;
  const real det = a * c - b * b;
  const real disc = std::sqrt(std::max(real(0), tr * tr * real(0.25) - det));
  k1 = real(0.5) * tr + disc;
  k2 = real(0.5) * tr - disc;
  real x = b;
  real y = k1 - a;
  const real n2 = x * x + y * y;
  if (n2 < real(1e-24)) {
    x = real(1);
    y = real(0);
  } else {
    const real inv = real(1) / std::sqrt(n2);
    x *= inv;
    y *= inv;
  }
  x1 = x;
  y1 = y;
  x2 = -y;
  y2 = x;
}

inline shape_w fit_shape_w(const std::vector<shape_sample> &samples) {
  shape_w W;
  int used = 0;
  mat3 ATA = mat3::Zero();
  vec3 Atb = vec3::Zero();
  for (const shape_sample &s : samples) {
    if (!std::isfinite(s.tx) || !std::isfinite(s.ty) || !std::isfinite(s.sx) ||
        !std::isfinite(s.sy))
      continue;
    const vec3 r1(s.tx, s.ty, real(0));
    const vec3 r2(real(0), s.tx, s.ty);
    ATA += r1 * r1.transpose();
    ATA += r2 * r2.transpose();
    Atb += s.sx * r1;
    Atb += s.sy * r2;
    ++used;
  }
  if (used < 2)
    return W;
  Eigen::JacobiSVD<mat3> svd(ATA, Eigen::ComputeThinU | Eigen::ComputeThinV);
  const vec3 sv = svd.singularValues();
  const vec3 coef = svd.solve(Atb);
  if (!coef.allFinite() || !(sv[0] > real(0)) || !(sv[2] > real(1e-8) * sv[0]))
    return W;
  W.k11 = coef[0];
  W.k12 = coef[1];
  W.k22 = coef[2];
  sym2_eig(W.k11, W.k12, W.k22, W.k1, W.k2, W.e1x, W.e1y, W.e2x, W.e2y);
  W.ok = std::isfinite(W.k1) && std::isfinite(W.k2);
  return W;
}

inline bool shape_umbilic(const shape_w &W) {
  const real scale = std::abs(W.k1) + std::abs(W.k2);
  if (!(scale > real(1e-8)))
    return true;
  return std::abs(W.k1 - W.k2) < real(1e-2) * scale;
}

inline real shape_quad(const shape_w &W, real tx, real ty) {
  return W.k11 * tx * tx + real(2) * W.k12 * tx * ty + W.k22 * ty * ty;
}

/// Curvature along world X in the tangent plane becomes K + sin(x). The
/// perpendicular curvature stays the one that was measured. x is the vertex
/// world coordinate.
inline real curvature_along(const shape_w &W, const vec3 &T, const vec3 &B,
                            const vec3 &N, const vec3 &dir_in) {
  if (!W.ok || !T.allFinite() || !B.allFinite() || !N.allFinite() ||
      !dir_in.allFinite())
    return std::numeric_limits<real>::quiet_NaN();
  vec3 D = dir_in - N * N.dot(dir_in);
  if (!(D.squaredNorm() > real(1e-16)))
    return std::numeric_limits<real>::quiet_NaN();
  D.normalize();
  return shape_quad(W, D.dot(T), D.dot(B));
}

/// Replace the normal curvature along `dir_in` with `k_new`. The perpendicular
/// curvature stays measured. A direction with no tangent part leaves `W`.
inline shape_w shape_with_direction_curvature(const shape_w &W, const vec3 &T,
                                              const vec3 &B, const vec3 &N,
                                              const vec3 &dir_in, real k_new) {
  shape_w out = W;
  if (!W.ok || !std::isfinite(k_new) || !T.allFinite() || !B.allFinite() ||
      !N.allFinite() || !dir_in.allFinite())
    return out;
  vec3 D = dir_in - N * N.dot(dir_in);
  if (!(D.squaredNorm() > real(1e-16)))
    return out;
  D.normalize();
  vec3 P = N.cross(D);
  if (!(P.squaredNorm() > real(1e-16)))
    return out;
  P.normalize();
  const real ex = D.dot(T);
  const real ey = D.dot(B);
  const real px = P.dot(T);
  const real py = P.dot(B);
  const real k_perp = shape_quad(W, px, py);
  if (!std::isfinite(k_perp))
    return out;
  out.k11 = k_new * ex * ex + k_perp * px * px;
  out.k22 = k_new * ey * ey + k_perp * py * py;
  out.k12 = k_new * ex * ey + k_perp * px * py;
  out.ok = std::isfinite(out.k11) && std::isfinite(out.k12) &&
           std::isfinite(out.k22);
  return out;
}

inline shape_w shape_with_x_axis_curvature(const shape_w &W, const vec3 &T,
                                           const vec3 &B, const vec3 &N,
                                           real x) {
  const real k = curvature_along(W, T, B, N, vec3::UnitX());
  if (!std::isfinite(k))
    return W;
  return shape_with_direction_curvature(W, T, B, N, vec3::UnitX(),
                                        hopf_curvature_plus_sin(k, x));
}

/// Curvature `k_new` along `dir_in` and 0 across it. A direction with no
/// tangent part leaves the shape unset.
inline shape_w shape_along_tangent(const vec3 &T, const vec3 &B, const vec3 &N,
                                   const vec3 &dir_in, real k_new) {
  shape_w out;
  if (!std::isfinite(k_new) || !T.allFinite() || !B.allFinite() ||
      !N.allFinite() || !dir_in.allFinite())
    return out;
  vec3 D = dir_in - N * N.dot(dir_in);
  if (!(D.squaredNorm() > real(1e-16)))
    return out;
  D.normalize();
  const real ex = D.dot(T);
  const real ey = D.dot(B);
  out.k11 = k_new * ex * ex;
  out.k22 = k_new * ey * ey;
  out.k12 = k_new * ex * ey;
  out.ok = std::isfinite(out.k11) && std::isfinite(out.k12) &&
           std::isfinite(out.k22);
  return out;
}

/// Minor-circle fiber. `tube` is that tangent after `fiber` rotates it.
/// `kappa` is the measured curvature there times cos(frequency * phase).
struct hopf_geodesic {
  vec3 major = vec3::Zero();
  vec3 tube = vec3::Zero();
  real kappa = 0;
  bool ok = false;
};

inline hopf_geodesic hopf_phase_geodesic(const shape_w &W, const vec3 &T,
                                         const vec3 &B, const vec3 &N,
                                         const vec3 &point, real frequency = 1,
                                         real fiber = 0) {
  hopf_geodesic g;
  vec3 major0 = vec3::Zero();
  vec3 tube0 = vec3::Zero();
  if (!hopf_revolution_axes(point, N, major0, tube0))
    return g;
  hopf_fiber_axes(major0, tube0, fiber, g.major, g.tube);
  const real k = curvature_along(W, T, B, N, g.tube);
  if (!std::isfinite(k))
    return g;
  const real freq = std::isfinite(frequency) ? frequency : real(1);
  g.kappa = k * std::cos(freq * hopf_normal_angle(N, g.major));
  g.ok = std::isfinite(g.kappa);
  return g;
}

/// Perpendicular curvature for Pinkall's shape operator. Measured fiber
/// curvature `k0` stays, and `k1 = -1/k0`, so `k0/k1 = -k0²`.
inline real pinkall_perpendicular_curvature(real k0) {
  if (!std::isfinite(k0) || !(std::abs(k0) > real(1e-8)))
    return std::numeric_limits<real>::quiet_NaN();
  const real k1 = -(k_hopf_sphere_curvature * k_hopf_sphere_curvature) / k0;
  return std::isfinite(k1) ? k1 : std::numeric_limits<real>::quiet_NaN();
}

/// Two fiber curvatures with mean `H` and cross ratio `k0 / k1 = skew`.
/// `skew` 1 gives the umbilic `H, H`.
inline void cross_ratio_curvatures(real H, real skew, real &k0, real &k1) {
  k0 = std::numeric_limits<real>::quiet_NaN();
  k1 = std::numeric_limits<real>::quiet_NaN();
  if (!std::isfinite(H) || !std::isfinite(skew))
    return;
  const real denom = real(1) + skew;
  if (!(std::abs(denom) > real(1e-8)))
    return;
  k0 = real(2) * H * skew / denom;
  k1 = real(2) * H / denom;
}

/// Hyperboloid pair in a chosen frame. Same-sign curvatures flip the smaller
/// one, so the axes passed in stay the saddle axes.
inline void hyperboloid_curvature_pair(real k0, real k1, real &ku, real &kv) {
  ku = k0;
  kv = k1;
  if (!std::isfinite(k0) || !std::isfinite(k1)) {
    ku = std::numeric_limits<real>::quiet_NaN();
    kv = std::numeric_limits<real>::quiet_NaN();
    return;
  }
  if (!(k0 * k1 > real(0)))
    return;
  if (std::abs(k0) >= std::abs(k1))
    kv = -k1;
  else
    ku = -k0;
}

/// Principal directions of `W` in the tangent frame, rotated by `twist` about
/// the normal. `u` is the first principal, `v` the second. `twist` 0 leaves
/// them on the measured axes.
inline bool principal_axes_twisted(const shape_w &W, const vec3 &T, const vec3 &B,
                                   const vec3 &N, real twist, vec3 &u, vec3 &v) {
  u = vec3::Zero();
  v = vec3::Zero();
  if (!W.ok || !T.allFinite() || !B.allFinite() || !N.allFinite() ||
      !std::isfinite(W.e1x) || !std::isfinite(W.e1y) || !std::isfinite(W.e2x) ||
      !std::isfinite(W.e2y))
    return false;
  vec3 u0 = W.e1x * T + W.e1y * B;
  vec3 v0 = W.e2x * T + W.e2y * B;
  u0 -= N * u0.dot(N);
  v0 -= N * v0.dot(N);
  if (!(u0.squaredNorm() > real(1e-16)) || !(v0.squaredNorm() > real(1e-16)))
    return false;
  u0.normalize();
  v0.normalize();
  const real ang = std::isfinite(twist) ? twist : real(0);
  const real c = std::cos(ang);
  const real s = std::sin(ang);
  u = c * u0 + s * v0;
  v = c * v0 - s * u0;
  return u.allFinite() && v.allFinite() && u.squaredNorm() > real(1e-16) &&
         v.squaredNorm() > real(1e-16);
}

/// `axial` on the stored shape operator. `angle` is signed radians from T
/// toward B in the coefficient basis. Normal curvature along that direction
/// is multiplied by `1 + scale`. The perpendicular curvature and the cross
/// term stay. `scale` 0 leaves the operator. Reversing the direction does
/// not change the result.
inline shape_w map_stored_axial(const shape_w &H, real angle, real scale) {
  shape_w out;
  if (!H.ok)
    return out;
  if (!std::isfinite(angle))
    angle = real(0);
  if (!std::isfinite(scale))
    scale = real(0);
  const real c = std::cos(angle);
  const real s = std::sin(angle);
  const real kappa = H.k11 * c * c + real(2) * H.k12 * c * s + H.k22 * s * s;
  if (!std::isfinite(kappa))
    return out;
  const real a = scale * kappa;
  out.k11 = H.k11 + a * c * c;
  out.k12 = H.k12 + a * c * s;
  out.k22 = H.k22 + a * s * s;
  if (!std::isfinite(out.k11) || !std::isfinite(out.k12) || !std::isfinite(out.k22))
    return shape_w{};
  sym2_eig(out.k11, out.k12, out.k22, out.k1, out.k2, out.e1x, out.e1y, out.e2x,
           out.e2y);
  out.ok = std::isfinite(out.k1) && std::isfinite(out.k2);
  return out.ok ? out : shape_w{};
}

/// `fiber_skew` on the stored shape operator. `skew` splits the mean onto the
/// principal axes, `twist` rotates those axes. The hyperboloid sign flip is
/// not part of this map. The chart is the coefficient basis, so no frame is
/// required.
inline shape_w map_stored_fiber_skew(const shape_w &H, real skew, real twist) {
  shape_w out;
  if (!H.ok)
    return out;
  const real mean = real(0.5) * (H.k1 + H.k2);
  real k0 = 0;
  real k1 = 0;
  cross_ratio_curvatures(mean, skew, k0, k1);
  if (!std::isfinite(k0) || !std::isfinite(k1))
    return out;
  const real ang = std::isfinite(twist) ? twist : real(0);
  const real c = std::cos(ang);
  const real s = std::sin(ang);
  const real e1x = c * H.e1x + s * H.e2x;
  const real e1y = c * H.e1y + s * H.e2y;
  const real e2x = c * H.e2x - s * H.e1x;
  const real e2y = c * H.e2y - s * H.e1y;
  out.k11 = k0 * e1x * e1x + k1 * e2x * e2x;
  out.k12 = k0 * e1x * e1y + k1 * e2x * e2y;
  out.k22 = k0 * e1y * e1y + k1 * e2y * e2y;
  if (!std::isfinite(out.k11) || !std::isfinite(out.k12) || !std::isfinite(out.k22))
    return shape_w{};
  sym2_eig(out.k11, out.k12, out.k22, out.k1, out.k2, out.e1x, out.e1y, out.e2x,
           out.e2y);
  out.ok = std::isfinite(out.k1) && std::isfinite(out.k2);
  return out.ok ? out : shape_w{};
}

/// One neighbor of a stored shape, in that vertex's principal frame.
/// `u`, `v` are lengths along e1 and e2. `psi` is the neighbor eigenframe's
/// turn relative to e1, in [-π/2, π/2]. `framed` is false when that turn
/// could not be read.
struct ring_shape_sample {
  real u = 0;
  real v = 0;
  real w = 0;
  real k1 = 0;
  real k2 = 0;
  real psi = 0;
  bool framed = false;
  /// Position-chart eigenvalues. Sharpen reads these, not the stored pair.
  real gk1 = 0;
  real gk2 = 0;
  bool geom = false;
  /// Signed angle of the neighbor position frame, in (-π, π]. No half-turn flip.
  real gpsi = 0;
  bool gframed = false;
};

inline shape_w shape_from_principals(real k1, real k2, real e1x, real e1y,
                                     real e2x, real e2y) {
  shape_w out;
  if (!std::isfinite(k1) || !std::isfinite(k2))
    return out;
  out.k11 = k1 * e1x * e1x + k2 * e2x * e2x;
  out.k12 = k1 * e1x * e1y + k2 * e2x * e2y;
  out.k22 = k1 * e1y * e1y + k2 * e2y * e2y;
  if (!std::isfinite(out.k11) || !std::isfinite(out.k12) || !std::isfinite(out.k22))
    return shape_w{};
  sym2_eig(out.k11, out.k12, out.k22, out.k1, out.k2, out.e1x, out.e1y, out.e2x,
           out.e2y);
  out.ok = std::isfinite(out.k1) && std::isfinite(out.k2);
  return out.ok ? out : shape_w{};
}

/// κ1 moves toward neighbors that lie on e1, κ2 toward those on e2.
inline shape_w map_ring_dupin(const shape_w &H,
                              const std::vector<ring_shape_sample> &samples) {
  if (!H.ok)
    return {};
  real sum1 = 0;
  real acc1 = 0;
  real sum2 = 0;
  real acc2 = 0;
  for (const ring_shape_sample &s : samples) {
    if (!(s.w > real(0)) || !std::isfinite(s.k1) || !std::isfinite(s.k2))
      continue;
    const real l2 = s.u * s.u + s.v * s.v;
    if (!(l2 > real(1e-16)))
      continue;
    const real a1 = s.w * (s.u * s.u) / l2;
    const real a2 = s.w * (s.v * s.v) / l2;
    acc1 += a1 * s.k1;
    sum1 += a1;
    acc2 += a2 * s.k2;
    sum2 += a2;
  }
  const real k1 = sum1 > real(1e-12) ? acc1 / sum1 : H.k1;
  const real k2 = sum2 > real(1e-12) ? acc2 / sum2 : H.k2;
  return shape_from_principals(k1, k2, H.e1x, H.e1y, H.e2x, H.e2y);
}

/// Hold κ1 κ2 and add `skew` to log|κ1/κ2|. Signs stay. `skew` 0 is the identity.
inline bool log_skew_pair(real k1, real k2, real skew, real &o1, real &o2) {
  o1 = k1;
  o2 = k2;
  const real a1 = std::abs(k1);
  const real a2 = std::abs(k2);
  if (!(a1 > real(1e-8)) || !(a2 > real(1e-8)) || !std::isfinite(skew) ||
      !std::isfinite(k1) || !std::isfinite(k2))
    return false;
  const real scale = std::exp(real(0.5) * (std::log(a1 / a2) + skew));
  const real root = std::sqrt(a1 * a2);
  if (!(scale > real(0)) || !std::isfinite(scale) || !std::isfinite(root))
    return false;
  const real m1 = root * scale;
  const real m2 = root / scale;
  if (!std::isfinite(m1) || !std::isfinite(m2))
    return false;
  o1 = std::copysign(m1, k1);
  o2 = std::copysign(m2, k2);
  return true;
}

/// Hold κ1 κ2 and set `|κ1/κ2|` to `ratio`. Signs stay. `ratio` 1 makes
/// equal principal magnitudes. A non-positive ratio fails.
inline bool absolute_ratio_pair(real k1, real k2, real ratio, real &o1,
                                real &o2) {
  o1 = k1;
  o2 = k2;
  if (!(ratio > real(0)) || !std::isfinite(ratio) || !std::isfinite(k1) ||
      !std::isfinite(k2))
    return false;
  const real a1 = std::abs(k1);
  const real a2 = std::abs(k2);
  if (!(a1 > real(1e-8)) || !(a2 > real(1e-8)))
    return false;
  const real gap = std::log(ratio) - std::log(a1 / a2);
  return log_skew_pair(k1, k2, gap, o1, o2);
}

/// log|κ| flows toward the neighbors on its own principal line. `skew` is then
/// added to log|κ1/κ2| and the product of that flowed pair stays.
inline shape_w map_ring_log_skew(const shape_w &H,
                                 const std::vector<ring_shape_sample> &samples,
                                 real skew) {
  if (!H.ok)
    return {};
  real sum1 = 0;
  real acc1 = 0;
  real sum2 = 0;
  real acc2 = 0;
  for (const ring_shape_sample &s : samples) {
    if (!(s.w > real(0)))
      continue;
    const real l2 = s.u * s.u + s.v * s.v;
    if (!(l2 > real(1e-16)))
      continue;
    const real w1 = s.w * (s.u * s.u) / l2;
    const real w2 = s.w * (s.v * s.v) / l2;
    const real a1 = std::abs(s.k1);
    const real a2 = std::abs(s.k2);
    if (a1 > real(1e-8) && std::isfinite(a1)) {
      acc1 += w1 * std::log(a1);
      sum1 += w1;
    }
    if (a2 > real(1e-8) && std::isfinite(a2)) {
      acc2 += w2 * std::log(a2);
      sum2 += w2;
    }
  }
  real k1 = H.k1;
  real k2 = H.k2;
  if (sum1 > real(1e-12) && std::abs(H.k1) > real(1e-8)) {
    const real mag = std::exp(acc1 / sum1);
    if (std::isfinite(mag))
      k1 = std::copysign(mag, H.k1);
  }
  if (sum2 > real(1e-12) && std::abs(H.k2) > real(1e-8)) {
    const real mag = std::exp(acc2 / sum2);
    if (std::isfinite(mag))
      k2 = std::copysign(mag, H.k2);
  }
  real s1 = k1;
  real s2 = k2;
  log_skew_pair(k1, k2, skew, s1, s2);
  return shape_from_principals(s1, s2, H.e1x, H.e1y, H.e2x, H.e2y);
}

/// Mean becomes the weighted neighbor mean. The trace-free part stays.
inline shape_w map_ring_cmc(const shape_w &H,
                            const std::vector<ring_shape_sample> &samples) {
  if (!H.ok)
    return {};
  real sw = 0;
  real sH = 0;
  for (const ring_shape_sample &s : samples) {
    if (!(s.w > real(0)) || !std::isfinite(s.k1) || !std::isfinite(s.k2))
      continue;
    sw += s.w;
    sH += s.w * real(0.5) * (s.k1 + s.k2);
  }
  if (!(sw > real(1e-12)))
    return H;
  const real shift = sH / sw - real(0.5) * (H.k1 + H.k2);
  shape_w out = H;
  out.k11 += shift;
  out.k22 += shift;
  sym2_eig(out.k11, out.k12, out.k22, out.k1, out.k2, out.e1x, out.e1y, out.e2x,
           out.e2y);
  out.ok = std::isfinite(out.k1) && std::isfinite(out.k2);
  return out.ok ? out : shape_w{};
}

/// Rotate the eigenframe by a signed angle in radians. `twist` and `fiber`
/// add. Positive turns e1 toward e2. The base is the position-chart frame
/// when `measured` is set, so the angle is the same on the next iteration.
/// A stored-only call rotates that stored frame by the neighbor angle plus
/// the same offset. `sharpen` is a length². 0 leaves κ1 and κ2.
/// Otherwise the mean of `measured` moves by `sharpen` times its second
/// difference on the position chart. The step is capped by that mean, and
/// κ1 κ2 of the chart stays. The stored pair is not the source: reading it
/// back would multiply the wave on every iteration.
inline shape_w map_ring_spin(const shape_w &H,
                             const std::vector<ring_shape_sample> &samples,
                             real twist, real fiber, real sharpen = 0,
                             const shape_w *measured = nullptr) {
  if (!H.ok)
    return {};
  const real tw = std::isfinite(twist) ? twist : real(0);
  const real fb = std::isfinite(fiber) ? fiber : real(0);
  const real sh = std::isfinite(sharpen) ? sharpen : real(0);
  const real offset = tw + fb;
  const real full_turn = real(2) * std::acos(real(-1));
  const shape_w *base = &H;
  real mean = 0;
  bool from_chart = false;
  if (measured && measured->ok && std::isfinite(measured->e1x) &&
      std::isfinite(measured->e1y) && std::isfinite(measured->e2x) &&
      std::isfinite(measured->e2y)) {
    real sw = 0;
    real sa = 0;
    for (const ring_shape_sample &s : samples) {
      if (!s.gframed || !(s.w > real(0)) || !std::isfinite(s.gpsi))
        continue;
      sw += s.w;
      sa += s.w * s.gpsi;
    }
    if (sw > real(1e-12))
      mean = sa / sw;
    base = measured;
    from_chart = true;
  }
  if (!from_chart) {
    real sw = 0;
    real sa = 0;
    for (const ring_shape_sample &s : samples) {
      if (!s.framed || !(s.w > real(0)) || !std::isfinite(s.psi))
        continue;
      sw += s.w;
      sa += s.w * s.psi;
    }
    if (sw > real(1e-12))
      mean = sa / sw;
  }
  real alpha = std::remainder(mean + offset, full_turn);
  if (!std::isfinite(alpha))
    alpha = real(0);
  if (alpha == real(0) && sh == real(0) && !from_chart)
    return H;
  const real c = std::cos(alpha);
  const real s = std::sin(alpha);
  const real e1x = c * base->e1x + s * base->e2x;
  const real e1y = c * base->e1y + s * base->e2y;
  const real e2x = c * base->e2x - s * base->e1x;
  const real e2y = c * base->e2y - s * base->e1y;
  real k1 = H.k1;
  real k2 = H.k2;
  if (sh != real(0) && measured && measured->ok && std::isfinite(measured->k1) &&
      std::isfinite(measured->k2)) {
    real nw = 0;
    real nH = 0;
    real h2 = 0;
    for (const ring_shape_sample &sm : samples) {
      if (!sm.geom || !(sm.w > real(0)) || !std::isfinite(sm.gk1) ||
          !std::isfinite(sm.gk2))
        continue;
      const real r2 = sm.u * sm.u + sm.v * sm.v;
      if (!(r2 > real(1e-16)))
        continue;
      nw += sm.w;
      nH += sm.w * real(0.5) * (sm.gk1 + sm.gk2);
      h2 += sm.w * r2;
    }
    if (nw > real(1e-12) && h2 > real(1e-16)) {
      const real mu = real(0.5) * (measured->k1 + measured->k2);
      const real mu_ring = nH / nw;
      const real h_avg = h2 / nw;
      real step = sh * (mu - mu_ring) / h_avg;
      const real cap = std::max(std::abs(mu), std::abs(mu_ring));
      if (!std::isfinite(step) || !(cap > real(0)))
        step = real(0);
      else if (std::abs(step) > cap)
        step = std::copysign(cap, step);
      real mu2 = mu + step;
      const real prod = measured->k1 * measured->k2;
      real disc2 = mu2 * mu2 - prod;
      if (!(disc2 > real(0))) {
        const real floor = prod > real(0) ? std::sqrt(prod) : real(0);
        mu2 = std::copysign(floor, mu2);
        disc2 = real(0);
      }
      if (std::isfinite(mu2) && std::isfinite(disc2)) {
        const real disc = std::sqrt(disc2);
        k1 = mu2 + disc;
        k2 = mu2 - disc;
      }
    }
  }
  return shape_from_principals(k1, k2, e1x, e1y, e2x, e2y);
}

inline bool fit_ring_gradient(const std::vector<ring_shape_sample> &samples,
                              bool use_frame, real &g1, real &g2, bool kappa1) {
  real uu = 0;
  real uv = 0;
  real vv = 0;
  real ub = 0;
  real vb = 0;
  for (const ring_shape_sample &s : samples) {
    if (!(s.w > real(0)))
      continue;
    real b = 0;
    if (use_frame) {
      if (!s.framed || !std::isfinite(s.psi))
        continue;
      b = s.psi;
    } else if (kappa1) {
      if (!std::isfinite(s.k1))
        continue;
      b = s.k1;
    } else {
      if (!std::isfinite(s.k2))
        continue;
      b = s.k2;
    }
    uu += s.w * s.u * s.u;
    uv += s.w * s.u * s.v;
    vv += s.w * s.v * s.v;
    ub += s.w * s.u * b;
    vb += s.w * s.v * b;
  }
  const real det = uu * vv - uv * uv;
  const real scale = uu + vv;
  if (!(det > real(1e-12) * scale * scale) && !(det > real(1e-12)))
    return false;
  g1 = (vv * ub - uv * vb) / det;
  g2 = (uu * vb - uv * ub) / det;
  return std::isfinite(g1) && std::isfinite(g2);
}

/// (κ1−κ2) ω(e1) = ∂2 κ1 and (κ1−κ2) ω(e2) = ∂1 κ2. Mean and frame stay.
/// A vanishing spin leaves the gap.
inline shape_w map_ring_codazzi(const shape_w &H,
                                const std::vector<ring_shape_sample> &samples) {
  if (!H.ok)
    return {};
  std::vector<ring_shape_sample> shifted = samples;
  for (ring_shape_sample &s : shifted) {
    s.k1 -= H.k1;
    s.k2 -= H.k2;
  }
  real d1k1 = 0;
  real d2k1 = 0;
  real d1k2 = 0;
  real d2k2 = 0;
  real w1 = 0;
  real w2 = 0;
  if (!fit_ring_gradient(shifted, false, d1k1, d2k1, true) ||
      !fit_ring_gradient(shifted, false, d1k2, d2k2, false) ||
      !fit_ring_gradient(samples, true, w1, w2, false))
    return H;
  (void)d1k1;
  (void)d2k2;
  const real wn = w1 * w1 + w2 * w2;
  if (!(wn > real(1e-16)))
    return H;
  const real delta = (w1 * d2k1 + w2 * d1k2) / wn;
  const real mean = real(0.5) * (H.k1 + H.k2);
  return shape_from_principals(mean + real(0.5) * delta, mean - real(0.5) * delta,
                               H.e1x, H.e1y, H.e2x, H.e2y);
}

/// `fiber_skew` on a measured operator: split the trace, rotate the axes, flip
/// the smaller curvature when the pair has one sign. Spoke landing, not the
/// stored-H row.
inline shape_w prescribe_fiber_skew(const shape_w &measured, const vec3 &T,
                                   const vec3 &B, const vec3 &N, real skew,
                                   real twist) {
  shape_w out;
  if (!measured.ok)
    return out;
  const real H = real(0.5) * (measured.k1 + measured.k2);
  real k0 = 0;
  real k1 = 0;
  real ku = 0;
  real kv = 0;
  cross_ratio_curvatures(H, skew, k0, k1);
  vec3 u = vec3::Zero();
  vec3 v = vec3::Zero();
  if (!principal_axes_twisted(measured, T, B, N, twist, u, v))
    return out;
  hyperboloid_curvature_pair(k0, k1, ku, kv);
  if (!std::isfinite(ku) || !std::isfinite(kv))
    return out;
  const shape_w line = shape_along_tangent(T, B, N, u, ku);
  return shape_with_direction_curvature(line, T, B, N, v, kv);
}

/// 1D jet: `kappa` along the minor-circle tangent and 0 across it.
inline shape_w shape_with_hopf_line_curvature(const shape_w &W, const vec3 &T,
                                              const vec3 &B, const vec3 &N,
                                              const vec3 &point,
                                              real frequency = 1,
                                              real fiber = 0) {
  const hopf_geodesic g =
      hopf_phase_geodesic(W, T, B, N, point, frequency, fiber);
  if (!g.ok)
    return {};
  return shape_along_tangent(T, B, N, g.tube, g.kappa);
}

/// Spoke `d` from the foot, onto the circle of curvature `g.kappa` in the
/// `(tube, N)` plane. The `major` component stays.
inline bool conformal_hopf_target(const hopf_geodesic &g, const vec3 &N,
                                  const vec3 &d, vec3 &target) {
  if (!g.ok || !d.allFinite() || !N.allFinite())
    return false;
  const real hn = d.dot(N);
  const vec3 t = d - hn * N;
  real uo = 0;
  real ho = 0;
  if (!hopf_curvature_land(g.kappa, t.dot(g.tube), hn, uo, ho))
    return false;
  target = uo * g.tube + t.dot(g.major) * g.major + ho * N;
  return target.allFinite();
}

inline bool frame_axes(const quat &q_in, vec3 &T, vec3 &B, vec3 &N) {
  if (!frame_quat_ok(q_in))
    return false;
  quat q = q_in;
  q.normalize();
  T = q * vec3::UnitX();
  N = q * vec3::UnitZ();
  if (N.squaredNorm() < 1e-24)
    return false;
  N.normalize();
  T -= T.dot(N) * N;
  if (T.squaredNorm() < 1e-24)
    return false;
  T.normalize();
  B = N.cross(T);
  return B.squaredNorm() > 1e-24;
}

inline bool tangent_basis(const vec3 &Nin, vec3 &T, vec3 &B) {
  if (Nin.squaredNorm() < 1e-24)
    return false;
  const vec3 N = Nin.normalized();
  const vec3 a = std::abs(N[0]) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
  T = N.cross(a);
  if (T.squaredNorm() < 1e-24)
    return false;
  T.normalize();
  B = N.cross(T);
  return B.squaredNorm() > 1e-24;
}

inline quat frame_quat_from_axes(const vec3 &T, const vec3 &B, const vec3 &N,
                                const quat &current) {
  mat3 R;
  R.col(0) = T;
  R.col(1) = B;
  R.col(2) = N;
  quat q(R);
  q.normalize();
  if (!frame_quat_ok(q))
    return current;
  if (q.coeffs().dot(current.coeffs()) < 0.0)
    q.coeffs() *= real(-1);
  return q;
}

struct curvature_frame {
  quat q = quat::Identity();
  bool data = false;
};

/// Local curvature frame from the one-ring face normals. `data` is set when
/// red is the principal axis of larger |κ|. An umbilic leaves `data` false.
inline curvature_frame local_curvature_frame(const asawa::shell::shell &M,
                                            const std::vector<vec3> &x, int i,
                                            const quat &current) {
  curvature_frame out;
  out.q = current;
  if (i < 0 || i >= static_cast<int>(x.size()))
    return out;
  const asawa::shell::VertId vi = asawa::shell::vert_id(i);
  if (static_cast<int>(M.vbegin(vi)) < 0)
    return out;
  const vec3 N = asawa::shell::vert_normal(M, vi, x);
  if (!N.allFinite() || N.squaredNorm() < 1e-12)
    return out;
  vec3 T0, B0;
  if (!tangent_basis(N, T0, B0))
    return out;

  vec3 Tcur, Bcur, Ncur;
  vec3 Tkeep = T0;
  if (frame_axes(current, Tcur, Bcur, Ncur)) {
    Tkeep = Tcur - Tcur.dot(N) * N;
    if (Tkeep.squaredNorm() < 1e-16)
      Tkeep = T0;
    else
      Tkeep.normalize();
  }

  std::vector<shape_sample> samples;
  M.const_for_each_vertex(vi, [&](asawa::shell::CornerId ci,
                                  const asawa::shell::shell &S) {
    const asawa::shell::FaceId f = S.face(ci);
    if (static_cast<int>(f) < 0 || static_cast<int>(S.fbegin(f)) < 0)
      return;
    const vec3 Nf = asawa::shell::face_normal(S, f, x);
    if (!Nf.allFinite())
      return;
    const vec3 c = asawa::shell::face_center(S, f, x);
    vec3 t = c - x[static_cast<size_t>(i)];
    t -= t.dot(N) * N;
    if (t.squaredNorm() < 1e-16)
      return;
    vec3 dn = Nf - N;
    dn -= dn.dot(N) * N;
    samples.push_back(
        shape_sample{t.dot(T0), t.dot(B0), -dn.dot(T0), -dn.dot(B0)});
  });

  const shape_w W = fit_shape_w(samples);
  if (!W.ok || shape_umbilic(W))
    return out;
  const vec3 e1 = (W.e1x * T0 + W.e1y * B0).normalized();
  const vec3 e2 = (W.e2x * T0 + W.e2y * B0).normalized();
  if (!e1.allFinite() || !e2.allFinite())
    return out;
  vec3 T = std::abs(W.k1) >= std::abs(W.k2) ? e1 : e2;
  if (T.dot(Tkeep) < 0.0)
    T = -T;
  if (T.squaredNorm() < 1e-16)
    return out;
  T.normalize();
  vec3 B = N.cross(T);
  if (B.squaredNorm() < 1e-24)
    return out;
  B.normalize();
  out.q = frame_quat_from_axes(T, B, N, current);
  out.data = frame_quat_ok(out.q);
  return out;
}

inline quat frame_hemisphere(const quat &ref, quat q) {
  if (!frame_quat_ok(q))
    return q;
  q.normalize();
  if (ref.coeffs().dot(q.coeffs()) < 0.0)
    q.coeffs() *= real(-1);
  return q;
}

/// Body log of q_i^{-1} q_j. z is the spin about the frame normal.
struct frame_body_omega {
  real x = 0;
  real y = 0;
  real z = 0;
  bool ok = false;
};

inline frame_body_omega frame_body_omega_of(const quat &qi_in, const quat &qj_in) {
  frame_body_omega w;
  if (!frame_quat_ok(qi_in) || !frame_quat_ok(qj_in))
    return w;
  const quat qi = qi_in.normalized();
  const quat qj = frame_hemisphere(qi, qj_in);
  const vec3 omega = frame_quat_log(qi.conjugate() * qj);
  if (!omega.allFinite())
    return w;
  w.x = omega.x();
  w.y = omega.y();
  w.z = omega.z();
  w.ok = true;
  return w;
}

inline quat frame_quat_exp(const vec3 &omega) {
  const real ang = omega.norm();
  if (!(ang > real(1e-14)) || !std::isfinite(ang))
    return quat::Identity();
  return quat(Eigen::AngleAxis<real>(ang, omega / ang));
}

inline shape_w transport_shape(const shape_w &Wj, const vec3 &Tj, const vec3 &Bj,
                              const vec3 &Ti, const vec3 &Bi, const vec3 &Ni) {
  auto apply = [&](const vec3 &v) {
    const real vx = v.dot(Tj);
    const real vy = v.dot(Bj);
    return (Wj.k11 * vx + Wj.k12 * vy) * Tj + (Wj.k12 * vx + Wj.k22 * vy) * Bj;
  };
  vec3 ST = apply(Ti);
  vec3 SB = apply(Bi);
  ST -= ST.dot(Ni) * Ni;
  SB -= SB.dot(Ni) * Ni;
  shape_w O;
  O.k11 = Ti.dot(ST);
  O.k22 = Bi.dot(SB);
  O.k12 = real(0.5) * (Ti.dot(SB) + Bi.dot(ST));
  O.ok = std::isfinite(O.k11) && std::isfinite(O.k12) && std::isfinite(O.k22);
  return O;
}

/// A 2×2 map, row-major. `m01` is the top-right entry.
struct linear_map2 {
  real m00 = 1;
  real m01 = 0;
  real m10 = 0;
  real m11 = 1;
  bool ok = false;
};

/// Scale, skew, and symmetric traceless shear of a 2×2.
struct linear_map_parts {
  real sigma = 0;
  real omega = 0;
  real h00 = 0;
  real h01 = 0;
  real h11 = 0;
};

inline linear_map_parts split_linear_map(const linear_map2 &M) {
  linear_map_parts p;
  p.sigma = real(0.5) * (M.m00 + M.m11);
  p.omega = real(0.5) * (M.m10 - M.m01);
  p.h00 = M.m00 - p.sigma;
  p.h11 = M.m11 - p.sigma;
  p.h01 = real(0.5) * (M.m01 + M.m10);
  return p;
}

inline vec2 apply_linear_map(const linear_map2 &M, const vec2 &e) {
  return vec2(M.m00 * e.x() + M.m01 * e.y(), M.m10 * e.x() + M.m11 * e.y());
}

/// Spoke target in the Dirichlet chart. `bend` is along û = (tx, ty). `twist`
/// is along the conormal ν = (−ty, tx).
inline vec2 darboux_spoke_target(real tx, real ty, real bend, real twist) {
  if (!std::isfinite(bend))
    bend = real(0);
  if (!std::isfinite(twist))
    twist = real(0);
  return vec2(bend * tx - twist * ty, bend * ty + twist * tx);
}

/// Measured map plus the gap to `skew` and `shear`. Scale stays. Shear is the
/// tangent component. The 45° shear is driven to 0. A failed fit stays failed.
inline linear_map2 cross_ratio_spoke_map(const linear_map2 &M, real skew, real shear) {
  linear_map2 O;
  O.ok = false;
  if (!M.ok)
    return O;
  if (!std::isfinite(skew))
    skew = real(0);
  if (!std::isfinite(shear))
    shear = real(0);
  const linear_map_parts p = split_linear_map(M);
  O.m00 = p.sigma + shear;
  O.m01 = -skew;
  O.m10 = skew;
  O.m11 = p.sigma - shear;
  O.ok = std::isfinite(O.m00) && std::isfinite(O.m01) && std::isfinite(O.m10) &&
         std::isfinite(O.m11);
  if (!O.ok)
    return linear_map2{};
  return O;
}

inline vec2 shape_dir_action(real k11, real k12, real k22, real tx, real ty) {
  return vec2(k11 * tx + k12 * ty, k12 * tx + k22 * ty);
}

/// Columns are the images of the coefficient basis `(k11, k12, k22)` under `(P W) t`.
inline void shape_transport_jacobian(const vec3 &Tj, const vec3 &Bj, const vec3 &Ti,
                                    const vec3 &Bi, const vec3 &Ni, real tx, real ty,
                                    real J[6]) {
  for (int c = 0; c < 3; ++c) {
    shape_w basis;
    if (c == 0)
      basis.k11 = 1;
    else if (c == 1)
      basis.k12 = 1;
    else
      basis.k22 = 1;
    const shape_w O = transport_shape(basis, Tj, Bj, Ti, Bi, Ni);
    const vec2 s = shape_dir_action(O.k11, O.k12, O.k22, tx, ty);
    J[c] = s.x();
    J[3 + c] = s.y();
  }
}

/// Least-squares `e ≈ M e0`. Empty or rank-deficient samples leave `ok` false.
inline linear_map2 fit_spoke_map(const std::vector<real> &w, const std::vector<vec2> &e,
                                 const std::vector<vec2> &e0) {
  linear_map2 M;
  M.ok = false;
  const int n = static_cast<int>(std::min(w.size(), std::min(e.size(), e0.size())));
  real a = 0;
  real b = 0;
  real c = 0;
  real n00 = 0;
  real n01 = 0;
  real n10 = 0;
  real n11 = 0;
  for (int i = 0; i < n; ++i) {
    if (!(w[static_cast<size_t>(i)] > real(0)) || !e[static_cast<size_t>(i)].allFinite() ||
        !e0[static_cast<size_t>(i)].allFinite())
      continue;
    const real wk = w[static_cast<size_t>(i)];
    const real rx = e0[static_cast<size_t>(i)].x();
    const real ry = e0[static_cast<size_t>(i)].y();
    const real ex = e[static_cast<size_t>(i)].x();
    const real ey = e[static_cast<size_t>(i)].y();
    a += wk * rx * rx;
    b += wk * rx * ry;
    c += wk * ry * ry;
    n00 += wk * ex * rx;
    n01 += wk * ex * ry;
    n10 += wk * ey * rx;
    n11 += wk * ey * ry;
  }
  const real det = a * c - b * b;
  if (!(std::abs(det) > real(1e-16)) || !std::isfinite(det))
    return M;
  const real inv = real(1) / det;
  M.m00 = (n00 * c - n01 * b) * inv;
  M.m01 = (-n00 * b + n01 * a) * inv;
  M.m10 = (n10 * c - n11 * b) * inv;
  M.m11 = (-n10 * b + n11 * a) * inv;
  M.ok = std::isfinite(M.m00) && std::isfinite(M.m01) && std::isfinite(M.m10) &&
         std::isfinite(M.m11);
  if (!M.ok)
    return linear_map2{};
  return M;
}

/// Least-squares scale `e ≈ σ e0`. A vanishing rest gradient returns 0.
inline real fit_isotropic_scale(const std::vector<real> &w, const std::vector<vec2> &e,
                                const std::vector<vec2> &e0) {
  const int n = static_cast<int>(std::min(w.size(), std::min(e.size(), e0.size())));
  real num = 0;
  real den = 0;
  for (int i = 0; i < n; ++i) {
    if (!(w[static_cast<size_t>(i)] > real(0)) || !e[static_cast<size_t>(i)].allFinite() ||
        !e0[static_cast<size_t>(i)].allFinite())
      continue;
    const real wk = w[static_cast<size_t>(i)];
    num += wk * e[static_cast<size_t>(i)].dot(e0[static_cast<size_t>(i)]);
    den += wk * e0[static_cast<size_t>(i)].squaredNorm();
  }
  if (!(den > real(1e-16)) || !std::isfinite(num))
    return real(0);
  const real sigma = num / den;
  return std::isfinite(sigma) ? sigma : real(0);
}


struct fair_spoke {
  index_t j = -1;
  real K = 0;
  char ok = 0;
};

struct fair_ring {
  index_t i = -1;
  std::vector<fair_spoke> spokes;
  /// Local face and edge-dof indices. Empty on the quaternion chart.
  std::vector<int> faces;
  std::vector<int> edges;
};

struct ring_frame {
  shape_w W;
  vec3 T = vec3::UnitX();
  vec3 B = vec3::UnitY();
  vec3 N = vec3::UnitZ();
  bool axes = false;
};

/// Shape operator of one vertex from neighbor quaternion normals in the frame
/// (T, B, N). Half the trace is H. The trace does not depend on which way T
/// points. The monge path reads the Crouzeix–Raviart operator instead.
inline ring_frame fit_ring_frame(const sim_block &pos, const sim_block &frame,
                                const vecX &q, const fair_ring &ring) {
  ring_frame out;
  const quat ui = frame.get_quat(ring.i, q);
  vec3 T, B, N;
  if (!frame_axes(ui, T, B, N))
    return out;
  out.T = T;
  out.B = B;
  out.N = N;
  out.axes = true;
  const vec3 xi = pos.get_vec3(ring.i, q);
  std::vector<shape_sample> samples;
  samples.reserve(ring.spokes.size());
  for (const fair_spoke &sp : ring.spokes) {
    if (!sp.ok)
      continue;
    quat uj = frame.get_quat(sp.j, q);
    if (!frame_quat_ok(uj))
      continue;
    uj.normalize();
    const vec3 Ns = uj * vec3::UnitZ();
    if (!Ns.allFinite() || Ns.squaredNorm() < real(1e-24))
      continue;
    const vec3 xj = pos.get_vec3(sp.j, q);
    vec3 t = xj - xi;
    t -= t.dot(N) * N;
    if (t.squaredNorm() < real(1e-16))
      continue;
    vec3 dn = Ns.normalized() - N;
    dn -= dn.dot(N) * N;
    if (!dn.allFinite())
      continue;
    samples.push_back(shape_sample{t.dot(T), t.dot(B), -dn.dot(T), -dn.dot(B)});
  }
  out.W = fit_shape_w(samples);
  return out;
}

/// 2×2 of a 3×3 shape operator in an orthonormal tangent frame.
inline shape_w shape_w_from_operator(const mat3 &S, const vec3 &T, const vec3 &B) {
  shape_w W;
  W.k11 = T.dot(S * T);
  W.k12 = T.dot(S * B);
  W.k22 = B.dot(S * B);
  if (!std::isfinite(W.k11) || !std::isfinite(W.k12) || !std::isfinite(W.k22))
    return W;
  sym2_eig(W.k11, W.k12, W.k22, W.k1, W.k2, W.e1x, W.e1y, W.e2x, W.e2y);
  W.ok = std::isfinite(W.k1) && std::isfinite(W.k2);
  return W;
}

inline shape_w shape_w_from_coeffs(const vec3 &c) {
  shape_w W;
  if (!c.allFinite())
    return W;
  W.k11 = c[0];
  W.k12 = c[1];
  W.k22 = c[2];
  sym2_eig(W.k11, W.k12, W.k22, W.k1, W.k2, W.e1x, W.e1y, W.e2x, W.e2y);
  W.ok = std::isfinite(W.k1) && std::isfinite(W.k2);
  return W;
}

inline vec3 coeffs_of_shape(const shape_w &W) {
  return vec3(W.k11, W.k12, W.k22);
}

/// Chart written by the vertex polar cache. Empty or unmarked vertices fail.
inline bool vertex_chart_at(const std::vector<vertex_chart> *charts, index_t i,
                            shape_w &W, vec3 &T, vec3 &B, vec3 &N) {
  if (!charts || i < 0 || static_cast<size_t>(i) >= charts->size())
    return false;
  const vertex_chart &c = (*charts)[static_cast<size_t>(i)];
  if (!c.ok || !c.h.allFinite() || !c.T.allFinite() || !c.B.allFinite() ||
      !c.N.allFinite() || !(c.T.squaredNorm() > real(1e-16)) ||
      !(c.B.squaredNorm() > real(1e-16)))
    return false;
  W = shape_w_from_coeffs(c.h);
  if (!W.ok)
    return false;
  T = c.T;
  B = c.B;
  N = c.N;
  return true;
}

/// Chart frame, with `S` read live from the edge store when `stars` is set:
/// `blocks[0]` positions, `blocks[1]` the edge store, `S` the local polar of
/// the vertex's incident edge tensors in the chart's `(T, B, N)`.
inline bool vertex_row_frame(const std::vector<sim_block::ptr> &blocks,
                             const std::vector<vertex_chart> *charts, index_t i,
                             const vecX &q, shape_w &W, vec3 &T, vec3 &B,
                             vec3 &N,
                             const kusama::weingarten::vertex_stars *stars = nullptr) {
  if (!vertex_chart_at(charts, i, W, T, B, N))
    return false;
  if (!stars || blocks.size() < 2 || !blocks[0] || !blocks[1] ||
      i >= static_cast<index_t>(stars->at.size()))
    return true;
  const vec3 h = kusama::weingarten::vertex_tensor(
      *stars, static_cast<int>(i), [&](int v) { return blocks[0]->get_vec3(v, q); },
      [&](int e) { return blocks[1]->get_vec3(e, q); }, T, B, N);
  const shape_w row = shape_w_from_coeffs(h);
  if (!row.ok)
    return true;
  W = row;
  return true;
}

/// Principal axis of a fitted ring, as a quaternion. An umbilic leaves `data`
/// false and keeps `current`.
inline curvature_frame curvature_frame_from_ring(const ring_frame &rf,
                                                const quat &current) {
  curvature_frame out;
  out.q = current;
  if (!rf.axes || !rf.W.ok || shape_umbilic(rf.W))
    return out;
  const vec3 e1 = rf.W.e1x * rf.T + rf.W.e1y * rf.B;
  const vec3 e2 = rf.W.e2x * rf.T + rf.W.e2y * rf.B;
  vec3 T = std::abs(rf.W.k1) >= std::abs(rf.W.k2) ? e1 : e2;
  if (!T.allFinite() || T.squaredNorm() < real(1e-16))
    return out;
  T.normalize();
  if (T.dot(rf.T) < real(0))
    T = -T;
  vec3 B = rf.N.cross(T);
  if (B.squaredNorm() < real(1e-24))
    return out;
  B.normalize();
  out.q = frame_quat_from_axes(T, B, rf.N, current);
  out.data = frame_quat_ok(out.q);
  return out;
}

inline fair_ring fair_ring_of(index_t i, const std::vector<index_t> &spokes,
                              const std::vector<real> &K,
                              const std::vector<char> &ok) {
  fair_ring ring;
  ring.i = i;
  const int n = static_cast<int>(spokes.size());
  ring.spokes.reserve(static_cast<size_t>(std::max(n, 0)));
  for (int k = 0; k < n; ++k) {
    fair_spoke sp;
    sp.j = spokes[static_cast<size_t>(k)];
    sp.K = k < static_cast<int>(K.size()) ? K[static_cast<size_t>(k)] : real(0);
    sp.ok = k < static_cast<int>(ok.size()) ? ok[static_cast<size_t>(k)] : 0;
    ring.spokes.push_back(sp);
  }
  return ring;
}

/// One spoke of a scalar on the one-ring. `(px, py)` is the spoke in the
/// chart tangent `(T, B)`, and `K` is its cotangent weight.
struct ring_scalar_sample {
  real s = 0;
  real px = 0;
  real py = 0;
  real K = 0;
};

/// Tangent gradient of one scalar. `div` is the cotangent sum
/// `Σ K (s_j − s)`, so a caller can step without a second pass.
struct ring_scalar_gradient {
  real gx = 0;
  real gy = 0;
  real div = 0;
  bool ok = false;
};

/// Least squares of `s_j − s` against the spoke tangent. A zero weight is
/// skipped. The fit does not know which scalar it was given.
inline ring_scalar_gradient
fit_ring_scalar_gradient(real s, const std::vector<ring_scalar_sample> &samples) {
  ring_scalar_gradient out;
  real axx = 0;
  real axy = 0;
  real ayy = 0;
  real bx = 0;
  real by = 0;
  for (const ring_scalar_sample &sp : samples) {
    if (!std::isfinite(sp.s) || !std::isfinite(sp.px) || !std::isfinite(sp.py) ||
        !std::isfinite(sp.K) || sp.K == real(0))
      continue;
    const real ds = sp.s - s;
    if (!std::isfinite(ds))
      continue;
    out.div += sp.K * ds;
    axx += sp.K * sp.px * sp.px;
    axy += sp.K * sp.px * sp.py;
    ayy += sp.K * sp.py * sp.py;
    bx += sp.K * sp.px * ds;
    by += sp.K * sp.py * ds;
  }
  if (!std::isfinite(out.div))
    out.div = 0;
  const real det = axx * ayy - axy * axy;
  const real scale = std::max(std::abs(axx), std::abs(ayy));
  if (!(scale > real(0)) || !(std::abs(det) > real(1e-16) * scale * scale))
    return out;
  const real gx = (ayy * bx - axy * by) / det;
  const real gy = (axx * by - axy * bx) / det;
  if (!std::isfinite(gx) || !std::isfinite(gy))
    return out;
  out.gx = gx;
  out.gy = gy;
  out.ok = true;
  return out;
}

/// `(k11, k12, k22)` of a diagonal operator in the chart tangent.
/// `k0` lies along `(e1x, e1y)` and `k1` along `(e2x, e2y)`.
inline vec3 principal_coeffs(real k0, real k1, real e1x, real e1y, real e2x,
                             real e2y) {
  return vec3(k0 * e1x * e1x + k1 * e2x * e2x,
              k0 * e1x * e1y + k1 * e2x * e2y,
              k0 * e1y * e1y + k1 * e2y * e2y);
}

/// Directional curvature of a chart along `axis`. The axis is projected into
/// the chart tangent. The quadratic form is even, so the sign of `axis` does
/// not matter and a neighbor umbilic cannot swap the two samples.
inline real chart_normal_curvature(const vertex_chart &c, const vec3 &axis) {
  if (!c.ok || !c.h.allFinite() || !c.T.allFinite() || !c.B.allFinite() ||
      !c.N.allFinite() || !axis.allFinite())
    return std::numeric_limits<real>::quiet_NaN();
  vec3 t = axis - c.N.dot(axis) * c.N;
  const real n2 = t.squaredNorm();
  if (!(n2 > real(1e-16)) || !t.allFinite())
    return std::numeric_limits<real>::quiet_NaN();
  t /= std::sqrt(n2);
  const real x = t.dot(c.T);
  const real y = t.dot(c.B);
  const real kappa =
      c.h[0] * x * x + real(2) * c.h[1] * x * y + c.h[2] * y * y;
  if (!std::isfinite(kappa))
    return std::numeric_limits<real>::quiet_NaN();
  return kappa;
}

/// Crouzeix–Raviart chart of one fan. Positions and coefficients are read by
/// id. The fan's faces are enough: an edge incident to the vertex has both of
/// its triangles in the fan, and those are the triangles that build its
/// derivative.
struct fan_chart {
  shape_w W;
  vec3 T = vec3::UnitX();
  vec3 B = vec3::UnitY();
  vec3 N = vec3::UnitZ();
  bool axes = false;
};

/// Called from a constraint constructor. The integer fan is not kept.
inline fan_local bake_fan(const std::vector<fan_face> &faces,
                          const std::vector<fan_face> &support,
                          const std::vector<fan_edge> &edges,
                          const std::vector<int> &tangent, int vertex) {
  fan_local out;
  if (faces.empty() || vertex < 0)
    return out;

  auto consider = [&](int v) {
    if (v >= 0)
      out.vids.push_back(v);
  };
  auto consider_faces = [&](const std::vector<fan_face> &src) {
    for (const fan_face &f : src) {
      for (int i = 0; i < 3; ++i)
        consider(f.v[static_cast<size_t>(i)]);
    }
  };
  consider_faces(faces);
  consider_faces(support);
  for (const fan_edge &fe : edges) {
    consider(fe.v0);
    consider(fe.v1);
  }
  std::sort(out.vids.begin(), out.vids.end());
  out.vids.erase(std::unique(out.vids.begin(), out.vids.end()), out.vids.end());
  std::unordered_map<int, int> vloc;
  for (int i = 0; i < static_cast<int>(out.vids.size()); ++i)
    vloc.emplace(out.vids[static_cast<size_t>(i)], i);
  auto intern = [&](int v) {
    if (v < 0)
      return -1;
    const auto it = vloc.find(v);
    return it == vloc.end() ? -1 : it->second;
  };
  const auto vit = vloc.find(vertex);
  if (vit == vloc.end())
    return fan_local{};
  out.center = vit->second;

  std::vector<int> dofs;
  auto consider_edge_faces = [&](const std::vector<fan_face> &src) {
    for (const fan_face &f : src) {
      for (int i = 0; i < 3; ++i) {
        if (f.e[static_cast<size_t>(i)] >= 0)
          dofs.push_back(f.e[static_cast<size_t>(i)]);
      }
    }
  };
  consider_edge_faces(faces);
  consider_edge_faces(support);
  for (const fan_edge &fe : edges) {
    if (fe.dof >= 0)
      dofs.push_back(fe.dof);
  }
  std::sort(dofs.begin(), dofs.end());
  dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
  std::unordered_map<int, int> eloc;
  out.ledges.resize(dofs.size());
  out.dof_of.assign(dofs.size(), -1);
  for (int i = 0; i < static_cast<int>(dofs.size()); ++i) {
    const int dof = dofs[static_cast<size_t>(i)];
    eloc.emplace(dof, i);
    out.dof_of[static_cast<size_t>(i)] = dof;
    cr_edge e;
    for (const fan_edge &fe : edges) {
      if (fe.dof == dof) {
        e.v0 = intern(fe.v0);
        e.v1 = intern(fe.v1);
        break;
      }
    }
    out.ledges[static_cast<size_t>(i)] = e;
  }
  auto intern_e = [&](int dof) {
    if (dof < 0)
      return -1;
    const auto it = eloc.find(dof);
    return it == eloc.end() ? -1 : it->second;
  };

  out.n_incident = static_cast<int>(faces.size());
  out.tris.reserve(faces.size() + support.size());
  auto push_tris = [&](const std::vector<fan_face> &src) {
    for (const fan_face &f : src) {
      cr_tri t;
      for (int i = 0; i < 3; ++i) {
        t.v[static_cast<size_t>(i)] = intern(f.v[static_cast<size_t>(i)]);
        t.e[static_cast<size_t>(i)] = intern_e(f.e[static_cast<size_t>(i)]);
      }
      out.tris.push_back(t);
    }
  };
  push_tris(faces);
  push_tris(support);

  out.tangent_e.reserve(tangent.size());
  for (const int dof : tangent) {
    const int local = intern_e(dof);
    if (local >= 0)
      out.tangent_e.push_back(local);
  }
  return out;
}

template <class PosAt, class CoeffAt>
inline fan_chart evaluate_fan(const fan_local &fan, int nE, PosAt pos_at,
                              CoeffAt coeff_at, vec3 *lock) {
  fan_chart out;
  if (fan.center < 0 || fan.vids.empty() || fan.n_incident <= 0)
    return out;

  std::vector<vec3> x(fan.vids.size(), vec3::Zero());
  for (int i = 0; i < static_cast<int>(fan.vids.size()); ++i)
    x[static_cast<size_t>(i)] = pos_at(fan.vids[static_cast<size_t>(i)]);

  const int nL = static_cast<int>(fan.ledges.size());
  std::vector<real> coeff(static_cast<size_t>(std::max(2 * nL, 0)), real(0));
  for (int e = 0; e < nL; ++e) {
    const int dof = e < static_cast<int>(fan.dof_of.size())
                        ? fan.dof_of[static_cast<size_t>(e)]
                        : -1;
    if (dof < 0)
      continue;
    coeff[static_cast<size_t>(e)] = coeff_at(dof);
    if (nE > 0)
      coeff[static_cast<size_t>(e + nL)] = coeff_at(dof + nE);
  }

  const std::vector<cr_face_operator> faceS =
      cr_face_operators(x, fan.tris, fan.ledges);
  std::vector<int> incident(static_cast<size_t>(std::max(fan.n_incident, 0)));
  std::iota(incident.begin(), incident.end(), 0);
  vec3 N;
  if (!cr_fan_normal(x, fan.tris, incident, fan.center, N))
    return out;
  vec3 T;
  if (!cr_dirichlet_tangent(x, fan.ledges, coeff, fan.tangent_e, N, T, lock)) {
    const vec3 a = std::abs(N.x()) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
    T = N.cross(a);
    if (!(T.squaredNorm() > real(1e-24)))
      return out;
    T.normalize();
  }
  vec3 B = N.cross(T);
  if (!(B.squaredNorm() > real(1e-24)))
    return out;
  B.normalize();
  T = B.cross(N);
  if (!T.allFinite() || !B.allFinite())
    return out;
  if (lock)
    *lock = T;
  out.T = T;
  out.B = B;
  out.N = N;
  out.axes = true;
  out.W = shape_w_from_operator(cr_vertex_operator(faceS, incident), T, B);
  return out;
}

/// One third of the incident face areas, spokes in ring order.
inline real ring_vertex_area(const sim_block &pos, const vecX &q,
                            const fair_ring &ring) {
  const vec3 xi = pos.get_vec3(ring.i, q);
  const int n = static_cast<int>(ring.spokes.size());
  real faces = 0;
  for (int k = 0; k < n; ++k) {
    const fair_spoke &a = ring.spokes[static_cast<size_t>(k)];
    const fair_spoke &b =
        ring.spokes[static_cast<size_t>((k + 1) % n)];
    const vec3 e0 = pos.get_vec3(a.j, q) - xi;
    const vec3 e1 = pos.get_vec3(b.j, q) - xi;
    const real cr = e0.cross(e1).norm();
    if (std::isfinite(cr))
      faces += real(0.5) * cr;
  }
  return faces / real(3);
}

/// Cotangent combination Σ K (H_j − H_i). Divide by vertex area for ΔH.
inline real cotan_diff_H(const fair_ring &ring, real Hi,
                        const std::vector<real> &H, const std::vector<char> &ok,
                        const std::vector<int> &which) {
  real acc = 0;
  for (const fair_spoke &sp : ring.spokes) {
    if (!sp.ok || !std::isfinite(sp.K))
      continue;
    if (sp.j < 0 || static_cast<size_t>(sp.j) >= which.size())
      continue;
    const int jr = which[static_cast<size_t>(sp.j)];
    if (jr < 0 || !ok[static_cast<size_t>(jr)])
      continue;
    acc += sp.K * (H[static_cast<size_t>(jr)] - Hi);
  }
  return acc;
}

/// Mean spoke length of a ring.
inline real ring_mean_edge(const sim_block &pos, const vecX &q,
                          const fair_ring &ring) {
  const vec3 xi = pos.get_vec3(ring.i, q);
  real sum = 0;
  int n = 0;
  for (const fair_spoke &sp : ring.spokes) {
    const real len = (pos.get_vec3(sp.j, q) - xi).norm();
    if (std::isfinite(len) && len > real(0)) {
      sum += len;
      ++n;
    }
  }
  return n > 0 ? sum / static_cast<real>(n) : 0;
}

/// Nearest real root of ΔZ + 2Z(Z² − K) = 0 to the measured H.
/// The Laplacian is frozen on the stored field: `sum_w` and `sum_wZ` are the
/// neighbor weights and weighted stored curvatures, and `Z_old` is this
/// vertex's stored value. K is the measured Gaussian curvature. A smooth
/// sphere returns H. A smooth cylinder returns 0. A root farther from H than
/// the current curvature is rejected: the frozen cotangent term is then
/// outside the cubic, and taking it writes a curvature the next jet overflows.
inline real project_willmore_curvature(real H, real Z_old, real Kcurv,
                                       real area, real sum_w, real sum_wZ) {
  if (!std::isfinite(H) || !std::isfinite(Z_old) || !std::isfinite(Kcurv) ||
      !std::isfinite(area) || !std::isfinite(sum_w) || !std::isfinite(sum_wZ))
    return H;
  if (!(area > real(1e-12)) || !(sum_w > real(1e-12)))
    return H;
  const real lap = sum_wZ - sum_w * Z_old;
  const real p = -Kcurv;
  const real q = lap / (real(4) * area);
  const real scale =
      std::max(std::abs(H), std::sqrt(std::max(real(0), std::abs(Kcurv))));
  auto accept = [&](real z) {
    if (!std::isfinite(z))
      return H;
    if (std::abs(z - H) > scale)
      return H;
    return z;
  };
  auto nearest = [&](std::initializer_list<real> roots) {
    real best = H;
    real best_d = std::numeric_limits<real>::infinity();
    for (const real z : roots) {
      if (!std::isfinite(z))
        continue;
      const real d = std::abs(z - H);
      if (d < best_d) {
        best = z;
        best_d = d;
      }
    }
    return accept(best);
  };
  const real disc = real(0.25) * q * q + (p * p * p) / real(27);
  if (p < real(0) && !(disc > real(0))) {
    const real r = std::sqrt(-p / real(3));
    const real r3 = r * r * r;
    if (!(r3 > real(0)))
      return H;
    const real arg = std::clamp((-real(0.5) * q) / r3, real(-1), real(1));
    const real theta = std::acos(arg);
    const real two_pi = real(2) * std::acos(real(-1));
    return nearest({real(2) * r * std::cos(theta / real(3)),
                    real(2) * r * std::cos((theta + two_pi) / real(3)),
                    real(2) * r * std::cos((theta + real(2) * two_pi) / real(3))});
  }
  if (disc < real(0))
    return H;
  const real s = std::sqrt(std::max(disc, real(0)));
  const real z = std::cbrt(-real(0.5) * q + s) + std::cbrt(-real(0.5) * q - s);
  return accept(z);
}

/// One descent step of ∫ H² along the vertex normal. g ~ 1/length³, so
/// g L⁴ is a length. The explicit step is limited to a fraction of L: the
/// cotangent residual is otherwise many edge lengths and the next iteration
/// refits an even larger H. Zero when the ring is spherical.
inline real willmore_normal_step(real H, real K, real lapH, real edge) {
  if (!std::isfinite(H) || !std::isfinite(K) || !std::isfinite(lapH))
    return 0;
  if (!(edge > real(1e-8)))
    return 0;
  const real g = lapH + real(2) * H * (H * H - K);
  if (!std::isfinite(g))
    return 0;
  const real L2 = edge * edge;
  const real raw = -g * real(0.05) * L2 * L2;
  // Each projection iteration enforces this target. At saturation, thirty
  // iterations move a vertex by about half an edge in one frame.
  const real lim = real(0.02) * edge;
  const real a = std::abs(raw);
  if (!(a > real(0)))
    return 0;
  return raw * lim / (lim + a);
}

struct monge_axes {
  real k_hi = 0;
  real k_lo = 0;
  real hx = 1;
  real hy = 0;
  real lx = 0;
  real ly = 1;
};

inline shape_w shape_with_axes(shape_w W) {
  if (!W.ok)
    return W;
  sym2_eig(W.k11, W.k12, W.k22, W.k1, W.k2, W.e1x, W.e1y, W.e2x, W.e2y);
  W.ok = std::isfinite(W.k1) && std::isfinite(W.k2);
  return W;
}

/// `Pᵀ` of the principal tangents. `k0` lies on `k0_tan` and `k1` on `k1_tan`.
/// `(e1x, e1y)` and `(e2x, e2y)` are those tangents in the polar basis.
inline bool invert_principal_frame(real k0, real k1, const vec3 &k0_tan,
                                   const vec3 &k1_tan, real e1x, real e1y,
                                   real e2x, real e2y, const vec3 &N, vec3 &T,
                                   vec3 &B, vec3 &h) {
  if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)) ||
      !k0_tan.allFinite() || !k1_tan.allFinite() || !std::isfinite(k0) ||
      !std::isfinite(k1))
    return false;
  const vec3 Nn = N.normalized();
  T = e1x * k0_tan + e2x * k1_tan;
  B = e1y * k0_tan + e2y * k1_tan;
  T -= Nn * T.dot(Nn);
  B -= Nn * B.dot(Nn);
  if (!(T.squaredNorm() > real(1e-16)) || !(B.squaredNorm() > real(1e-16)) ||
      !T.allFinite() || !B.allFinite())
    return false;
  T.normalize();
  B = Nn.cross(T);
  if (!(B.squaredNorm() > real(1e-16)))
    return false;
  B.normalize();
  const vec3 t0 = k0_tan - Nn * k0_tan.dot(Nn);
  const vec3 t1 = k1_tan - Nn * k1_tan.dot(Nn);
  if (!(t0.squaredNorm() > real(1e-16)) || !(t1.squaredNorm() > real(1e-16)))
    return false;
  const vec3 u0 = t0.normalized();
  const vec3 u1 = t1.normalized();
  h = principal_coeffs(k0, k1, u0.dot(T), u0.dot(B), u1.dot(T), u1.dot(B));
  return h.allFinite();
}


inline monge_axes monge_axes_of(const shape_w &W) {
  monge_axes a;
  if (std::abs(W.k1) >= std::abs(W.k2)) {
    a.k_hi = W.k1;
    a.k_lo = W.k2;
    a.hx = W.e1x;
    a.hy = W.e1y;
    a.lx = W.e2x;
    a.ly = W.e2y;
  } else {
    a.k_hi = W.k2;
    a.k_lo = W.k1;
    a.hx = W.e2x;
    a.hy = W.e2y;
    a.lx = W.e1x;
    a.ly = W.e1y;
  }
  return a;
}

inline monge_axes specialize_monge_axes(monge_axes a, monge_shape shape) {
  if (shape == monge_shape::hyperboloid && a.k_hi * a.k_lo > real(0))
    a.k_lo = -a.k_lo;
  return a;
}

inline void directed_stretch::project(const vecX &q, vecX &p) {
  if (_blocks.empty() || _ids.empty() || _d0.empty())
    return;
  shape_w goal;
  vec3 frame_T = vec3::Zero();
  vec3 frame_B = vec3::Zero();
  vec3 frame_N = vec3::Zero();
  const bool charted =
      vertex_chart_at(_charts, _ids[0], goal, frame_T, frame_B, frame_N);
  const shape_w axes = shape_with_axes(goal);
  vec3 u = vec3::UnitX();
  vec3 v = vec3::UnitY();
  real sx = 0;
  real sy = 0;
  bool scale = false;
  if (charted && axes.ok) {
    const monge_axes ax = monge_axes_of(axes);
    const vec3 minor = ax.hx * frame_T + ax.hy * frame_B;
    const vec3 major = ax.lx * frame_T + ax.ly * frame_B;
    if (minor.allFinite() && major.allFinite() &&
        minor.squaredNorm() > real(1e-16) &&
        major.squaredNorm() > real(1e-16)) {
      stretch_frame(minor, major, _angle, u, v);
      const real c = std::cos(_angle);
      const real s = std::sin(_angle);
      real ux = c * ax.hx + s * ax.lx;
      real uy = c * ax.hy + s * ax.ly;
      real vx = -s * ax.hx + c * ax.lx;
      real vy = -s * ax.hy + c * ax.ly;
      if (const real l = std::hypot(_dir.dot(frame_T), _dir.dot(frame_B)); l > real(1e-8)) {
        ux = _dir.dot(frame_T) / l;
        uy = _dir.dot(frame_B) / l;
        vx = -uy;
        vy = ux;
        u = ux * frame_T + uy * frame_B;
        v = vx * frame_T + vy * frame_B;
      }
      sx = stretch_axis_scale(axes, ux, uy, _sx, _power);
      sy = stretch_axis_scale(axes, vx, vy, _sy, _power);
      scale = u.allFinite() && v.allFinite();
    }
  }
  // Stretch the rest spokes in the frame-start chart, then fit the rotation
  // to those: an exact projection onto {R S d0}, so the row is variational.
  std::vector<vec3> c(_d0.size());
  for (size_t k = 0; k < _d0.size(); ++k)
    c[k] = scale ? scale_spoke(scale_spoke(_d0[k], u, sx), v, sy) : _d0[k];
  // Rows are scaled by w_k, so the energy weighs spoke k by w_k^2.
  const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
  mat3 M = mat3::Zero();
  for (size_t k = 0; k < _edge_weights.size() && k < _d0.size(); ++k) {
    const real wk = _edge_weights[k];
    const vec3 dq = _blocks[0]->get_vec3(_ids[k + 1], q) - q0;
    M += wk * wk * dq * c[k].transpose();
  }
  const mat3 R = polar_rotation(M);
  _d_rest.resize(_d0.size());
  for (size_t k = 0; k < _edge_weights.size() && k < _d0.size(); ++k) {
    const real wk = _edge_weights[k];
    const vec3 spoke = R * c[k];
    _d_rest[k] = spoke;
    p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = _w * wk * spoke;
  }
}

inline vec3 tilt_orientation::chart_axis(const std::vector<vertex_chart> *charts, int i,
                                         principal_axis axis) {
  shape_w W;
  vec3 T, B, N;
  if (!vertex_chart_at(charts, i, W, T, B, N))
    return vec3::Zero();
  const bool on_k0 = axis == principal_axis::k0;
  const vec3 a = (on_k0 ? W.e1x : W.e2x) * T + (on_k0 ? W.e1y : W.e2y) * B;
  if (!a.allFinite() || !(a.squaredNorm() > real(1e-16)))
    return vec3::Zero();
  return a.normalized();
}

inline void tilt_orientation::orient(const std::vector<vertex_chart> *charts) {
  ready = true;
  const int n = static_cast<int>(nbr.size());
  sign.assign(static_cast<size_t>(n), real(1));
  std::vector<vec3> a(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    a[static_cast<size_t>(i)] = chart_axis(charts, i, axis);
  std::vector<char> seen(static_cast<size_t>(n), 0);
  std::vector<int> queue;
  for (int s = 0; s < n; ++s) {
    if (seen[static_cast<size_t>(s)] || a[static_cast<size_t>(s)].isZero())
      continue;
    seen[static_cast<size_t>(s)] = 1;
    queue.assign(1, s);
    for (size_t h = 0; h < queue.size(); ++h) {
      const size_t i = static_cast<size_t>(queue[h]);
      for (int j : nbr[i]) {
        if (j < 0 || j >= n || seen[static_cast<size_t>(j)] || a[static_cast<size_t>(j)].isZero())
          continue;
        seen[static_cast<size_t>(j)] = 1;
        sign[static_cast<size_t>(j)] =
            (sign[i] * a[i]).dot(a[static_cast<size_t>(j)]) < real(0) ? real(-1) : real(1);
        queue.push_back(j);
      }
    }
  }
}

inline void directed_tilt::hold_turn() {
  _target = _d0;
  shape_w W;
  vec3 T, B, N;
  if (_ids.empty() || !vertex_chart_at(_charts, _ids[0], W, T, B, N))
    return;
  const vec3 Nn = N.normalized();
  const bool on_k0 = _spec.axis == principal_axis::k0;
  const real ax = on_k0 ? W.e1x : W.e2x, ay = on_k0 ? W.e1y : W.e2y;
  vec3 a = ax * T + ay * B;
  if (!a.allFinite() || !(a.squaredNorm() > real(1e-16)))
    return;
  a.normalize();
  real l = 0;
  for (const vec3 &d : _d0)
    l += d.norm();
  l /= real(std::max<size_t>(_d0.size(), 1));
  if (!(l > real(0)))
    return;
  const real m = curvature_measure(_spec.metric, W.k1, W.k2);
  real theta = _spec.gain * m;
  vec3 axis = a;
  if (_spec.mode == tilt_mode::spin) {
    axis = Nn;
  } else if (_spec.mode == tilt_mode::rigid && _spec.sign == tilt_sign::propagate) {
    if (_field) {
      if (!_field->ready)
        _field->orient(_charts);
      const size_t i = static_cast<size_t>(_ids[0]);
      if (i < _field->sign.size())
        axis = _field->sign[i] * a;
    } else if (a.dot(T) < real(0)) {
      axis = -a;
    }
  } else if (_spec.mode == tilt_mode::rigid) {
    // Least-squares gradient of the metric over the ring in (T, B).
    real suu = 0, suv = 0, svv = 0, sum = 0, svm = 0;
    for (size_t k = 0; k < _d0.size(); ++k) {
      shape_w Wj;
      vec3 Tj, Bj, Nj;
      if (!vertex_chart_at(_charts, _ids[k + 1], Wj, Tj, Bj, Nj))
        continue;
      const real dm = curvature_measure(_spec.metric, Wj.k1, Wj.k2) - m;
      const real u = _d0[k].dot(T), v = _d0[k].dot(B);
      suu += u * u;
      suv += u * v;
      svv += v * v;
      sum += u * dm;
      svm += v * dm;
    }
    const real det = suu * svv - suv * suv;
    if (!(det > real(1e-12) * (suu + svv) * (suu + svv)))
      return;
    const vec3 grad = ((svv * sum - suv * svm) / det) * T + ((suu * svm - suv * sum) / det) * B;
    theta = _spec.gain * l * grad.dot(Nn.cross(a));
  }
  theta = stretch_axis_scale(W, ax, ay, theta, _spec.power);
  if (!std::isfinite(theta) || theta == real(0))
    return;
  const vec3 b = Nn.cross(a);
  for (size_t k = 0; k < _d0.size(); ++k) {
    const vec3 &d = _d0[k];
    real angle = theta;
    if (_spec.mode == tilt_mode::twist)
      angle = theta * d.dot(a) / l;
    else if (_spec.mode == tilt_mode::bend)
      angle = theta * d.dot(b) / l;
    _target[k] = Eigen::AngleAxis<real>(angle, axis) * d;
  }
}

inline void directed_tilt::project(const vecX &q, vecX &p) {
  (void)q;
  if (!_held) {
    hold_turn();
    _held = true;
  }
  for (size_t k = 0; k < _edge_weights.size() && k < _target.size(); ++k)
    p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = _w * _edge_weights[k] * _target[k];
}

inline void axis_ratio::project(const vecX &q, vecX &p) {
  if (_blocks.empty() || !_blocks[0] || _ids.empty() || _d0.empty())
    return;
  const size_t n = std::min(_edge_weights.size(), _d0.size());
  const vec3 q0 = _blocks[0]->get_vec3(_ids[0], q);
  std::vector<vec3> dq(n);
  for (size_t k = 0; k < n; ++k)
    dq[k] = _blocks[0]->get_vec3(_ids[k + 1], q) - q0;
  // Rows are scaled by w_k, so the energy weighs spoke k by w_k^2.
  auto w2 = [&](size_t k) { return _edge_weights[k] * _edge_weights[k]; };
  auto fit_rotation = [&](const std::vector<vec3> &c) {
    mat3 M = mat3::Zero();
    for (size_t k = 0; k < n; ++k)
      M += w2(k) * dq[k] * c[k].transpose();
    return polar_rotation(M);
  };
  std::vector<vec3> c(_d0.begin(), _d0.begin() + static_cast<long>(n));
  shape_w goal;
  vec3 T, B, N;
  vec3 u = vec3::Zero();
  if (vertex_chart_at(_charts, _ids[0], goal, T, B, N)) {
    N.normalize();
    if (const shape_w axes = shape_with_axes(goal); axes.ok) {
      const monge_axes ax = monge_axes_of(axes);
      u = std::cos(_angle) * (ax.hx * T + ax.hy * B) + std::sin(_angle) * (ax.lx * T + ax.ly * B);
    }
    u -= N * u.dot(N);
  }
  if (u.squaredNorm() > real(1e-16) && u.allFinite() && _ratio > real(0) && std::isfinite(_ratio)) {
    u.normalize();
    const vec3 v = N.cross(u);
    std::vector<real> al(n), be(n);
    std::vector<vec3> nrm(n);
    real den_u = 0, den_v = 0;
    for (size_t k = 0; k < n; ++k) {
      al[k] = _d0[k].dot(u);
      be[k] = _d0[k].dot(v);
      nrm[k] = _d0[k] - al[k] * u - be[k] * v;
      den_u += _edge_weights[k] * al[k] * al[k];
      den_v += _edge_weights[k] * be[k] * be[k];
    }
    real t = _ratio;
    if (_absolute)
      t = den_u > real(1e-16) && den_v > real(1e-16) ? _ratio * std::sqrt(den_v / den_u) : real(0);
    if (t > real(0) && std::isfinite(t)) {
      // c_k = cs P_k + b Q_k + n_k: P_k = U(cs = 1, b = 0) t_k, Q_k the skew.
      std::vector<vec3> P(n), Q(n);
      for (size_t k = 0; k < n; ++k) {
        P[k] = t * al[k] * u + be[k] * v;
        Q[k] = be[k] * u + al[k] * v;
      }
      const real b_max = real(0.95) * std::sqrt(t);
      real cs = 1, b = 0;
      mat3 R = fit_rotation(c);
      for (int it = 0; it < 16; ++it) {
        real pp = 0, pq = 0, qq = 0, pr = 0, qr = 0;
        for (size_t k = 0; k < n; ++k) {
          const vec3 r = R.transpose() * dq[k] - nrm[k];
          pp += w2(k) * P[k].squaredNorm();
          pq += w2(k) * P[k].dot(Q[k]);
          qq += w2(k) * Q[k].squaredNorm();
          pr += w2(k) * P[k].dot(r);
          qr += w2(k) * Q[k].dot(r);
        }
        const real det = pp * qq - pq * pq;
        if (!(det > real(1e-30) * std::max(pp * qq, real(1e-300))))
          break;
        real cs_n = (qq * pr - pq * qr) / det;
        real b_n = (pp * qr - pq * pr) / det;
        cs_n = std::max(cs_n, real(1e-3));
        b_n = std::clamp(b_n, -b_max * cs_n, b_max * cs_n);
        // Area outside the band: the fit on its edge `t c^2 - b^2 = D`,
        // `c = sqrt((D + b^2) / t)`, by golden section over `b` (the energy
        // is the same quadratic, cheap from the sums).
        const real area = t * cs_n * cs_n - b_n * b_n;
        if (!_area.contains(area)) {
          const real D = _area.clamp(area);
          auto c_of = [&](real bb) { return std::sqrt((D + bb * bb) / t); };
          auto energy = [&](real bb) {
            const real cc = c_of(bb);
            return cc * cc * pp + real(2) * cc * bb * pq + bb * bb * qq - real(2) * cc * pr -
                   real(2) * bb * qr;
          };
          // |b| <= b_max c keeps U definite: b^2 <= b_max^2 (D + b^2) / t.
          const real bm2 = b_max * b_max / t;
          real hi = bm2 < real(1) ? std::sqrt(bm2 * D / (real(1) - bm2)) : real(10) * std::sqrt(D);
          real lo = -hi;
          const real phi = real(0.5) * (std::sqrt(real(5)) - real(1));
          real x1 = hi - phi * (hi - lo), x2 = lo + phi * (hi - lo);
          real f1 = energy(x1), f2 = energy(x2);
          for (int g = 0; g < 60; ++g) {
            if (f1 < f2) {
              hi = x2;
              x2 = x1;
              f2 = f1;
              x1 = hi - phi * (hi - lo);
              f1 = energy(x1);
            } else {
              lo = x1;
              x1 = x2;
              f1 = f2;
              x2 = lo + phi * (hi - lo);
              f2 = energy(x2);
            }
          }
          b_n = real(0.5) * (lo + hi);
          cs_n = c_of(b_n);
        }
        const bool done = std::abs(cs_n - cs) + std::abs(b_n - b) < real(1e-12) * cs_n;
        cs = cs_n;
        b = b_n;
        for (size_t k = 0; k < n; ++k)
          c[k] = cs * P[k] + b * Q[k] + nrm[k];
        R = fit_rotation(c);
        if (done)
          break;
      }
    }
  }
  const mat3 R = fit_rotation(c);
  _d_rest.resize(_d0.size());
  for (size_t k = 0; k < n; ++k) {
    const vec3 spoke = R * c[k];
    _d_rest[k] = spoke;
    p.block(_id0 + static_cast<index_t>(3 * k), 0, 3, 1) = _w * _edge_weights[k] * spoke;
  }
}

inline bool monge_is_cylinder(monge_shape shape) {
  return shape == monge_shape::cylinder_minor ||
         shape == monge_shape::cylinder_major;
}

/// Circle of one principal curvature, tangent at the foot. `major` selects the
/// smaller-magnitude curvature. The other principal direction is not moved.
inline bool monge_cylinder_target(monge_shape shape, const shape_w &Wi,
                                 const vec3 &T, const vec3 &B, const vec3 &N,
                                 const vec3 &d, vec3 &target) {
  if (!monge_is_cylinder(shape) || !Wi.ok || !d.allFinite() || !T.allFinite() ||
      !B.allFinite() || !N.allFinite())
    return false;
  const shape_w A = shape_with_axes(Wi);
  if (!A.ok)
    return false;
  const monge_axes axes = monge_axes_of(A);
  const bool major = shape == monge_shape::cylinder_major;
  const real k = major ? axes.k_lo : axes.k_hi;
  const real ex = major ? axes.lx : axes.hx;
  const real ey = major ? axes.ly : axes.hy;
  if (!std::isfinite(k) || !(std::abs(k) > real(1e-8)))
    return false;
  vec3 tube = ex * T + ey * B;
  tube -= N * tube.dot(N);
  if (!(tube.squaredNorm() > real(1e-16)))
    return false;
  tube.normalize();
  vec3 axis = N.cross(tube);
  if (!(axis.squaredNorm() > real(1e-16)))
    return false;
  axis.normalize();
  const real hn = d.dot(N);
  const vec3 tang = d - hn * N;
  real uo = 0;
  real ho = 0;
  if (!hopf_curvature_land(k, tang.dot(tube), hn, uo, ho))
    return false;
  target = tang.dot(axis) * axis + uo * tube + ho * N;
  return target.allFinite();
}

inline real monge_prin_quad(const monge_axes &a, real tx, real ty) {
  const real u = tx * a.hx + ty * a.hy;
  const real v = tx * a.lx + ty * a.ly;
  return a.k_hi * u * u + a.k_lo * v * v;
}

inline real monge_jet_height(monge_shape shape, const shape_w &Wi,
                             const shape_w *Wji, real tx, real ty) {
  const monge_axes ai =
      specialize_monge_axes(monge_axes_of(shape_with_axes(Wi)), shape);
  const real h2 = real(0.5) * monge_prin_quad(ai, tx, ty);
  real h3 = 0;
  if (Wji != nullptr && Wji->ok) {
    const monge_axes aj =
        specialize_monge_axes(monge_axes_of(shape_with_axes(*Wji)), shape);
    h3 = real(1.0 / 6.0) *
         (monge_prin_quad(aj, tx, ty) - monge_prin_quad(ai, tx, ty));
  }
  real h = h2;
  if (std::isfinite(h3))
    h += h3;
  return std::isfinite(h) ? h : std::numeric_limits<real>::quiet_NaN();
}

/// `inner` is the equator at distance `|R − r|`, section `S = R − s`.
/// The outer equator is `S = R + s`. NaN past the silhouette.
inline real monge_torus_section_height(const albers::osculating_torus &tor,
                                       real tx, real ty, const vec3 &T,
                                       const vec3 &B, const vec3 &N,
                                       bool inner) {
  if (!tor.valid || !(tor.r > real(1e-12)) || !tor.axis.allFinite() ||
      !T.allFinite() || !B.allFinite() || !N.allFinite())
    return std::numeric_limits<real>::quiet_NaN();
  if (inner && !(tor.R > tor.r + real(1e-8)))
    return std::numeric_limits<real>::quiet_NaN();
  const vec3 revolution = tor.axis.normalized();
  const vec3 p = tx * T + ty * B;
  // `tor.axis` is the revolution axis \(n\times e_{\mathrm{lo}}\). The tube
  // direction is that axis. The major tangent is axis × N.
  vec3 e_major = revolution.cross(N);
  if (e_major.squaredNorm() < real(1e-16))
    return std::numeric_limits<real>::quiet_NaN();
  e_major.normalize();
  const real t_tube = p.dot(revolution);
  const real t_major = p.dot(e_major);
  if (!(std::abs(t_tube) < tor.r))
    return std::numeric_limits<real>::quiet_NaN();
  const real s =
      std::sqrt(std::max(real(0), tor.r * tor.r - t_tube * t_tube));
  const real S = inner ? (tor.R - s) : (tor.R + s);
  if (inner && !(S > real(0)))
    return std::numeric_limits<real>::quiet_NaN();
  const real under = S * S - t_major * t_major;
  if (!(under > real(0)))
    return std::numeric_limits<real>::quiet_NaN();
  const real side = static_cast<real>(tor.sign);
  const real dist = inner ? (tor.R - tor.r) : (tor.R + tor.r);
  const real root = std::sqrt(under);
  const real h = inner ? side * (root - dist) : side * (dist - root);
  return std::isfinite(h) ? h : std::numeric_limits<real>::quiet_NaN();
}

inline albers::osculating_torus monge_torus_of(const shape_w &Wi, const vec3 &T,
                                              const vec3 &B, const vec3 &N) {
  const mat3 Wm = Wi.k11 * T * T.transpose() +
                  Wi.k12 * (T * B.transpose() + B * T.transpose()) +
                  Wi.k22 * B * B.transpose();
  return albers::osculating_torus_from_shape_operator(vec3::Zero(), N, Wm);
}

/// Circular section of the osculating torus, as a height over the tangent
/// plane. The 2-jet matches ½ tᵀ W t. Same-sign principal curvatures are the
/// outer equator (`|κ_lo| = 1/(R+r)`); opposite signs are the inner equator
/// (`|κ_lo| = 1/(R−r)`). Past the silhouette, fall back to that jet.
inline real monge_torus_height(const shape_w &Wi, real tx, real ty,
                               const vec3 &T, const vec3 &B, const vec3 &N) {
  const real quad = real(0.5) * shape_quad(Wi, tx, ty);
  if (!Wi.ok || !T.allFinite() || !B.allFinite() || !N.allFinite())
    return std::numeric_limits<real>::quiet_NaN();
  if (!(std::abs(Wi.k11) + std::abs(Wi.k12) + std::abs(Wi.k22) > real(1e-8)))
    return quad;
  const shape_w A = shape_with_axes(Wi);
  const bool inner = A.ok && (A.k1 * A.k2 < real(0));
  const albers::osculating_torus tor = monge_torus_of(Wi, T, B, N);
  const real h = monge_torus_section_height(tor, tx, ty, T, B, N, inner);
  return std::isfinite(h) ? h : quad;
}

inline real monge_fairing_height(monge_shape shape, const shape_w &Wi, real tx,
                                 real ty, const shape_w *Wji, const vec3 &T,
                                 const vec3 &B, const vec3 &N) {
  if (!Wi.ok)
    return std::numeric_limits<real>::quiet_NaN();
  if (shape == monge_shape::quadric) {
    const real h2 = real(0.5) * shape_quad(Wi, tx, ty);
    real h3 = 0;
    if (Wji != nullptr && Wji->ok)
      h3 = real(1.0 / 6.0) * (shape_quad(*Wji, tx, ty) - shape_quad(Wi, tx, ty));
    real h = h2;
    if (std::isfinite(h3))
      h += h3;
    return std::isfinite(h) ? h : std::numeric_limits<real>::quiet_NaN();
  }
  if (shape == monge_shape::sphere) {
    const shape_w A = shape_with_axes(Wi);
    const real H = real(0.5) * (A.k1 + A.k2);
    const real h = real(0.5) * H * (tx * tx + ty * ty);
    return std::isfinite(h) ? h : std::numeric_limits<real>::quiet_NaN();
  }
  if (shape == monge_shape::osculating_torus)
    return monge_torus_height(Wi, tx, ty, T, B, N);
  if (monge_is_cylinder(shape)) {
    if (!std::isfinite(tx) || !std::isfinite(ty))
      return std::numeric_limits<real>::quiet_NaN();
    vec3 hit;
    if (!monge_cylinder_target(shape, Wi, T, B, N, tx * T + ty * B, hit))
      return real(0);
    const real h = hit.dot(N);
    return std::isfinite(h) ? h : std::numeric_limits<real>::quiet_NaN();
  }
  return monge_jet_height(shape, Wi, Wji, tx, ty);
}

/// The graph is the osculating chart. It stops once the spoke is longer than
/// the tighter principal radius: past that, ½ κ t² is not a surface point.
inline bool monge_graph_in_chart(const shape_w &W, real tlen) {
  if (!std::isfinite(tlen) || tlen < real(0))
    return false;
  const shape_w A = shape_with_axes(W);
  if (!A.ok)
    return false;
  const real k = std::max(std::abs(A.k1), std::abs(A.k2));
  if (!std::isfinite(k))
    return false;
  return tlen * k <= real(1);
}

/// Nearest listed landing. Cylinder objectives are the circular tube. The
/// others are graphs, so the spoke keeps its tangential part, and only while
/// the spoke sits in the osculating chart. An empty list is the unrestricted
/// quadric.
inline bool select_monge_spoke(const std::vector<objective> &objectives,
                              const vec3 &d, const shape_w &Wi,
                              const shape_w *Wji, const vec3 &T, const vec3 &B,
                              const vec3 &N, vec3 &target) {
  if (!d.allFinite() || !N.allFinite() || !T.allFinite() || !B.allFinite())
    return false;
  const real hn = d.dot(N);
  const vec3 tang = d - hn * N;
  std::vector<monge_shape> shapes;
  shapes.reserve(objectives.size());
  for (const objective &o : objectives) {
    if (const auto *m = std::get_if<monge_fairing>(&o))
      shapes.push_back(m->shape);
  }
  if (shapes.empty())
    shapes.push_back(monge_shape::quadric);
  vec3 best = d;
  real best_abs = std::numeric_limits<real>::infinity();
  bool any = false;
  for (const monge_shape shape : shapes) {
    vec3 hit;
    if (monge_is_cylinder(shape)) {
      if (!monge_cylinder_target(shape, Wi, T, B, N, d, hit))
        continue;
    } else {
      const real tlen = tang.norm();
      if (!monge_graph_in_chart(Wi, tlen))
        continue;
      if (Wji != nullptr && Wji->ok && !monge_graph_in_chart(*Wji, tlen))
        continue;
      const real h = monge_fairing_height(shape, Wi, tang.dot(T), tang.dot(B),
                                          Wji, T, B, N);
      if (!std::isfinite(h))
        continue;
      hit = tang + h * N;
    }
    if (!hit.allFinite())
      continue;
    const real err = (hit - d).norm();
    if (!any || err < best_abs) {
      best = hit;
      best_abs = err;
      any = true;
    }
  }
  if (!any)
    return false;
  target = best;
  return true;
}

/// Among the enabled graphs, the one closest to the current spoke height.
/// An empty list is the unrestricted quadric.
inline real select_monge_height(const std::vector<objective> &objectives,
                                real h_now, const shape_w &Wi, real tx, real ty,
                                const shape_w *Wji, const vec3 &T, const vec3 &B,
                                const vec3 &N) {
  std::vector<monge_shape> shapes;
  shapes.reserve(objectives.size());
  for (const objective &o : objectives) {
    if (const auto *m = std::get_if<monge_fairing>(&o))
      shapes.push_back(m->shape);
  }
  if (shapes.empty())
    shapes.push_back(monge_shape::quadric);
  real best = 0;
  real best_abs = std::numeric_limits<real>::infinity();
  bool any = false;
  for (const monge_shape shape : shapes) {
    const real h = monge_fairing_height(shape, Wi, tx, ty, Wji, T, B, N);
    if (!std::isfinite(h))
      continue;
    const real err = std::abs(h - h_now);
    if (!any || err < best_abs) {
      best = h;
      best_abs = err;
      any = true;
    }
  }
  return any ? best : std::numeric_limits<real>::quiet_NaN();
}

/// One-ring nodes for the connection eigenproblem, circulation order.
/// `spoke_z` is the foot spin. `link_z` is the spin to the next kept neighbor.
inline std::vector<hopf_ring_node>
hopf_ring_nodes_of(const quat &ui_in, const std::vector<vec3> &ts,
                   const std::vector<quat> &qs) {
  std::vector<hopf_ring_node> nodes;
  const int n = static_cast<int>(std::min(ts.size(), qs.size()));
  if (!frame_quat_ok(ui_in) || n < 3)
    return nodes;
  const quat ui = ui_in.normalized();
  struct kept {
    vec3 t = vec3::Zero();
    quat q = quat::Identity();
  };
  std::vector<kept> rows;
  rows.reserve(static_cast<size_t>(n));
  for (int k = 0; k < n; ++k) {
    if (!ts[static_cast<size_t>(k)].allFinite() ||
        !(ts[static_cast<size_t>(k)].squaredNorm() > real(1e-16)))
      continue;
    quat qj = qs[static_cast<size_t>(k)];
    if (!frame_quat_ok(qj))
      continue;
    qj.normalize();
    rows.push_back(kept{ts[static_cast<size_t>(k)], qj});
  }
  if (rows.size() < 3)
    return nodes;
  nodes.resize(rows.size());
  for (size_t k = 0; k < rows.size(); ++k) {
    const frame_body_omega spoke = frame_body_omega_of(ui, rows[k].q);
    const size_t next = (k + 1) % rows.size();
    const frame_body_omega link = frame_body_omega_of(rows[k].q, rows[next].q);
    nodes[k].t = rows[k].t;
    nodes[k].spoke_z = spoke.ok ? spoke.z : real(0);
    nodes[k].link_z = link.ok ? link.z : real(0);
  }
  return nodes;
}

/// Neighbor fan passed into `frame_fair::create`. The constructor bakes the
/// Stein fan and keeps the spoke ids. These face lists are not members.
struct nbr_fan {
  index_t j = -1;
  std::vector<fan_face> faces;
  std::vector<fan_face> support;
  std::vector<fan_edge> edges;
  std::vector<int> tangent;
  std::vector<index_t> spokes;
  std::vector<real> spoke_k;
  std::vector<char> spoke_ok;
};

/// Neighbor ids owned by one fairing row.
struct nbr_ids {
  index_t j = -1;
  fan_local local;
  std::vector<index_t> spokes;
  std::vector<real> spoke_k;
  std::vector<char> spoke_ok;
};

/// One cotangent row per vertex, the same stencil as `laplacian` and `bending`.
/// The 3-vector is Σ K * target. A Monge target is the graph edge t + h N.
/// An Albers target is the spoke's closest point on that surface. `world_x`
/// lands on the hyperboloid. `hopf_line` is the 1D parabola along the phase
/// geodesic. `conformal_hopf` is the osculating circle in that same plane.
/// Positions are block 0. A Stein row uses block 1 for edge coefficients and
/// block 2 for the shape operator. A quaternion row uses block 1 for frames.
class frame_fair : public block_constraint {
public:
  DEFINE_CREATE_FUNC(frame_fair)

  frame_fair(std::vector<index_t> ids, std::vector<real> weights,
             std::vector<fan_face> faces, std::vector<fan_face> support,
             std::vector<fan_edge> edges, std::vector<int> tangent,
             std::vector<nbr_fan> nbrs, std::vector<objective> objectives,
             bool stein, int nE, int debug_vert, const real &w,
             std::vector<sim_block::ptr> blocks,
             const std::vector<vertex_chart> *charts = nullptr,
             std::vector<index_t> fit_ids = {}, real spine_step = real(1),
             const kusama::weingarten::vertex_stars *stars = nullptr)
      : block_constraint(std::move(ids), w, std::move(blocks)),
        _edge_weights(std::move(weights)),
        _objectives(std::move(objectives)), _stein(stein), _nE(nE),
        _debug_vert(debug_vert), _charts(charts), _stars(stars),
        _fit_ids(std::move(fit_ids)), _spine_step(spine_step) {
    const int vi = _ids.empty() ? -1 : _ids[0];
    if (_stein)
      _fan = bake_fan(faces, support, edges, tangent, vi);
    _nbrs.reserve(nbrs.size());
    for (nbr_fan &nb : nbrs) {
      nbr_ids row;
      row.j = nb.j;
      if (_stein)
        row.local = bake_fan(nb.faces, nb.support, nb.edges, nb.tangent, nb.j);
      row.spokes = std::move(nb.spokes);
      row.spoke_k = std::move(nb.spoke_k);
      row.spoke_ok = std::move(nb.spoke_ok);
      _nbrs.push_back(std::move(row));
    }
    _rest_initialized = true;
  }

  virtual std::string name() { return "frame_fair"; }

  virtual void project(const vecX &q, vecX &p) {
    if (_blocks.empty() || _ids.empty())
      return;
    const index_t vi = _ids[0];
    const int n = static_cast<int>(_edge_weights.size());
    std::vector<index_t> spokes;
    spokes.reserve(static_cast<size_t>(std::max(n, 0)));
    for (int k = 0; k < n && k + 1 < static_cast<int>(_ids.size()); ++k)
      spokes.push_back(_ids[static_cast<size_t>(k + 1)]);
    const std::vector<char> spoke_ok(spokes.size(), 1);
    const fair_ring ring =
        fair_ring_of(vi, spokes, _edge_weights, spoke_ok);

    auto pos_at = [&](int i) { return _blocks[0]->get_vec3(i, q); };
    auto coeff_at = [&](int i) -> real {
      if (_blocks.size() < 2)
        return real(0);
      const index_t off = _blocks[1]->get_offset_idx(0);
      const index_t k = off + i;
      if (k < 0 || k >= q.size())
        return real(0);
      return q(k);
    };

    const bool quat_chart = !_stein && _blocks.size() >= 2;
    const spine_set spines = spine_set_of(_objectives);

    ring_frame rf;
    if (_charts) {
      rf.axes =
          vertex_row_frame(_blocks, _charts, vi, q, rf.W, rf.T, rf.B, rf.N, _stars);
    } else if (_stein) {
      const fan_chart ev = evaluate_fan(_fan, _nE, pos_at, coeff_at, &_lock);
      rf.W = ev.W;
      rf.T = ev.T;
      rf.B = ev.B;
      rf.N = ev.N;
      rf.axes = ev.axes;
    } else if (quat_chart) {
      rf = fit_ring_frame(*_blocks[0], *_blocks[1], q, ring);
    }
    const bool chart = rf.axes && rf.W.ok;

    struct held {
      index_t j = -1;
      shape_w stored;
      vec3 T = vec3::Zero();
      vec3 B = vec3::Zero();
      vec3 N = vec3::Zero();
      bool axes = false;
    };
    std::vector<held> helds;
    const bool have_shape = _stein && _blocks.size() >= 3 && !_charts;
    if (_charts || have_shape || quat_chart) {
      helds.reserve(_nbrs.size());
      for (const nbr_ids &nb : _nbrs) {
        held h;
        h.j = nb.j;
        if (_charts) {
          h.axes = vertex_row_frame(_blocks, _charts, nb.j, q, h.stored, h.T, h.B,
                                    h.N, _stars);
        } else if (_stein) {
          const fan_chart ev =
              evaluate_fan(nb.local, _nE, pos_at, coeff_at, nullptr);
          h.T = ev.T;
          h.B = ev.B;
          h.N = ev.N;
          h.axes = ev.axes;
        } else if (quat_chart) {
          const fair_ring nr =
              fair_ring_of(nb.j, nb.spokes, nb.spoke_k, nb.spoke_ok);
          const ring_frame nrf =
              fit_ring_frame(*_blocks[0], *_blocks[1], q, nr);
          h.T = nrf.T;
          h.B = nrf.B;
          h.N = nrf.N;
          h.axes = nrf.axes;
        }
        if (have_shape)
          h.stored = shape_w_from_coeffs(_blocks[2]->get_vec3(nb.j, q));
        helds.push_back(h);
      }
    }
    auto held_of = [&](index_t j) -> const held * {
      for (const held &h : helds) {
        if (h.j == j)
          return &h;
      }
      return nullptr;
    };

    const vec3 xi = _blocks[0]->get_vec3(vi, q);

    spine_landing spine_land;
    if (spines.any() && chart) {
      shape_w goal = rf.W;
      const shape_w axes = shape_with_axes(goal);
      if (axes.ok) {
        const monge_axes ax = monge_axes_of(axes);
        spine_curvature sc;
        sc.tube = ax.hx * rf.T + ax.hy * rf.B;
        sc.ruling = ax.lx * rf.T + ax.ly * rf.B;
        sc.normal = rf.N;
        sc.k_hi = ax.k_hi;
        sc.k_lo = ax.k_lo;
        // Goal S → `sc`. Rest centers from canal fit on the 2-ring; land
        // the 1-ring with goal κ / axis.
        std::vector<vec3> dproject;
        dproject.reserve(spokes.size());
        for (const index_t j : spokes)
          dproject.push_back(_blocks[0]->get_vec3(j, q) - xi);
        std::vector<vec3> dfit;
        if (!_fit_ids.empty()) {
          dfit.reserve(_fit_ids.size());
          for (const index_t j : _fit_ids)
            dfit.push_back(_blocks[0]->get_vec3(j, q) - xi);
        }
        spine_land =
            choose_spine_shapes(sc, dproject, spines, dfit, _spine_step);
      }
    }

    const bool dump = _debug_vert >= 0 && vi == _debug_vert && !_debug_printed;

    // Per-spoke rows: a single cotan sum is blind to isotropic canal seating
    // (centered ring ⇒ Σ K d = Σ K hit = 0 even when every spoke is off).
    index_t row = _id0;
    for (int i = 0; i < n; ++i) {
      const index_t j = spokes[static_cast<size_t>(i)];
      const real K = _edge_weights[static_cast<size_t>(i)];
      const vec3 d = _blocks[0]->get_vec3(j, q) - xi;
      vec3 target = d;
      if (d.allFinite() && spines.any() &&
          spine_land.fit != spine_kind::none) {
        if (static_cast<size_t>(i) < spine_land.ok.size() &&
            spine_land.ok[static_cast<size_t>(i)])
          target = spine_land.hit[static_cast<size_t>(i)];
      } else if (d.allFinite() && chart) {
        shape_w Wi = rf.W;
        shape_w carried;
        const shape_w *Wj = nullptr;
        if (const held *nb = held_of(j)) {
          if (nb->axes && nb->stored.ok) {
            carried = transport_shape(nb->stored, nb->T, nb->B, rf.T, rf.B, rf.N);
            if (carried.ok)
              Wj = &carried;
          }
        }
        vec3 hit;
        if (select_monge_spoke(_objectives, d, Wi, Wj, rf.T, rf.B, rf.N, hit))
          target = hit;
      }
      vec3 rhs = vec3::Zero();
      if (d.allFinite())
        rhs = d;
      if (target.allFinite())
        rhs = target;
      vec3 rowp = _w * K * rhs;
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
    if (dump)
      _debug_printed = true;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    const int n = static_cast<int>(_edge_weights.size());
    for (int k = 0; k < n && k + 1 < static_cast<int>(_ids.size()); ++k) {
      const index_t jj =
          _blocks[0]->get_offset_idx(_ids[static_cast<size_t>(k + 1)]);
      const real K = _edge_weights[static_cast<size_t>(k)];
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, ii + ax, -_w * K));
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, jj + ax, _w * K));
      id0 += 3;
    }
  }

  std::vector<real> _edge_weights;
  fan_local _fan;
  std::vector<nbr_ids> _nbrs;
  std::vector<objective> _objectives;
  vec3 _lock = vec3::Zero();
  bool _stein = false;
  int _nE = 0;
  int _debug_vert = -1;
  bool _debug_printed = false;
  const std::vector<vertex_chart> *_charts = nullptr;
  /// Live `S` from the edge store (`blocks[1]`); null keeps the chart's.
  const kusama::weingarten::vertex_stars *_stars = nullptr;
  /// 2-ring vertex ids for osculating canal fit (exclude the center). Empty
  /// falls back to the 1-ring project spokes inside `choose_spine_shapes`.
  std::vector<index_t> _fit_ids;
  /// Canal-radial residual bite in [0, 1]. `1` lands fully each project.
  real _spine_step = real(1);
};

/// Closest point on a Calder fit of one vertex ring. The samples are the
/// spoke positions. Each sample normal is the fan of the two triangles that
/// meet that spoke. The only block is positions.
class albers_ring_fair : public block_constraint {
public:
  DEFINE_CREATE_FUNC(albers_ring_fair)

  albers_ring_fair(std::vector<index_t> ids, std::vector<real> weights,
                   std::vector<objective> objectives, const real &w,
                   std::vector<sim_block::ptr> blocks)
      : block_constraint(std::move(ids), w, std::move(blocks)),
        _edge_weights(std::move(weights)), _objectives(std::move(objectives)) {
    _rest_initialized = true;
  }

  virtual std::string name() { return "albers_ring_fair"; }

  virtual void project(const vecX &q, vecX &p) {
    if (_blocks.empty() || !_blocks[0] || _ids.size() < 2)
      return;
    const albers_mask am = albers_mask_of(_objectives);
    const index_t vi = _ids[0];
    const int n = static_cast<int>(_edge_weights.size());
    std::vector<index_t> spokes;
    spokes.reserve(static_cast<size_t>(std::max(n, 0)));
    for (int k = 0; k < n && k + 1 < static_cast<int>(_ids.size()); ++k)
      spokes.push_back(_ids[static_cast<size_t>(k + 1)]);
    const vec3 xi = _blocks[0]->get_vec3(vi, q);
    const vec3 foot_N = fan_normal(q, xi, spokes);

    real length_sum = 0;
    int length_n = 0;
    std::vector<albers_fit_sample> samples;
    samples.reserve(spokes.size());
    for (int k = 0; k < static_cast<int>(spokes.size()); ++k) {
      const vec3 dp = _blocks[0]->get_vec3(spokes[static_cast<size_t>(k)], q) - xi;
      if (!dp.allFinite() || !(dp.norm() > real(1e-8)))
        continue;
      length_sum += dp.norm();
      ++length_n;
      const vec3 Ns = spoke_normal(q, xi, spokes, k);
      if (!(Ns.squaredNorm() > real(1e-16)))
        continue;
      albers_fit_sample s;
      s.dp = dp;
      s.N = Ns;
      s.w = 1;
      samples.push_back(s);
    }
    const real l0 =
        length_n > 0 ? length_sum / static_cast<real>(length_n) : real(0);

    albers::vec10 quadric_Q;
    bool quadric_ok = false;
    vec3 sphere_c = vec3::Zero();
    real sphere_r = 0;
    bool sphere_ok = false;
    albers::vec14 cyclide_Q;
    bool cyclide_ok = false;
    if (am.any() && foot_N.squaredNorm() > real(1e-16) &&
        static_cast<int>(samples.size()) >= 4) {
      if (am.quadric)
        quadric_ok = fit_albers_quadric(samples, l0, foot_N, quadric_Q);
      if (am.sphere)
        sphere_ok = fit_albers_sphere(samples, l0, foot_N, sphere_c, sphere_r);
      if (am.cyclide)
        cyclide_ok = fit_albers_cyclide(samples, l0, foot_N, cyclide_Q);
    }

    // Same as frame_fair: closest-point landings need one row per spoke.
    index_t row = _id0;
    for (int i = 0; i < n; ++i) {
      const index_t j = spokes[static_cast<size_t>(i)];
      const real K = _edge_weights[static_cast<size_t>(i)];
      const vec3 d = _blocks[0]->get_vec3(j, q) - xi;
      vec3 target = d;
      if (d.allFinite() && am.any()) {
        vec3 best = d;
        real best_err = std::numeric_limits<real>::infinity();
        bool any = false;
        auto consider = [&](const vec3 &hit) {
          if (!hit.allFinite())
            return;
          const real err = (hit - d).norm();
          if (!any || err < best_err) {
            best = hit;
            best_err = err;
            any = true;
          }
        };
        vec3 hit = d;
        if (quadric_ok && land_albers_quadric(quadric_Q, d, hit))
          consider(hit);
        if (sphere_ok && land_albers_sphere(sphere_c, sphere_r, d, hit))
          consider(hit);
        if (cyclide_ok && land_albers_cyclide(cyclide_Q, d, hit))
          consider(hit);
        if (any)
          target = best;
      }
      vec3 rowp = vec3::Zero();
      if (d.allFinite() && target.allFinite())
        rowp = _w * K * target;
      if (!rowp.allFinite())
        rowp = vec3::Zero();
      p.block(row, 0, 3, 1) = rowp;
      row += 3;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    if (_blocks.empty() || !_blocks[0] || _ids.empty())
      return;
    const index_t ii = _blocks[0]->get_offset_idx(_ids[0]);
    const int n = static_cast<int>(_edge_weights.size());
    for (int k = 0; k < n && k + 1 < static_cast<int>(_ids.size()); ++k) {
      const index_t jj =
          _blocks[0]->get_offset_idx(_ids[static_cast<size_t>(k + 1)]);
      const real K = _edge_weights[static_cast<size_t>(k)];
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, ii + ax, -_w * K));
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(id0 + ax, jj + ax, _w * K));
      id0 += 3;
    }
  }

private:
  vec3 fan_normal(const vecX &q, const vec3 &xi,
                  const std::vector<index_t> &spokes) const {
    vec3 N = vec3::Zero();
    const int n = static_cast<int>(spokes.size());
    for (int k = 0; k < n; ++k) {
      const vec3 a = _blocks[0]->get_vec3(spokes[static_cast<size_t>(k)], q) - xi;
      const vec3 b =
          _blocks[0]->get_vec3(spokes[static_cast<size_t>((k + 1) % n)], q) - xi;
      const vec3 cr = a.cross(b);
      if (cr.allFinite())
        N += cr;
    }
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
      return vec3::Zero();
    return N.normalized();
  }

  vec3 spoke_normal(const vecX &q, const vec3 &xi,
                    const std::vector<index_t> &spokes, int k) const {
    const int n = static_cast<int>(spokes.size());
    if (n < 2)
      return vec3::Zero();
    auto at = [&](int i) {
      const int s = (i % n + n) % n;
      return _blocks[0]->get_vec3(spokes[static_cast<size_t>(s)], q);
    };
    const vec3 a = at(k - 1) - xi;
    const vec3 b = at(k) - xi;
    const vec3 c = at(k + 1) - xi;
    const vec3 N = a.cross(b) + b.cross(c);
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
      return vec3::Zero();
    return N.normalized();
  }

  std::vector<real> _edge_weights;
  std::vector<objective> _objectives;
};

/// Same spokes as frame_fair. Each row is the height of one spoke along the
/// one-ring normal. That height moves by one Willmore step. The tangent of
/// the spoke is not a row. `_rings` is this vertex plus the rings its
/// Laplacian reads. Only this vertex's spokes are written.
class frame_minimizer : public block_constraint {
public:
  DEFINE_CREATE_FUNC(frame_minimizer)

  frame_minimizer(std::vector<fair_ring> rings, index_t self, const real &w,
                  std::vector<sim_block::ptr> blocks)
      : block_constraint({}, w, std::move(blocks)), _rings(std::move(rings)),
        _self(self) {
    _rest_initialized = true;
  }

  virtual std::string name() { return "frame_minimizer"; }

  virtual void project(const vecX &q, vecX &p) {
    if (_blocks.size() < 2)
      return;
    const int nr = static_cast<int>(_rings.size());
    std::vector<real> H(static_cast<size_t>(nr), 0);
    std::vector<real> Kc(static_cast<size_t>(nr), 0);
    std::vector<real> A(static_cast<size_t>(nr), 0);
    std::vector<char> ok(static_cast<size_t>(nr), 0);
    int max_id = 0;
    for (const fair_ring &ring : _rings) {
      max_id = std::max(max_id, ring.i);
      for (const fair_spoke &sp : ring.spokes)
        max_id = std::max(max_id, sp.j);
    }
    std::vector<int> which(static_cast<size_t>(max_id + 1), -1);
    for (int r = 0; r < nr; ++r)
      which[static_cast<size_t>(_rings[static_cast<size_t>(r)].i)] = r;

    for (int r = 0; r < nr; ++r) {
      const fair_ring &ring = _rings[static_cast<size_t>(r)];
      const ring_frame rf = fit_ring_frame(*_blocks[0], *_blocks[1], q, ring);
      A[static_cast<size_t>(r)] = ring_vertex_area(*_blocks[0], q, ring);
      if (!rf.axes || !rf.W.ok)
        continue;
      const real Hr = real(0.5) * (rf.W.k1 + rf.W.k2);
      const real Kr = rf.W.k1 * rf.W.k2;
      if (!std::isfinite(Hr) || !std::isfinite(Kr))
        continue;
      H[static_cast<size_t>(r)] = Hr;
      Kc[static_cast<size_t>(r)] = Kr;
      ok[static_cast<size_t>(r)] = 1;
    }

    std::vector<real> delta(static_cast<size_t>(nr), 0);
    for (int r = 0; r < nr; ++r) {
      if (!ok[static_cast<size_t>(r)])
        continue;
      const real lap = cotan_diff_H(_rings[static_cast<size_t>(r)],
                                   H[static_cast<size_t>(r)], H, ok, which);
      const real area = A[static_cast<size_t>(r)];
      const real lapH = area > real(1e-12) ? lap / (real(2) * area) : 0;
      const real edge = ring_mean_edge(*_blocks[0], q, _rings[static_cast<size_t>(r)]);
      delta[static_cast<size_t>(r)] = willmore_normal_step(
          H[static_cast<size_t>(r)], Kc[static_cast<size_t>(r)], lapH, edge);
    }

    index_t row = _id0;
    for (int r = 0; r < nr; ++r) {
      const fair_ring &ring = _rings[static_cast<size_t>(r)];
      if (ring.i != _self)
        continue;
      const vec3 xi = _blocks[0]->get_vec3(ring.i, q);
      const real step = delta[static_cast<size_t>(r)];
      int s = 0;
      for (const fair_spoke &sp : ring.spokes) {
        if (s >= static_cast<int>(_spoke_normals.size()))
          break;
        const vec3 N = _spoke_normals[static_cast<size_t>(s)];
        ++s;
        if (!(N.squaredNorm() > real(1e-16)))
          continue;
        const vec3 d = _blocks[0]->get_vec3(sp.j, q) - xi;
        real h = 0;
        if (d.allFinite())
          h = d.dot(N) - step;
        if (!std::isfinite(h))
          h = 0;
        const real rhs = _w * sp.K * h;
        p(row) = std::isfinite(rhs) ? rhs : real(0);
        ++row;
      }
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    _spoke_normals.clear();
    if (_blocks.empty())
      return;
    for (const fair_ring &ring : _rings) {
      if (ring.i != _self)
        continue;
      const vec3 N = one_ring_normal(ring);
      const index_t ii = _blocks[0]->get_offset_idx(ring.i);
      for (const fair_spoke &sp : ring.spokes) {
        if (!std::isfinite(sp.K) || !(N.squaredNorm() > real(1e-16))) {
          _spoke_normals.push_back(vec3::Zero());
          continue;
        }
        const index_t jj = _blocks[0]->get_offset_idx(sp.j);
        for (int ax = 0; ax < 3; ++ax) {
          triplets.push_back(trip(id0, ii + ax, -_w * sp.K * N[ax]));
          triplets.push_back(trip(id0, jj + ax, _w * sp.K * N[ax]));
        }
        _spoke_normals.push_back(N);
        ++id0;
      }
    }
  }

  vec3 one_ring_normal(const fair_ring &ring) const {
    const auto *pos = dynamic_cast<const vec3_block *>(_blocks[0].get());
    if (pos == nullptr)
      return vec3::Zero();
    const std::vector<vec3> &x = pos->_x;
    if (ring.i < 0 || static_cast<size_t>(ring.i) >= x.size())
      return vec3::Zero();
    const vec3 xi = x[static_cast<size_t>(ring.i)];
    vec3 N = vec3::Zero();
    const int n = static_cast<int>(ring.spokes.size());
    for (int k = 0; k < n; ++k) {
      const int j = ring.spokes[static_cast<size_t>(k)].j;
      const int j1 = ring.spokes[static_cast<size_t>((k + 1) % n)].j;
      if (j < 0 || j1 < 0 || static_cast<size_t>(j) >= x.size() ||
          static_cast<size_t>(j1) >= x.size())
        continue;
      const vec3 cr =
          (x[static_cast<size_t>(j)] - xi).cross(x[static_cast<size_t>(j1)] - xi);
      if (cr.allFinite())
        N += cr;
    }
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-16)))
      return vec3::Zero();
    return N.normalized();
  }

  std::vector<fair_ring> _rings;
  std::vector<vec3> _spoke_normals;
  index_t _self = -1;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif