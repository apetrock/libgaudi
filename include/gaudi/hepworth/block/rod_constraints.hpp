
#ifndef __HEP_ROD_BLOCK_CONSTRAINTS__
#define __HEP_ROD_BLOCK_CONSTRAINTS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>

#include <iostream>
#include <limits>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"

#include "../projection_constraint.hpp"
#include "block_constraint.hpp"
#include "sim_block.hpp"
#include "gaudi/geometry_logger.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

class smooth : public block_constraint {
public:
  typedef std::shared_ptr<smooth> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<smooth>(ids, w, blocks);
  }

  smooth(const std::vector<index_t> &ids, const real &w,
         std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}
  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    index_t i0 = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t in = this->_ids[2];

    vec3 qm = _blocks[0]->get_vec3(ip, q);
    vec3 q0 = _blocks[0]->get_vec3(i0, q);
    vec3 qp = _blocks[0]->get_vec3(in, q);
    vec3 dqm = qm - q0;
    vec3 dqp = qp - q0;
    vec3 N = -(dqm + dqp).normalized();
    dqm = va::reject(N, dqm);
    dqp = va::reject(N, dqp);
    p.block(_id0, 0, 3, 1) = _w * (dqm + dqp);
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i0 = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t ip = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t in = _blocks[0]->get_offset_idx(this->_ids[2]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i0 + ax, -2.0 * _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, ip + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, in + ax, _w));
    id0 += 3;
  }
};

// Centered 3rd difference → 0 (min-kink regularizer).
// ids: {i-2, i-1, i, i+1, i+2}; A uses integer stencil (-1, 2, 0, -2, 1).
class min_kink : public block_constraint {
public:
  typedef std::shared_ptr<min_kink> ptr;
  static constexpr int kStencil = 5;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<min_kink>(ids, w, blocks);
  }

  min_kink(const std::vector<index_t> &ids, const real &w,
           std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    assert(ids.size() >= kStencil);
  }
  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX & /*q*/, vecX &p) {
    p.block(_id0, 0, 3, 1).setZero();
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    // δ³ ~ -x_{-2} + 2 x_{-1} - 2 x_{+1} + x_{+2}  (center coeff 0)
    static const real kCoeff[kStencil] = {-1.0, 2.0, 0.0, -2.0, 1.0};
    for (int k = 0; k < kStencil; ++k) {
      if (kCoeff[k] == 0.0)
        continue;
      const index_t ik = _blocks[0]->get_offset_idx(this->_ids[k]);
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(_id0 + ax, ik + ax, kCoeff[k] * _w));
    }
    id0 += 3;
  }
};

// Spherical cubic Hermite / Bezier on S³ from endpoint jets; value-pull middle.
// ids: {i-2, i-1, i, i+1, i+2}. Curve uses ends + log-relatives; mid is target only.
class squad_smooth : public block_constraint {
public:
  typedef std::shared_ptr<squad_smooth> ptr;
  static constexpr int kStencil = 5;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<squad_smooth>(ids, w, blocks);
  }

  squad_smooth(const std::vector<index_t> &ids, const real &w,
               std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    assert(ids.size() >= kStencil);
  }
  virtual std::string name() { return typeid(*this).name(); }

  static quat hemisphere_align(const quat &ref, quat q) {
    if (ref.coeffs().dot(q.coeffs()) < 0.0)
      q.coeffs() = -q.coeffs();
    return q.normalized();
  }

  // Unit quat → so(3) vector θû (θ ∈ (-π, π]).
  static vec3 quat_log_vec(quat q) {
    q.normalize();
    if (q.w() < 0.0)
      q.coeffs() = -q.coeffs();
    const vec3 v = q.vec();
    const real vn = v.norm();
    if (vn < 1.0e-14)
      return vec3::Zero();
    const real w = std::max(real(-1.0), std::min(real(1.0), q.w()));
    const real angle = real(2.0) * std::atan2(vn, w);
    return (angle / vn) * v;
  }

  static quat quat_exp_vec(const vec3 &omega) {
    const real angle = omega.norm();
    if (angle < 1.0e-14)
      return quat::Identity();
    return quat(Eigen::AngleAxis<real>(angle, omega / angle));
  }

  static vec3 relative_omega(const quat &a, const quat &b) {
    const quat bn = hemisphere_align(a, b);
    return quat_log_vec(a.conjugate() * bn);
  }

  // Cubic Slerp De Casteljau on S³.
  static quat slerp_bezier(const quat &A, const quat &B, const quat &C,
                           const quat &D, real t) {
    const float tf = static_cast<float>(t);
    const quat ab = va::slerp(A, B, tf);
    const quat bc = va::slerp(B, C, tf);
    const quat cd = va::slerp(C, D, tf);
    const quat abc = va::slerp(ab, bc, tf);
    const quat bcd = va::slerp(bc, cd, tf);
    return va::slerp(abc, bcd, tf).normalized();
  }

  // Ends Q0,Q4; jets from one-step relatives, scaled to full [0,1] window (4 steps).
  static quat eval_mid(const quat Q[kStencil]) {
    quat A = Q[0].normalized();
    quat Q1 = hemisphere_align(A, Q[1]);
    quat Q3 = hemisphere_align(A, Q[3]);
    quat D = hemisphere_align(A, Q[4]);
    Q3 = hemisphere_align(D, Q3);

    const vec3 wa = real(4.0) * relative_omega(A, Q1);
    const vec3 wb = real(4.0) * relative_omega(Q3, D);
    const quat B = (A * quat_exp_vec(wa / real(3.0))).normalized();
    const quat C = (D * quat_exp_vec(-wb / real(3.0))).normalized();
    return slerp_bezier(A, B, C, D, real(0.5));
  }

  virtual void project(const vecX &q, vecX &p) {
    quat Q[kStencil];
    for (int i = 0; i < kStencil; ++i)
      Q[i] = _blocks[0]->get_quat(this->_ids[i], q).normalized();

    quat target = eval_mid(Q);
    const quat mid = Q[2];
    target = hemisphere_align(mid, target);
    if (target.coeffs().hasNaN()) {
      p.block(_id0, 0, 4, 1) = _w * vec4(mid.coeffs().data());
      return;
    }
    p.block(_id0, 0, 4, 1) = _w * vec4(target.coeffs().data());
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t i = _blocks[0]->get_offset_idx(this->_ids[2]);
    for (int ax = 0; ax < 4; ++ax)
      triplets.push_back(trip(_id0 + ax, i + ax, _w));
    id0 += 4;
  }
};

/// Position-only rest-length spring (1D edge stretch). No quaternions.
/// PD form matching shell edge_strain / legacy edge_stretch energy (ℓ−ℓ₀)².
class edge_stretch : public block_constraint {
public:
  typedef std::shared_ptr<edge_stretch> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l0, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<edge_stretch>(ids, w, l0, blocks);
  }

  edge_stretch(const std::vector<index_t> &ids, const real &w, const real &l0,
               std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _l0(l0) {}

  virtual std::string name() { return "edge_stretch"; }

  virtual void project(const vecX &q, vecX &p) {
    const index_t i = this->_ids[0];
    const index_t j = this->_ids[1];
    const vec3 q0 = _blocks[0]->get_vec3(i, q);
    const vec3 q1 = _blocks[0]->get_vec3(j, q);
    vec3 dq = q1 - q0;
    real l = dq.norm();
    l = std::max(l, real(1e-8));
    dq /= l;
    p.block(_id0, 0, 3, 1) = _w * dq;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    const index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);
    const real s = _w / std::max(_l0, real(1e-8));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, -s));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, j + ax, s));
    id0 += 3;
  }

  real _l0 = 1.0;
};

class stretch_shear : public block_constraint {
public:
  typedef std::shared_ptr<stretch_shear> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l0, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<stretch_shear>(ids, w, l0, blocks);
  }

  // New dual-weight constructor
  static ptr create(const std::vector<index_t> &ids, const real &w1, const real &w2,
                    const real &l0, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<stretch_shear>(ids, w1, w2, l0, blocks);
  }

  stretch_shear(const std::vector<index_t> &ids, const real &w, const real &l0,
                std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _l0(l0), _w1(w), _w2(w) {}
  
  // New dual-weight constructor
  stretch_shear(const std::vector<index_t> &ids, const real &w1, const real &w2, const real &l0,
                std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w1, blocks), _l0(l0), _w1(w1), _w2(w2) {}
  
  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];
    index_t k = this->_ids[2];

    vec3 q0 = _blocks[0]->get_vec3(i, q);
    vec3 q1 = _blocks[0]->get_vec3(j, q);

    real l0 = _l0;
    real l = (q1 - q0).norm();
    vec3 dq = (q1 - q0).normalized();

    quat u = _blocks[1]->get_quat(k, q);

    vec3 d2 = u * vec3(0, 0, 1);
    // d2.normalize();
    //geometry_logger::line(q0, q0 + 0.1 * d2, vec4(1.0, 0.0, 0.0, 1.0));
    quat du = quat::FromTwoVectors(d2, dq);

    u = du * u;
    // u.normalize();
    if (d2.hasNaN()) {
      std::cout << __PRETTY_FUNCTION__ << " d2 is nan" << std::endl;
      exit(0);
    }
    if (u.coeffs().hasNaN()) {
      std::cout << __PRETTY_FUNCTION__ << " u is nan" << std::endl;
      exit(0);
    }

    // l = std::clamp(l / _l0, 0.1, 1.5);
    // p.block(_id0, 0, 3, 1) = _w * l * d2;
    p.block(_id0, 0, 3, 1) = _w1 * d2;  // Use _w1 for stretch/shear component

    // p.block(k, 0, 4, 1) += _w * q.block(k, 0, 4, 1);
    p.block(_id0 + 3, 0, 4, 1) = _w2 * vec4(u.coeffs().data());  // Use _w2 for rotation component
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t ii = _blocks[1]->get_offset_idx(this->_ids[2]);

    // index_t Nv = this->_ids[3];
    // index_t k = 3 * Nv + 4 * this->_ids[2];

    // std::cout << ii << " " << k << std::endl;
    // std::cout << "ii: " << ii << std::endl;
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, -_w1 / _l0));  // Use _w1
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, j + ax, _w1 / _l0));   // Use _w1

    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, ii + ax, _w2));  // Use _w2
    id0 += 7;
  }
  real _l0;
  real _w1;  // Weight for stretch/shear component
  real _w2;  // Weight for rotation component
};

class straight : public block_constraint {
public:
  typedef std::shared_ptr<straight> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<straight>(ids, w, blocks);
  }

  straight(const std::vector<index_t> &ids, const real &w,
             std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}
  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {

    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];

    quat ui = _blocks[0]->get_quat(ii, q);
    quat uj = _blocks[0]->get_quat(jj, q);
    quat uij = va::slerp(ui, uj, 0.5);

    uij.normalize();

    ui = uij;
    uj = uij;

    p.block(_id0 + 0, 0, 4, 1) = _w * vec4(ui.coeffs().data());
    p.block(_id0 + 4, 0, 4, 1) = _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);

    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 4 + ax, j + ax, _w));
    id0 += 8;
  }
};


namespace detail {
inline quat rest_relative(const std::vector<quat> &u, index_t i, index_t j) {
  return u[i].inverse() * u[j];
}

// omega from relative-error quat dO = O^{-1} O0; split bend (⊥ e2) / twist (∥ e2).
inline void split_bend_twist_omega(const quat &dO, vec3 &omega_b, vec3 &omega_t) {
  Eigen::AngleAxis<real> aa(dO);
  const vec3 omega = aa.angle() * aa.axis();
  const vec3 e2 = vec3::UnitZ();
  omega_t = omega.dot(e2) * e2;
  omega_b = omega - omega_t;
}

inline quat omega_to_quat(const vec3 &omega) {
  const real ang = omega.norm();
  if (ang < 1.0e-16)
    return quat::Identity();
  return quat(Eigen::AngleAxis<real>(ang, omega / ang));
}
} // namespace detail

// Bend-only Cosserat hinge: corrects ω_b (⊥ material e2 / tangent).
class bend : public block_constraint {
public:
  typedef std::shared_ptr<bend> ptr;

  static ptr create(const std::vector<index_t> &ids, const std::vector<quat> &u,
                    const real &w, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<bend>(ids, u, w, blocks);
  }

  bend(const std::vector<index_t> &ids, const std::vector<quat> &u, const real &w,
       std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    _O = detail::rest_relative(u, ids[0], ids[1]);
  }

  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    quat ui = _blocks[0]->get_quat(ii, q);
    quat uj = _blocks[0]->get_quat(jj, q);
    quat O = ui.inverse() * uj;
    quat dO = O.inverse() * _O;

    vec3 omega_b, omega_t;
    detail::split_bend_twist_omega(dO, omega_b, omega_t);
    dO = detail::omega_to_quat(omega_b);

    ui = (ui * dO.inverse()).normalized();
    uj = (uj * dO).normalized();

    p.block(_id0 + 0, 0, 4, 1) = _w * vec4(ui.coeffs().data());
    p.block(_id0 + 4, 0, 4, 1) = _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 4 + ax, j + ax, _w));
    id0 += 8;
  }

  quat _O;
};

// Twist-only Cosserat hinge: corrects ω_t (∥ material e2 / tangent).
class twist : public block_constraint {
public:
  typedef std::shared_ptr<twist> ptr;

  static ptr create(const std::vector<index_t> &ids, const std::vector<quat> &u,
                    const real &w, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<twist>(ids, u, w, blocks);
  }

  twist(const std::vector<index_t> &ids, const std::vector<quat> &u, const real &w,
        std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    _O = detail::rest_relative(u, ids[0], ids[1]);
  }

  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    quat ui = _blocks[0]->get_quat(ii, q);
    quat uj = _blocks[0]->get_quat(jj, q);
    quat O = ui.inverse() * uj;
    quat dO = O.inverse() * _O;

    vec3 omega_b, omega_t;
    detail::split_bend_twist_omega(dO, omega_b, omega_t);
    dO = detail::omega_to_quat(omega_t);

    ui = (ui * dO.inverse()).normalized();
    uj = (uj * dO).normalized();

    p.block(_id0 + 0, 0, 4, 1) = _w * vec4(ui.coeffs().data());
    p.block(_id0 + 4, 0, 4, 1) = _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 4 + ax, j + ax, _w));
    id0 += 8;
  }

  quat _O;
};


class angle : public block_constraint {
public:
  typedef std::shared_ptr<angle> ptr;

  static ptr create(const std::vector<index_t> &ids, vec3 z, real phi,
                    const real &w, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<angle>(ids, z, phi, w, blocks);
  }

  angle(const std::vector<index_t> &ids, vec3 z, real phi, const real &w,
        std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _z(z), _phi(phi) {}
  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {

    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];

    quat ui = _blocks[0]->get_quat(ii, q).normalized();
    quat uj = _blocks[0]->get_quat(jj, q).normalized();

    /*
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    index_t Nv = this->_ids[2];
    index_t i = 3 * Nv + 4 * ii;
    index_t j = 3 * Nv + 4 * jj;
    quat ui = quat(q.block(i, 0, 4, 1).data()).normalized();
    quat uj = quat(q.block(j, 0, 4, 1).data()).normalized();
*/
    vec3 zi = ui * _z;
    vec3 zj = uj * _z;

    vec3 N = zi.cross(zj).normalized();
    real thet = atan2(zi.cross(zj).dot(N), zi.dot(zj));
    real dthet = thet - _phi;
    quat ugi(Eigen::AngleAxisd(0.5 * dthet, N));
    quat ugj(Eigen::AngleAxisd(-0.5 * dthet, N));
    ui = ugi * ui;
    uj = ugj * uj;

    // zi = ui * _z;
    // zj = uj * _z;
    // real thetp = atan2(zi.cross(zj).dot(N), zi.dot(zj));
    //  std::cout << thet << " " << thetp << std::endl;
    //  gg::geometry_geometry_logger::line(q0, q0 + 0.05 * zi, vec4(1.5, 0.0, 0.0, 1.0));
    //  gg::geometry_geometry_logger::line(q0, q0 + 0.05 * zj, vec4(0.0, 1.5, 0.0, 1.0));
    if (ui.coeffs().hasNaN()) {
      std::cout << __PRETTY_FUNCTION__ << " ui is nan" << std::endl;
      exit(0);
    }
    if (uj.coeffs().hasNaN()) {
      std::cout << __PRETTY_FUNCTION__ << " uj is nan" << std::endl;
      exit(0);
    }

    p.block(_id0 + 0, 0, 4, 1) = _w * vec4(ui.coeffs().data());
    p.block(_id0 + 4, 0, 4, 1) = _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 0 + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 4 + ax, j + ax, _w));

    id0 += 8;
  }
  real _phi = 0.0;

  vec3 _z = vec3(1.0, 0.0, 0.0);
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif
