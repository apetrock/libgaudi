
#ifndef __HEP_ROD_BLOCK_CONSTRAINTS__
#define __HEP_ROD_BLOCK_CONSTRAINTS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <limits>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#include "gaudi/common.h"

#include "../projection_constraint.hpp"
#include "block_constraint.hpp"
#include "sim_block.hpp"

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

class stretch_shear : public block_constraint {
public:
  typedef std::shared_ptr<stretch_shear> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l0, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<stretch_shear>(ids, w, l0, blocks);
  }

  stretch_shear(const std::vector<index_t> &ids, const real &w, const real &l0,
                std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _l0(l0) {}
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
    //logger::line(q0, q0 + 0.1 * d2, vec4(1.0, 0.0, 0.0, 1.0));
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
    p.block(_id0, 0, 3, 1) = _w * d2;

    // p.block(k, 0, 4, 1) += _w * q.block(k, 0, 4, 1);
    p.block(_id0 + 3, 0, 4, 1) = _w * vec4(u.coeffs().data());
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
      triplets.push_back(trip(_id0 + ax, i + ax, -_w / _l0));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, j + ax, _w / _l0));

    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(_id0 + 3 + ax, ii + ax, _w));
    id0 += 7;
  }
  real _l0;
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


class bend_twist : public block_constraint {
public:

  //DEFINE_CREATE_FUNC(bend_twist)
  static ptr create(const std::vector<index_t> &ids, 
             const std::vector<quat> &u, const real &w,
             std::vector<sim_block::ptr> blocks) {
    return std::make_shared<bend_twist>(ids, u, w, blocks);
  }
  bend_twist(const std::vector<index_t> &ids, 
             const std::vector<quat> &u, const real &w,
             std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {

        quat ui = u[ids[0]];
        quat uj = u[ids[1]];
        _O = ui.inverse() * uj;
      }

  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {

    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];

    quat ui = _blocks[0]->get_quat(ii, q);
    quat uj = _blocks[0]->get_quat(jj, q);
    quat O = ui.inverse() * uj;
    quat dO = O.inverse() * _O;
    //dO = uj * ui.inverse() * _
    ui = ui * dO.inverse();
    uj = uj * dO;

    ui.normalize();
    uj.normalize();

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
    //  gg::geometry_logger::line(q0, q0 + 0.05 * zi, vec4(1.5, 0.0, 0.0, 1.0));
    //  gg::geometry_logger::line(q0, q0 + 0.05 * zj, vec4(0.0, 1.5, 0.0, 1.0));
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

void Preprocess(const std::vector<vec3> &points, std::vector<vec3> &X,
                vec3 &average, vec6 &mu, mat3 &F0, mat36 &F1, mat6 &F2) {

  real n = real(points.size());
  average = vec3::Zero();
  for (int i = 0; i < points.size(); ++i) {
    average += points[i];
  }

  average /= n;

  for (int i = 0; i < points.size(); ++i) {
    X[i] = points[i] - average;
  }

  vec6 zero = vec6::Zero();
  std::vector<vec6> products(points.size(), vec6::Zero());
  mu = zero;

  for (int i = 0; i < n; ++i) {
    products[i][0] = X[i][0] * X[i][0];
    products[i][1] = X[i][0] * X[i][1];
    products[i][2] = X[i][0] * X[i][2];
    products[i][3] = X[i][1] * X[i][1];
    products[i][4] = X[i][1] * X[i][2];
    products[i][5] = X[i][2] * X[i][2];
    mu[0] += products[i][0];
    mu[1] += 2 * products[i][1];
    mu[2] += 2 * products[i][2];
    mu[3] += products[i][3];
    mu[4] += 2 * products[i][4];
    mu[5] += products[i][5];
  }
  mu /= n;
  F0 = mat3::Zero();
  F1 = mat36::Zero();
  F2 = mat6::Zero();
  for (int i = 0; i < n; ++i) {
    vec6 delta;
    delta[0] = products[i][0] - mu[0];
    delta[1] = 2.0 * products[i][1] - mu[1];
    delta[2] = 2.0 * products[i][2] - mu[2];
    delta[3] = products[i][3] - mu[3];
    delta[4] = 2.0 * products[i][4] - mu[4];
    delta[5] = products[i][5] - mu[5];
    F0(0, 0) += products[i][0];
    F0(0, 1) += products[i][1];
    F0(0, 2) += products[i][2];
    F0(1, 1) += products[i][3];
    F0(1, 2) += products[i][4];
    F0(2, 2) += products[i][5];
    F1 += X[i] * delta.transpose();
    F2 += delta * delta.transpose();
  }
  F0 /= n;
  F0(1, 0) = F0(0, 1);
  F0(2, 0) = F0(0, 2);
  F0(2, 1) = F0(1, 2);
  F1 /= n;
  F2 /= n;
}

real G(std::vector<vec3> &X, vec6 mu, mat3 F0, mat36 F1, mat6 F2, vec3 W,
       vec3 &PC, real &rSqr) {
  real n = real(X.size());
  mat3 P = mat3::Identity() - W * W.transpose(); // P = I = W * WˆT
  // S = {{0, =w2, w1}, {w2, 0, =w0}, {=w1, w0, 0}}, inner braces are rows
  mat3 S;
  S.block(0, 0, 3, 1) = vec3(0.0, -W[2], W[1]);
  S.block(0, 1, 3, 1) = vec3(W[2], 0.0, -W[0]);
  S.block(0, 2, 3, 1) = vec3(-W[1], W[0], 0.0);

  mat3 A = P * F0 * P;
  mat3 hatA = -(S * A * S);
  mat3 hatAA = hatA * A;
  real trace = hatAA.trace();
  mat3 Q = hatA / trace;
  vec6 p = {P(0, 0), P(0, 1), P(0, 2), P(1, 1), P(1, 2), P(2, 2)};
  vec3 alpha = F1 * p;
  vec3 beta = Q * alpha;
  real error =
      (p.dot(F2 * p) - 4.0 * alpha.dot(beta) + 4.0 * beta.dot(F0 * beta)) / n;
  PC = beta;
  rSqr = p.dot(mu) + beta.dot(beta);
  return error;
}

real FitCylinder(const std::vector<vec3> &points, const vec3 &W, real &rSqr,
                 vec3 &C) {

  std::vector<vec3> X(points.size());
  vec6 mu;
  vec3 avg;
  mat3 F0;
  mat36 F1;
  mat6 F2;
  C = vec3::Zero();
  Preprocess(points, X, avg, mu, F0, F1, F2);
  real error = G(X, mu, F0, F1, F2, W, C, rSqr);
  // Choose imax and jmax as desired for the level of granularity you
  // want for sampling W vectors on the hemisphere.
  real minError = std::numeric_limits<real>::infinity();

  C += avg;
  return minError;
}

#if 1
class helicitiy : public block_constraint {
public:
  typedef std::shared_ptr<helicitiy> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<helicitiy>(ids, w, blocks);
  }

  helicitiy(const std::vector<index_t> &ids, const real &w,
            std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}
  virtual std::string name() { return typeid(*this).name(); }

  vec3 find_axis(const vecX &q) {

    index_t i0 = this->_ids[0];
    vec3 q0 = _blocks[0]->get_vec3(i0, q);

    int k = 0;
    vec3 cen = vec3::Zero();
    mat3 U = mat3::Zero();
    vec3 dq = vec3::Zero();
    for (int i = 1; i < this->_ids.size() - 1; i++) {
      int ii0 = this->_ids[i + 0];
      int ii1 = this->_ids[i + 1];
      if (i0 == ii1 || i0 == ii1)
        continue;

      vec3 qi0 = q.block(3 * ii0, 0, 3, 1);
      vec3 qi1 = q.block(3 * ii1, 0, 3, 1);
      vec3 ceni;
      real ri;
      va::estimate_3D_circle(q0, qi0, qi1, ceni, ri);
      if (std::isnan(ri))
        continue;
      vec3 q_avg = 0.5 * (qi0 + qi1);
      vec3 dq0 = q0 - ceni;
      vec3 dqi0 = qi0 - ceni;
      vec3 dqi1 = qi1 - ceni;
      dq += (dqi0 + dqi1).normalized();
      U += dqi0 * dqi0.transpose();
      U += dqi1 * dqi1.transpose();

      cen += ceni;
      k++;
    }
    dq /= real(k);
    U /= real(k);
    cen /= real(k);

    Eigen::JacobiSVD<mat3> svd(U, Eigen::ComputeFullU);
    U = svd.matrixU();
    vec3 s = svd.singularValues();
    vec3 T0 = U.col(0).transpose();
    vec3 T1 = U.col(1).transpose();
    vec3 T2 = U.col(2).transpose();
    vec3 Ti = q0 - cen;
    real t0 = abs(Ti.dot(T0));
    real t1 = abs(Ti.dot(T1));
    real t2 = abs(Ti.dot(T2));
    T0 = va::sgn(Ti.dot(T0)) * T0;
    T1 = va::sgn(dq.dot(T1)) * T1;
    
    // eigen::angleaxis rotation about T0, '
    //std::random_device rd;
    //s/td::mt19937 gen(rd());
    //std::uniform_real_distribution<> dis(-1.0, 1.0);

    //Eigen::AngleAxis<real> R(0.01 * dis(gen) * M_PI, T1);

    //T2 = R * T2;
    //  logger::line(q0, q0 + 0.1 * T2, vec4(0.5, 0.0, 1.0, 1.0));
    int i = 0;
    i = t0 < t1 && t0 < t2 ? 0 : i;
    i = t1 < t0 && t1 < t2 ? 1 : i;
    i = t2 < t0 && t2 < t1 ? 2 : i;
    return T2;
  }

  virtual void project(const vecX &q, vecX &p) {
    index_t i0 = this->_ids[0];
    vec3 q0 = _blocks[0]->get_vec3(i0, q);

    vec3 dq = vec3::Zero();
    int k = 0;
    std::vector<vec3> qs(this->_ids.size() - 1);
    for (int i = 1; i < this->_ids.size(); i++) {
      int ii = this->_ids[i + 0];
      vec3 qi = _blocks[0]->get_vec3(ii, q);
      qs[i - 1] = qi;
    }

    vec3 W = find_axis(q);
    real r2;
    vec3 C;
    real e = FitCylinder(qs, W, r2, C);
    vec3 qc = va::project_on_line(C, vec3(C + W), q0);
    // random number between 0 and 1 using std::mersenne_twister_engine
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    vec3 qp = qc + (1.0 + dis(gen) * 0.001) * sqrt(r2) * (q0 - qc).normalized();
    // logger::line(qc, qp, vec4(0.5, 0.0, 1.0, 1.0));
    p.block(_id0, 0, 3, 1) += _w * qp;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    // return;
    _id0 = id0;
    index_t ii = _blocks[0]->get_offset_idx(this->_ids[0]);
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, ii + ax, 1.0 * _w));
    id0 += 3;
  }
};
#endif
} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif