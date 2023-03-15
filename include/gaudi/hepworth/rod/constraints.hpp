
#ifndef __HEP_ROD_CONSTRAINTS__
#define __HEP_ROD_CONSTRAINTS__

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

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace rod {
class projection_constraint {
public:
  typedef std::shared_ptr<projection_constraint> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<projection_constraint>(ids, w);
  }

  projection_constraint(const std::vector<index_t> &ids, const real &w)
      : _ids(ids), _w(w) {}
  virtual void project(const vecX &q, vecX &p){};
  virtual void fill_A(std::vector<trip> &triplets){};

  std::vector<index_t> _ids;

  real _w;
};

class growth : public projection_constraint {
public:
  typedef std::shared_ptr<growth> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l) {
    return std::make_shared<growth>(ids, w, l);
  }

  growth(const std::vector<index_t> &ids, const real &w, const real &l)
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

    p.block(3 * i, 0, 3, 1) = _w * _l * dq;
    // p.block(3 * j, 0, 3, 1) = -_w * dq / l;
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * i + ax, -1.0 * _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * j + ax, 1.0 * _w));
    // for (int ax = 0; ax < 3; ax++)
    //   triplets.push_back(trip(3 * j + ax, 3 * j + ax, -_w / _l));
    // for (int ax = 0; ax < 3; ax++)
    //   triplets.push_back(trip(3 * j + ax, 3 * i + ax, _w / _l));
  }
  real _l = 1.0;
};

class smooth : public projection_constraint {
public:
  typedef std::shared_ptr<smooth> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<smooth>(ids, w);
  }

  smooth(const std::vector<index_t> &ids, const real &w)
      : projection_constraint(ids, w) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i0 = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t in = this->_ids[2];

    vec3 qm = q.block(3 * ip, 0, 3, 1);
    vec3 q0 = q.block(3 * i0, 0, 3, 1);
    vec3 qp = q.block(3 * in, 0, 3, 1);
    vec3 dqm = qm - q0;
    vec3 dqp = qp - q0;
    vec3 N = -(dqm + dqp).normalized();
    dqm = va::reject(N, dqm);
    dqp = va::reject(N, dqp);
    p.block(3 * i0, 0, 3, 1) += _w * (dqm + dqp);
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    index_t i0 = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t in = this->_ids[2];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i0 + ax, 3 * i0 + ax, -2.0 * _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i0 + ax, 3 * ip + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i0 + ax, 3 * in + ax, _w));
  }
};

class cylinder : public projection_constraint {
public:
  typedef std::shared_ptr<cylinder> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<cylinder>(ids, w);
  }

  cylinder(const std::vector<index_t> &ids, const real &w)
      : projection_constraint(ids, w) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i0 = this->_ids[0];

    real r = 0.0;
    vec3 cen = vec3::Zero();
    int k = 0;
    for (int i = 1; i < this->_ids.size() - 2; i++) {
      int i0 = this->_ids[i + 0];
      int i1 = this->_ids[i + 1];
      int i2 = this->_ids[i + 2];

      vec3 q0 = q.block(3 * i0, 0, 3, 1);
      vec3 q1 = q.block(3 * i1, 0, 3, 1);
      vec3 q2 = q.block(3 * i2, 0, 3, 1);
      vec3 ceni;
      real ri;
      va::estimate_3D_circle(q0, q1, q2, ceni, ri);
      cen += ceni;
      r += ri;
      k++;
    }
    cen /= real(k);
    r /= real(k);
    mat3 U = mat3::Zero();
    for (int i = 1; i < this->_ids.size(); i++) {
      int ii = this->_ids[i];
      vec3 qi = q.block(3 * ii, 0, 3, 1);
      vec3 dq = qi - cen;
      U += dq * dq.transpose();
    }

    Eigen::JacobiSVD<mat3> svd(U, Eigen::ComputeFullU);
    U = svd.matrixU();
    vec3 s = svd.singularValues();
    vec3 N = U.col(2).transpose();
    r = sqrt(s[0]);
    //  r = 0.5;
    //    r = max(r, 0.3);
    vec3 qi = q.block(3 * i0, 0, 3, 1);
    vec3 qc = va::project_on_line(cen, vec3(cen + N), qi);
    vec3 dq = r * (qi - qc).normalized();
    real t = s[1] / 5.0 / s[0];
    // std::cout << t << std::endl;
    vec3 qp = qc + dq;
    qp = va::mix(t, qi, qp);
    //  qp = qi;
    //   gg::geometry_logger::line(qi, qc, vec4(1.0, 1.0, 0.0, 1.0));
    //   gg::geometry_logger::line(qi, qp, vec4(1.0, 0.0, 0.0, 1.0));

    p.block(3 * i0, 0, 3, 1) += _w * (qp);
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    // return;

    int ii = this->_ids[0];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * ii + ax, 3 * ii + ax, 1.0 * _w));
  }
};

class stretch_shear : public projection_constraint {
public:
  typedef std::shared_ptr<stretch_shear> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l0) {
    return std::make_shared<stretch_shear>(ids, w, l0);
  }

  stretch_shear(const std::vector<index_t> &ids, const real &w, const real &l0)
      : projection_constraint(ids, w), _l0(l0) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];
    index_t ii = this->_ids[2];
    index_t Nv = this->_ids[3];
    index_t k = 3 * Nv + 4 * ii;

    vec3 q0 = q.block(3 * i, 0, 3, 1);
    vec3 q1 = q.block(3 * j, 0, 3, 1);
    real l0 = _l0;
    real l = (q1 - q0).norm();
    vec3 dq = (q1 - q0).normalized();

    quat u = quat(q.block(k, 0, 4, 1).data()).normalized();
    vec3 d2 = u * vec3(0, 0, 1);
    // d2.normalize();

    quat du = quat::FromTwoVectors(d2, dq);

    u = du * u;
    // u.normalize();
    p.block(3 * i, 0, 3, 1) += _w * d2;
    // p.block(k, 0, 4, 1) += _w * q.block(k, 0, 4, 1);
    p.block(k, 0, 4, 1) += _w * vec4(u.coeffs().data());
  }

  virtual void fill_A(std::vector<trip> &triplets) {
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];
    index_t ii = this->_ids[2];
    index_t Nv = this->_ids[3];
    index_t k = 3 * Nv + 4 * ii;

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * i + ax, -_w / _l0));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * j + ax, _w / _l0));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(k + ax, k + ax, _w));
  }
  real _l0;
};

class bend_twist : public projection_constraint {
public:
  typedef std::shared_ptr<bend_twist> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<bend_twist>(ids, w);
  }

  bend_twist(const std::vector<index_t> &ids, const real &w)
      : projection_constraint(ids, w) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    index_t Nv = this->_ids[2];
    index_t i = 3 * Nv + 4 * ii;
    index_t j = 3 * Nv + 4 * jj;
    quat ui = quat(q.block(i, 0, 4, 1).data());
    quat uj = quat(q.block(j, 0, 4, 1).data());

#if 0
    vec3 v1 = ui.vec().normalized();
    vec3 v2 = uj.vec().normalized();

    // Step 2: Compute the midpoint on the unit sphere
    vec3 vmid = (v1 + v2).normalized();

    // Step 3: Compute the quaternion corresponding to the midpoint
    quat q_min = quat::FromTwoVectors(v1, vmid);

    // Step 4: Rotate both quaternions to point in the same direction
    ui = q_min * ui;
    uj = q_min.conjugate() * uj;
#elif 1
    quat uij = ui.slerp(0.5, uj);
    /*
    vec3 q0 = q.block(3 * ii, 0, 3, 1);
    vec3 d20 = ui * vec3(0, 0, 0.1);
    gg::geometry_logger::line(q0, q0 + d20, vec4(1.0, 0.0, 0.0, 1.0));
    vec3 d21 = uj * vec3(0, 0, 0.1);
    gg::geometry_logger::line(q0, q0 + d21, vec4(0.0, 0.0, 1.0, 1.0));
*/
    ui = uij;
    uj = uij;
    /*
    vec3 d22 = ui * vec3(0, 0, 0.1);
    gg::geometry_logger::line(q0, q0 + d22, vec4(1.0, 0.0, 1.0, 1.0));
    vec3 d23 = uj * vec3(0, 0, 0.1);
    gg::geometry_logger::line(q0, q0 + d23, vec4(1.0, 1.0, 0.0, 1.0));
*/
    // ui.normalize();
    // uj.normalize();

#endif

    p.block(i, 0, 4, 1) += _w * vec4(ui.coeffs().data());
    p.block(j, 0, 4, 1) += _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(std::vector<trip> &triplets) {
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    index_t Nv = this->_ids[2];
    index_t i = 3 * Nv + 4 * ii;
    index_t j = 3 * Nv + 4 * jj;
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(i + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(j + ax, j + ax, _w));
  }
};

class angle : public projection_constraint {
public:
  typedef std::shared_ptr<angle> ptr;

  static ptr create(const std::vector<index_t> &ids, vec3 z, real phi,
                    const real &w) {
    return std::make_shared<angle>(ids, z, phi, w);
  }

  angle(const std::vector<index_t> &ids, vec3 z, real phi, const real &w)
      : projection_constraint(ids, w), _z(z), _phi(phi) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    index_t Nv = this->_ids[2];
    index_t i = 3 * Nv + 4 * ii;
    index_t j = 3 * Nv + 4 * jj;
    quat ui = quat(q.block(i, 0, 4, 1).data()).normalized();
    quat uj = quat(q.block(j, 0, 4, 1).data()).normalized();

    vec3 zi = ui * _z;
    vec3 zj = uj * _z;

    vec3 q0 = q.block(3 * ii, 0, 3, 1);
    // gg::geometry_logger::line(q0, q0 + 0.05 * zi, vec4(0.75, 0.0, 0.0, 1.0));
    // gg::geometry_logger::line(q0, q0 + 0.05 * zj, vec4(0.0, 0.75, 0.0, 1.0));

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

    p.block(i, 0, 4, 1) += _w * vec4(ui.coeffs().data());
    p.block(j, 0, 4, 1) += _w * vec4(uj.coeffs().data());
  }

  virtual void fill_A(std::vector<trip> &triplets) {

    index_t ii = this->_ids[0];
    index_t jj = this->_ids[1];
    index_t Nv = this->_ids[2];
    index_t i = 3 * Nv + 4 * ii;
    index_t j = 3 * Nv + 4 * jj;
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(i + ax, i + ax, _w));
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(j + ax, j + ax, _w));
  }
  real _phi = 0.0;
  vec3 _z = vec3(1.0, 0.0, 0.0);
};

class smooth_bend : public projection_constraint {
public:
  typedef std::shared_ptr<smooth_bend> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<smooth_bend>(ids, w);
  }

  smooth_bend(const std::vector<index_t> &ids, const real &w)
      : projection_constraint(ids, w) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t iip = this->_ids[0];
    index_t ii0 = this->_ids[1];
    index_t iin = this->_ids[2];
    index_t Nv = this->_ids[3];
    index_t ip = 3 * Nv + 4 * iip;
    index_t i0 = 3 * Nv + 4 * ii0;
    index_t in = 3 * Nv + 4 * iin;
    quat uip = quat(q.block(ip, 0, 4, 1).data());
    quat ui0 = quat(q.block(i0, 0, 4, 1).data());
    quat uin = quat(q.block(in, 0, 4, 1).data());

    ui0 = uip.slerp(0.5, uin).normalized();

    p.block(ii0, 0, 4, 1) += _w * vec4(ui0.coeffs().data());
  }

  virtual void fill_A(std::vector<trip> &triplets) {
    index_t iip = this->_ids[0];
    index_t ii0 = this->_ids[1];
    index_t iin = this->_ids[2];
    index_t Nv = this->_ids[3];
    index_t ip = 3 * Nv + 4 * iip;
    index_t i0 = 3 * Nv + 4 * ii0;
    index_t in = 3 * Nv + 4 * iin;
    for (int ax = 0; ax < 4; ax++)
      triplets.push_back(trip(i0 + ax, i0 + ax, _w));
  }
};
} // namespace rod
} // namespace hepworth
} // namespace gaudi
#endif