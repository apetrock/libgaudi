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
#include "primitive_objects.hpp"

#ifndef __ASAWA_ROD__
#define __ASAWA_ROD__

namespace gaudi {
namespace asawa {
typedef int index_t;

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

class length : public projection_constraint {
public:
  typedef std::shared_ptr<length> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l) {
    return std::make_shared<length>(ids, w, l);
  }

  length(const std::vector<index_t> &ids, const real &w, const real &l)
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

    real l0 = _l0;
    vec3 q0 = q.block(3 * i, 0, 3, 1);
    vec3 q1 = q.block(3 * j, 0, 3, 1);
    vec3 dq = 1.0 / l0 * (q1 - q0);

    quat u = quat(q.block(k, 0, 4, 1).data());
    u.normalize();
    vec3 d2 = u * vec3(0, 0, 1);
    quat du = quat::FromTwoVectors(dq, d2);
    u = u * du;
    //  std::cout << "a: " << dq.transpose() << std::endl;
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
    quat O = ui.conjugate() * uj;
#if 1
    vec3 v1 = ui.vec().normalized();
    vec3 v2 = uj.vec().normalized();

    // Step 2: Compute the midpoint on the unit sphere
    vec3 vmid = (v1 + v2).normalized();

    // Step 3: Compute the quaternion corresponding to the midpoint
    quat q_min = quat::FromTwoVectors(v1, vmid);

    // Step 4: Rotate both quaternions to point in the same direction
    ui = q_min * ui * q_min.inverse();
    uj = q_min * uj * q_min.inverse();
#endif
#if 0
    real k = 2.0 / M_PI * acos(fabs((ui.inverse() * uj).w()));

    // Step 2: Compute the average quaternion
    quat q_avg =
        ui * exp((1.0 / 2.0) * log((ui.inverse() * uj).inverse())) * uj;

    // Step 3: Rotate both quaternions to point in the same direction
    ui = q_avg * ui * q_avg.inverse();
    uj = q_avg * uj * q_avg.inverse();
#endif
    p.block(i, 0, 4, 1) += _w * vec4(uj.coeffs().data());
    p.block(j, 0, 4, 1) += _w * vec4(ui.coeffs().data());
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

class rod {
public:
  typedef std::shared_ptr<rod> ptr;

  static ptr create() {
    ptr R = std::make_shared<rod>();
    return R;
  }
  static ptr create(const std::vector<vec3> &verts, bool loop = true) {
    ptr R = std::make_shared<rod>();
    R->insert_strand(verts, loop);
    return R;
  }

  rod() {}

  void insert_strand(const std::vector<vec3> &verts, bool loop) {

    __x.insert(__x.end(), verts.begin(), verts.end());

    int vN = verts.size();
    int cbegin = __corners_prev.size();

    for (int i = 0; i < vN; i++) {
      if (i == 0) {
        if (loop) {
          __corners_prev.push_back(cbegin + vN - 1);
          __corners_next.push_back(cbegin + i + 1);
        } else {
          __corners_prev.push_back(-1);
          __corners_next.push_back(cbegin + 1);
        }
      } else if (i == vN - 1) {
        if (loop) {
          __corners_prev.push_back(cbegin + i - 1);
          __corners_next.push_back(cbegin + 0);
        } else {
          __corners_prev.push_back(cbegin + i - 1);
          __corners_next.push_back(-1);
        }
      } else {
        __corners_prev.push_back(cbegin + i - 1);
        __corners_next.push_back(cbegin + i + 1);
      }
    }
    _init_params();
  }

  void _init_params() {
    __l0.resize(__corners_next.size());
    __v.resize(__corners_next.size());
    __M.resize(__corners_next.size());
    __J.resize(__corners_next.size());
    __o.resize(__corners_next.size());

    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_next[i] == -1)
        continue;
      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      real l = (q1 - q0).norm();
      __l0[i] = l;
      __v[i] = vec3::Zero();
      real M = M_PI * _r * _r * l;
      real J = M * _r * _r;

      __M[i][0] = M;
      __M[i][1] = M;
      __M[i][2] = M;

      __J[i][0] = 0.0;
      __J[i][1] = J;
      __J[i][2] = J;
      __J[i][3] = __J[i][1] + __J[i][2];
    }
    _update_frames();
  }

  quat _calc_frame(const index_t &ip, const index_t &i, const index_t &in) {
    vec3 cp = __x[ip];
    vec3 c0 = __x[i];
    vec3 cn = __x[in];

    vec3 dn0 = cn - c0;
    vec3 d0p = c0 - cp;

    vec3 d0 = dn0.cross(d0p).normalized();
    vec3 d1 = dn0.cross(d0).normalized();
    vec3 d2 = d0.cross(d1).normalized();

    mat3 F;
    F.col(0) = d0;
    F.col(1) = d1;
    F.col(2) = d2;

    quat q(F);
    q.normalize();

    return q;
  }

  void _update_frames() {
    __u.resize(__corners_next.size());
    for (int i = 0; i < __corners_next.size(); i++) {

      if (next(i) == -1)
        __u[i] = _calc_frame(prev(prev(i)), prev(i), i);
      else if (prev(i) == -1)
        __u[i] = _calc_frame(i, next(i), next(next(i)));
      else
        __u[i] = _calc_frame(prev(i), i, next(i));

      index_t j = __corners_next[i];
      vec3 q0 = __x[i];
      vec3 q1 = __x[j];
      vec3 dq = q1 - q0;
      dq.normalize();
      //__o[i] = quat(0.0, dq[0], dq[1], dq[2]);
      __o[i] = quat(0.0, 0.0, 0.0, 0.0);
    }
  }
  // accessors
  index_t other(index_t id) const { return 2 * (id / 2) + (id + 1) % 2; };

  index_t next(index_t i) { return __corners_next[i]; }
  index_t prev(index_t i) { return __corners_prev[i]; }

  std::vector<vec3> &corner_verts() { return __x; }
  const std::vector<vec3> &corner_verts() const { return __x; }

  vecX solve(matS &A, vecX &b) {

    // Eigen::ConjugateGradient<matS, Eigen::Upper> solver;
    Eigen::SimplicialLDLT<matS> solver;
    solver;

    solver.compute(A);

    if (solver.info() != Eigen::Success) {
      // decomposition failed
      std::cout << ".....decomposition error! " << std::endl;
    }
    vecX x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      // solving failed
      std::cout << ".....solve error! " << std::endl;
    }

    return x;
  }

  matS mass() {

    vecX M = to(__M);
    vecX J = to(__J);

    vecX MJ = concat(M, J);
    std::vector<trip> triplets;
    for (int i = 0; i < MJ.rows(); i++) {
      triplets.push_back(trip(i, i, MJ[i]));
    }
    matS Mx(MJ.size(), MJ.size());
    Mx.setFromTriplets(triplets.begin(), triplets.end());
    return Mx;
  }

  void step() {
    std::vector<projection_constraint::ptr> constraints;
    int Ni = __corners_next.size();
    int Nm = 3 * Ni + 4 * Ni;
    real h = 0.1;

    std::vector<quat> _u0(__u);
    for (size_t i = 0; i < __corners_next.size(); i++) {
      quat O = __o[i];
      quat u = __u[i];

      real J = J;
      quat sO = O; // torque terms + h / J
      quat su = u;
      su.coeffs() += 0.5 * h * (u * sO).coeffs();
      __u[i] = su;
      __u[i].normalize();
    }
    for (int i = 0; i < __corners_next.size(); i++) {
      index_t j = __corners_next[i];
      __l0[i] *= 1.01;
      constraints.push_back(length::create({i, j, i, Ni}, 4.0, 0.80 * __l0[i]));
    }
#if 1
    for (int i = 0; i < __corners_next.size(); i++) {
      index_t j = __corners_next[i];
      constraints.push_back(stretch_shear::create({i, j, i, Ni}, 0.1, __l0[i]));
    }
#endif
#if 1
    for (int i = 0; i < __corners_next.size(); i++) {
      index_t j = __corners_next[i];
      constraints.push_back(bend_twist::create({i, j, Ni}, 0.75));
    }
#endif
    matS A(Nm, Nm);
    matS M = mass();

    M.setIdentity();
    M *= 1.0 / h / h;

    std::vector<trip> triplets;
    for (auto &constraint : constraints) {
      constraint->fill_A(triplets);
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    matS AtA = A.transpose() * A;
    std::cout << "A   sum: " << A.sum() << std::endl;
    std::cout << "AtA sum: " << AtA.sum() << std::endl;

    vecX f = vecX::Zero(3 * Ni);
    // f.block(0, 0, 3, 1) = 0.1 * vec3(0.0, 1.0, 0.0);

    vecX x = to(__x);
    vecX v = to(__v);
    vecX u = to(__u);
    vecX o = to(__o);

    vecX x0 = x;
    vecX u0 = u;

    vecX q = concat(x, u);

    std::cout << x.size() << " " << v.size() << " " << f.size() << std::endl;
    x += h * v + h * h * f + 0e-12 * vecX::Random(3 * Ni);
    vecX s = concat(x, u);

    vecX p = vecX::Zero(Nm);

    for (int k = 0; k < 10; k++) {
      p.setZero();
      for (auto &constraint : constraints) {
        constraint->project(q, p);
      }

      std::cout << "pnorm: " << p.norm() << std::endl;
      vecX b = M * s + A.transpose() * p;
      // std::cout << p << std::endl;
      matS MAtA = M + AtA;
      q = solve(MAtA, b);
    }

    split(q, x, u);
    std::cout << " x norm: " << (x - x0).norm() << std::endl;
    std::cout << " u norm: " << (u - u0).norm() << std::endl;
    real damp = 0.9;
    v = (1.0 - damp) / h * (x - x0);
    from(__v, v);

    from(__x, x);
    from(__u, u);
    for (int i = 0; i < __u.size(); i++) {
      __o[0] = _u0[i].conjugate() * __u[i];
      __o[0].coeffs() *= 2.0 * (1.0 - damp) / h;
    }
  }

  // std::vector<vec3> &ghost_verts() { return __ghost_verts; }
  // const std::vector<vec3> &ghost_verts() const { return __ghost_verts; }
  /*
 std::array<vec3, 2> length_gradient(index_t i) {
   vec3 p0 = __x[__corners_prev[i]];
   vec3 p1 = __x[i];
   vec3 dp = p1 - p0;
   real li = dp.norm();
   real dl = li - l0;
   real wi0 = w0 / (w0 + w1);
   real wi1 = w1 / (w0 + w1);
   vec3 dp0 = wi0 * dl * dp / li;
   vec3 dp1 = -wi1 * dl * dp / li;
   return {dp0, dp1};
 }

 std::array<vec3, 3> perpendicular_bisector_gradient() {
   vec3 p0 = __x[__corners_prev[i]];
   vec3 p1 = __x[i];
   vec3 pg = __ghost_verts[i];

   vec3 pm = 0.5 * (p1 + p0);
   vec3 dp10 = p1 - p0;
   vec3 dp0g = p0 - pg;
   vec3 dpg1 = pg - p1;
   vec3 dpgm - pg - pm;
   real li = dp.norm();
   real dl = li - l0;
   real lam =           //
       dpgm.dot(dp10) / //
       (w0 * dp0g.norm() + w1 * dpg1.norm() + w2 * dp10.norm());
   vec3 dp0 = -w0 * lam * dp0g;
   vec3 dp1 = -w1 * lam * dpg1;
   vec3 dpg = -w2 * lam * dp10;
   return {dp0, dp1, dpg};
 }

 std::array<vec3, 3> diameter_gradient() {
   vec3 p0 = __x[__corners_prev[i]];
   vec3 p1 = __x[i];
   vec3 pg = __ghost_verts[i];

   vec3 pm = 0.5 * (p1 + p0);
   vec3 dpgm - pg - pm;
   vec3 ldpgm = dpgm.norm();
   real li = dp.norm();
   real dl = li - l0;
   real lam =         //
       (ldpgm - lg) / //
       (0.25 * w0 + 0.25 * w1 + w2);
   vec3 dp0 = 0.5 * w0 * lam * dpgm;
   vec3 dp1 = 0.5 * w1 * lam * dpgm;
   vec3 dpg = -w2 * lam * dpgm;
   return {dp0, dp1, dpg};
 }

 std::array<mat3, 9> frame_gradient() {
   auto X[](vec3 x)->mat3 {
     mat3 X = mat3::Zero();
     X.col(0) = vec3(0, x[2], -x[1]);
X.col(1) = vec3(-x[2],  0, x[0]);
X.col(2) = vec3(x[1], -x[0],  0);
return X;
} vec3 p0 = __x[__corners_prev[i]];
vec3 p1 = __x[i];
vec3 pg = __ghost_verts[i];

vec3 d10 = p1 - p0;
vec3 d0g = p0 - pg;
vec3 dg1 = pg - p1;

real lx = d10.cross(dp0g).norm();
real l10 = dp.norm();

vec3 d0 = F.col(0);
vec3 d1 = F.col(1);
vec3 d2 = F.col(2);
mat3 I = mat3::Identity();
mat3 Id2 = I - d2 * d2.transpose();
mat3 Id1 = I - d1 * d1.transpose();
mat3 d2dp0 = -1.0 / l10 * Id2;
mat3 d2dp1 = 1.0 / l10 * Id2;
mat3 d2dp2 = mat3::Zero();

mat3 d1dp0 = 1.0 / lx * Id1 * X(d10);
mat3 d1dp1 = 1.0 / lx * Id1 * X(d0g);
mat3 d1dp2 = 1.0 / lx * Id1 * X(dg1);

mat3 d0dp0 = -X(F.col(2)) * d1p0 + X(F.col(1)) * d2p0;
mat3 d0dp1 = -X(F.col(2)) * d1p1 + X(F.col(1)) * d2p1;
mat3 d0dp2 = -X(F.col(2)) * d1p2;
return {d0dp0, d0dp1, d0dp2, //
       d1dp0, d1dp1, d1dp2, //
       d2dp0, d2dp1, d2dp2};
}

void darboux_gradient() {
 vec3 p0 = __x[__corners_prev[i]];
 vec3 p1 = __x[i];
 vec3 p2 = __x[__corners_next[i]];
 vec3 pgA = __ghost_verts[i];
 vec3 pgB = __ghost_verts[__corners_next[i]];
 mat3 DA = __frames[i];
 mat3 DB = __frames[__corners_next[i]];
}
*/

  void debug() {
    for (int i = 0; i < __corners_next.size(); i++) {
      if (__corners_prev[i] == -1)
        continue;
      vec3 c1 = __x[__corners_next[i]];
      vec3 c0 = __x[i];
      vec3 v = __v[i];

      gg::geometry_logger::line(c0, c1, vec4(0.0, 0.7, 1.0, 1.0));
      gg::geometry_logger::line(c0, c0 + v, vec4(0.8, 0.3, 0.0, 1.0));
      quat u = __u[i];
#if 0
      gg::geometry_logger::line(c0, c0 + u * vec3(0.05, 0.0, 0.0),
                                vec4(1.0, 0.0, 0.0, 1.0));
      gg::geometry_logger::line(c0, c0 + u * vec3(0.0, 0.05, 0.0),
                                vec4(0.0, 1.0, 0.0, 1.0));
      gg::geometry_logger::line(c0, c0 + u * vec3(0.0, 0.0, 0.05),
                                vec4(0.0, 0.0, 1.0, 1.0));
#endif
    }
  }

  real _r = 0.1;

  std::vector<vec3> __M;
  std::vector<vec4> __J;

  std::vector<index_t> __corners_next;
  std::vector<index_t> __corners_prev;

  std::vector<vec3> __x;
  std::vector<vec3> __v; // velocity;
  std::vector<quat> __u;
  std::vector<quat> __o; // omega

  std::vector<real> __l0;
};

} // namespace asawa
} // namespace gaudi
#endif