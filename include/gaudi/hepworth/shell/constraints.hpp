
#ifndef __HEP_SHELL_CONSTRAINTS__
#define __HEP_SHELL_CONSTRAINTS__

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

#include "../projection_constraint.hpp"
namespace gaudi {
namespace hepworth {
namespace shell {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class pinned : public projection_constraint {
public:
  typedef std::shared_ptr<pinned> ptr;

  static ptr create(const std::vector<index_t> &ids, const vec3 &p,
                    const real &w) {
    return std::make_shared<pinned>(ids, p, w);
  }

  pinned(const std::vector<index_t> &ids, const vec3 &p, const real &w)
      : projection_constraint(ids, w), _p(p) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    vec3 q0 = from3(i, q);
    p.block(_id0, 0, 3, 1) = _w * _p;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = this->_ids[0];
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, 3 * i + ax, _w));
    id0 += 3;
  }
  vec3 _p;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class edge_strain : public projection_constraint {
public:
  typedef std::shared_ptr<edge_strain> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l, const real &g) {
    return std::make_shared<edge_strain>(ids, w, l, g);
  }

  edge_strain(const std::vector<index_t> &ids, const real &w, const real &l,
              const real &g)
      : projection_constraint(ids, w), _l(l), _g(g) {}

  virtual void project(const vecX &q, vecX &p) {

    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    vec3 q0 = from3(i, q);
    vec3 q1 = from3(j, q);
    vec3 dq = q1 - q0;

    real l = dq.norm();
    dq /= l;
    // std::cout << l << " " << _l << " " << l / _l << std::endl;
    //   gg::geometry_logger::line(q0, q0 + 0.01 * l / _l * dq,
    //                             vec4(1.0, 0.0, 0.0, 1.0));
    l = std::clamp(l / _l, 0.01, 2.0);

    p.block(_id0, 0, 3, 1) = _w * _g * l * dq;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {

    _id0 = id0;
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, 3 * i + ax, -_w / _l));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, 3 * j + ax, _w / _l));
    id0 += 3;
  }
  real _l = -1.0;
  real _g = 1.0;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class bending : public projection_constraint {
public:
  typedef std::shared_ptr<bending> ptr;

  static ptr create(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights,
                    const std::vector<vec3> &x, const real &w) {
    return std::make_shared<bending>(ids, edge_weights, x, w);
  }

  bending(const std::vector<index_t> &ids,
          const std::vector<real> &edge_weights, const std::vector<vec3> &x,
          const real &w)
      : projection_constraint(ids, w), _edge_weights(edge_weights) {
    _vg = vec3::Zero();
    // todo:: calculate area/cotans here.
    index_t i = this->_ids[0];
    vec3 qi0 = x[i];
    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];
      vec3 qi = x[j];
      vec3 dqi = qi - qi0;
      _vg += K * dqi;
    }
  }

  virtual void project(const vecX &q, vecX &p) {

    vec3 vf = vec3::Zero();
    index_t i = this->_ids[0];
    vec3 qi0 = from3(i, q);
    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];

      vec3 qi = from3(j, q);
      vec3 dqi = qi - qi0;
      // K = 1.0;
      vf += K * dqi;
    }

    vec3 Rvg = vf * _vg.squaredNorm() / vf.squaredNorm();
    if (vf.squaredNorm() > 1e-12)
      p.block(_id0, 0, 3, 1) = _w * Rvg;
    else
      p.block(_id0, 0, 3, 1) = qi0;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = this->_ids[0];
    real Km = 0.0;

    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];
      if (std::isnan(K))
        abort();

      // K = 1.0;
      Km += K;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + ax, 3 * j + ax, _w * K));
    }

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, 3 * i + ax, -_w * Km));

    id0 += 3;
  }

  std::vector<real> _edge_weights;
  vec3 _vg;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class laplacian : public projection_constraint {
public:
  typedef std::shared_ptr<laplacian> ptr;

  static ptr create(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights,
                    const std::vector<vec3> &x, const real &w) {
    return std::make_shared<laplacian>(ids, edge_weights, x, w);
  }

  laplacian(const std::vector<index_t> &ids,
            const std::vector<real> &edge_weights, const std::vector<vec3> &x,
            const real &w)
      : projection_constraint(ids, w), _edge_weights(edge_weights) {
    _vg = vec3::Zero();
    index_t i = this->_ids[0];
    vec3 qi0 = x[i];
    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];
      vec3 qi = x[j];
      vec3 dqi = qi - qi0;
      _vg += K * dqi;
    }
  }

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    // to3(i, p, _w * _vg);
    vec3 vf = vec3::Zero();
    vec3 qi0 = from3(i, q);
    mat3 M = mat3::Zero();

    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];

      vec3 qi = from3(j, q);
      vec3 dqi = qi - qi0;
      // K = 1.0;
      vf += K * dqi;
      M += K * dqi * dqi.transpose();
    }

    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU);
    mat3 U = svd.matrixU();
    vec3 s = svd.singularValues();

    // p.block(_id0, 0, 3, 1) = _w * vf;
    p.block(_id0, 0, 3, 1) = _w * s[0] * U.col(0);

    // p.block(_id0, 0, 3, 1) += _w * _vg;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = this->_ids[0];
    real Km = 0.0;

    for (int k = 1; k < this->_ids.size(); k++) {
      index_t j = this->_ids[k];
      real K = _edge_weights[k - 1];
      // K = 1.0;
      Km += K;
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + ax, 3 * j + ax, _w * K));
    }

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, 3 * i + ax, -_w * Km));

    id0 += 3;
  }
  std::vector<real> _edge_weights;
  vec3 _vg;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class area : public projection_constraint {
public:
  typedef std::shared_ptr<area> ptr;

  static ptr create(const std::vector<index_t> &ids, const std::vector<vec3> &x,
                    real rangeMin, real rangeMax, const real &w) {
    return std::make_shared<area>(ids, x, rangeMin, rangeMax, w);
  }

  area(const std::vector<index_t> &ids, const std::vector<vec3> &x,
       real rangeMin, real rangeMax, real w)
      : projection_constraint(ids, w), _rangeMin(rangeMin),
        _rangeMax(rangeMax) {
    assert(ids.size() == 3);
    mat32 edges, P;
    edges.col(0) = x[_ids[1]] - x[_ids[0]];
    edges.col(1) = x[_ids[2]] - x[_ids[0]];
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    _rest = (P.transpose() * edges).inverse();
    real A = (P.transpose() * edges).determinant() / 2.0f;
    _w *= std::sqrt(std::abs(A));
  }

  virtual void project(const vecX &q, vecX &p) const {
    mat32 edges, P;
    edges.col(0) = from3(_ids[1], q) - from3(_ids[0], q);
    edges.col(1) = from3(_ids[2], q) - from3(_ids[0], q);
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

    std::cout << " --- " << std::endl;
    std::cout << P << std::endl;
    std::cout << F << std::endl;

    p.block(3 * _id0 + 0, 0, 3, 1) = _w * PF.col(0);
    p.block(3 * _id0 + 3, 0, 3, 1) = _w * PF.col(1);
  }

  ///////////////////////////////////////////////////////////////////////////////
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    int n = 2;
    for (int i = 0; i < n; ++i) {
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, 3 * _ids[0] + ax,
                                -_w * (_rest(0, i) + _rest(1, i))));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(
            trip(_id0 + 3 * i + ax, 3 * _ids[1] + ax, _w * _rest(0, i)));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(
            trip(_id0 + 3 * i + ax, 3 * _ids[2] + ax, _w * _rest(1, i)));
    }
    id0 += 3 * n;
  }
  real _rangeMax = 1.0;
  real _rangeMin = 0.0;
  mat2 _rest;
};
#endif

class cross : public projection_constraint {
public:
  typedef std::shared_ptr<cross> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    const real &l) {
    return std::make_shared<cross>(ids, w, l);
  }

  cross(const std::vector<index_t> &ids, const real &w, const real &l)
      : projection_constraint(ids, w), _lambda(l) {}

  virtual void project(const vecX &q, vecX &p) {
    //    i0 B jp
    //  A /| /
    //   / |/ C
    // ip D j0
    index_t i0 = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t j0 = this->_ids[2];
    index_t jp = this->_ids[3];

    vec3 qi0 = q.block(3 * i0, 0, 3, 1);
    vec3 qip = q.block(3 * ip, 0, 3, 1);
    vec3 qj0 = q.block(3 * j0, 0, 3, 1);
    vec3 qjp = q.block(3 * jp, 0, 3, 1);
    vec3 qm = 0.5 * (qi0 + qj0);
    vec3 A = qip - qi0;
    vec3 B = qi0 - qjp;
    vec3 C = qjp - qj0;
    vec3 D = qj0 - qip;
    mat3 AC = A * C.transpose();
    mat3 DB = D * B.transpose();
    mat3 ACDB = DB.inverse() * AC;

    Eigen::JacobiSVD<mat3> svd(AC, Eigen::ComputeFullU);
    mat3 U = svd.matrixU();
    vec3 s = svd.singularValues();
    vec3 Ts = U.col(0).transpose();
    vec3 Ns = U.col(1).transpose();
    vec3 Bs = U.col(2).transpose();

    real lA = A.norm();
    real lB = B.norm();
    real lC = C.norm();
    real lD = D.norm();

    real lp = sqrt(_lambda * lD * lB);

    gg::geometry_logger::line(qm, qm + 0.1 * s[0] * Ts,
                              vec4(1.0, 0.0, 0.0, 1.0));
    gg::geometry_logger::line(qm, qm + 0.1 * s[1] * Ns,
                              vec4(0.0, 1.0, 0.0, 1.0));
    gg::geometry_logger::line(qm, qm + 0.1 * s[2] * Bs,
                              vec4(0.0, 0.0, 1.0, 1.0));
    std::cout << s.transpose() << std::endl;
    // gg::geometry_logger::line(qj0, qj0 - 0.5 * lp * AC,
    //                          vec4(1.0, 0.0, 0.0, 1.0));
    // p.block(3 * i0, 0, 3, 1) += _w * s[0] * Ts;
    // p.block(3 * j0, 0, 3, 1) -= _w * s[0] * Ts;
  }
  virtual void fill_A(std::vector<trip> &triplets) {
    return;
    index_t i = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t j = this->_ids[2];
    index_t jp = this->_ids[3];

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * i + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * i + ax, 3 * ip + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * j + ax, 3 * j + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(3 * j + ax, 3 * jp + ax, _w));
  }
  real _lambda = 1.0;
};

} // namespace shell
} // namespace hepworth
} // namespace gaudi
#endif