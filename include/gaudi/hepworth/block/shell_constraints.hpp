
#ifndef __HEP_SHELL_BLOCK_CONSTRAINTS__
#define __HEP_SHELL_BLOCK_CONSTRAINTS__

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
#include "sim_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class edge_strain : public block_constraint {
public:
  DEFINE_CREATE_FUNC(edge_strain)

  edge_strain(const std::vector<index_t> &ids, const real &l, const real &g,
              const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _l(l), _g(g) {}

  virtual void project(const vecX &q, vecX &p) {

    index_t i = this->_ids[0];
    index_t j = this->_ids[1];

    vec3 q0 = _blocks[0]->get_vec3(i, q);
    vec3 q1 = _blocks[0]->get_vec3(j, q);
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
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[1]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, -_w / _l));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, j + ax, _w / _l));
    id0 += 3;
  }
  real _l = -1.0;
  real _g = 1.0;
};

///////////////////////////////////////////////////////////////////////////////

typedef std::function<void(const index_t &j, const real &k, const vec3 &dq)>
    edge_fcn;

void one_ring_iter(const std::vector<index_t> &ids,
                   const std::vector<real> &weights, //
                   std::vector<vec3> x,              //
                   edge_fcn fcn) {
  index_t i = ids[0];
  vec3 qi0 = x[i];
  for (int k = 1; k < ids.size(); k++) {
    index_t j = ids[k];
    real K = weights[k - 1];
    vec3 qi = x[j];
    vec3 dqi = qi - qi0;
    fcn(j, K, dqi);
  }
}

void one_ring_iter(const std::vector<index_t> &ids,
                   const std::vector<real> &weights, //
                   const vecX &q, sim_block &block,  //
                   edge_fcn fcn) {
  index_t i = ids[0];
  vec3 qi0 = block.get_vec3(i, q);

  for (int k = 1; k < ids.size(); k++) {
    index_t j = ids[k];
    real K = weights[k - 1];
    vec3 qi = block.get_vec3(j, q);
    vec3 dqi = qi - qi0;
    fcn(j, K, dqi);
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

class bending : public block_constraint {
public:
  DEFINE_CREATE_FUNC(bending)

  bending(const std::vector<index_t> &ids,
          const std::vector<real> &edge_weights, const std::vector<vec3> &x,
          const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    _vg = vec3::Zero();
    // todo:: calculate area/cotans here.
    one_ring_iter(ids, edge_weights, x,
                  [&](const index_t &j, const real &k, const vec3 &dq) {
                    _vg += k * dq;
                  });
  }

  virtual void project(const vecX &q, vecX &p) {

    vec3 vf = vec3::Zero();
    one_ring_iter(
        _ids, _edge_weights, q, *_blocks[0],
        [&](const index_t &j, const real &k, const vec3 &dq) { vf += k * dq; });

    vec3 qi0 = _blocks[0]->get_vec3(_ids[0], q);
    vec3 Rvg = vf * _vg.squaredNorm() / vf.squaredNorm();

    if (vf.squaredNorm() > 1e-12)
      p.block(_id0, 0, 3, 1) = _w * Rvg;
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

class laplacian : public block_constraint {
public:
  DEFINE_CREATE_FUNC(laplacian)

  static ptr create(const std::vector<index_t> &ids,
                    const std::vector<real> &edge_weights,
                    const std::vector<vec3> &x, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<laplacian>(ids, edge_weights, x, w, blocks);
  }

  laplacian(const std::vector<index_t> &ids,
            const std::vector<real> &edge_weights, const std::vector<vec3> &x,
            const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {}

  virtual void project(const vecX &q, vecX &p) {
    index_t i = this->_ids[0];
    // to3(i, p, _w * _vg);
    mat3 M = mat3::Zero();
    one_ring_iter(_ids, _edge_weights, q, *_blocks[0],
                  [&](const index_t &j, const real &k, const vec3 &dq) {
                    M += k * dq * dq.transpose();
                  });

    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU);
    mat3 U = svd.matrixU();
    vec3 s = svd.singularValues();

    // p.block(_id0, 0, 3, 1) = _w * s[0] * U.col(0);
    p.block(_id0, 0, 3, 1) = vec3::Zero();
    // p.block(_id0, 0, 3, 1) += _w * _vg;
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    one_ring_fill_A(_id0, _w, _ids, _edge_weights, *_blocks[0], triplets);
    id0 += 3;
  }
  std::vector<real> _edge_weights;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class vec_align : public block_constraint {
public:
  DEFINE_CREATE_FUNC(vec_align)

  vec_align(const std::vector<index_t> &ids,
            const std::vector<real> &edge_weights, const std::vector<vec3> &x,
            const real &w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _edge_weights(edge_weights) {
    _dq0.clear();
    one_ring_iter(ids, edge_weights, x,
                  [&](const index_t &j, const real &k, const vec3 &dq) {
                    _dq0.push_back(dq);
                  });
  }

  virtual void project(const vecX &q, vecX &p) {
    vec3 qi = _blocks[0]->get_vec3(_ids[0], q);

    mat3 M = mat3::Zero();
    int i = 0;
    vec3 dqi = vec3::Zero();
    real K = 0;
    one_ring_iter(_ids, _edge_weights, q, *_blocks[0],
                  [&](const index_t &j, const real &k, const vec3 &dq) {
                    dqi += k * dq;
                    M += k * _dq0[i++] * dq.transpose();
                    // M += k * dq * dq.transpose();
                    K += k;
                  });
    // M /= K;
    // dqi /= K;

    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeThinV | Eigen::ComputeThinU);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    vec3 s = svd.singularValues();
    vec3 N = U.col(2);

    Eigen::AngleAxis<real> aa(0.5, N);
    mat3 R = aa.toRotationMatrix();
    index_t id = 0;
    one_ring_iter(_ids, _edge_weights, q, *_blocks[0],
                  [&](const index_t &j, const real &k, const vec3 &dq) {
                    vec3 rdq = R * dq;
                    p.block(_id0 + id, 0, 3, 1) = _w * k * rdq;
                    // gg::geometry_logger::line(qi, qi + 0.5 * rdq, vec4(1.0,
                    // 0.0, 0.0, 1.0));
                    id += 3;
                  });

    // gg::geometry_logger::line(qi, qi + 0.01 * v0, vec4(1.0, 0.0, 0.0, 1.0));
    // gg::geometry_logger::line(qi, qi + dqr, vec4(0.7, 0.0, 1.0, 1.0));
    // gg::geometry_logger::line(qi, qi + dqi, vec4(1.0, 0.0, 0.0, 1.0));
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t ii = _blocks[0]->get_offset_idx(_ids[0]);

    index_t id = 0;
    for (int j = 1; j < _ids.size(); j++) {
      real K = _edge_weights[j - 1];
      index_t jj = _blocks[0]->get_offset_idx(_ids[j]);
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + ax, ii + ax, -_w * K));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + ax, jj + ax, _w * K));

      id += 3;
      id0 += 3;
    }
  }

  std::vector<vec3> _dq0;
  std::vector<real> _edge_weights;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#if 1
class area : public block_constraint {
public:
  DEFINE_CREATE_FUNC(area)

  area(const std::vector<index_t> &ids, const std::vector<vec3> &x,
       real rangeMin, real rangeMax, real w, std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _rangeMin(rangeMin),
        _rangeMax(rangeMax) {
    assert(ids.size() == 3);
    mat32 edges, P;
    edges.col(0) = x[_ids[1]] - x[_ids[0]];
    edges.col(1) = x[_ids[2]] - x[_ids[0]];
    P.col(0) = edges.col(0).normalized();
    P.col(1) =
        (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    _rest = 0.98 * (P.transpose() * edges).inverse();
    real A = (P.transpose() * edges).determinant() / 2.0f;
    _w *= std::sqrt(std::abs(A));
  }

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

    // p.block(_id0 + 0, 0, 3, 1) = _w * PF.col(0);
    // p.block(_id0 + 3, 0, 3, 1) = _w * PF.col(1);

    p.block(_id0 + 0, 0, 3, 1) = vec3::Zero();
    p.block(_id0 + 3, 0, 3, 1) = vec3::Zero();
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
                                -_w * (_rest(0, i) + _rest(1, i))));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i1 + ax, _w * _rest(0, i)));
      for (int ax = 0; ax < 3; ax++)
        triplets.push_back(trip(_id0 + 3 * i + ax, i2 + ax, _w * _rest(1, i)));
    }
    id0 += 3 * n;
  }
  real _rangeMax = 1.0;
  real _rangeMin = 0.0;
  mat2 _rest;
};
#endif

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

    // gg::geometry_logger::line(q0, q0 + edges.col(0), vec4(0.0, 1.0,
    // 0.0, 1.0)); gg::geometry_logger::line(q0, q0 + edges.col(1),
    // vec4(0.0, 1.0, 0.0, 1.0));

    vec3 N = edges.col(0).cross(edges.col(1)).normalized();
    real d = N.dot(_N0);
    quat qN = quat::FromTwoVectors(N, _N0);
    edges.col(0) = qN * edges.col(0);
    edges.col(1) = qN * edges.col(1);

    // gg::geometry_logger::line(q0, q0 + edges.col(0), vec4(1.0, 0.0,
    // 0.0, 1.0)); gg::geometry_logger::line(q0, q0 + edges.col(1), vec4(1.0,
    // 0.0, 0.0, 1.0));

    // gg::geometry_logger::line(qc, qc + 0.025 * N, vec4(0.0, 1.0, 0.5, 1.0));
    // gg::geometry_logger::line(qc, qc + 0.025 * _N0, vec4(1.0, 0.0,
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

class cross : public block_constraint {
public:
  typedef std::shared_ptr<cross> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &l,
                    const real &w, std::vector<sim_block::ptr> blocks) {
    return std::make_shared<cross>(ids, w, l, blocks);
  }

  cross(const std::vector<index_t> &ids, const real &l, const real &w,
        std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _lambda(l) {}

  virtual void project(const vecX &q, vecX &p) {
    //    i0 B jp
    //  A /| /
    //   / |/ C
    // ip D j0
    index_t i0 = this->_ids[0];
    index_t ip = this->_ids[1];
    index_t j0 = this->_ids[2];
    index_t jp = this->_ids[3];

    vec3 qi0 = _blocks[0]->get_vec3(i0, q);
    vec3 qi1 = _blocks[0]->get_vec3(ip, q);
    vec3 qj0 = _blocks[0]->get_vec3(j0, q);
    vec3 qj1 = _blocks[0]->get_vec3(jp, q);

    vec3 qm = 0.5 * (qi0 + qj0);
    vec3 A = qi1 - qi0;
    vec3 B = qi0 - qj1;
    vec3 C = qj1 - qj0;
    vec3 D = qj0 - qi1;
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
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    index_t ip = _blocks[0]->get_offset_idx(this->_ids[1]);
    index_t j = _blocks[0]->get_offset_idx(this->_ids[2]);
    index_t jp = _blocks[0]->get_offset_idx(this->_ids[3]);

    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(i + ax, i + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(i + ax, ip + ax, _w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(j + ax, j + ax, -_w));
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(j + ax, jp + ax, _w));
  }
  real _lambda = 1.0;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif