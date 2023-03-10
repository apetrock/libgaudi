#ifndef __PENKO_CONSTRAINTS__
#define __PENKO_CONSTRAINTS__

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "Eigen/src/SparseCore/SparseUtil.h"

#include <GaudiGraphics/geometry_logger.h>

#include <gaudi/asawa/coordinate_interface.hpp>
#include <gaudi/asawa/shell.hpp>
#include <gaudi/vec_addendum.h>

#include <gaudi/calder/harmonic_integrators.hpp>

using namespace asawa;

namespace hepworth {
typedef Eigen::Triplet<double> triplet;

template <typename SPACE> class vec_interface {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;

  static coordinate_type from(const vecX &vals, size_t i) {
    return coordinate_type(vals[3 * i + 0], vals[3 * i + 1], vals[3 * i + 2]);
  };

  static void to_plus(const coordinate_type &c, vecX &vals, size_t i) {
    vals[3 * i + 0] += c[0];
    vals[3 * i + 1] += c[1];
    vals[3 * i + 2] += c[2];
  }

  static void to(const coordinate_type &c, vecX &vals, size_t i) {
    vals[3 * i + 0] = c[0];
    vals[3 * i + 1] = c[1];
    vals[3 * i + 2] = c[2];
  }

  // these will be function pointers with a default
  static vecX to(const coordinate_array &positions) {
    vecX out(3 * positions.size());
    int i = 0;

    for (int i = 0; i < positions.size(); i++) {
      coordinate_type p = positions[i];
      to(p, out, i);
    }

    return out;
  }

  static void from(coordinate_array &positions, const vecX &x) {
    for (int i = 0; i < positions.size(); i++) {
      positions[i] = from(x, i);
    }
  }

  static vec3 coord(const int &i, const vecX &vals) { return from(vals, i); }

  static vec3 coord(vertex_ptr v, const vecX &vals) {
    size_t i = v->position_in_set();
    return from(vals, i);
  }

  static vec3 coord(face_vertex_ptr fv, const vecX &vals) {
    vertex_ptr v = fv->vertex();

    return coord(v, vals);
    // size_t i = v->position_in_set();
    // return from(vals, i);
  }

  static vec3 normal(const int &i0, const int &i1, const int &i2,
                     const vecX &vals) {
    vec3 v0 = coord(i0, vals);
    vec3 v1 = coord(i1, vals);
    vec3 v2 = coord(i2, vals);
    return va::calculate_normal(v0, v1, v2);
  }

  static vec3 normal(face_ptr f, const vecX &vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();
    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return va::calculate_normal(v0, v1, v2);
  }

  static real area(face_ptr f, const vecX &vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();

    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return va::calculate_area(v0, v1, v2);
  }

  static vec3 center(face_ptr f, const vecX &vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();

    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return 1.0 / 3.0 * (v0 + v1 + v2);
  }

  static vec3 center(edge_ptr e, const vecX &vals) {

    face_vertex_ptr fv0 = e->v1();
    face_vertex_ptr fv1 = e->v2();

    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    return 1.0 / 2.0 * (v0 + v1);
  }

  template <typename InputIterator>
  static void setFromTriplets(const InputIterator &begin,
                              const InputIterator &end, vecX &x) {
    x.setZero();
    for (InputIterator it(begin); it != end; ++it)
      x(it->col()) += it->value();
  }
};

template <typename SPACE> class constraint {
public:
  M2_TYPEDEFS;
  typedef std::shared_ptr<constraint<SPACE>> ptr;
  static ptr create() { return std::make_shared<constraint<SPACE>>(); }

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  virtual void preprocess(const vecX &v) {}

  virtual void fill_gradient(vecX &G) {}
  virtual void get_grad_triplets(std::vector<triplet> &triplets) {}
  virtual void get_hess_triplets(std::vector<triplet> &triplets) {}
  virtual void update(const vecX &vals) {}
  virtual void init_rest() {}

  virtual real evaluate_constraint() { return 0.0; }
  virtual void set_weight(real w) { _mu = w; };
  real _mu = 1e-4;
};

template <typename SPACE> class point : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<point<SPACE>> ptr;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  static ptr create(size_t i, const vec3 &t, const real &mu = 0.1) {
    return std::make_shared<point<SPACE>>(i, t, mu);
  }

  typedef Eigen::SparseMatrix<real> sparmat;

  point(size_t i, const vec3 &t, const real &mu = 0.1) : _i(i), _t(t) {
    this->set_weight(mu);
  }

  virtual void set_target(const vec3 &t) { _t = t; }

  virtual void update(const vecX &vals) {
    _p = vec_interface<SPACE>::from(vals, _i);
  }

  virtual real evaluate_constraint() {
    coordinate_type dp = _p - _t;
    return this->_mu * dp.transpose() * dp;
  }

  vec3 local_gradient() {
    coordinate_type dp = 2.0 * (_p - _t);
    vec3 G = this->_mu * dp;
    return G;
  }

  mat33 local_hessien() { return 2.0 * this->_mu * mat33::Identity(3, 3); }

  virtual void fill_gradient(vecX &G) {
    vec3 g = local_gradient();
    size_t ii[] = {3 * _i + 0, 3 * _i + 1, 3 * _i + 2};
    for (int i = 0; i < 3; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    mat33 H = local_hessien();
    // std::cout << H << std::endl;
    size_t ii[] = {3 * _i + 0, 3 * _i + 1, 3 * _i + 2};
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        triplets.push_back(triplet(ii[i], ii[j], H(i, j)));
      }
    }
    // std::cout << endl;
  }

  size_t _i = -1;
  vec3 _p, _t;
}; // namespace asawa

template <typename SPACE> class edge_length : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<edge_length<SPACE>> ptr;

  static ptr create(size_t i0, size_t i1) {
    return std::make_shared<edge_length<SPACE>>(i0, i1);
  }

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 6> mat36;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  typedef Eigen::SparseMatrix<real> sparmat;

  real _w = 1e-4;
  edge_length(size_t i0, size_t i1) : _i0(i0), _i1(i1) { this->set_weight(_w); }

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
  }

  virtual real evaluate_constraint() {
    coordinate_type dp = _p1 - _p0;
    return this->_mu * 0.5 * dp.transpose() * dp;
  }

  vec6 local_gradient() {
    coordinate_type dp = this->_mu * (_p1 - _p0);
    vec6 G;
    G.segment(0, 3) = -dp;
    G.segment(3, 3) = dp;
    return G;
  }

  mat66 local_hessien() {
    mat3 I = this->_mu * mat3::Identity(3, 3);
    mat66 H;
    H.block(0, 0, 3, 3) = I;
    H.block(0, 3, 3, 3) = -I;
    H.block(3, 3, 3, 3) = I;
    // symmetric matrix
    H.block(3, 0, 3, 3) = -I;
    return H;
  }

  virtual void fill_gradient(vecX &G) {
    vec6 g = local_gradient();
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    mat66 H = local_hessien();
    // std::cout << H << std::endl;
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        triplets.push_back(triplet(ii[i], ii[j], H(i, j)));
      }
    }
    // std::cout << endl;
  }

  size_t _i0 = -1, _i1 = -1;
  coordinate_type _p0, _p1;

}; // namespace asawa

template <typename SPACE> class edge_stretch : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<edge_stretch<SPACE>> ptr;

  static ptr create(size_t i0, size_t i1) {
    return std::make_shared<edge_stretch<SPACE>>(i0, i1);
  }

  static ptr create(size_t i0, size_t i1, real l) {
    return std::make_shared<edge_stretch<SPACE>>(i0, i1, l);
  }

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 6> mat36;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  typedef Eigen::SparseMatrix<real> sparmat;

  real _w = 1e-2; // convenience local variable
  edge_stretch(size_t i0, size_t i1) : _i0(i0), _i1(i1), _init_rest(true) {
    this->set_weight(_w);
  }

  edge_stretch(size_t i0, size_t i1, real l)
      : _i0(i0), _i1(i1), _lc(l), _init_rest(false) {
    this->set_weight(_w);
  }

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
  }

  virtual void init_rest() {
    if (_init_rest) {
      _lc = (_p0 - _p1).norm();
    }
  }

  virtual bool safe() {
    if ((_p1 - _p0).norm() < 1e-4)
      return false;
    return true;
  }

  virtual real evaluate_constraint() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    real de = l - _lc;
    return 0.5 * this->_mu * de * de;
  }

  vec6 local_gradient() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    coordinate_type dldx = dp * invl;

    if (isnan(invl))
      std::cout << __FUNCTION__ << " NAN invl: " << std::endl;

    real de = l - _lc;
    coordinate_type dCdxi = this->_mu * de * dldx;

    vec6 G;
    G.segment(0, 3) = -dCdxi;
    G.segment(3, 3) = dCdxi;
    return G;
  }
#if 1
  mat66 local_hessien() {

    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    if (isnan(invl))
      std::cout << __FUNCTION__ << " NAN invl: " << std::endl;

    // coordinate_type dldx = dp * invl;
    // real de = l - _lc;
    mat3 I = mat3::Identity(3, 3);
    vec3 dldx = invl * dp;
    mat3 dldlT = dldx * dldx.transpose();
    mat3 d2CDX2 = this->_mu * (-I + _lc * invl * (I - dldlT));

    mat66 H;
    H.block(0, 0, 3, 3) = d2CDX2;
    H.block(0, 3, 3, 3) = -d2CDX2;
    H.block(3, 3, 3, 3) = d2CDX2;
    // symmetric matrix
    H.block(3, 0, 3, 3) = -d2CDX2;
    //    std::cout << H << std::endl;
    //    std::cout << "=================" << std::endl;

    return H;
  }

#else
  mat66 local_hessien() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real il = 1.0 / l;
    coordinate_type dldx = dp * il;

    double deltaRatio = 1.0 - _lc * il;
    mat3 llT = dldx * dldx.transpose();
    mat3 I = mat3::Identity(3, 3);
    mat3 d2CDX2 = _mu * (llT + deltaRatio * (I - llT.transpose()));

    mat66 H;
    H.block(0, 0, 3, 3) = d2CDX2;
    H.block(0, 3, 3, 3) = -d2CDX2;
    H.block(3, 3, 3, 3) = d2CDX2;
    // symmetric matrix
    H.block(3, 0, 3, 3) = -d2CDX2;
    //    std::cout << H << std::endl;
    //    std::cout << "=================" << std::endl;

    return -H;
  }
#endif

  virtual void fill_gradient(vecX &G) {
    if (!safe())
      return;
    vec6 g = local_gradient();
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    if (!safe())
      return;
    mat66 H = local_hessien();
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        triplets.push_back(triplet(ii[i], ii[j], H(i, j)));
      }
    }
    // std::cout << endl;
  }

  size_t _i0 = -1, _i1 = -1;
  coordinate_type _p0, _p1;
  real _lc;
  // real _E = 0.1;
  // real _T = 1;
  bool _init_rest = false;
}; // namespace asawa

template <typename SPACE> class shear : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<shear<SPACE>> ptr;

  static ptr create(size_t i0, size_t i1) {
    return std::make_shared<shear<SPACE>>(i0, i1);
  }

  static ptr create(size_t i0, size_t i1, real l) {
    return std::make_shared<shear<SPACE>>(i0, i1, l);
  }

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 6> mat36;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  typedef Eigen::SparseMatrix<real> sparmat;

  real _w = 1e-2; // convenience local variable
  shear(size_t i0, size_t i1) : _i0(i0), _i1(i1), _init_rest(true) {
    this->set_weight(_w);
  }

  shear(size_t i0, size_t i1, real l)
      : _i0(i0), _i1(i1), _lc(l), _init_rest(false) {
    this->set_weight(_w);
  }

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
    _p2 = vec_interface<SPACE>::from(vals, _i1);
  }

  virtual void init_rest() {
    if (_init_rest) {
      _lc = (_p0 - _p1).norm();
    }
  }

  virtual bool safe() {
    if ((_p1 - _p0).norm() < 1e-4)
      return false;
    return true;
  }

  /*
   *                x3
   *            /|
   *           / |
   *          /  |
   *         /   |  e2
   *        /    |
   *       /     |
   *      /      |
   *   x1 <-------  x2
   *          e1
   *
   */

  virtual real evaluate_constraint() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    real de = l - _lc;
    return 0.5 * this->_mu * de * de;
  }

  vec6 local_gradient() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    coordinate_type dldx = dp * invl;
    real de = l - _lc;
    coordinate_type dCdxi = this->_mu * de * dldx;

    vec6 G;
    G.segment(0, 3) = -dCdxi;
    G.segment(3, 3) = dCdxi;
    return G;
  }
#if 1
  mat66 local_hessien() {

    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    // coordinate_type dldx = dp * invl;
    // real de = l - _lc;
    mat3 I = mat3::Identity(3, 3);
    vec3 dldx = invl * dp;
    mat3 dldlT = dldx * dldx.transpose();
    mat3 d2CDX2 = this->_mu * (-I + _lc * invl * (I - dldlT));

    mat66 H;
    H.block(0, 0, 3, 3) = d2CDX2;
    H.block(0, 3, 3, 3) = -d2CDX2;
    H.block(3, 3, 3, 3) = d2CDX2;
    // symmetric matrix
    H.block(3, 0, 3, 3) = -d2CDX2;
    //    std::cout << H << std::endl;
    //    std::cout << "=================" << std::endl;

    return H;
  }

#else
  mat66 local_hessien() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real il = 1.0 / l;
    coordinate_type dldx = dp * il;

    double deltaRatio = 1.0 - _lc * il;
    mat3 llT = dldx * dldx.transpose();
    mat3 I = mat3::Identity(3, 3);
    mat3 d2CDX2 = _mu * (llT + deltaRatio * (I - llT.transpose()));

    mat66 H;
    H.block(0, 0, 3, 3) = d2CDX2;
    H.block(0, 3, 3, 3) = -d2CDX2;
    H.block(3, 3, 3, 3) = d2CDX2;
    // symmetric matrix
    H.block(3, 0, 3, 3) = -d2CDX2;
    //    std::cout << H << std::endl;
    //    std::cout << "=================" << std::endl;

    return -H;
  }
#endif

  virtual void fill_gradient(vecX &G) {
    if (!safe())
      return;
    vec6 g = local_gradient();
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    if (!safe())
      return;
    mat66 H = local_hessien();
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2,
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2};
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 6; j++) {
        triplets.push_back(triplet(ii[i], ii[j], H(i, j)));
      }
    }
    // std::cout << endl;
  }

  size_t _i0 = -1, _i1 = -1;
  coordinate_type _p0, _p1, _p2;
  real _lc;
  // real _E = 0.1;
  // real _T = 1;
  bool _init_rest = false;
}; // namespace asawa

// -----------------------------------------------------------------------

template <typename SPACE> class cross_ratio : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<cross_ratio<SPACE>> ptr;
  static ptr create(size_t ii, size_t ij, size_t ik, size_t il) {
    return std::make_shared<cross_ratio<SPACE>>(ii, ij, ik, il);
  }

  static ptr create(size_t ii, size_t ij, size_t ik, size_t il, real k) {
    return std::make_shared<cross_ratio<SPACE>>(ii, ij, ik, il, k);
  }

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  typedef Eigen::SparseMatrix<real> sparmat;
  real _w = 1e-6;
  cross_ratio(size_t ii, size_t ij, size_t ik, size_t il)
      : _ii(ii), _ij(ij), _ik(ik), _il(il), _init_rest(true) {
    this->set_weight(_w);
  }
  cross_ratio(size_t ii, size_t ij, size_t ik, size_t il, real k)
      : _ii(ii), _ij(ij), _ik(ik), _il(il), _k(k), _init_rest(false) {
    this->set_weight(_w);
  }

  virtual void init_rest() {
    if (_init_rest)
      _k = cross();
  }
  bool safe() { return true; }
  virtual void update(const vecX &vals) {

    _pi = vec_interface<SPACE>::from(vals, _ii);
    _pj = vec_interface<SPACE>::from(vals, _ij);
    _pk = vec_interface<SPACE>::from(vals, _ik);
    _pl = vec_interface<SPACE>::from(vals, _il);
    //_k = cross();
    // real k = cross();
    // real C = 0.999;
    //_k = C * k + (1.0 - C) * 1.0;
  }

  real cross() {
    coordinate_type dpil = _pi - _pl;
    coordinate_type dpjk = _pj - _pk;
    coordinate_type dplj = _pl - _pj;
    coordinate_type dpki = _pk - _pi;

    real lil = sqrt(va::norm2(dpil) + 1e-2);
    real ljk = sqrt(va::norm2(dpjk) + 1e-2);
    real llj = sqrt(va::norm2(dplj) + 1e-2);
    real lki = sqrt(va::norm2(dpki) + 1e-2);
    real cross = lil * ljk / (llj * lki);
    return cross;
  }

  vec3 gradlnl(vec3 dx) {
    real xTx = va::norm2(dx) + 1e-2;
    return dx / xTx;
  };

  mat33 hesslnl(vec3 dx) {
    real xTx = va::norm2(dx) + 1e-2;
    real ixTx = 1.0 / xTx;
    mat33 I;
    I.setIdentity();
    mat33 H = I * ixTx - dx * dx.transpose() * ixTx * ixTx;

    return H;
  };

  virtual real evaluate_constraint() {
    real C = cross() - _k;
    return 0.5 * this->_mu * C * C;
  }

  vec12 local_gradient() {
    vec3 dpil = _pi - _pl;
    vec3 dpjk = _pj - _pk;
    vec3 dplj = _pl - _pj;
    vec3 dpki = _pk - _pi;
    vec3 gil = gradlnl(dpil);
    vec3 gjk = gradlnl(dpjk);
    vec3 glj = gradlnl(dplj);
    vec3 gki = gradlnl(dpki);

    real X = cross();
    real C = X - _k;

    vec3 gi = (gil + gki);
    vec3 gj = (gjk + glj);
    vec3 gk = -(gjk + gki);
    vec3 gl = -(gil + glj);
    vec12 G;

    G.block(0, 0, 3, 1) = gi;
    G.block(3, 0, 3, 1) = gj;
    G.block(6, 0, 3, 1) = gk;
    G.block(9, 0, 3, 1) = gl;

    // std::cout << X << std::endl;
#if 0
    std::cout << G.transpose() << std::endl;
    std::cout << X << " " << C << std::endl;
    std::cout << _ii << ": " << _pi.transpose() << std::endl;
    std::cout << _ij << ": " << _pj.transpose() << std::endl;
    std::cout << _ik << ": " << _pk.transpose() << std::endl;
    std::cout << _il << ": " << _pl.transpose() << std::endl;

    std::cout << std::endl;
#endif
    return this->_mu * C * X * G;
    // return -_mu * X * G;
  }

  mat1212 local_hessien() {
    vec3 dpil = _pi - _pl;
    vec3 dpjk = _pj - _pk;
    vec3 dplj = _pl - _pj;
    vec3 dpki = _pk - _pi;

    vec3 gLNil = gradlnl(dpil);
    vec3 gLNjk = gradlnl(dpjk);
    vec3 gLNlj = gradlnl(dplj);
    vec3 gLNki = gradlnl(dpki);

    // lnX = Lil - Llj + Ljk - Lki
    real X = cross();
    real C = X - _k;

    vec3 gLNfi = gLNil + gLNki;
    vec3 gLNfj = gLNjk + gLNlj;
    vec3 gLNfk = -(gLNjk + gLNki);
    vec3 gLNfl = -(gLNil + gLNlj);

    vec12 D2;
    D2.block(0, 0, 3, 1) = gLNfi;
    D2.block(3, 0, 3, 1) = gLNfj;
    D2.block(6, 0, 3, 1) = gLNfk;
    D2.block(9, 0, 3, 1) = gLNfl;

    mat33 hLNil = hesslnl(dpil);
    mat33 hLNjk = hesslnl(dpjk);
    mat33 hLNlj = hesslnl(dplj);
    mat33 hLNki = hesslnl(dpki);

    mat33 hLNfii = hLNil + hLNki;
    mat33 hLNfil = hLNil;
    mat33 hLNfik = hLNki;

    mat33 hLNfjj = hLNjk + hLNlj;
    mat33 hLNfjl = hLNlj;
    mat33 hLNfjk = hLNjk;

    mat33 hLNfkk = -(hLNjk + hLNki);
    mat33 hLNfkj = -hLNjk;
    mat33 hLNfki = -hLNki;

    mat33 hLNfll = -(hLNil + hLNlj);
    mat33 hLNfli = -hLNil;
    mat33 hLNflj = -hLNlj;

    mat1212 H;
    H.setZero();
    H.block(0, 0, 3, 3) = hLNfii;
    H.block(0, 6, 3, 3) = hLNfil;
    H.block(0, 9, 3, 3) = hLNfik;

    H.block(3, 3, 3, 3) = hLNfjj;
    H.block(3, 6, 3, 3) = hLNfjk;
    H.block(3, 9, 3, 3) = hLNfjl;

    H.block(6, 6, 3, 3) = hLNfkk;
    H.block(6, 0, 3, 3) = hLNfki;
    H.block(6, 3, 3, 3) = hLNfkj;

    H.block(9, 9, 3, 3) = hLNfll;
    H.block(9, 0, 3, 3) = hLNfli;
    H.block(9, 3, 3, 3) = hLNflj;

    H = this->_mu * (C * X * H + (C * X + X * X) * D2 * D2.transpose());
    return H;
  }

  virtual void fill_gradient(vecX &G) {
    if (!safe())
      return;

    vec12 g = local_gradient();
    // std::cout << " lGNorm: " << g.norm() << std::endl;
    /*
        size_t ii[] = {3 * _ii + 0, 3 * _ii + 1, 3 * _ii + 2, //
                       3 * _ij + 0, 3 * _ij + 1, 3 * _ij + 2, //
                       3 * _ik + 0, 3 * _ik + 1, 3 * _ik + 2, //
                       3 * _il + 0, 3 * _il + 1, 3 * _il + 2};
        for (int i = 0; i < 12; i++) {
          G(ii[i]) += g(i);
        }
    */
    //    j
    // k <|> l
    //    i
    vec3 N0 = va::calculate_normal(_pi, _pk, _pj);
    vec3 N1 = va::calculate_normal(_pj, _pl, _pi);
    vec3 N = (N0 + N1).normalized();
    size_t ii[] = {_ik, _ii, _ij, _il};
    vec3 Ns[] = {N0, N, N, N1};
    vec3 ps[] = {_pk, _pi, _pj, _pl};

    for (int i = 0; i < 4; i++) {
      vec3 gi = g.segment(3 * i, 3);
      // gg::geometry_logger::line(ps[i], ps[i] + 1000 * gi,
      //                           vec4(0.0, 1.0, 0.0, 1.0));

      // gi = va::reject(N, gi);
      // gg::geometry_logger::line(ps[i], ps[i] + 1000 * gi,
      //                           vec4(1.0, 0.0, 0.0, 1.0));
      vec_interface<SPACE>::to_plus(gi, G, ii[i]);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    if (!safe())
      return;

    mat1212 H = local_hessien();

    size_t ii[] = {3 * _ii + 0, 3 * _ii + 1, 3 * _ii + 2, //
                   3 * _ij + 0, 3 * _ij + 1, 3 * _ij + 2, //
                   3 * _ik + 0, 3 * _ik + 1, 3 * _ik + 2, //
                   3 * _il + 0, 3 * _il + 1, 3 * _il + 2};
    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 12; j++) {
        triplets.push_back(triplet(ii[i], ii[j], H(i, j)));
      }
    }
    // std::cout << endl;
  }

  size_t _ii = -1, _ij = -1;
  size_t _ik = -1, _il = -1;

  coordinate_type _pi, _pj, _pk, _pl;
  real _k;

  bool _init_rest = false;
}; // namespace asawa

///////////////////////////////////////////////////////////

template <typename SPACE> class psi {
public:
  M2_TYPEDEFS
  typedef std::shared_ptr<psi<SPACE>> ptr;
  static ptr create() { return std::make_shared<psi<SPACE>>(); }

  virtual void calc(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e,  //
                    real &psi,                                 //
                    real &psi_p,                               //
                    real &psi_pp) = 0;
  virtual void init(const real &p) = 0;
  virtual real rest() { return 0.0; }
};

template <typename SPACE> class tan_psi : public psi<SPACE> {
  M2_TYPEDEFS
public:
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef std::shared_ptr<tan_psi<SPACE>> ptr;
  static ptr create() { return std::make_shared<tan_psi<SPACE>>(); }

  static ptr create(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e) {
    return std::make_shared<tan_psi<SPACE>>(N0, N1, e);
  }

  tan_psi() {}
  tan_psi(const typename SPACE::coordinate_type &N0, //
          const typename SPACE::coordinate_type &N1, //
          const typename SPACE::coordinate_type &e) {
    init(calc_phi(N0, N1, e));
  }

  static real static_calc_phi(const coordinate_type &N0, //
                              const coordinate_type &N1, //
                              const coordinate_type &e) {
    real sint = va::sin(N0, N1, e);
    real cost = va::cos(N0, N1, e);
    if (cost == 0.0)
      return 0.0;
    real tant = sint / cost;
    return tant;
  }

  real calc_phi(const coordinate_type &N0, //
                const coordinate_type &N1, //
                const coordinate_type &e) {
    return static_calc_phi(N0, N1, e);
  };

  void calc_phi_derivatives(const real &phi, //
                            real &phi_p,     //
                            real &phi_pp) {
    real phi2_1 = phi * phi + 1;
    phi_p = phi2_1;
    phi_pp = 0.5 * phi_p * phi;

    return;
  };

  void calc_psi(const real &phi,     //
                const real &phi_hat, //
                const real &phi_p,   //
                const real &phi_pp,  //
                real &psi,           //
                real &psi_p,         //
                real &psi_pp) {
    real dph = phi - phi_hat;
    psi = dph * dph;
    psi_p = 2.0 * dph * phi_p;
    psi_pp = 2.0 * (dph * phi_pp + phi_p * phi_p);
  };

  virtual void calc(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e,  //
                    real &psi,                                 //
                    real &psi_p,                               //
                    real &psi_pp) {
    real phi1 = calc_phi(N0, N1, e);
    real phi_p = 0, phi_pp = 0;
    calc_phi_derivatives(phi1, phi_p, phi_pp);
    calc_psi(phi1, phi0, phi_p, phi_pp, psi, psi_p, psi_pp);
  };

  void init(const real &p0) { phi0 = p0; };

  virtual real rest() { return phi0; }

  real phi0 = 0.0;
};

template <typename SPACE> class sin_psi : public psi<SPACE> {
  M2_TYPEDEFS
public:
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef std::shared_ptr<sin_psi<SPACE>> ptr;
  static ptr create() { return std::make_shared<sin_psi<SPACE>>(); }

  static ptr create(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e) {
    return std::make_shared<sin_psi<SPACE>>(N0, N1, e);
  }

  sin_psi() {}
  sin_psi(const typename SPACE::coordinate_type &N0, //
          const typename SPACE::coordinate_type &N1, //
          const typename SPACE::coordinate_type &e) {}

  real calc_phi(const coordinate_type &N0, //
                const coordinate_type &N1, //
                const coordinate_type &e) {
    // 2.0 tan(thet/2)
    real sint = va::sin(N0, N1, e);
    return sint;
  };

  void calc_phi_derivatives(const coordinate_type &N0, //
                            const coordinate_type &N1, //
                            const coordinate_type &e,  //
                            real &phi_p,               //
                            real &phi_pp) {

    real sint = va::sin(N0, N1, e);
    real cost = va::cos(N0, N1, e);

    phi_p = 0.5 * cost;
    phi_pp = -0.25 * sint;

    return;
  };

  void calc_psi(const real &phi,     //
                const real &phi_hat, //
                const real &phi_p,   //
                const real &phi_pp,  //
                real &psi,           //
                real &psi_p,         //
                real &psi_pp) {
    psi = 0.5 * phi * phi;
    psi_p = phi * phi_p;
    psi_pp = phi * phi_pp + phi_p * phi_p;
  };

  virtual void calc(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e,  //
                    real &psi,                                 //
                    real &psi_p,                               //
                    real &psi_pp) {
    real phi1 = calc_phi(N0, N1, e);
    real phi_p = 0, phi_pp = 0;
    calc_phi_derivatives(N0, N1, e, phi_p, phi_pp);
    calc_psi(phi1, 0.0, phi_p, phi_pp, psi, psi_p, psi_pp);
  };

  virtual void init(const real &p){};
};

template <typename SPACE> class cos_psi : public psi<SPACE> {
  M2_TYPEDEFS
public:
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef std::shared_ptr<cos_psi<SPACE>> ptr;
  static ptr create() { return std::make_shared<cos_psi<SPACE>>(); }

  static ptr create(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e) {
    return std::make_shared<cos_psi<SPACE>>(N0, N1, e);
  }

  cos_psi() {}
  cos_psi(const typename SPACE::coordinate_type &N0, //
          const typename SPACE::coordinate_type &N1, //
          const typename SPACE::coordinate_type &e) {}

  real calc_phi(const coordinate_type &N0, //
                const coordinate_type &N1, //
                const coordinate_type &e) {
    // 2.0 tan(thet/2)
    return va::cos(N0, N1, e);
  };

  void calc_phi_derivatives(const coordinate_type &N0, //
                            const coordinate_type &N1, //
                            const coordinate_type &e,  //
                            real &phi_p,               //
                            real &phi_pp) {

    real sint = va::sin(N0, N1, e);
    real cost = va::cos(N0, N1, e);

    phi_p = -0.5 * sint;
    phi_pp = -0.25 * cost;

    return;
  };

  void calc_psi(const real &phi,     //
                const real &phi_hat, //
                const real &phi_p,   //
                const real &phi_pp,  //
                real &psi,           //
                real &psi_p,         //
                real &psi_pp) {
    psi = 0.5 * phi * phi;
    psi_p = phi * phi_p;
    psi_pp = phi * phi_pp + phi_p * phi_p;
  };

  virtual void calc(const typename SPACE::coordinate_type &N0, //
                    const typename SPACE::coordinate_type &N1, //
                    const typename SPACE::coordinate_type &e,  //
                    real &psi,                                 //
                    real &psi_p,                               //
                    real &psi_pp) {
    real phi1 = calc_phi(N0, N1, e);
    real phi_p = 0, phi_pp = 0;
    calc_phi_derivatives(N0, N1, e, phi_p, phi_pp);
    calc_psi(phi1, 0.0, phi_p, phi_pp, psi, psi_p, psi_pp);
  };

  virtual void init(const real &p){};
};

//////////////////////////////////////////////////////////
// edge_bend
//////////////////////////////////////////////////////////

template <typename SPACE> class edge_bend : public constraint<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<edge_bend<SPACE>> ptr;

  static ptr create(const int &i0, const int &i1, //
                    const int &i2, const int &i3) {
    return std::make_shared<edge_bend<SPACE>>(i0, i1, i2, i3);
  }

  real _w = 1e-3;
  edge_bend(const int &i0, const int &i1, //
            const int &i2, const int &i3)
      : _i0(i0), _i1(i1), _i2(i2), _i3(i3) {
    this->constraint<SPACE>::set_weight(_w);
  }

  // need to cache areas, normals from faces

  //         x00
  //         /\
  //        /  \
  //    e02/    \e01
  //      /  f0  \
  //     /        \
  //    /          \
  //    ------------
  // x1     /e0/      x2
  //    ------------
  //    \          /
  //     \   f1   /
  //      \      /
  //   e12p\    /e11
  //        \  /
  //         \/
  //         x10
  //
  // Edge orientation: e0,e1,e2 point away from x0
  //                      e3,e4 point away from x1

  virtual bool safe() {
    coordinate_type x00 = _p0;
    coordinate_type x1 = _p1;
    coordinate_type x2 = _p2;
    coordinate_type x10 = _p3;
    A0 = va::calculate_area(x00, x1, x2);
    A1 = va::calculate_area(x10, x2, x1);

    if (A0 < 1e-10)
      return false;
    if (A1 < 1e-10)
      return false;
    if (A0 / A1 > 4.0)
      return false;
    if (A1 / A0 > 4.0)
      return false;

    if (std::isinf(_psi_fcn->rest())) {
      return false;
    }

    if (std::isnan(_psi_fcn->rest()))
      return false;

    if (fabs(_psi_fcn->rest()) > 4.0)
      return false;

    return true;
  }

  void update_intermediates() {
    coordinate_type x00 = _p0;
    coordinate_type x1 = _p1;
    coordinate_type x2 = _p2;
    coordinate_type x10 = _p3;

    N0 = va::calculate_normal(x00, x1, x2);
    N1 = va::calculate_normal(x10, x2, x1);
    A0 = va::calculate_area(x00, x1, x2);
    A1 = va::calculate_area(x10, x2, x1);

    e0 = x2 - x1;
    e02 = x00 - x1;
    e01 = x00 - x2;
    e12 = x10 - x1;
    e11 = x10 - x2;

    l0 = e0.norm();
    l01 = e01.norm();
    l02 = e02.norm();
    l11 = e11.norm();
    l12 = e12.norm();

    e0.normalize();
    e02.normalize();
    e01.normalize();
    e12.normalize();
    e11.normalize();

    ih00 = 0.5 * l0 / A0;
    ih01 = 0.5 * l01 / A0;
    ih02 = 0.5 * l02 / A0;

    ih10 = 0.5 * l0 / A1;
    ih11 = 0.5 * l11 / A1;
    ih12 = 0.5 * l12 / A1;

    cos01 = va::dot(e0, e02);
    cos02 = -va::dot(e0, e01);
    cos11 = va::dot(e0, e12);
    cos12 = -va::dot(e0, e11);
  }

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
    _p2 = vec_interface<SPACE>::from(vals, _i2);
    _p3 = vec_interface<SPACE>::from(vals, _i3);
  }

  virtual void set_weight(real k) { constraint<SPACE>::set_weight(k); }

  void set_weight(real h, real v, real Y) {
    real k = _k = Y * h * h * h / 24.0 / (1.0 - v);
    this->set_weight(k);
  }

  real calc_a(real k, real a0, real a1, real l) {
    real denom = a0 + a1;
    real C = 3.0 * k * l * l / denom;
    C = denom = a0 < 1e-5 ? 0 : C;
    C = denom = a1 < 1e-5 ? 0 : C;
    C = a0 / a1 > 4.0 ? 0 : C;
    C = a1 / a0 > 4.0 ? 0 : C;

    return C;
  }

  virtual void init_rest() {
    update_intermediates();
    _ai = this->calc_a(this->_mu, A0, A1, l0);
    _psi_fcn = tan_psi<SPACE>::create(N0, N1, e0);
    //_psi_fcn = cos_psi<SPACE>::create(N0, N1, e0);
    //_psi_fcn = sin_psi<SPACE>::create(N0, N1, e0);
  }

  virtual real evaluate_constraint() { return _psi; }

  virtual void preprocess(const vecX &v) {
    _psi_fcn->calc(N0, N1, e0, _psi, _psi_p, _psi_pp);
    _psi = _ai * _psi;
    _psi_p = _ai * _psi_p;
    _psi_pp = _ai * _psi_pp;

    update_intermediates();
    _grad = vec12::Zero();
    _hess = mat1212::Zero();
  }

  vec12 d_theta() {
    _grad.segment(0, 3) = -ih00 * N0;
    _grad.segment(3, 3) = cos02 * ih01 * N0 + cos12 * ih11 * N1;
    _grad.segment(6, 3) = cos01 * ih02 * N0 + cos11 * ih12 * N1;
    _grad.segment(9, 3) = -ih10 * N1;
    return _grad;
  }

  vec12 local_gradient() {
    vec12 dtheta = d_theta();
    if (isnan(dtheta.norm()))
      std::cout << __FUNCTION__ << ":" << __LINE__
                << " NAN dtheta: " << dtheta.transpose() << std::endl;

    if (isnan(_psi_p))
      std::cout << __FUNCTION__ << ":" << __LINE__ << " NAN _psi_p "
                << std::endl;
    assert(!isnan(_psi_p * dtheta.norm()));
    return _psi_p * dtheta;
  }

  mat1212 local_hessien() {
    coordinate_type m00 = va::cross(e0, N0);
    coordinate_type m10 = va::cross(e0, N1);

    coordinate_type m01 = va::cross(e01, N0);
    coordinate_type m02 = va::cross(e02, N0);
    coordinate_type m11 = va::cross(e11, N1);
    coordinate_type m12 = va::cross(e12, N1);

    mat33 M00 = N0 * m00.transpose();
    mat33 M01 = N0 * m01.transpose();
    mat33 M02 = N0 * m02.transpose();

    mat33 M10 = N1 * m10.transpose();
    mat33 M11 = N1 * m11.transpose();
    mat33 M12 = N1 * m12.transpose();

    auto S = [](mat33 A) { return A + A.transpose(); };

    int i0 = 0;
    int i1 = 3;
    int i2 = 6;
    int i3 = 9;

    mat33 M00t = M00.transpose();
    mat33 M01t = M01.transpose();
    mat33 M02t = M02.transpose();

    mat33 M10t = M10.transpose();
    mat33 M11t = M11.transpose();
    mat33 M12t = M12.transpose();

    mat33 N00 = M00 / l0 / l0;
    mat33 N10 = M10 / l0 / l0;
    mat33 Q00 = ih00 * ih00 * M00;
    mat33 Q01 = ih00 * ih01 * M01;
    mat33 Q02 = ih00 * ih02 * M02;

    mat33 Q10 = ih10 * ih10 * M10;
    mat33 Q11 = ih10 * ih11 * M11;
    mat33 Q12 = ih10 * ih12 * M12;

    mat33 P011 = ih01 * ih01 * cos01 * M01t;
    mat33 P022 = ih02 * ih02 * cos02 * M02t;
    mat33 P010 = ih01 * ih00 * cos01 * M00t;
    mat33 P020 = ih02 * ih00 * cos02 * M00t;
    mat33 P012 = ih01 * ih02 * cos01 * M02t;
    mat33 P021 = ih02 * ih02 * cos02 * M01t;

    mat33 P111 = ih11 * ih11 * cos11 * M11t;
    mat33 P122 = ih12 * ih12 * cos12 * M12t;
    mat33 P110 = ih11 * ih10 * cos11 * M10t;
    mat33 P120 = ih12 * ih10 * cos12 * M10t;
    mat33 P112 = ih11 * ih12 * cos11 * M12t;
    mat33 P121 = ih12 * ih12 * cos12 * M11t;

    mat1212 &H = _hess;
    H.block(i0, i0, 3, 3) = -S(Q00);
    H.block(i1, i1, 3, 3) = S(P011) - N00 + //
                            S(P111) - N10;
    H.block(i2, i2, 3, 3) = S(P022) - N00 + //
                            S(P122) - N10;
    H.block(i3, i3, 3, 3) = -S(Q10);

    mat33 H10 = P010 - Q01;
    mat33 H20 = P020 - Q02;
    mat33 H13 = P110 - Q11;
    mat33 H23 = P120 - Q12;

    mat33 H12 = P012 + P021.transpose() + N00 + //
                P112 + P121.transpose() + N10;

    H.block(i1, i0, 3, 3) = H10;
    H.block(i0, i1, 3, 3) = H10.transpose();

    H.block(i2, i0, 3, 3) = H20;
    H.block(i0, i2, 3, 3) = H20.transpose();

    H.block(i1, i3, 3, 3) = H13;
    H.block(i3, i1, 3, 3) = H13.transpose();

    H.block(i2, i3, 3, 3) = H23;
    H.block(i3, i2, 3, 3) = H23.transpose();

    H.block(i1, i2, 3, 3) = H12;
    H.block(i2, i1, 3, 3) = H12.transpose();
    vec12 dtheta = d_theta();

    H = _psi_p * H + _psi_pp * dtheta * dtheta.transpose();
    if (isnan(H.norm()))
      std::cout << __FUNCTION__ << ":" << __LINE__
                << " NAN H: " << H.transpose() << std::endl;
    //    if (_psi_p > 0)
    //      std::cout << _ai << " " << _psi_p << " " << _psi_pp << std::endl;
    assert(!isnan(H.norm()));
    return -H;
  }

  virtual void fill_gradient(vecX &G) {
    if (!safe())
      return;
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2, //
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2, //
                   3 * _i2 + 0, 3 * _i2 + 1, 3 * _i2 + 2, //
                   3 * _i3 + 0, 3 * _i3 + 1, 3 * _i3 + 2};
    vec12 g = local_gradient();
    for (int i = 0; i < 12; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    if (!safe())
      return;
#if 1
    size_t ii[] = {3 * _i0 + 0, 3 * _i0 + 1, 3 * _i0 + 2, //
                   3 * _i1 + 0, 3 * _i1 + 1, 3 * _i1 + 2, //
                   3 * _i2 + 0, 3 * _i2 + 1, 3 * _i2 + 2, //
                   3 * _i3 + 0, 3 * _i3 + 1, 3 * _i3 + 2};

    const mat1212 Hi = local_hessien();
    for (int i = 0; i < 12; i++) {
      for (int j = 0; j < 12; j++) {
        triplets.push_back(triplet(ii[i], ii[j], Hi(i, j)));
      }
    }
#endif
  }

  real _k = 0.0;
  real _ai;

  typename psi<SPACE>::ptr _psi_fcn;
  real _psi;    // psi derivative
  real _psi_p;  // psi derivative
  real _psi_pp; // psi double derivative

  vec12 _grad;
  mat1212 _hess;

  coordinate_type e0, e02, e01, e12, e11;
  real l0, l01, l02, l11, l12; // inverse_edge_lengths

  real A0, A1;
  coordinate_type N0, N1;

  real cos01, cos02, cos11, cos12;
  real ih00, ih01, ih02;
  real ih10, ih11, ih12;

  size_t _i0 = -1, _i1 = -1, _i2 = -1, _i3 = -1;
  coordinate_type _p0, _p1, _p2, _p3;
};
#if 1
//////////////////////////////////////////////////////////
// edge_mem_bend
//////////////////////////////////////////////////////////

template <typename SPACE> class willmore_bend : public edge_bend<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<willmore_bend<SPACE>> ptr;

  static ptr create(const int &i0, const int &i1, //
                    const int &i2, const int &i3) {
    return std::make_shared<willmore_bend<SPACE>>(i0, i1, i2, i3);
  }

  real _w = 1e-3;

  willmore_bend(const int &i0, const int &i1, //
                const int &i2, const int &i3)
      : edge_bend<SPACE>(i0, i1, i2, i3) {
    this->set_weight(_w);
  }

  virtual void init_rest() {
    this->update_intermediates();
    this->_ai = this->calc_a(this->_mu, this->A0, this->A1, this->l0);
    // this->_psi_fcn = tan_psi<SPACE>::create(N0, N1, e0);
    //_psi_fcn = cos_psi<SPACE>::create(N0, N1, e0);
    this->_psi_fcn = sin_psi<SPACE>::create(this->N0, this->N1, this->e0);
  }
};

#endif

#if 1
//////////////////////////////////////////////////////////
// edge_mem_bend
//////////////////////////////////////////////////////////

template <typename SPACE> class edge_mem_bend : public edge_bend<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<edge_mem_bend<SPACE>> ptr;

  static ptr create(const int &i0, const int &i1, //
                    const int &i2, const int &i3, const real &phi0) {
    return std::make_shared<edge_mem_bend<SPACE>>(i0, i1, i2, i3, phi0);
  }

  real _w0 = 1e-3;
  real _w1 = 1e-3;

  edge_mem_bend(const int &i0, const int &i1, //
                const int &i2, const int &i3, const real &phi0)
      : edge_bend<SPACE>(i0, i1, i2, i3), _phi0(phi0) {
    this->set_weight(_w0, _w1);
  }

  void set_weight(real w0, real w1) {
    _k0 = w0;
    _k1 = w1;
  }

  virtual bool safe() {
    coordinate_type x00 = this->_p0;
    coordinate_type x1 = this->_p1;
    coordinate_type x2 = this->_p2;
    coordinate_type x10 = this->_p3;
    this->A0 = va::calculate_area(x00, x1, x2);
    this->A1 = va::calculate_area(x10, x2, x1);
    /*
        if (A0 / A1 > 4.0)
          return false;
        if (A1 / A0 > 4.0)
          return false;

        if (std::isinf(_psi_fcn->rest())) {
          return false;
        }

        if (std::isnan(_psi_fcn->rest()))
          return false;

        if (fabs(_psi_fcn->rest()) > 4.0)
          return false;
    */
    return true;
  }

  virtual void init_rest() {
    this->update_intermediates();

    _ai0 = this->calc_a(_k0, this->A0, this->A1, this->l0);
    _ai1 = this->calc_a(_k1, this->A0, this->A1, this->l0);

    _psi_fcn0 = tan_psi<SPACE>::create();
    _psi_fcn0->init(_phi0);

    _psi_fcn1 = tan_psi<SPACE>::create(this->N0, this->N1, this->e0);
  }

  virtual real evaluate_constraint() { return this->_psi; }

  virtual void preprocess(const vecX &v) {

    coordinate_type N0 = this->N0;
    coordinate_type N1 = this->N1;
    coordinate_type e0 = this->e0;
    real psi0, psi0_p, psi0_pp;
    _psi_fcn0->calc(N0, N1, e0, psi0, psi0_p, psi0_pp);
    real psi1, psi1_p, psi1_pp;
    _psi_fcn1->calc(N0, N1, e0, psi1, psi1_p, psi1_pp);

    this->_psi = _ai0 * psi0 + _ai1 * psi1;
    this->_psi_p = _ai0 * psi0_p + _ai1 * psi1_p;
    this->_psi_pp = _ai0 * psi0_pp + _ai1 * psi1_pp;
    this->update_intermediates();
    this->_grad = vec12::Zero();
    this->_hess = mat1212::Zero();
  }

  real _k0 = 0.0;
  real _k1 = 0.0;
  real _ai0;
  real _ai1;

  real _phi0 = 0.0;

  typename psi<SPACE>::ptr _psi_fcn0;
  typename psi<SPACE>::ptr _psi_fcn1;
};

#endif

//////////////////////////////////////////////////////////
// total_bend
//////////////////////////////////////////////////////////

template <typename SPACE> class total_bend : public constraint<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<total_bend<SPACE>> ptr;

  static ptr create(surf_ptr surf) {
    return std::make_shared<total_bend<SPACE>>(surf);
  }

  total_bend(surf_ptr surf) : _surf(surf) {
    coordinate_array positions = asawa::ci::get_coordinates<SPACE>(_surf);
    _positions = vec_interface<SPACE>::to(positions);
    init();
  }

  vec3 coord(face_vertex_ptr fv, const vecX &vals) {
    return vec_interface<SPACE>::coord(fv, vals);
  }

  vec3 normal(face_ptr f, const vecX &vals) {
    return vec_interface<SPACE>::normal(f, vals);
  }

  real area(face_ptr f, const vecX &vals) {
    real A = vec_interface<SPACE>::area(f, vals);
    return A;
  }

  std::function<real(const coordinate_type &, //
                     const coordinate_type &, //
                     const coordinate_type)>

      calc_phi = [](const coordinate_type &N0, //
                    const coordinate_type &N1, //
                    const coordinate_type &e) {
        real sint = 0.5 * va::norm(vec3(N1 - N0));
        real cost = 0.5 * va::norm(vec3(N1 + N0));
        real tant = va::sgn(va::determinant(N0, N1, e)) * sint / cost;
        return tant;
      };

  std::function<void(const real &, //
                     real &, real &)>

      calc_phi_derivatives = [](const real &phi, //
                                real &phi_p,     //
                                real &phi_pp) {
        real phi2_1 = phi * phi + 1;
        phi_p = 0.5 * phi2_1;
        phi_pp = phi_p * phi;

        return;
      };

  std::function<void(const real &, //
                     const real &, //
                     const real &, //
                     const real &, //
                     const real &, //
                     real &, real &, real &, real &)>
      calc_psi = [](const real &phi,     //
                    const real &phi_hat, //
                    const real &phi_p,   //
                    const real &phi_pp,  //
                    const real &k,       //
                    const real &ai,      //
                    real &psi,           //
                    real &psi_p,         //
                    real &psi_pp) {
        real dph = phi - phi_hat;
        real C = 2.0 * k * ai;

        psi = 0.5 * C * dph * dph;

        psi_p = C * dph * phi_p;

        psi_pp = C * (dph * phi_pp + phi_p * phi_p);

        return psi_pp;
      };

  void init() {
    std::vector<edge_ptr> edges = _surf->get_edges();
    std::vector<face_ptr> faces = _surf->get_faces();

    _phi0 = std::vector<real>(edges.size(), 0);
    _ai = std::vector<real>(edges.size(), 0);
    real h = 0.01;
    real v = 0.5;
    real Y = 0.1;
    //_k = Y * h * h * h / 24.0 / (1.0 - v);
    _k = 0.0001;
    int i = 0;

    for (auto f : faces) {
      real A = area(f, _positions);
      _Abar += A;
    }

    _Abar /= float(faces.size());
    for (auto e : edges) {
      face_vertex_ptr fv0 = e->v1();
      face_vertex_ptr fv1 = e->v2();
      face_ptr f0 = fv0->face();
      face_ptr f1 = fv1->face();

      real A0 = area(f0, _positions);
      real A1 = area(f1, _positions);

      coordinate_type c0 = coord(fv0, _positions);
      coordinate_type c1 = coord(fv1, _positions);

      coordinate_type N0 = normal(f0, _positions);
      coordinate_type N1 = normal(f1, _positions);

      coordinate_type e_vec = c1 - c0;
      _phi0[i] = calc_phi(N0, N1, e_vec);
#if 0
      if (i == 0)
        _phi0[i] *= 1.5;
#endif
      real denom = A0 + A1;
      denom = max(denom, 2e-5);
      _ai[i] = 3.0 * va::norm2(e_vec) / denom;
      i++;
    }

    _e_len = std::vector<real>(edges.size());  // inverse_edge_lengths
    _phi1 = std::vector<real>(edges.size());   // phi1
    _psi = std::vector<real>(edges.size());    // psi derivative
    _psi_p = std::vector<real>(edges.size());  // psi derivative
    _psi_pp = std::vector<real>(edges.size()); // psi double derivative
    _grad = std::vector<vec12>(edges.size(), vec12::Zero());
    _hess0 = std::vector<mat99>(faces.size(),
                                mat99::Zero()); // hess
    _hess1 = std::vector<mat1212>(edges.size(),
                                  mat1212::Zero()); // hess
  }

  virtual real evaluate_constraint() {
    std::vector<edge_ptr> edges = _surf->get_edges();
    int i = 0;
    real C = 0.0;
    for (auto e : edges) {
      C += _psi[i];
      i++;
    }
    return C;
  }

  void init_edges() {
    edge_array edges = _surf->get_edges();

    int i = 0;
    // need to cache areas, normals from faces
    for (auto e : edges) {

      face_vertex_ptr fv0 = e->v1();
      face_vertex_ptr fv1 = e->v2();
      face_ptr f0 = fv0->face();
      face_ptr f1 = fv1->face();

      coordinate_type c0 = coord(fv0, _positions);
      coordinate_type c1 = coord(fv1, _positions);

      coordinate_type N0 = normal(f0, _positions);
      coordinate_type N1 = normal(f1, _positions);

      coordinate_type e_vec = c1 - c0;

      _phi1[i] = calc_phi(N0, N1, e_vec);
      _e_len[i] = va::norm(e_vec);
      real phi_p = 0.0, phi_pp = 0.0;

      calc_phi_derivatives(_phi1[i], phi_p, phi_pp);
      calc_psi(_phi1[i], _phi0[i], phi_p, phi_pp, _k, _ai[i], _psi[i],
               _psi_p[i], _psi_pp[i]);

      _grad[i] = vec12::Zero();
      _hess1[i] = mat1212::Zero(); // inverse_edge_lengths

      i++;
    }

    i = 0;
    face_array faces = _surf->get_faces();
    _Abar = 0;
    for (auto f : faces) {
      _hess0[i++] = mat99::Zero(); // inverse_edge_lengths
      real A = area(f, _positions);
      _Abar += A;
    }
    //_Abar /= float(faces.size());
  }

  void finish_grad_hess() {
    int i = 0;
    edge_array edges = _surf->get_edges();
    for (auto e : edges) {
      const vec12 &dt = _grad[i];
      const real &p = _psi_p[i];
      const real &pp = _psi_pp[i];
      _grad[i] = p * dt;
      _hess1[i] += pp * dt * dt.transpose();
      //_hess1[i] = p * mat1212::Identity() + pp * dt * dt.transpose();

      //_hess1[i] *= -1;
      //_hess1[i] = dt * dt.transpose();

      i++;
    }
  }

  virtual void preprocess(const vecX &v) { preprocess_edges(); }

  virtual void preprocess_edges() {
    init_edges();
    edge_array edges = _surf->get_edges();

    int i = 0;
    // need to cache areas, normals from faces

    //         x00
    //         /\
    //        /  \
    //    e02/    \e01
    //      /  f0  \
    //     /        \
    //    /          \
    //    ------------
    // x1     /e0/      x2
    //    ------------
    //    \          /
    //     \   f1   /
    //      \      /
    //   e12p\    /e11
    //        \  /
    //         \/
    //         x10
    //
    // Edge orientation: e0,e1,e2 point away from x0
    //                      e3,e4 point away from x1
    int ii = 0;
    for (auto e : edges) {
      int i = e->position_in_set();

      const real &p = _psi_p[i];
      face_vertex_ptr fv1 = e->v1();
      face_vertex_ptr fv2 = e->v2();
      face_vertex_ptr fv00 = e->v1()->prev();
      face_vertex_ptr fv10 = e->v2()->prev();
      /*
            if (ii == 1)
              break;
            if (fv00->vertex()->position_in_set() != 0)
              continue;
      */
      ii++;
      face_ptr f0 = e->v1()->face();
      face_ptr f1 = e->v2()->face();

      coordinate_type x00 = coord(fv00, _positions);
      coordinate_type x1 = coord(fv1, _positions);
      coordinate_type x2 = coord(fv2, _positions);
      coordinate_type x10 = coord(fv10, _positions);

      coordinate_type N0 = normal(f0, _positions);
      coordinate_type N1 = normal(f1, _positions);
      real A0 = area(f0, _positions);
      real A1 = area(f1, _positions);

      coordinate_type e0 = x2 - x1;
      coordinate_type e02 = x00 - x1;
      coordinate_type e01 = x00 - x2;
      coordinate_type e12 = x10 - x1;
      coordinate_type e11 = x10 - x2;

      real l0 = e0.norm();
      real l01 = e01.norm();
      real l02 = e02.norm();
      real l11 = e11.norm();
      real l12 = e12.norm();

      e0.normalize();
      e02.normalize();
      e01.normalize();
      e12.normalize();
      e11.normalize();

      real ih00 = 0.5 * l0 / A0;
      real ih01 = 0.5 * l01 / A0;
      real ih02 = 0.5 * l02 / A0;

      real ih10 = 0.5 * l0 / A1;
      real ih11 = 0.5 * l11 / A1;
      real ih12 = 0.5 * l12 / A1;

      real cos01 = va::dot(e0, e02);
      real cos02 = -va::dot(e0, e01);
      real cos11 = va::dot(e0, e12);
      real cos12 = -va::dot(e0, e11);

      vec12 &G = _grad[i];

      G.segment(0, 3) = -ih00 * N0;
      G.segment(3, 3) = cos02 * ih01 * N0 + cos12 * ih11 * N1;
      G.segment(6, 3) = cos01 * ih02 * N0 + cos11 * ih12 * N1;
      G.segment(9, 3) = -ih10 * N1;
      real il0 = 1.0 / l0;

      coordinate_type m00 = va::cross(e0, N0);
      coordinate_type m10 = va::cross(e0, N1);

      coordinate_type m01 = va::cross(e01, N0);
      coordinate_type m02 = va::cross(e02, N0);
      coordinate_type m11 = va::cross(e11, N1);
      coordinate_type m12 = va::cross(e12, N1);

      mat33 M00 = N0 * m00.transpose();
      mat33 M01 = N0 * m01.transpose();
      mat33 M02 = N0 * m02.transpose();

      mat33 M10 = N1 * m10.transpose();
      mat33 M11 = N1 * m11.transpose();
      mat33 M12 = N1 * m12.transpose();

      auto S = [](mat33 A) { return A + A.transpose(); };

      mat1212 &H = _hess1[i];

      int i0 = 0;
      int i1 = 3;
      int i2 = 6;
      int i3 = 9;

      mat33 M00t = M00.transpose();
      mat33 M01t = M01.transpose();
      mat33 M02t = M02.transpose();

      mat33 M10t = M10.transpose();
      mat33 M11t = M11.transpose();
      mat33 M12t = M12.transpose();

      mat33 N00 = il0 * il0 * M00;
      mat33 N10 = il0 * il0 * M10;
      mat33 Q00 = ih00 * ih00 * M00;
      mat33 Q01 = ih00 * ih01 * M01;
      mat33 Q02 = ih00 * ih02 * M02;

      mat33 Q10 = ih10 * ih10 * M10;
      mat33 Q11 = ih10 * ih11 * M11;
      mat33 Q12 = ih10 * ih12 * M12;

      mat33 P011 = ih01 * ih01 * cos01 * M01t;
      mat33 P022 = ih02 * ih02 * cos02 * M02t;
      mat33 P010 = ih01 * ih00 * cos01 * M00t;
      mat33 P020 = ih02 * ih00 * cos02 * M00t;
      mat33 P012 = ih01 * ih02 * cos01 * M02t;
      mat33 P021 = ih02 * ih02 * cos02 * M01t;

      mat33 P111 = ih11 * ih11 * cos11 * M11t;
      mat33 P122 = ih12 * ih12 * cos12 * M12t;
      mat33 P110 = ih11 * ih10 * cos11 * M10t;
      mat33 P120 = ih12 * ih10 * cos12 * M10t;
      mat33 P112 = ih11 * ih12 * cos11 * M12t;
      mat33 P121 = ih12 * ih12 * cos12 * M11t;

      H.block(i0, i0, 3, 3) = -S(Q00);
      H.block(i1, i1, 3, 3) = S(P011) - N00 + //
                              S(P111) - N10;
      H.block(i2, i2, 3, 3) = S(P022) - N00 + //
                              S(P122) - N10;
      H.block(i3, i3, 3, 3) = -S(Q10);

      mat33 H10 = P010 - Q01;
      mat33 H20 = P020 - Q02;
      mat33 H13 = P110 - Q11;
      mat33 H23 = P120 - Q12;

      mat33 H12 = P012 + P021.transpose() + N00 + //
                  P112 + P121.transpose() + N10;

      H.block(i1, i0, 3, 3) = H10;
      H.block(i0, i1, 3, 3) = H10.transpose();

      H.block(i2, i0, 3, 3) = H20;
      H.block(i0, i2, 3, 3) = H20.transpose();

      H.block(i1, i3, 3, 3) = H13;
      H.block(i3, i1, 3, 3) = H13.transpose();

      H.block(i2, i3, 3, 3) = H23;
      H.block(i3, i2, 3, 3) = H23.transpose();

      H.block(i1, i2, 3, 3) = H12;
      H.block(i2, i1, 3, 3) = H12.transpose();

      _hess1[i] = p * H;
    }

    finish_grad_hess();
  }

  virtual void preprocess_tries() {
    face_array faces = _surf->get_faces();
    init_edges();

    int i = 0;
    for (auto f : faces) {

      real A = area(f, _positions);
      coordinate_type N = normal(f, _positions);

      face_vertex_ptr fv0 = f->fbegin();
      face_vertex_ptr fv1 = fv0->next();
      face_vertex_ptr fv2 = fv1->next();

      // coordinates are taken from point opposite the edge

      coordinate_type p0 = coord(fv0, _positions);
      coordinate_type p1 = coord(fv1, _positions);
      coordinate_type p2 = coord(fv2, _positions);

      edge_ptr ep0 = fv1->edge();
      edge_ptr ep1 = fv2->edge();
      edge_ptr ep2 = fv0->edge();

      int o0 = ep0->side(fv1);
      int o1 = ep1->side(fv2);
      int o2 = ep2->side(fv0);

      int ie0 = ep0->position_in_set();
      int ie1 = ep1->position_in_set();
      int ie2 = ep2->position_in_set();

      real psip0 = _psi_p[ie0];
      real psip1 = _psi_p[ie1];
      real psip2 = _psi_p[ie2];

      real l0 = _e_len[ie0];
      real l1 = _e_len[ie1];
      real l2 = _e_len[ie2];

      real il0 = 1.0 / l0;
      real il1 = 1.0 / l1;
      real il2 = 1.0 / l2;

      coordinate_type e0 = il0 * (p2 - p1);
      coordinate_type e1 = il1 * (p0 - p2);
      coordinate_type e2 = il2 * (p1 - p0);

      real ih0 = 0.5 * l0 / A;
      real ih1 = 0.5 * l1 / A;
      real ih2 = 0.5 * l2 / A;

      real cos0 = -va::dot(e2, e1);
      real cos1 = -va::dot(e0, e2);
      real cos2 = -va::dot(e1, e0);

      auto grad = [](vec12 &dtheta, const int &orient, const vec3 &N,      //
                     const real &cos0, const real &cos1, const real &cos2, //
                     const real &ih0, const real &ih1, const real &ih2     //
                  ) {
        vec3 dt0 = -ih0 * N;
        vec3 dt1 = cos2 * ih1 * N;
        vec3 dt2 = cos1 * ih2 * N;
        if (orient == 0) {
          dtheta.block(0, 0, 3, 1) += dt0;
          dtheta.block(3, 0, 3, 1) += dt1;
          dtheta.block(6, 0, 3, 1) += dt2;
        } else {
          dtheta.block(9, 0, 3, 1) += dt0;
          dtheta.block(6, 0, 3, 1) += dt1;
          dtheta.block(3, 0, 3, 1) += dt2;
        }
      };

      grad(_grad[ie0], o0, N, //
           cos0, cos1, cos2,  //
           ih0, ih1, ih2);
      grad(_grad[ie1], o1, N, //
           cos1, cos2, cos0,  //
           ih1, ih2, ih0);
      grad(_grad[ie2], o2, N, //
           cos2, cos0, cos1,  //
           ih2, ih0, ih1);

      coordinate_type m0 = va::cross(e0, N);
      coordinate_type m1 = va::cross(e1, N);
      coordinate_type m2 = va::cross(e2, N);

      mat33 M0 = N * m0.transpose();
      mat33 M1 = N * m1.transpose();
      mat33 M2 = N * m2.transpose();

      int iv0 = fv0->vertex()->position_in_set();
      int iv1 = fv1->vertex()->position_in_set();
      int iv2 = fv2->vertex()->position_in_set();

      auto Htri = [iv0, iv1,
                   iv2](mat99 &H,                                             //
                        const int &o0, const int &o1, const int &o2,          //
                        const real &psi0, const real &psi1, const real &psi2, //
                        const real &il0, const real &il1, const real &il2,    //
                        const real &ih0, const real &ih1, const real &ih2,    //
                        const real &cos0, const real &cos1, const real &cos2, //
                        const mat33 &M0, const mat33 &M1, const mat33 &M2     //
                  ) {
        auto dag = [](const mat33 &M, int orient) {
          mat33 Mt = M.transpose();
          return orient ? Mt : M;
          // return orient ? M.transpose() : M;
        };

        real d0 = psi2 * cos1 + psi1 * cos2 - psi0;
        real d1 = psi0 * cos2 + psi2 * cos0 - psi1;
        real d2 = psi1 * cos0 + psi0 * cos1 - psi2;

        mat33 M0t = M0.transpose();
        mat33 M1t = M1.transpose();
        mat33 M2t = M2.transpose();
        mat33 R0 = psi0 * il0 * il0 * M0;
        mat33 R1 = psi1 * il1 * il1 * M1;
        mat33 R2 = psi2 * il2 * il2 * M2;

        mat33 H00 = ih0 * ih0 * (d0 * M0t + d0 * M0) - R2 - R1;
        mat33 H11 = ih1 * ih1 * (d1 * M1t + d1 * M1) - R0 - R2;
        mat33 H22 = ih2 * ih2 * (d2 * M2t + d2 * M2) - R1 - R0;
        mat33 H01 = ih0 * ih1 * (d0 * M1t + d1 * M0) + dag(R2, o2);
        mat33 H12 = ih1 * ih2 * (d1 * M2t + d2 * M1) + dag(R0, o0);
        mat33 H20 = ih2 * ih0 * (d2 * M0t + d0 * M2) + dag(R1, o1);

        int i0 = 0;
        int i1 = 3;
        int i2 = 6;

        H.block(i0, i0, 3, 3) = H00;
        H.block(i1, i1, 3, 3) = H11;
        H.block(i2, i2, 3, 3) = H22;

        H.block(i0, i1, 3, 3) = H01;
        H.block(i1, i0, 3, 3) = H01.transpose();

        H.block(i1, i2, 3, 3) = H12;
        H.block(i2, i1, 3, 3) = H12.transpose();

        H.block(i2, i0, 3, 3) = H20;
        H.block(i0, i2, 3, 3) = H20.transpose();
      };

      Htri(_hess0[i],           //
           o0, o1, o2,          //
           psip0, psip1, psip2, //
           il0, il1, il2,       //
           ih0, ih1, ih2,       //
           cos0, cos1, cos2,    //
           M0, M1, M2);
      //_hess0[i] /= _Abar;

      i++;
    }

    finish_grad_hess();
  }

  virtual void fill_gradient(vecX &G) {
    edge_array edges = _surf->get_edges();
    for (auto &e : edges) {
      vertex_ptr v1 = e->v1()->vertex();
      vertex_ptr v2 = e->v2()->vertex();
      vertex_ptr v0 = e->v1()->prev()->vertex();
      vertex_ptr v3 = e->v2()->prev()->vertex();

      size_t i0 = v0->position_in_set();
      size_t i1 = v1->position_in_set();
      size_t i2 = v2->position_in_set();
      size_t i3 = v3->position_in_set();

      size_t ii[] = {3 * i0 + 0, 3 * i0 + 1, 3 * i0 + 2, //
                     3 * i1 + 0, 3 * i1 + 1, 3 * i1 + 2, //
                     3 * i2 + 0, 3 * i2 + 1, 3 * i2 + 2, //
                     3 * i3 + 0, 3 * i3 + 1, 3 * i3 + 2};
      vec12 g = _grad[e->position_in_set()];
      for (int i = 0; i < 12; i++) {
        G(ii[i]) += g(i);
      }
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
    face_array faces = _surf->get_faces();
#if 0
    for (auto &f : faces) {

      face_vertex_ptr fv0 = f->fbegin();
      face_vertex_ptr fv1 = fv0->next();
      face_vertex_ptr fv2 = fv1->next();

      vertex_ptr v0 = fv0->vertex();
      vertex_ptr v1 = fv1->vertex();
      vertex_ptr v2 = fv2->vertex();

      size_t i0 = v0->position_in_set();
      size_t i1 = v1->position_in_set();
      size_t i2 = v2->position_in_set();

      size_t ii[] = {3 * i0 + 0, 3 * i0 + 1, 3 * i0 + 2, //
                     3 * i1 + 0, 3 * i1 + 1, 3 * i1 + 2, //
                     3 * i2 + 0, 3 * i2 + 1, 3 * i2 + 2};
      const mat99 &Hi = _hess0[f->position_in_set()];

      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          triplets.push_back(triplet(ii[i], ii[j], Hi(i, j)));
        }
      }
    }
#endif
#if 0
    edge_array edges = _surf->get_edges();
    for (auto &e : edges) {
      vertex_ptr v1 = e->v1()->vertex();
      vertex_ptr v2 = e->v2()->vertex();
      vertex_ptr v0 = e->v1()->prev()->vertex();
      vertex_ptr v3 = e->v2()->prev()->vertex();

      size_t i0 = v0->position_in_set();
      size_t i1 = v1->position_in_set();
      size_t i2 = v2->position_in_set();
      size_t i3 = v3->position_in_set();

      size_t ii[] = {3 * i0 + 0, 3 * i0 + 1, 3 * i0 + 2, //
                     3 * i1 + 0, 3 * i1 + 1, 3 * i1 + 2, //
                     3 * i2 + 0, 3 * i2 + 1, 3 * i2 + 2, //
                     3 * i3 + 0, 3 * i3 + 1, 3 * i3 + 2};

      const mat1212 &Hi = _hess1[e->position_in_set()];
      for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
          triplets.push_back(triplet(ii[i], ii[j], Hi(i, j)));
        }
      }
    }
#endif
  }

  virtual void update(const vecX &vals) { _positions = vals; }
  real _Abar = 0.0;
  real _k = 0.0;
  surf_ptr _surf;
  vecX _positions;

  std::vector<real> _phi0;
  std::vector<real> _ai;

  std::vector<vec12> _grad;
  std::vector<mat99> _hess0;
  std::vector<mat1212> _hess1;

  std::vector<real> _e_len; // inverse_edge_lengths
  std::vector<real> _phi1;  // phi1

  std::vector<real> _psi;    // psi derivative
  std::vector<real> _psi_p;  // psi derivative
  std::vector<real> _psi_pp; // psi double derivative
};

#if 1
template <typename SPACE> class internal_collisions : public constraint<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<internal_collisions<SPACE>> ptr;

  static ptr create(surf_ptr surf, real reg) {
    return std::make_shared<internal_collisions<SPACE>>(surf, reg);
  }

  internal_collisions(surf_ptr surf, real reg) : _surf(surf), _reg(reg) {}

  virtual real evaluate_constraint() { return 0.0; }

  virtual void preprocess(const vecX &vel) {
    coordinate_array p = ci::get_coordinates<SPACE>(_surf);
    coordinate_array N = ci::get_vertex_normals<SPACE>(_surf);

    vector<vec3> forces = calc_f(_surf, 0.25 * _reg);
    std::cout << " checking internal collisions" << std::endl;
    //    gg::geometry_logger::field(p, forces, 1.0, gg::PresetColor::rainbow);

    for (int i = 0; i < forces.size(); i++) {
      vec3 &f = forces[i];
      vec3 v = vel.segment(3 * i, 3);
      real fnorm = f.norm();
      real mx = 0.5 * fnorm;
      real denom = f.dot(f);
      denom = denom < 1e-16 ? 1 : denom;
      // f = 0.5 * v.dot(f) / denom * f;
      f = 0.75 * va::sgn(v.dot(f)) * v.dot(f) / denom * f;
      real fn = f.norm();

      f = fn > mx ? f * mx / fn : f;
      f = fn == 0.0 ? 0.0 * f : f;
    }

    _f = vec_interface<SPACE>::to(forces);
  }

  vector<typename SPACE::vec3> calc_f(asawa::surf<SPACE> *mesh,
                                      typename SPACE::real regLength = 0.5) {
    M2_TYPEDEFS;

    typedef Eigen::Matrix<real, 3, 1> vec3;
    typedef Eigen::Matrix<real, 4, 1> vec4;

    using Avg_Integrator =
        calder::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    std::vector<vec3> face_normals = ci::get_normals<SPACE>(mesh);
    coordinate_array vertex_normals = ci::get_vertex_normals<SPACE>(mesh);
    coordinate_array evalPoints = ci::get_coordinates<SPACE>(mesh);

    if (evalPoints.empty())
      return std::vector<vec3>();

    auto pre = [face_normals](const vector<triangle_type> &tris, ANode &node,
                              ATree &tree, vec3 &netCharge,
                              coordinate_type &avgPoint,
                              coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      netCharge = z::zero<vec3>();

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        netCharge += w * face_normals[ii];
        netWeight += w;
      }

      avgPoint /= netWeight;
    };

    auto compute = [&vertex_normals, regLength](
                       int i_c, const vec3 &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> vec3 {
      vec3 out = z::zero<vec3>();

      coordinate_type Nv = vertex_normals[i_c];
      coordinate_type dp = pc - pe;
      T dist = va::norm(dp);
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          auto tri = tris[ii];
          auto Nf = tri.normal();
          bool itx = tri.rayIntersect(pe, pe + Nv, dist);

          if (itx) {
            // std::cout << itx << " " << dist << std::endl;
            gg::geometry_logger::line(pe, pe + dist * Nv,
                                      vec4(1.0, 0.0, 1.0, 1.0));
            out += dist * Nv;
          }
        }
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }
    vector<vec3> u(evalPoints.size(), z::zero<vec3>());
    Avg_Integrator integrator;
    integrator.integrate(face_normals, triangles, evalPoints, u, pre, compute);
    int i = 0;

    return u;
  }

  virtual void fill_gradient(vecX &G) { G += _f; }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {}
  vecX _f;
  surf_ptr _surf;
  real _reg;
};

#else
template <typename SPACE> class internal_collisions : public constraint<SPACE> {

public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 9, 1> vec9;
  typedef Eigen::Matrix<real, 9, 9> mat99;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<internal_collisions<SPACE>> ptr;

  static ptr create(surf_ptr surf, real reg) {
    return std::make_shared<internal_collisions<SPACE>>(surf, reg);
  }

  internal_collisions(surf_ptr surf, real reg) : _surf(surf), _reg(reg) {}

  virtual real evaluate_constraint() { return 0.0; }

  virtual void preprocess(const vecX &vel) {
    coordinate_array p = ci::get_coordinates<SPACE>(_surf);
    coordinate_array N = ci::get_vertex_normals<SPACE>(_surf);

    std::cout << "vevl" << std::endl;
    vector<vec3> velv(p);
    vec_interface<SPACE>::from(velv, vel);
    std::cout << "calc: " << velv.size() << std::endl;

    vector<vec3> forces = this->calc_f(_surf, N, 1.0 * _reg);
    std::cout << " checking internal collisions" << std::endl;
    std::cout << " forces: " << forces.size() << std::endl;
    std::cout << " forces: " << p.size() << std::endl;
    real C = _reg;
    // gg::geometry_logger::field(p, forces, 2.0, gg::PresetColor::red);
    // gg::geometry_logger::field(p, N, 1.0, gg::PresetColor::blue);

    for (int i = 0; i < forces.size(); i++) {
      vec3 &f = forces[i];
      f -= N[i];
      // f = vec3::Zero();
      vec3 v = vel.segment(3 * i, 3);
      real fnorm = f.norm();
      real mx = 0.5 * fnorm;
      real denom = f.dot(f);
      denom = denom < 1e-16 ? 1 : denom;
      // f = 0.5 * v.dot(f) / denom * f;
      f = -1.0 * va::sgn(v.dot(f)) * v.dot(f) / denom * f;
      real fn = f.norm();
      f = fn > mx ? f * mx / fn : f;
    }
    gg::geometry_logger::field(p, forces, 1.0, gg::PresetColor::blue);

    _f = vec_interface<SPACE>::to(forces);
  }

  vector<typename SPACE::vec3>
  calc_f(asawa::surf<SPACE> *mesh,
         const std::vector<typename SPACE::vec3> &vertVals,
         typename SPACE::real regLength = 0.5) {
    M2_TYPEDEFS;

    typedef Eigen::Matrix<real, 3, 1> vec3;
    typedef Eigen::Matrix<real, 4, 1> vec4;

    using Avg_Integrator =
        asawa::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

    using ATree = typename Avg_Integrator::Tree;
    using ANode = typename Avg_Integrator::Node;

    std::vector<vec3> face_normals = ci::get_normals<SPACE>(mesh);
    coordinate_array vertex_normals = ci::get_vertex_normals<SPACE>(mesh);
    coordinate_array evalPoints = ci::get_coordinates<SPACE>(mesh);

    std::vector<vec3> faceVals = ci::verts_to_faces<SPACE>(vertVals, mesh);

    if (evalPoints.empty())
      return std::vector<vec3>();

    auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                          ATree &tree, vec3 &netCharge,
                          coordinate_type &avgPoint,
                          coordinate_type &avgNormal) -> void {
      avgPoint = coordinate_type(0, 0, 0);
      avgNormal = coordinate_type(0, 0, 0);
      netCharge = coordinate_type(0, 0, 0);

      T netWeight = 0;

      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        triangle_type tri = tris[ii];
        T w = tri.area();
        avgPoint += w * tri.center();
        avgNormal += w * tri.normal();
        netCharge += w * faceVals[ii];
        netWeight += w;
      }

      avgPoint /= netWeight;
    };

    auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 0
      T dist3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = 1.0 / (dist3);
      return kappa / 4.0 / M_PI;
#elif 1
      // T dx = dist - C;
      T dx = dist;

      T C2 = C * C;
      T X = dx * dx / 2.0 / C2;
      T kappa = exp(-X) / pow(2.0 * M_PI * C2, 0.5);
      return kappa;
#endif
    };

    std::vector<real> K(evalPoints.size());
    auto compute = [&vertex_normals, &faceVals, regLength, &K, &computeK](
                       int i_c, const vec3 &wq, const coordinate_type &pc,
                       const coordinate_type &pe, const coordinate_type &N,
                       const vector<triangle_type> &tris, ANode &node,
                       ATree &tree) -> vec3 {
      vec3 out = z::zero<vec3>();

      coordinate_type Nv = vertex_normals[i_c];
      coordinate_type dp = pc - pe;
      T dist = va::norm(dp);
      if (node.isLeaf()) {
        for (int i = node.begin; i < node.begin + node.size; i++) {
          int ii = tree.permutation[i];
          vec3 qi = faceVals[ii];
          auto tri = tris[ii];
          auto w = tri.area();
          auto c = tri.center();
          coordinate_type dp = c - pe;
          T dist = va::norm(dp);

          T k = computeK(dist, regLength);
          out += w * k * qi;
          K[i_c] += w * k;
        }
      } else {
        coordinate_type dp = pc - pe;
        T dist = dp.norm();
        T k = computeK(dist, regLength);
        out += k * wq;
        K[i_c] += k * N.norm();
      }
      return out;
    };

    vector<face_ptr> &faces = mesh->get_faces();
    vector<coordinate_type> normals;
    vector<triangle_type> triangles;
    for (int i = 0; i < faces.size(); i++) {
      if (!mesh->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    for (auto t : triangles) {
      normals.push_back(t.normal());
    }

    vector<vec3> u(evalPoints.size(), z::zero<vec3>());
    Avg_Integrator integrator;
    integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);
    int i = 0;
    /*
        for (auto p : evalPoints) {
          gg::geometry_logger::line(p, p + u[i++], vec4(1.0, 0.0, 1.0, 1.0));
        }
    */
    i = 0;

    for (auto &ui : u) {
      ui /= K[i++];
      // i++;
    }

    return u;
  }

  virtual void fill_gradient(vecX &G) { G += _f; }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {}
  vecX _f;
  surf_ptr _surf;
  real _reg;
};
#endif

} // namespace hepworth
#endif