
#ifndef __TWOMANIFOLD_POSITION_OPTIMIZATION__
#define __TWOMANIFOLD_POSITION_OPTIMIZATION__
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "Eigen/src/SparseCore/SparseUtil.h"
#include "Eigen/src/SparseLU/SparseLU.h"
#include "Eigen/src/SparseQR/SparseQR.h"
#include "coordinate_interface.hpp"

#include <GaudiGraphics/geometry_logger.h>

#include "manifold/harmonic_integrators.hpp"
#include "manifold/m2.hpp"

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace m2 {

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

  static void to(coordinate_type c, vecX &vals, size_t i) {
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

  static vec3 normal(face_ptr f, const vecX &vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();
    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return m2::va::calculate_normal(v0, v1, v2);
  }

  static real area(face_ptr f, const vecX &vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();

    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return m2::va::calculate_area(v0, v1, v2);
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

  virtual void preprocess() {}
  virtual void fill_gradient(vecX &G) {}
  virtual void get_grad_triplets(std::vector<triplet> &triplets) {}
  virtual void get_hess_triplets(std::vector<triplet> &triplets) {}
  virtual void update(const vecX &vals) {}
  virtual real evaluate_constraint() { return 0.0; }

  mat32 calc_f(const coordinate_type &c0, //
               const coordinate_type &c1, //
               const coordinate_type &c2, //
               const coordinate_type &U,  //
               const coordinate_type &V) {

    coordinate_type dx0 = c1 - c0;
    coordinate_type dx1 = c2 - c0;
    vec2 uv0(0, 0);
    vec2 uv1(va::dot(dx0, U), va::dot(dx0, V));
    vec2 uv2(va::dot(dx1, U), va::dot(dx1, V));
    mat2 duv;
    duv << //
        uv1[0] - uv0[0],
        uv2[0] - uv0[0], //
        uv1[1] - uv0[1], //
        uv2[1] - uv0[1];
    mat32 dx;
    dx.block(0, 0, 1, 3) << dx0.transpose();
    dx.block(0, 1, 1, 3) << dx1.transpose();
    mat32 F = dx * duv.inverse();
    return F;
  }

  mat32 calc_f(const coordinate_type &c0, //
               const coordinate_type &c1, //
               const coordinate_type &c2) {
    coordinate_type N = va::calculate_normal<SPACE>(c0, c1, c2);
    coordinate_type U = va::cross(N, coordinate_type(0.0, 0.0, 1.0));
    coordinate_type V = va::cross(N, V);
    return calc_uv(c0, c1, c2, U, V);
  }
};

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

  edge_length(size_t i0, size_t i1) : _i0(i0), _i1(i1) {}

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
  }

  mat36 jacobian() {
    mat36 J;
    J.setZero();
    J.block(0, 0, 3, 3) = -mat3::Identity(3, 3);
    J.block(0, 3, 3, 3) = mat3::Identity(3, 3);
    return J;
  }

  virtual real evaluate_constraint() {
    coordinate_type dp = _p1 - _p0;
    return _mu * dp.transpose() * dp;
  }

  vec6 local_gradient() {
    coordinate_type dp = 2.0 * (_p1 - _p0);
    mat36 J = jacobian();
    vec6 G = _mu * J.transpose() * dp;
    // out.setZero();
    // out.block(0, 0, 3, 1) = dCdxi;
    // out.block(3, 0, 3, 1) = -dCdxi;

    return G;
  }

  mat66 local_hessien() {
    mat3 I = 2.0 * mat3::Identity(3, 3);
    mat36 J = jacobian();
    mat66 H = _mu * J.transpose() * I * J;
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
  real _mu = 0.1;
}; // namespace m2

template <typename SPACE> class edge_stretch : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<edge_stretch<SPACE>> ptr;

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

  edge_stretch(size_t i0, size_t i1, real l) : _i0(i0), _i1(i1), _lc(l) {}

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
  }

  mat36 jacobian() {
    mat36 J;
    J.setZero();
    J.block(0, 0, 3, 3) = -mat3::Identity(3, 3);
    J.block(0, 3, 3, 3) = mat3::Identity(3, 3);
    return J;
  }

  virtual real evaluate_constraint() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    real de = l - _lc;
    return _mu * de * de;
  }

  vec6 local_gradient() {
    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    coordinate_type dldx = dp * invl;
    real de = l - _lc;
    coordinate_type dCdxi = 2.0 * _mu * de * dldx;
    // coordinate_type dCdxi = 2.0 * _mu * (dp - _lc * dldx);

    mat36 J = jacobian();
    vec6 G = J.transpose() * dCdxi;
    // out.setZero();
    // out.block(0, 0, 3, 1) = dCdxi;
    // out.block(3, 0, 3, 1) = -dCdxi;

    return G;
  }

  mat66 local_hessien() {

    coordinate_type dp = _p1 - _p0;
    real l = va::norm(dp);
    real invl = 1.0 / l;
    // coordinate_type dldx = dp * invl;
    // real de = l - _lc;
    mat3 I = mat3::Identity(3, 3);
    mat3 dldlT = dp * dp.transpose() * invl * invl;
    mat3 d2ldx2 = I - dldlT;
    // mat3 d2CDX2 = 2.0 * (de * d2ldx2 + 2.0 * dldlT);
    mat3 d2CDX2 = 2.0 * _mu * (-I + _lc * invl * d2ldx2);
    // std::cout << "norm: " << l << " " << _lc << " " << d2CDX2.norm()
    //           << std::endl;

    mat36 J = jacobian();
    // vec6 G = local_gradient();
    // mat66 Hd = G * G.transpose();

    // return Hd;
    mat66 H = J.transpose() * d2CDX2 * J;
    /*
        std::cout << d2CDX2 << std::endl;
        std::cout << std::endl;
        std::cout << d2CDX2p << std::endl;
        std::cout << " --------------- " << std::endl;
        std::cout << std::endl;
    */
    if (0) {
      std::cout << "l0/l: " << _lc * invl << std::endl;
      std::cout << dldlT << std::endl;
      std::cout << std::endl;

      std::cout << d2CDX2 << std::endl;
      std::cout << std::endl;

      std::cout << H << std::endl;
      std::cout << std::endl;
    }
    // std::cout << H << std::endl;
    // std::cout << H.inverse() << std::endl;
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
  real _lc;
  real _mu = 1.0;
}; // namespace m2

template <typename SPACE> class cross_ratio : public constraint<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<cross_ratio<SPACE>> ptr;

  static ptr create(size_t ii, size_t ij, size_t ik, size_t il, real k) {
    return std::make_shared<cross_ratio<SPACE>>(ii, ij, ik, il, k);
  }

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 3> mat33;
  typedef Eigen::Matrix<real, 12, 1> vec12;
  typedef Eigen::Matrix<real, 12, 12> mat1212;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  typedef Eigen::SparseMatrix<real> sparmat;

  cross_ratio(size_t ii, size_t ij, size_t ik, size_t il, real k)
      : _ii(ii), _ij(ij), _ik(ik), _il(il), _k(k) {}

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
    return _mu * C * C;
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
    return 2.0 * _mu * C * X * G;
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

    H = 2.0 * _mu * (C * X * H + (C * X + X * X) * D2 * D2.transpose());
    return H;
  }

  virtual void fill_gradient(vecX &G) {
    vec12 g = local_gradient();
    // std::cout << " lGNorm: " << g.norm() << std::endl;

    size_t ii[] = {3 * _ii + 0, 3 * _ii + 1, 3 * _ii + 2, //
                   3 * _ij + 0, 3 * _ij + 1, 3 * _ij + 2, //
                   3 * _ik + 0, 3 * _ik + 1, 3 * _ik + 2, //
                   3 * _il + 0, 3 * _il + 1, 3 * _il + 2};
    for (int i = 0; i < 12; i++) {
      G(ii[i]) += g(i);
    }
  }

  virtual void get_hess_triplets(std::vector<triplet> &triplets) {
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
  real _mu = 1.0;
}; // namespace m2

#if 0
template <typename SPACE> class shear : public constraint<SPACE> {

  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 6> mat36;
  typedef Eigen::Matrix<real, 6, 6> mat66;

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
    _p2 = vec_interface<SPACE>::from(vals, _i2);
  }

  mat36 jacobian() {
    mat36 J;
    J.setZero();
    J.block(0, 0, 3, 3) = mat3::Identity(3, 3);
    J.block(0, 3, 3, 3) = -mat3::Identity(3, 3);
    return J;
  }

  virtual real evaluate_constraint() {
    const mat32 F = calc_f(_p0, _p1, _p1, _u, _v);
    return (_a * F).transpose() * F * _b;
  }

  mat32 local_gradient() {
    const mat32 F = calc_f(_p0, _p1, _p1, _u, _v);
    return F * (_a * _b.transpose() + _b * _a.transpose());
  };

  vec6 flatten(const mat32 &A) const {
    vec6 column;
    unsigned int index = 0;
    for (unsigned int j = 0; j < 2; j++)
      for (unsigned int i = 0; i < 3; i++, index++)
        column[index] = A(i, j);
    return column;
  }

  void local_hessien(const mat32 &F, mat66 &pPpF) const {
    const vec2 u(1.0, 0.0);
    const vec2 v(0.0, 1.0);
    const real I6 = (F * u).transpose() * (F * v);
    const real signI6 = (I6 >= 0) ? 1.0 : -1.0;
    mat66 H = mat66::Zero();
    H(3, 0) = H(4, 1) = H(5, 2) = H(0, 3) = H(1, 4) = H(2, 5) = 1.0;
    const vec6 g = flatten(F * (u * v.transpose() + v * u.transpose()));
    // get the novel eigenvalue
    const real I2 = F.squaredNorm();
    const real lambda0 = 0.5 * (I2 + sqrt(I2 * I2 + 12.0 * I6 * I6));
    // get the novel eigenvector
    // the H multiply is a column swap; could be optimized more
    const vec6 q0 = (I6 * H * g + lambda0 * g).normalized();
    mat66 T = mat66::Identity();
    T = 0.5 * (T + signI6 * H);
    const vec6 Tq = T * q0;
    const real normTq = Tq.squaredNorm();
    pPpF = fabs(I6) * (T - (Tq * Tq.transpose()) / normTq) +
           lambda0 * (q0 * q0.transpose());
    // half from mu and leading 2 on Hessian cancel
    pPpF *= _mu;
  }

  real _mu = 1.0;
  size_t _i0 = -1, _i1 = -1, _i2 = -1;
  coordinate_type _p0, _p1, _p2, _u, _v, _a, _b;
  real _la = 1.0, _lb = 1.0;
};

template <typename SPACE> class stretch : public constraint<SPACE> {

  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  virtual void update(const vecX &vals) {
    _p0 = vec_interface<SPACE>::from(vals, _i0);
    _p1 = vec_interface<SPACE>::from(vals, _i1);
    _p2 = vec_interface<SPACE>::from(vals, _i2);
  }

  virtual real evaluate_constraint() {
    const mat32 F = calc_f(_p0, _p1, _p1, _u, _v);
    real aFtFa = (_a * F).transpose() * F * _a;
    real bFtFb = (_b * F).transpose() * F * _b;
    real nFa = sqrt(aFtFa);
    real nFb = sqrt(bFtFb);
    real dla = nFa - _la;
    real dlb = nFb - _lb;

    return dla * dla + dlb + dlb;
  }

  mat32 local_gradient() {
    const mat32 F = calc_f(_p0, _p1, _p1, _u, _v);
    real aFtFa = (_a * F).transpose() * F * _a;
    real bFtFb = (_b * F).transpose() * F * _b;
    real lFa = sqrt(aFtFa);
    real lFb = sqrt(bFtFb);
    real Faa = F * _a * _a.transpose();
    real Fbb = F * _b * _b.transpose();

    mat32 dFa = Faa / lFa;
    mat32 dFb = Fbb / lFb;
    mat32 dCdx = 2.0 * (lFa - _la) * dFa + (lFb - _lb) * dFb;
    return dCdx;
  };

  void local_hessian(const mat32 &F, mat66 &H) const {
    H.setZero();
    const vec2 u(1.0, 0.0);
    const vec2 v(0.0, 1.0);
    const real I5u = (F * u).transpose() * (F * u);
    const real I5v = (F * v).transpose() * (F * v);
    const real invSqrtI5u = 1.0 / sqrt(I5u);
    const real invSqrtI5v = 1.0 / sqrt(I5v);
    // set the block diagonals, build the rank-three
    // subspace with all-(1 / invSqrtI5) eigenvalues
    H(0, 0) = H(1, 1) = H(2, 2) = std::max((1.0 - invSqrtI5u), 0.0);
    H(3, 3) = H(4, 4) = H(5, 5) = std::max((1.0 - invSqrtI5v), 0.0);
    // modify the upper block diagonal, bump the single
    // outer-product eigenvalue back to just 1, unless it
    // was clamped, then just set it directly to 1
    const vec3 fu = F.col(0).normalized();
    const real uCoeff = (1.0 - invSqrtI5u >= 0.0) ? invSqrtI5u : 1.0;
    H.block<3, 3>(0, 0) += uCoeff * (fu * fu.transpose());
    // modify the lower block diagonal similarly
    const vec3 fv = F.col(1).normalized();
    const real vCoeff = (1.0 - invSqrtI5v >= 0.0) ? invSqrtI5v : 1.0;
    H.block<3, 3>(3, 3) += vCoeff * (fv * fv.transpose());
    // the leading 2 is absorbed by the mu / 2 coefficient
    H *= _mu;
  }

  real _mu = 1.0;
  size_t _i0 = -1, _i1 = -1, _i2 = -1;
  coordinate_type _p0, _p1, _p2, _u, _v, _a, _b;
  real _la = 1.0, _lb = 1.0;
};
#endif

template <typename SPACE> class bend : public constraint<SPACE> {

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

  typedef std::shared_ptr<bend<SPACE>> ptr;

  static ptr create(surf_ptr surf) {
    return std::make_shared<bend<SPACE>>(surf);
  }

  bend(surf_ptr surf) : _surf(surf) {
    coordinate_array positions = m2::ci::get_coordinates<SPACE>(_surf);
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
    return vec_interface<SPACE>::area(f, vals);
  }

  std::function<real(const coordinate_type &, //
                     const coordinate_type &, //
                     const coordinate_type)>

      calc_phi = [](const coordinate_type &N0, //
                    const coordinate_type &N1, //
                    const coordinate_type &e) {
        real sint = 0.5 * va::norm(vec3(N1 - N0));
        real cost = 0.5 * va::norm(vec3(N1 + N0));
        real tant = m2::va::sgn(va::determinant(N0, N1, e)) * sint / cost;
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
    _k = 0.001;
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
      _ai[i] = 3.0 * va::norm2(e_vec) / (A0 + A1);
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

  void preprocess_edges() {
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
      _e_len[i] = m2::va::norm(e_vec);
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
      _hess1[i] = pp * dt * dt.transpose();
      //_hess1[i] = dt * dt.transpose();

      i++;
    }
  }

  virtual void preprocess() {
    face_array faces = _surf->get_faces();
    preprocess_edges();

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

      e0 *= 100.0;
      e1 *= 100.0;
      e2 *= 100.0;
      N = e0.cross(-e1);
      A = N.norm();
      N /= A;
      A /= 2;

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

      // l0 = 1.0;
      // l1 = 1.0;
      // l2 = 1.0;

      // psip0 = 1.0;
      // psip1 = 1.0;
      // psip2 = 1.0;

      // A = 1.0;

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
      /*
            coordinate_type cen = 0.333 * (p0 + p1 + p2);
            gg::geometry_logger::get_instance().line(p0, p0 + 0.02 * m0,
                                                     vec4(1.0, 0.0, 0.0, 1.0));
            gg::geometry_logger::get_instance().line(p1, p1 + 0.02 * m1,
                                                     vec4(0.0, 1.0, 0.0, 1.0));
            gg::geometry_logger::get_instance().line(p2, p2 + 0.02 * m2,
                                                     vec4(0.0, 0.0, 1.0, 1.0));
            gg::geometry_logger::get_instance().line(cen, cen + 0.02 * N,
                                                     vec4(1.0, 0.0, 0.0, 1.0));
      */

#if 1
      mat33 M0 = N * m0.transpose();
      mat33 M1 = N * m1.transpose();
      mat33 M2 = N * m2.transpose();
#else
      mat33 M0 = m0 * N.transpose();
      mat33 M1 = m1 * N.transpose();
      mat33 M2 = m2 * N.transpose();
#endif

#if 1
      if (abs(psip0 + psip1 + psip2)) {
        std::cout << " A: " << A << std::endl;
      }
#endif
#if 0
      if (abs(psip0 + psip1 + psip2)) {
        std::cout << " psi: " << psip0 << " " << psip1 << " " << psip2
                  << std::endl;
        std::cout << " ih: " << ih0 << " " << ih1 << " " << ih2 << std::endl;
        std::cout << " il: " << il0 << " " << il1 << " " << il2 << std::endl;
        std::cout << " cos: " << cos0 << " " << cos1 << " " << cos2
                  << std::endl;
        std::cout << " A: " << A << std::endl;
        std::cout << cos0 << " " << cos1 << " " << cos2 << std::endl;
        std::cout << va::abs_cos(p0, p1, p2) << " " << va::abs_cos(p1, p2, p0)
                  << " " << va::abs_cos(p2, p0, p1) << std::endl;
      }
#endif

      auto Htri = [](mat99 &H,                                             //
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

        //  mat33 H01 = ih0 * ih1 * (d0 * M1t + d1 * M0);
        //  mat33 H12 = ih1 * ih2 * (d1 * M2t + d2 * M1);
        //  mat33 H20 = ih2 * ih0 * (d2 * M0t + d0 * M2);

        // mat33 H00 = -R2 - R1;
        // mat33 H11 = -R0 - R2;
        // mat33 H22 = -R1 - R0;
        // mat33 H01 = dag(R2, o2);
        // mat33 H12 = dag(R0, o0);
        // mat33 H20 = dag(R1, o1);

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
#if 0
      if (abs(psip0 + psip1 + psip2)) {
        /*
        std::cout << _grad[ie0] << std::endl;
        std::cout << " " << std::endl;
        std::cout << _grad[ie1] << std::endl;
        std::cout << " " << std::endl;
        std::cout << _grad[ie2] << std::endl;
        */
        std::cout << " H1: " << _hess0[i] << std::endl;
      }

#endif
      /*
            Htri(_hess0[i], 1, oriente1, //
                 psip1, psip2, psip0,    //
                 il1, il2, il0,          //
                 ih1, ih2, ih0,          //
                 cos1, cos2, cos0,       //
                 M1, M2, M0);

            Htri(_hess0[i], 2, oriente2, //
                 psip2, psip0, psip1,    //
                 il2, il0, il1,          //
                 ih2, ih0, ih1,          //
                 cos2, cos0, cos1,       //
                 M2, M0, M1);
      */
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
#if 1
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
#if 1
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

template <typename SPACE> class objective_function {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  typedef std::shared_ptr<objective_function<SPACE>> ptr;

  static ptr create() { return std::make_shared<objective_function<SPACE>>(); }
  objective_function() {}

  virtual vecX get_x() const = 0;
  virtual void set_x(const vecX &x) = 0;
  virtual void preprocess() {}
  virtual real calc_objective() = 0;

  virtual void visualize(const sparmat &M, const vecX &x){};
  virtual void visualize(const vecX &v, const vecX &x, vec4 color){};

  virtual vecX project_gradient(const vecX &x) = 0;
  virtual vecX calc_gradient(const vecX &x) = 0;
  virtual sparmat calc_hessian(const vecX &x) = 0;
};

template <typename SPACE>
class rosenbrock_function : public objective_function<SPACE> {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  typedef std::shared_ptr<rosenbrock_function<SPACE>> ptr;

  static ptr create() { return std::make_shared<rosenbrock_function<SPACE>>(); }
  rosenbrock_function() { X = vecX(2); }

  virtual vecX get_x() const { return X; };
  virtual void set_x(const vecX &x) { X = x; };

  virtual real calc_objective() {
    real ax = a - X[0];
    real yxx = X[1] - X[0] * X[0];

    return ax * ax + b * yxx * yxx;
  };

  virtual vecX project_gradient(const vecX &x) { return x; };
  virtual vecX calc_gradient(const vecX &x) {
    double xV = X[0];
    double yV = X[1];
    double dfdx = 2 * (-a + xV + 2 * b * xV * (xV * xV - yV));
    double dfdy = 2 * b * (yV - xV * xV);

    assert(std::isfinite(dfdx));
    assert(std::isfinite(dfdy));

    // Store it
    vecX grad(2);
    grad[0] += dfdx; //+ grad[0];
    grad[1] += dfdy; //+ grad[1]; /// ? Adds ? Or replace ? Or appends ?
    return grad;
  };

  virtual sparmat calc_hessian(const vecX &x) {
    double xV = X[0];
    double yV = X[1];

    // Empty the Tripletd ?
    std::vector<triplet> hess_triplets;

    hess_triplets.push_back(
        triplet(0, 0, -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2));
    hess_triplets.push_back(triplet(0, 1, -4 * b * xV));
    hess_triplets.push_back(triplet(1, 0, -4 * b * xV));
    hess_triplets.push_back(triplet(1, 1, 2 * b));

    sparmat H(2, 2);

    H.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return H;
  };

  real a = 1.0, b = 100.0;
  vecX X;
};

template <typename SPACE>
class constraint_set : public objective_function<SPACE> {

public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 3> mat33;

  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<constraint_set<SPACE>> ptr;
  static ptr create(const surf_ptr surf) {
    return std::make_shared<constraint_set<SPACE>>(surf);
  }

  constraint_set(const surf_ptr surf) : _surf(surf) {
    _reg = m2::ci::geometric_mean_length<SPACE>(_surf);
  };

  template <typename MAT> mat33 get33Block(MAT M, int i, int j) {
    real c00 = M.coeff(3 * i + 0, 3 * j + 0);
    real c01 = M.coeff(3 * i + 0, 3 * j + 1);
    real c02 = M.coeff(3 * i + 0, 3 * j + 2);

    real c10 = M.coeff(3 * i + 1, 3 * j + 0);
    real c11 = M.coeff(3 * i + 1, 3 * j + 1);
    real c12 = M.coeff(3 * i + 1, 3 * j + 2);

    real c20 = M.coeff(3 * i + 2, 3 * j + 0);
    real c21 = M.coeff(3 * i + 2, 3 * j + 1);
    real c22 = M.coeff(3 * i + 2, 3 * j + 2);
    mat33 Mb;
    Mb.row(0) << c00, c01, c02;
    Mb.row(1) << c10, c11, c12;
    Mb.row(2) << c20, c21, c22;

    return Mb;
  }

  void minmax(mat33 M, real &min, real &max) {
    coordinate_type t0 = M.block(0, 0, 3, 1);
    coordinate_type t1 = M.block(0, 1, 3, 1);
    coordinate_type t2 = M.block(0, 2, 3, 1);
    real nt0 = t0.norm();
    real nt1 = t1.norm();
    real nt2 = t2.norm();
    min = std::min(min, std::min(nt0, std::min(nt1, nt2)));
    max = std::max(max, std::max(nt0, std::max(nt1, nt2)));
  }

  template <typename MAT>
  void minmax(MAT M, int i, int j, real &min, real &max) {
    mat33 Mb = get33Block(M, i, j);
    minmax(Mb, min, max);
  }

  template <typename MAT>
  void visualizeBlock(MAT M, coordinate_type c, int i, int j, real C = 1.0) {
    mat33 Mb = get33Block(M, i, j);
    gg::geometry_logger::get_instance().frame(Mb, c, C);
  }

  virtual void visualize(const sparmat &Ms, const vecX &x) {
    vertex_array verts = _surf->get_vertices();
    edge_array edges = _surf->get_edges();
    real min = std::numeric_limits<real>::max();
    real max = 0;

    for (int k = 0; k < Ms.outerSize(); ++k)
      for (typename sparmat::InnerIterator it(Ms, k); it; ++it) {
        it.value();
        it.row();   // row index
        it.col();   // col index (here it is equal to k)
        it.index(); // inner index, here it is equal to it.row()
      }

    for (auto e : edges) {
      vec3 c = vec_interface<SPACE>::center(e, x);
      int ei = e->v1()->vertex()->position_in_set();
      int ej = e->v2()->vertex()->position_in_set();
      minmax(Ms, ei, ej, min, max);
    }
    for (auto v : verts) {
      vec3 c = vec_interface<SPACE>::coord(v, x);
      int vi = v->position_in_set();
      minmax(Ms, vi, vi, min, max);
    }
    max = std::min(100.0, max);
    real scale = 1.0 / (max - min);
    std::cout << " min max scale: " << min << " " << max << " " << scale
              << std::endl;
#if 1
    for (auto e : edges) {

      vec3 c = vec_interface<SPACE>::center(e, x);
      int ei = e->v1()->vertex()->position_in_set();
      int ej = e->v2()->vertex()->position_in_set();
      visualizeBlock(Ms, c, ei, ej, 0.1 * scale);
    }
#endif
#if 1
    for (auto v : verts) {

      vec3 c = vec_interface<SPACE>::coord(v, x);
      int vi = v->position_in_set();
      visualizeBlock(Ms, c, vi, vi, 0.1 * scale);
    }
#endif
  };

  virtual void visualize(const vecX &g, const vecX &x,
                         vec4 color = vec4(0.0, 0.6, 0.7, 1.0)) {
    vertex_array verts = _surf->get_vertices();
#if 1
    for (auto v : verts) {
      int vi = v->position_in_set();
      vec3 c = vec_interface<SPACE>::coord(v, x);
      vec3 gi = vec_interface<SPACE>::coord(v, g);
      real C = 2.0;
      gg::geometry_logger::get_instance().line(c, c + C * gi, color);
    }
#endif
  };

  virtual void update_constraints(const vecX &x) {
    for (auto c : constraints) {
      c->update(x);
    }
  }

  virtual void set_positions(const coordinate_array &positions) {
    _N = positions.size();
    vecX x = vec_interface<SPACE>::to(positions);
    set_x(x);
  }

  virtual coordinate_array get_positions() const {
    coordinate_array positions(_N, coordinate_type(0, 0, 0));
    vec_interface<SPACE>::from(positions, _x);
    return positions;
  }

  virtual vecX get_x() const {
    std::cout << " x norm: " << _x.norm() << std::endl;
    return _x;
  };

  virtual void set_x(const vecX &x) {
    update_constraints(x);
    _x = x;
  };

  virtual real calc_objective() {
    real C = 0.0;
    for (auto c : constraints) {
      C += c->evaluate_constraint();
    }
    return C;
  }

#if 1
  virtual vecX project_gradient(const vecX &g) {
    coordinate_array positions = get_positions();
    coordinate_array gc(positions);
    vec_interface<SPACE>::from(gc, g);
    coordinate_array gf =
        m2::ci::verts_to_faces<SPACE, coordinate_type>(gc, _surf);
    m2::mesh_calculator<SPACE> calc;
    // this needs to be reworked to take a surface topology and a position set
    coordinate_array gn = calc.harmonicAvg(_surf, gf, positions, _reg);

    return vec_interface<SPACE>::to(gn);
  };
#endif

  virtual void preprocess() {
    for (auto c : constraints) {
      c->preprocess();
    }
  }

  virtual vecX calc_gradient(const vecX &x) {
    int N = _N;
    std::vector<triplet> grad_triplets;
    vecX G(3 * N);

    G.setZero();
    for (auto c : constraints) {
      c->fill_gradient(G);
    }
    // visualize(G, x);
    std::cout << " grad sum: " << G.sum() << std::endl;
    std::cout << " grad norm: " << G.norm() << std::endl;

    return G;
  }

  virtual sparmat calc_hessian(const vecX &x) {
    int N = _N;
    std::vector<triplet> hess_triplets;
    for (auto c : constraints) {
      c->get_hess_triplets(hess_triplets);
    }
    sparmat H(3 * N, 3 * N);

    H.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

    // visualize(H, x);
    size_t Nr = H.rows();
    size_t Nc = H.cols();
    sparmat I(Nr, Nc);

    I.setIdentity();

    // H.setIdentity();
    {
      real sum = 0.0;
      for (auto t : hess_triplets) {
        sum += t.value();
      }

      std::cout << " trip sum: " << sum << std::endl;
      std::cout << "    H sum: " << H.sum() << std::endl;
      std::cout << "    H norm: " << H.norm() << std::endl;
    }
    H = I - 1.0 * H;
    std::cout << "    I-H sum: " << H.sum() << std::endl;
    std::cout << "    I-H norm: " << H.norm() << std::endl;
    return H;
  }

  void add_constraint(typename constraint<SPACE>::ptr c) {
    this->constraints.push_back(c);
  }

  size_t _N = 0;
  std::vector<typename constraint<SPACE>::ptr> constraints;
  surf_ptr _surf;
  vecX _x;
  real _reg = 0.1;
};

template <typename SPACE> class solver {

  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  virtual void solve(typename objective_function<SPACE>::ptr fcn) = 0;
};

template <typename SPACE> class gradient_descent_solver : public solver<SPACE> {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  gradient_descent_solver(){

  };

  real calc_lambda(const vecX &g, typename objective_function<SPACE>::ptr fcn) {
    real C = fcn->calc_objective();
    real denom = g.transpose() * g;
    return C / denom;
  }

  virtual void solve(typename objective_function<SPACE>::ptr fcn) {
    vecX xk0 = fcn->get_x();
    real tol = 1;
    real alpha = 1.0;
    real beta = 0.75;
    int k = 0;
    while (tol > 1e-10 && k++ < 20) {
      fcn->preprocess();
      vecX g = fcn->calc_gradient(xk0);

      // g = fcn->project_gradient(g);
      std::cout << "gnorm: " << g.norm() << std::endl;

      if (g.norm() < 1e-10)
        return;

      real lambda = calc_lambda(g, fcn);
      std::cout << "lambda: " << lambda << std::endl;

      vecX xk1 = xk0 - lambda * g;
      tol = (xk1 - xk0).norm();
      // std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << "====>>tol: " << tol << std::endl;

      fcn->set_x(xk1);
      _g = g;
      std::swap(xk0, xk1);
    }

    fcn->set_x(xk0);
  }
  virtual vecX get_gradient() { return _g; };
  vecX _g;
};

template <typename SPACE> class newton_raphson_solver : public solver<SPACE> {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  newton_raphson_solver(){

  };

  virtual vecX solve(sparmat &A, vecX &b) {
#if 1
    // Eigen::SimplicialLLT<sparmat> solver;
    Eigen::SimplicialLDLT<sparmat> solver;
    // Eigen::ConjugateGradient<sparmat> solver;
    solver.compute(A);
    std::cout << " solver det: " << solver.determinant() << std::endl;
#else
    // Eigen::SparseQR<sparmat, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SparseLU<sparmat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
#endif
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

  real calc_lambda(const vecX &g, typename objective_function<SPACE>::ptr fcn) {
    real C = fcn->calc_objective();
    real denom = g.transpose() * g;
    return C / denom;
  }

  virtual void solve(typename objective_function<SPACE>::ptr fcn) {

    vecX xk0 = fcn->get_x();

    real tol = 1;
    real alpha = 1.0;
    real beta = 1.0;
    int k = 0;
    while (tol > 1e-10 && k < 20) {
      fcn->preprocess();

      std::vector<triplet> hess_triplets;
      std::vector<triplet> grad_triplets;
      vecX g = fcn->calc_gradient(xk0);
      sparmat H = fcn->calc_hessian(xk0);

#if 1
      vecX h = solve(H, g);
      std::cout << "h min/max : " << h.minCoeff() << " " << h.maxCoeff()
                << std::endl;
      std::cout << "h norm: " << h.norm() << std::endl;
      std::cout << "h sum: " << h.sum() << std::endl;
      std::cout << "g-h: " << (g - h).norm() << std::endl;

      // fcn->visualize(h, xk0, vec4(0.75, 0.0, 0.6, 1.0));

      vecX xk1 = xk0 - h;

      // return;
      // break;
#else
      real lambda = calc_lambda(g, fcn);
      vecX xk1 = H * xk0 - g;
      xk1 = solve(H, xk1);
      _g = g;
#endif
      // h = fcn->project_gradient(h, positions);

      // real lambda = calc_lambda(h, fcn);
      alpha *= beta;
      tol = (xk1 - xk0).norm();
      // if (k > 1)
      //   break;
      //  std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << "====>>tol: " << tol << std::endl; //<< '\r'
                                                      //      if (k > 1)
                                                      //        break;
      fcn->set_x(xk1);
      std::swap(xk0, xk1);
      k++;
    }
  }
};

template <typename SPACE> class optimizer {

  M2_TYPEDEFS;

public:
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  optimizer() {}

  void test_rosenbrock() {
    typename rosenbrock_function<SPACE>::ptr fcn =
        rosenbrock_function<SPACE>::create();
    vecX x(2);
    x[0] = 12.0;
    x[1] = 0.2125;

    fcn->set_x(x);
    // gradient_descent_solver<SPACE> solver;
    newton_raphson_solver<SPACE> solver;
    solver.solve(fcn);
    std::cout << "--- rosenbrock test ---" << std::endl;
    std::cout << fcn->get_x() << std::endl;
    std::cout << "--- --------------- ---" << std::endl;
  }

  void update(typename objective_function<SPACE>::ptr objective_function) {
    test_rosenbrock();
#if 1
    newton_raphson_solver<SPACE> solver;
#else
    gradient_descent_solver<SPACE> solver;
#endif
    solver.solve(objective_function);
  }

  std::vector<mat3> block_diag;
}; // class cloth

} // namespace m2
#endif