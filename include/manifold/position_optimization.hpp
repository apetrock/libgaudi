
#ifndef __TWOMANIFOLD_POSITION_OPTIMIZATION__
#define __TWOMANIFOLD_POSITION_OPTIMIZATION__
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "Eigen/src/SparseCore/SparseUtil.h"
#include "Eigen/src/SparseLU/SparseLU.h"
#include "Eigen/src/SparseQR/SparseQR.h"
#include "coordinate_interface.hpp"

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

  vec3 coord(face_vertex_ptr fv, vecX vals) {
    vertex_ptr v = fv->vertex();
    size_t i = v->position_in_set();
    return vec_interface<SPACE>::from(vals, i);
  }

  vec3 normal(face_ptr f, vecX vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();
    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return m2::va::calculate_normal(v0, v1, v2);
  }

  real area(face_ptr f, vecX vals) {
    face_vertex_ptr fv0 = f->fbegin();
    face_vertex_ptr fv1 = fv0->next();
    face_vertex_ptr fv2 = fv1->next();

    vec3 v0 = coord(fv0, vals);
    vec3 v1 = coord(fv1, vals);
    vec3 v2 = coord(fv2, vals);
    return m2::va::calculate_area(v0, v1, v2);
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
        phi_p = 0.5 * phi * phi + 0.5;
        phi_pp = 0.5 * (phi * phi + 1.0) * phi;

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
        real C = k * ai;
        psi = C * dph * dph;

        psi_p = 2.0 * C * dph * phi_p;

        psi_pp = 2.0 * C * (dph * phi_pp + phi_p * phi_p);

        return psi_pp;
      };

  virtual real evaluate_constraint() {
    std::vector<edge_ptr> edges = _surf->get_edges();

    int i = 0;
    real C = 0.0;
    for (auto e : edges) {
      real ai = _ai[i];
      real phi0 = _phi0[i];
      real phi1 = _phi1[i];
      real dphi = phi1 - phi0;
      C += ai * _k * dphi * dphi;
      i++;
    }
    return C;
  }

  void init() {
    std::vector<edge_ptr> edges = _surf->get_edges();

    _phi0 = std::vector<real>(edges.size(), 0);
    _ai = std::vector<real>(edges.size(), 0);
    real h = 0.01;
    real v = 0.5;
    real Y = 0.1;
    _k = Y * h * h * h / 24.0 / (1.0 - v);
    //_k = 0.0001;
    int i = 0;
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
      _ai[i] = 3.0 * va::norm2(e_vec) / (A0 + A1);
      i++;
    }

    _inv_e_len = std::vector<real>(edges.size()); // inverse_edge_lengths
    _unit_e_vec =
        std::vector<coordinate_type>(edges.size()); // unit edge vector
    _phi1 = std::vector<real>(edges.size());        // phi1
    _psi_p = std::vector<real>(edges.size());       // psi derivative
    _psi_pp = std::vector<real>(edges.size());      // psi double derivative
    _grad = std::vector<vec12>(edges.size(), vec12::Zero());
    _hess = std::vector<mat1212>(edges.size(),
                                 mat1212::Zero()); // inverse_edge_lengths
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

      real A0 = area(f0, _positions);
      real A1 = area(f1, _positions);

      coordinate_type c0 = coord(fv0, _positions);
      coordinate_type c1 = coord(fv1, _positions);

      coordinate_type N0 = normal(f0, _positions);
      coordinate_type N1 = normal(f1, _positions);

      coordinate_type e_vec = c1 - c0;
      real il = 1.0 / m2::va::norm(e_vec);

      _phi1[i] = calc_phi(N0, N1, e_vec);
      _inv_e_len[i] = il;
      _unit_e_vec[i] = e_vec * il;
      real phi_p, phi_pp, psi;
      calc_phi_derivatives(_phi1[i], phi_p, phi_pp);
      calc_psi(_phi1[i], _phi0[i], phi_p, phi_pp, _k, _ai[i], psi, _psi_p[i],
               _psi_pp[i]);
      // std::cout << _k << " " << _ai[i] << " dphi: " << phi1[i] - _phi0[i] <<
      // " "
      //           << psi << " " << psi_p[i] << " " << psi_pp[i] << std::endl;
      _grad[i] = vec12::Zero();
      _hess[i] = mat1212::Zero(); // inverse_edge_lengths

      i++;
    }
  }

  void finish_grad_hess() {
    int i = 0;
    edge_array edges = _surf->get_edges();

    for (auto e : edges) {
      const vec12 &dt = _grad[i];
      const mat1212 &Hi = _hess[i];
      const real &p = _psi_p[i];
      const real &pp = _psi_pp[i];
#if 0
      if (p > 1e-10) {
        std::cout << "--------------" << std::endl;
        std::cout << p << " " << pp << std::endl;
        std::cout << dt << std::endl;
        std::cout << Hi << std::endl;
        std::cout << _grad[i] << std::endl;
        std::cout << _hess[i] << std::endl;
    }
#endif
      _grad[i] = p * dt;
      _hess[i] = Hi + pp * dt * dt.transpose();

      // He[i] = -He[i];
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
      coordinate_type p0 = coord(fv2, _positions);
      coordinate_type p1 = coord(fv0, _positions);
      coordinate_type p2 = coord(fv1, _positions);

      int ifv0 = fv2->position_in_set();
      int ifv1 = fv0->position_in_set();
      int ifv2 = fv1->position_in_set();

      edge_ptr ep0 = fv0->edge();
      edge_ptr ep1 = fv1->edge();
      edge_ptr ep2 = fv2->edge();

      int oriente0 = ep0->side(fv0);
      int oriente1 = ep1->side(fv1);
      int oriente2 = ep2->side(fv2);

      int ie0 = fv0->edge()->position_in_set();
      int ie1 = fv1->edge()->position_in_set();
      int ie2 = fv2->edge()->position_in_set();

      real il0 = _inv_e_len[ie0];
      real il1 = _inv_e_len[ie1];
      real il2 = _inv_e_len[ie2];

      coordinate_type e0 = il0 * (p2 - p1);
      coordinate_type e1 = il1 * (p0 - p2);
      coordinate_type e2 = il2 * (p1 - p0);

      real ih0 = 0.5 / il0 / A;
      real ih1 = 0.5 / il1 / A;
      real ih2 = 0.5 / il2 / A;

      // checkpoint

      real cos0 = -va::dot(e2, e1);
      real cos1 = -va::dot(e0, e2);
      real cos2 = -va::dot(e1, e0);

      auto calc_dtheta = [](vec12 &dtheta, const int &orient, //
                            const real cos0, const real ih0,  //
                            const real cos1, const real ih1,  //
                            const real cos2, const real ih2,  //
                            const vec3 &N) {
        vec3 dt0 = -1.0 * ih0 * N;
        vec3 dt1 = cos2 * ih1 * N;
        vec3 dt2 = cos1 * ih2 * N;
        if (orient == 0) {
          dtheta.block(0, 0, 3, 1) += dt0;
          dtheta.block(3, 0, 3, 1) += dt1;
          dtheta.block(6, 0, 3, 1) += dt2;
        } else {
          dtheta.block(3, 0, 3, 1) += dt2;
          dtheta.block(6, 0, 3, 1) += dt1;
          dtheta.block(9, 0, 3, 1) += dt0;
        }
      };

      calc_dtheta(_grad[ie0], oriente0, //
                  cos0, ih0, cos1, ih1, cos2, ih2, N);
      calc_dtheta(_grad[ie1], oriente1, //
                  cos1, ih1, cos2, ih2, cos0, ih0, N);
      calc_dtheta(_grad[ie2], oriente2, //
                  cos2, ih2, cos0, ih0, cos1, ih1, N);

      coordinate_type m0 = va::cross(e0, N);
      coordinate_type m1 = va::cross(e1, N);
      coordinate_type m2 = va::cross(e2, N);

      mat3 M0 = N * m0.transpose();
      mat3 M1 = N * m1.transpose();
      mat3 M2 = N * m2.transpose();

      mat3 N0 = M0 * il0 * il0;
      mat3 N1 = M1 * il1 * il1;
      mat3 N2 = M2 * il2 * il2;

      real psi_p0 = _psi_p[ie0];
      real psi_p1 = _psi_p[ie1];
      real psi_p2 = _psi_p[ie2];

      real sig0 = 1; // no boundaries
      real sig1 = 1; // no boundaries
      real sig2 = 1; // no boundaries

      real c0 = sig0 * psi_p0;
      real c1 = sig1 * psi_p1;
      real c2 = sig2 * psi_p2;

      real d0 = c2 * cos1 + c1 * cos2 - c0;
      real d1 = c0 * cos2 + c2 * cos0 - c1;
      real d2 = c1 * cos0 + c0 * cos1 - c2;

      mat3 R0 = c0 * N0;
      mat3 R1 = c1 * N1;
      mat3 R2 = c2 * N2;

      mat3 R0d = oriente0 == 0 ? R0 : R0.transpose();
      mat3 R1d = oriente1 == 0 ? R1 : R1.transpose();
      mat3 R2d = oriente2 == 0 ? R2 : R2.transpose();

      auto calc_Htri = [](mat1212 &H, const int &orient,                     //
                          const real &d0, const real &d1, const real &d2,    //
                          const real &ih0, const real &ih1, const real &ih2, //
                          const mat33 &M0, const mat33 &M1, const mat33 &M2, //
                          const mat33 &R0, const mat33 &R1, const mat33 &R2, //
                          const mat33 &R0d, const mat33 &R1d,
                          const mat33 &R2d //
                       ) {
        mat33 H00 = ih0 * ih0 * (d0 * M0.transpose() + d0 * M0) - R1 - R2;
        mat33 H01 = ih0 * ih1 * (d0 * M1.transpose() + d1 * M0) + R2d;
        mat33 H11 = ih1 * ih1 * (d1 * M1.transpose() + d1 * M1) - R2 - R0;
        mat33 H12 = ih1 * ih2 * (d1 * M2.transpose() + d2 * M1) + R0d;
        mat33 H22 = ih2 * ih2 * (d2 * M2.transpose() + d2 * M2) - R1 - R0;

        if (orient == 0) {
          H.block(0, 0, 3, 3) += H00;
          H.block(0, 3, 3, 3) += H01;
          // H.block(0, 226, 3, 3) += Zero;
          H.block(3, 0, 3, 3) += H01.transpose();
          H.block(3, 3, 3, 3) += H11;
          H.block(3, 6, 3, 3) += H12;
          // H.block(6, 0, 3, 3) += Zero;
          H.block(6, 3, 3, 3) += H12.transpose();
          H.block(6, 6, 3, 3) += H22;
        } else {
          H.block(9, 9, 3, 3) += H00;
          H.block(9, 6, 3, 3) += H01;
          // H.block(9, 3, 3, 3) += Zero;
          H.block(6, 9, 3, 3) += H01.transpose();
          H.block(6, 6, 3, 3) += H11;
          H.block(6, 3, 3, 3) += H12;
          // H.block(3, 9, 3, 3) += Zero;
          H.block(3, 6, 3, 3) += H12.transpose();
          H.block(3, 3, 3, 3) += H22;
        }
      };

      calc_Htri(_hess[ie0], oriente0, //
                d0, d1, d2,           //
                ih0, ih1, ih2,        //
                M0, M1, M2,           //
                R0, R1, R2,           //
                R0d, R1d, R2d);

      calc_Htri(_hess[ie1], oriente1, //
                d1, d2, d0,           //
                ih1, ih2, ih0,        //
                M1, M2, M0,           //
                R1, R2, R0,           //
                R1d, R2d, R0d);

      calc_Htri(_hess[ie2], oriente2, //
                d2, d0, d1,           //
                ih2, ih0, ih1,        //
                M2, M0, M1,           //
                R2, R0, R1,           //
                R2d, R0d, R1d);

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
      const mat1212 &Hi = _hess[e->position_in_set()];

      for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
          triplets.push_back(triplet(ii[i], ii[j], Hi(i, j)));
        }
      }
    }
  }

  virtual void update(const vecX &vals) { _positions = vals; }

  real _k = 0.0;
  surf_ptr _surf;
  vecX _positions;

  std::vector<real> _phi0;
  std::vector<real> _ai;

  std::vector<vec12> _grad;
  std::vector<mat1212> _hess;

  std::vector<real> _inv_e_len;             // inverse_edge_lengths
  std::vector<coordinate_type> _unit_e_vec; // unit edge vector
  std::vector<real> _phi1;                  // phi1
  std::vector<real> _psi_p;                 // psi derivative
  std::vector<real> _psi_pp;                // psi double derivative
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

  virtual vecX project_gradient(const vecX &x) = 0;
  virtual vecX calc_gradient() = 0;
  virtual sparmat calc_hessian() = 0;
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
  virtual vecX calc_gradient() {
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

  virtual sparmat calc_hessian() {
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
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<constraint_set<SPACE>> ptr;
  static ptr create(const surf_ptr surf) {
    return std::make_shared<constraint_set<SPACE>>(surf);
  }

  constraint_set(const surf_ptr surf) : _surf(surf) {
    _reg = m2::ci::geometric_mean_length<SPACE>(_surf);
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

  virtual vecX calc_gradient() {
    int N = _N;
    std::vector<triplet> grad_triplets;
    vecX G(3 * N);
    G.setZero();
    for (auto c : constraints) {
      c->fill_gradient(G);
    }
    std::cout << " grad sum: " << G.sum() << std::endl;
    return G;
  }

  virtual sparmat calc_hessian() {
    int N = _N;
    std::vector<triplet> hess_triplets;
    for (auto c : constraints) {
      c->get_hess_triplets(hess_triplets);
    }
    sparmat H(3 * N, 3 * N);

    H.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

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
    }
    H = I - H;

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
      vecX g = fcn->calc_gradient();

      // g = fcn->project_gradient(g);
      std::cout << "gnorm: " << g.norm() << std::endl;

      if (g.norm() < 1e-10)
        return;

      real lambda = calc_lambda(g, fcn);
      std::cout << "lambda: " << lambda << std::endl;

      vecX xk1 = xk0 - lambda * g;
      tol = (xk1 - xk0).norm();
      // std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << "      tol: " << tol << std::endl;

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
    while (tol > 1e-10 && k++ < 10) {
      fcn->preprocess();

      std::vector<triplet> hess_triplets;
      std::vector<triplet> grad_triplets;
      vecX g = fcn->calc_gradient();
      sparmat H = fcn->calc_hessian();

      _H = H;
      // break;
#if 1
      vecX h = solve(H, g);
      vecX xk1 = xk0 - h;
      _g = h;
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

      // std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << "      tol: " << tol << std::endl; //<< '\r'

      fcn->set_x(xk1);
      std::swap(xk0, xk1);
    }
  }

  virtual vecX extract_gradient() { return _g; };
  virtual vector<mat3> extract_hessian_block_diag() {
    int block_size = 3;
    int N = _H.rows();
    int Nb = N / block_size;

    vector<mat3> blocks(Nb);
    for (auto &block : blocks)
      block.setZero();

#if 0
    sparmat E(N, block_size);
    std::vector<triplet> triplets;
    std::cout << N << " " << Nb << std::endl;
    for (int i = 0; i < Nb; i++) {
      for (int j = 0; j < block_size; j++) {
        triplets.push_back(triplet(block_size * i + j, j, 1.0));
      }
    }

    E.setFromTriplets(triplets.begin(), triplets.end());
    E = _H * E;
    std::cout << E.rows() << " " << E.cols() << std::endl;

    std::cout << blocks.size() << std::endl;
    std::cout << E.rows() << " " << E.cols() << std::endl;

    for (int k = 0; k < E.outerSize(); ++k)
      for (typename sparmat::InnerIterator it(E, k); it; ++it) {
        real val = it.value();
        int row = it.row(); // row index
        int col = it.col(); // col index (here it is equal to k)
        int i = row / block_size;
        int ii = row % block_size;
        //assert(i < blocks.size());
        mat3 &block = blocks[i];
        if(ii < 3 && col < 3)
          block(ii, col) = val;
        // it.index(); // inner index, here it is equal to it.row()
      }
#else
    for (int k = 0; k < _H.outerSize(); ++k)
      for (typename sparmat::InnerIterator it(_H, k); it; ++it) {
        real val = it.value();
        int row = it.row(); // row index
        int col = it.col(); // col index (here it is equal to k)
        int i = row / block_size;
        int ii = row % block_size;
        int jj = col % block_size;
        // assert(i < blocks.size());
        mat3 &block = blocks[i];
        if (abs(row - col) < 3)
          block(ii, jj) = val;
        // it.index(); // inner index, here it is equal to it.row()
      }
#endif
    std::cout << "4" << std::endl;
    return blocks;
  };

  vecX _g;
  sparmat _H;
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