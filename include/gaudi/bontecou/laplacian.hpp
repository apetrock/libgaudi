/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __KUSAMA_LAPLACE_REF_MAT__
#define __KUSAMA_LAPLACE_REF_MAT__

// #include "Eigen/src/SparseCore/SparseMatrix.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gaudi/common.h"
#include "gaudi/sparse_solver.h"

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <math.h>
#include "gaudi/geometry_logger.hpp"

namespace gaudi {
namespace bontecou {
#ifndef GAUDI_BONTECOU_SCALAR_ALIASES
#define GAUDI_BONTECOU_SCALAR_ALIASES
using real = double;
using index_t = int;
#endif

using triplet = Eigen::Triplet<double>;
using vec3 = Eigen::Matrix<real, 3, 1>;
using vec4 = Eigen::Matrix<real, 4, 1>;
using corner1 = std::array<index_t, 1>;
using corner2 = std::array<index_t, 2>;
using corner4 = std::array<index_t, 4>;

template <int S>
Eigen::SparseMatrix<real>
// refactor to iterate over edges/faces, etc
// build_edge_vert
// build_face_vert
//..etc
build_lap(asawa::shell::shell &M,     //
          const std::vector<vec3> &x, //
          std::function<real(asawa::shell::shell &M, asawa::shell::CornerId c,
                             const std::vector<vec3> &x)>
              func_ij,
          bool set_ij = true) {

  index_t vert_count = M.vert_count();
  index_t edge_count = M.edge_count();

  if (!M.verts_are_dense_packed()) {
    std::cerr << "[bontecou::build_lap] Shell vertices must be dense-packed: every "
                 "slot 0..vert_count-1 active (vbegin >= 0). "
                 "Call asawa::shell::pack(M) and keep vertex data permuted with it.\n";
    throw std::runtime_error("bontecou::build_lap: shell vertices not dense-packed");
  }
  if (static_cast<size_t>(x.size()) != static_cast<size_t>(vert_count)) {
    std::cerr << "[bontecou::build_lap] Position count " << x.size()
              << " != vert_count " << vert_count << "\n";
    throw std::runtime_error("bontecou::build_lap: position vector size mismatch");
  }

  std::vector<triplet> tripletList;
  tripletList.reserve(S * vert_count + S * edge_count);

  real Kmin = 9999;
  real Kmax = -9999;

  // Rows/columns use vertex id == row of x[vid]; valid only when verts_are_dense_packed().
  for (index_t vid = 0; vid < vert_count; ++vid) {
    const asawa::shell::VertId v = asawa::shell::vert_id(static_cast<int>(vid));

    real Km = 0.0;
    M.for_each_vertex(
        v, [&](asawa::shell::CornerId c, asawa::shell::shell &Ms) {
          const index_t j = static_cast<index_t>(Ms.vert(Ms.next(c)));
          real K = func_ij(Ms, c, x);
          Km += K;
          if (!set_ij)
            return;
          for (int kk = 0; kk < S; kk++)
            tripletList.push_back(triplet(S * vid + kk, S * j + kk, K));
        });

    Kmin = std::min(Kmin, Km);
    Kmax = std::max(Kmax, Km);
    for (int k = 0; k < S; k++)
      tripletList.push_back(triplet(S * vid + k, S * vid + k, -Km));
  }
  std::cout << " min/max K: " << Kmin << "/" << Kmax << std::endl;

  Eigen::SparseMatrix<real> mat(S * vert_count, S * vert_count);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());

  return mat;
}

void print_lap(asawa::shell::shell &M,     //
               const std::vector<vec3> &x, //
               std::function<real(asawa::shell::shell &M, asawa::shell::CornerId c,
                                  const std::vector<vec3> &x)>
                   func_ij,
               bool set_ij = true) {

  auto verts = M.get_vert_range();
  std::vector<index_t> edges = M.get_edge_vert_ids();

  int i = 0;
  real Kmin = 9999;
  real Kmax = -9999;
  for (auto v : M.get_vert_range()) {
    real Km = 0.0;
    M.for_each_vertex(
        asawa::shell::vert_id(v),
        [&Km, &x, &set_ij, i, func_ij](asawa::shell::CornerId c,
                                        asawa::shell::shell &M) {
          index_t j = M.vert(M.next(c));
          real K = func_ij(M, c, x);
          Km -= K;
          if (!set_ij)
            return;
          std::cout << " " << K;
        });

    std::cout << " - " << Km << std::endl;
    M.vprintv(asawa::shell::vert_id(v));
    Kmin = std::min(Kmin, Km);
    Kmax = std::max(Kmax, Km);
    i++;
  }
  std::cout << " min/max K: " << Kmin << "/" << Kmax << std::endl;
}

class laplacian {
public:
  typedef Eigen::SparseMatrix<real> sparmat;

  /// When \p defer_init_C is true, only \ref initM runs; stiffness \ref _matC is left
  /// zero-sized until \ref set_stiffness or \ref initC is called.
  explicit laplacian(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
                     bool defer_init_C = false)
      : __M(M), __x(x) {
    initM();
    if (!defer_init_C)
      initC();
  }

  ~laplacian() {}

  /// Replace cotan stiffness with an externally-built matrix (e.g. anisotropic cotan).
  /// Caller must supply \c nv×\c nv sparse matrix matching \ref build_lap conventions.
  void set_stiffness(const sparmat &C) {
    _matC = C;
    if (_matM.rows() == 0)
      initM();
  }

  void printC() {
    print_lap(
        *__M, __x, //
        [](asawa::shell::shell &M, asawa::shell::CornerId c,
           const std::vector<vec3> &x) {
          asawa::shell::CornerId c0p = M.prev(c);
          asawa::shell::CornerId c1p = M.prev(M.other(c));

          return asawa::shell::weak_cotan(M, c0p, x) +
                 asawa::shell::weak_cotan(M, c1p, x);
        });
  }

  void initC() {
    _matC = build_lap<1>(
        *__M, __x, //
        [](asawa::shell::shell &M, asawa::shell::CornerId c,
           const std::vector<vec3> &x) {
          asawa::shell::CornerId c0p = M.prev(c);
          asawa::shell::CornerId c1p = M.prev(M.other(c));
          return asawa::shell::weak_cotan(M, c0p, x) +
                 asawa::shell::weak_cotan(M, c1p, x);
        });
  }

  void initM() {
    _matM = build_lap<1>(
        *__M, __x,
        [](asawa::shell::shell &M, asawa::shell::CornerId c,
           const std::vector<vec3> &x) {
          real aj = asawa::shell::face_area(M, M.face(c), x);
          aj = -max(aj, 1e-6);
          return aj;
        },
        false);
  }

  void init() {
    this->initM();
    this->initC();
  }

  std::vector<real> multC(const std::vector<real> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(U.data(), U.size());

    Eigen::VectorXd Le = _matC * Ue;
    return std::vector<real>(Le.data(), Le.data() + Le.rows() * Le.cols());
  }

  std::vector<real> multM(const std::vector<real> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(U.data(), U.size());

    Eigen::VectorXd Le = _matM * Ue;
    return std::vector<real>(Le.data(), Le.data() + Le.rows() * Le.cols());
  }

  virtual vecX solve_2(sparmat &A, vecX &b) {
    // this is correct way, but does not work at moment
    m_solver S(A);
    if (!S.success()) {
      std::cout << "Solve failed" << std::endl;
      return vecX();
    }
    return S.solve(b);
  }

  virtual vecX solve(sparmat &A, vecX &b) {
#if 1
    // Eigen::SimplicialLLT<sparmat> solver;
    Eigen::SimplicialLDLT<sparmat> solver;
    //   Eigen::ConjugateGradient<sparmat> solver;
    //   Eigen::ConjugateGradient<sparmat, Eigen::Lower | Eigen::Upper,
    //                          Eigen::DiagonalPreconditioner<double>> solver;
    // Eigen::ConjugateGradient<sparmat, Eigen::Lower | Eigen::Upper,
    //                         Eigen::IncompleteCholesky<double>>
    solver;

    solver.compute(A);
//    std::cout << " solver det: " << solver.determinant() << std::endl;
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
#if 0
    real mn = 0.0;
    std::cout << mn << std::endl;
    for (int i = 0; i < x.rows(); i++) {
      mn += x[i];
    }

    mn /= real(x.rows());
    for (int i = 0; i < x.rows(); i++) {
      if (x[i] > 1000.0 * mn)
        std::cout << "fooey: " << i << " " << mn << " " << x[i] << std::endl;
    }
#endif
    return x;
  }

  std::vector<real> solve(const std::vector<real> &f) {

    using namespace asawa;
    // sparmat L = initML();

    Eigen::SparseMatrix<real> I(_matC.rows(), _matC.cols());
    I.setIdentity();
    sparmat L = _matC + 1e-6 * I;

    vecX fe = to(f);
    // std::cout << " L sum: " << L.sum() << std::endl;
    // std::cout << " f norm: " << fe.norm() << std::endl;
    vecX x = solve(L, fe);
    // std::cout << " x norm: " << x.norm() << std::endl;
    return from(x);
  }

  std::vector<real> diffuse(const std::vector<real> &f, const real &dt) {

    using namespace asawa;

    sparmat M = _matM;
    sparmat L = _matC;
    Eigen::SparseMatrix<real> I(_matC.rows(), _matC.cols());
    I.setIdentity();
#if 1
    sparmat A = M - dt * L;
    vecX fe = M * to(f);
#else
    sparmat A = I - dt * L;
    vecX fe = to(f);
#endif
    return from(solve(A, fe));
  }

  std::vector<real> diffuse2(const std::vector<real> &u0, const real &h) {

    using namespace asawa;

    sparmat M = _matM;
    sparmat L = _matC;
    Eigen::SparseMatrix<real> I(_matC.rows(), _matC.cols());
    I.setIdentity();
#if 1
    sparmat A = 2.0 * M - h * L;
    vecX u0v = to(u0);
    vecX b = 2.0 * M * u0v + h * L * u0v;
#else
    sparmat A = I - dt * L;
    vecX fe = to(f);
#endif
    return from(solve(A, b));
  }

  std::vector<real> heatDist(const std::vector<real> &f, real dt) {

    std::vector<real> u = diffuse(f, dt);
    int i = 0;

#if 0
    for (int i = 0; i < __M->vert_count(); i++) {
      vec3 N = vert_normal(*__M, i, __x);
      vec3 xi = __x[i];
      // gg::geometry_logger::line(xi, xi + f[i] * N, vec4(0.9, 0.2,
      // 0.25, 1.0));
      gg::geometry_logger::line(xi, xi + u[i] * N, vec4(0.9, 0.2, 0.25, 1.0));
    }
#endif

    std::vector<vec3> gradU = asawa::shell::gradient(*__M, u, __x);

    i = 0;

    std::for_each(gradU.begin(), gradU.end(), [](auto &v) {
      v.normalize();
      v *= -1;
    });
#if 0
    for (int i = 0; i < __M->face_count(); i++) {
      vec3 xi = asawa::shell::face_center(*__M, i, __x);
      vec3 gu = 0.02 * gradU[i];
      std::cout << gu.transpose() << std::endl;
      gg::geometry_logger::line(xi, xi + gu, vec4(0.8, 1.0, 0.35, 1.0));
    }
#endif

    std::vector<real> divu = asawa::shell::divergence(*__M, gradU, __x);

    std::vector<real> x = solve(divu);
    auto res = std::minmax_element(begin(x), end(x));
    real mn = *res.first;
    real mx = *res.second;

    if (mx - mn > 1e-16) {
      std::for_each(begin(x), end(x),
                    [mn, mx](real &v) { v = (v - mn) / (mx - mn); });
      // std::for_each(begin(x), end(x), [mn, mx](real &v) { v = (v - mn); });
    }
    // std::vector<real> x = lap.diffuse(divu, dt);

    return x;
  }

  std::vector<mat3> heatFrame(const std::vector<real> &f, real dt) {
    // frame centered on the faces
    std::vector<real> u = diffuse(f, dt);
    int i = 0;
    std::vector<vec3> N = asawa::shell::face_normals(*__M, __x);
    std::vector<vec3> gradU = asawa::shell::gradient(*__M, u, __x);

    std::vector<mat3> frame(__M->face_count(), mat3::Identity());
    i = 0;
    for (auto gi : gradU) {
      gi.normalize();
      gi *= -1;
      mat3 F;
      F.col(0) = gi;
      F.col(2) = N[i];
      F.col(1) = gi.cross(N[i]).normalized();
      frame[i] = F;
      i++;
    }
    return frame;
  }

  /// Cotangent stiffness matrix (weak Laplacian / graph Laplacian sign used here).
  const sparmat &stiffness() const { return _matC; }

  bool inited = false;

private:
  asawa::shell::shell::ptr __M;
  const std::vector<vec3> &__x;
  Eigen::SparseMatrix<real> _matC;
  Eigen::SparseMatrix<real> _matM;
};

class laplacian3 {

public:
  laplacian3(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
             bool unitary = false)
      : __M(M), __x(x), __unitary(unitary) {
    this->init();
  }

  ~laplacian3() {}
  /*L = MinvC*/
  // need one that projects into null space
  void initC() {
    _matC = build_lap<3>(
        *__M, __x, //
        [this](asawa::shell::shell &M, asawa::shell::CornerId c,
               const std::vector<vec3> &x) {
          asawa::shell::CornerId c0p = M.prev(c);
          asawa::shell::CornerId c1p = M.prev(M.other(c));
          if (__unitary)
            return real(1.0);
          return asawa::shell::weak_cotan(M, c0p, x) +
                 asawa::shell::weak_cotan(M, c1p, x);
        });
  }

  void init() { this->initC(); }

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  std::vector<vec3> mult(const std::vector<vec3> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue = to(U);
    Eigen::VectorXd Le = _matC * Ue;
    std::vector<vec3> Uc(U);
    from(Uc, Le);
    return Uc;
  }

  // need simple null space smoothing
  std::vector<vec3> smooth(const std::vector<vec3> &U, real c0, real c1) {

    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd U0 = to(U);
    Eigen::SparseMatrix<real> I(_matC.rows(), _matC.cols());
    I.setIdentity();
    Eigen::VectorXd U1 = (I + c0 * _matC) * U0;
    Eigen::VectorXd U2 = (I - c1 * _matC) * U1;
    std::vector<vec3> Uc(U);
    from(Uc, U2);
    return Uc;
  }
  bool inited = false;

private:
  Eigen::SparseMatrix<real> _matC;

  bool __unitary;
  asawa::shell::shell::ptr __M;
  const std::vector<vec3> &__x;
};

std::array<real, 2> grey_scott(std::array<real, 2> uv0a, const real &f,
                               const real &k, const real &h, const int N = 40) {
  // first order formulation
  vec2 uv0 = vec2(uv0a[0], uv0a[1]);
  vec2 uvi = uv0;
  for (int j = 0; j < N; j++) {
    real u = uvi[0];
    real v = uvi[1];
    real uv2 = u * v * v;
    real v2 = v * v;
    real uv = u * v;
    vec2 G;
    G(0) = (-uv2 + f * (1.0 - u));
    G(1) = (uv2 - (f + k) * v);

    vec2 F = uvi - uv0 - h * G;
    mat2 dG;
    dG(0, 0) = -(v2 + f);
    dG(0, 1) = -2.0 * uv;
    dG(1, 0) = v2;
    dG(1, 1) = 2.0 * uv - (f + k);

    mat2 dF = mat2::Identity() - h * dG;
    // Compute LU decomposition of dF
    Eigen::PartialPivLU<mat2> lu = dF.partialPivLu();
    // Solve the linear system using LU decomposition
    uvi += lu.solve(-F);

    if (F.norm() < 1.0e-12)
      break;
  }
  return {uvi[0], uvi[1]};
}

std::array<real, 2> grey_scott_2(std::array<real, 2> uv0a, const real &f,
                                 const real &k, const real &h,
                                 const int N = 40) {
  // second order formulation
  vec2 uv0 = vec2(uv0a[0], uv0a[1]);
  real u0 = uv0[0];
  real v0 = uv0[1];
  real uv20 = u0 * v0 * v0;
  vec2 G0;
  G0(0) = (-uv20 + f * (1.0 - u0));
  G0(1) = (uv20 - (f + k) * v0);

  vec2 uvi = uv0;

  for (int j = 0; j < N; j++) {
    real u = uvi[0];
    real v = uvi[1];
    real uv2 = u * v * v;

    real v2 = v * v;
    real uv = u * v;
    vec2 G;
    G(0) = (-uv2 + f * (1.0 - u));
    G(1) = (uv2 - (f + k) * v);

    vec2 F = uvi - uv0 - 0.5 * h * (G + G0);
    mat2 dG;
    dG(0, 0) = -(v2 + f);
    dG(0, 1) = -2.0 * uv;
    dG(1, 0) = v2;
    dG(1, 1) = 2.0 * uv - (f + k);

    mat2 dF = mat2::Identity() - 0.5 * h * dG;
    // Compute LU decomposition of dF
    Eigen::PartialPivLU<mat2> lu = dF.partialPivLu();
    // Solve the linear system using LU decomposition
    uvi += lu.solve(-F);

    if (F.norm() < 1.0e-12)
      break;
  }
  return {uvi[0], uvi[1]};
}

} // namespace bontecou
} // namespace gaudi
#endif
