/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __ASAWA_LAPLACE_REF_MAT__
#define __ASAWA_LAPLACE_REF_MAT__

//#include "Eigen/src/SparseCore/SparseMatrix.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gaudi/asawa/datum_x.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/manifold.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <math.h>

namespace gaudi {
namespace bontecou {
using real = double;
using index_t = int;

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
build_lap(asawa::manifold &M,         //
          const std::vector<vec3> &x, //
          std::function<real(asawa::manifold &M, index_t c,
                             const std::vector<vec3> &x)>
              func_ij,
          bool set_ij = true) {

  std::vector<index_t> verts = M.get_vert_range();
  std::vector<index_t> edges = M.get_edge_vert_ids();

  std::vector<triplet> tripletList;
  tripletList.reserve(S * verts.size() + S * edges.size());

  int i = 0;
  real Kmin = 9999;
  real Kmax = -9999;
  for (auto v : M.get_vert_range()) {
    real Km = 0.0;
    M.for_each_vertex(v, [&Km, &x, &set_ij, i, &tripletList,
                          func_ij](index_t c, asawa::manifold &M) {
      index_t j = M.vert(M.next(c));
      real K = func_ij(M, c, x);
      Km -= K;
      if (!set_ij)
        return;
      for (int k = 0; k < S; k++)
        tripletList.push_back(triplet(S * i + k, S * j + k, K));
    });

    Kmin = std::min(Kmin, Km);
    Kmax = std::max(Kmax, Km);

    for (int k = 0; k < S; k++)
      tripletList.push_back(triplet(S * i + k, S * i + k, Km));
    i++;
  }
  std::cout << " min/max K: " << Kmin << "/" << Kmax << std::endl;

  Eigen::SparseMatrix<real> mat(S * verts.size(), S * verts.size());
  mat.setFromTriplets(tripletList.begin(), tripletList.end());

  return mat;
}

void print_lap(asawa::manifold &M,         //
               const std::vector<vec3> &x, //
               std::function<real(asawa::manifold &M, index_t c,
                                  const std::vector<vec3> &x)>
                   func_ij,
               bool set_ij = true) {

  std::vector<index_t> verts = M.get_vert_range();
  std::vector<index_t> edges = M.get_edge_vert_ids();

  int i = 0;
  real Kmin = 9999;
  real Kmax = -9999;
  for (auto v : M.get_vert_range()) {
    real Km = 0.0;
    M.for_each_vertex(
        v, [&Km, &x, &set_ij, i, func_ij](index_t c, asawa::manifold &M) {
          index_t j = M.vert(M.next(c));
          real K = func_ij(M, c, x);
          Km -= K;
          if (!set_ij)
            return;
          std::cout << " " << K;
        });

    std::cout << " - " << Km << std::endl;
    M.vprintv(v);
    Kmin = std::min(Kmin, Km);
    Kmax = std::max(Kmax, Km);
    i++;
  }
  std::cout << " min/max K: " << Kmin << "/" << Kmax << std::endl;
}

class laplacian {
public:
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  laplacian(asawa::manifold::ptr M, const std::vector<vec3> &x)
      : __M(M), __x(x) {
    this->init();
  }

  ~laplacian() {}

  void printC() {
    print_lap(*__M, __x, //
              [](asawa::manifold &M, index_t c, const std::vector<vec3> &x) {
                index_t c0p = M.prev(c);
                index_t c1p = M.prev(M.other(c));
                
                return cotan(M, c0p, x) + cotan(M, c1p, x);
              });
  }

  void initC() {
    _matC = build_lap<1>(
        *__M, __x, //
        [](asawa::manifold &M, index_t c, const std::vector<vec3> &x) {
          index_t c0p = M.prev(c);
          index_t c1p = M.prev(M.other(c));

          return cotan(M, c0p, x) + cotan(M, c1p, x);
        });
  }

  void initM() {
    _matM = build_lap<1>(
        *__M, __x,
        [](asawa::manifold &M, index_t c, const std::vector<vec3> &x) {
          real aj = asawa::face_area(M, M.face(c), x);
          aj = aj < 1e-4 ? 1.0 : aj;
          return -aj / 3.0;
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

  vecX to(const std::vector<real> &U) {
    vecX Ue = Eigen::Map<const vecX, Eigen::Unaligned>(U.data(), U.size());
    return Ue;
  }

  std::vector<real> from(vecX U) {
    return std::vector<real>(U.data(), U.data() + U.rows() * U.cols());
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
    sparmat A = M - dt * L;
    // vecX fe = _matM * to(f);
    vecX fe = to(f);

    return from(solve(A, fe));
  }

  std::vector<real> heatDist(const std::vector<real> &f, real dt) {

    std::vector<real> u = diffuse(f, dt);
    int i = 0;

#if 1
    for (int i = 0; i < __M->vert_count(); i++) {
      vec3 N = vert_normal(*__M, i, __x);
      vec3 xi = __x[i];
      gg::geometry_logger::line(xi, xi + f[i] * N, vec4(0.9, 0.2, 0.25, 1.0));
    }
#endif

    std::vector<vec3> gradU = asawa::gradient(*__M, u, __x);

    i = 0;

    std::for_each(gradU.begin(), gradU.end(), [](auto &v) {
      v.normalize();
      v *= -1;
    });
#if 0
    for (int i = 0; i < __M->face_count(); i++) {
      vec3 xi = asawa::face_center(*__M, i, __x);
      vec3 gu = 0.02 * gradU[i];
      gg::geometry_logger::line(xi, xi + gu, vec4(0.8, 1.0, 0.35, 1.0));
    }
#endif

    std::vector<real> divu = asawa::divergence(*__M, gradU, __x);

    std::vector<real> x = solve(divu);
    auto res = std::minmax_element(begin(x), end(x));
    real mn = *res.first;
    real mx = *res.second;

    if (mx - mn > 1e-16) {
      std::for_each(begin(x), end(x),
                    [mn, mx](real &v) { v = (v - mn) / (mx - mn); });
    }
    // std::vector<real> x = lap.diffuse(divu, dt);

    return x;
  }

  bool inited = false;

private:
  asawa::manifold::ptr __M;
  const std::vector<vec3> &__x;
  Eigen::SparseMatrix<real> _matC;
  Eigen::SparseMatrix<real> _matM;
};

class laplacian3 {

public:
  laplacian3(asawa::manifold::ptr M, const std::vector<vec3> &x)
      : __M(M), __x(x) {
    this->init();
  }

  ~laplacian3() {}
  /*L = MinvC*/

  void initC() {
    _matC = build_lap<3>(
        *__M, __x, //
        [](asawa::manifold &M, index_t c, const std::vector<vec3> &x) {
          index_t c0p = M.prev(c);
          index_t c1p = M.prev(M.other(c));

          return cotan(M, c0p, x) + cotan(M, c1p, x);
        });
  }

  void init() { this->initC(); }

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;

  static vec3 from(const vecX &vals, size_t i) {
    return vec3(vals[3 * i + 0], vals[3 * i + 1], vals[3 * i + 2]);
  };

  static void to(vec3 c, vecX &vals, size_t i) {
    vals[3 * i + 0] = c[0];
    vals[3 * i + 1] = c[1];
    vals[3 * i + 2] = c[2];
  }

  // these will be function pointers with a default
  static vecX to(const std::vector<vec3> &positions) {
    vecX out(3 * positions.size());
    int i = 0;

    for (int i = 0; i < positions.size(); i++) {
      vec3 p = positions[i];
      to(p, out, i);
    }

    return out;
  }

  static void from(std::vector<vec3> &positions, const vecX &x) {
    for (int i = 0; i < positions.size(); i++) {
      positions[i] = from(x, i);
    }
  }

  std::vector<vec3> mult(const std::vector<vec3> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue = to(U);
    Eigen::VectorXd Le = _matC * Ue;
    std::vector<vec3> Uc(U);
    from(Uc, Le);
    return Uc;
  }

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

  asawa::manifold::ptr __M;
  const std::vector<vec3> &__x;
};
} // namespace bontecou
} // namespace gaudi
#endif
