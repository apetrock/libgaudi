/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __LAPLACE_MAT__
#define __LAPLACE_MAT__

//#include "Eigen/src/SparseCore/SparseMatrix.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "manifold/asawa/m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <math.h>

namespace bontecou {
  
using namespace asawa;
template <typename SPACE, typename TYPE> class laplacian_base {
  M2_TYPEDEFS;

public:
  laplacian_base(asawa::surf<SPACE> *surf) : _surf(surf) { this->init(); }

  ~laplacian_base() {}

  virtual void init() { std::cout << "in base init" << std::endl; }

  Eigen::SparseMatrix<real>
  build(std::function<real(typename surf<SPACE>::face_vertex_ptr fv, real &Km)>
            func) {
    auto vertices = _surf->get_vertices();
    auto edges = _surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;
    for (auto v : _surf->get_vertices()) {

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList, func](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            real K = func(fv, Km);
            if (K > 0)
              tripletList.push_back(triplet(i, j, K));
          });
      tripletList.push_back(triplet(i, i, Km));
      i++;
    }

    Eigen::SparseMatrix<real> mat(vertices.size(), vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
  }

  Eigen::SparseMatrix<real>
  build(std::function<real(face_vertex_ptr fv, real &Km)> funcij,
        std::function<real(real &Km)> funcii) {
    auto vertices = _surf->get_vertices();
    auto edges = _surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;
    for (auto v : _surf->get_vertices()) {

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList, funcij](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            real K = funcij(fv, Km);
            if (K > 0)
              tripletList.push_back(triplet(i, j, K));
          });

      Km = funcii(Km);
      tripletList.push_back(triplet(i, i, Km));
      i++;
    }

    Eigen::SparseMatrix<real> mat(vertices.size(), vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
  }

  asawa::surf<SPACE> *_surf;
};

template <typename SPACE, typename TYPE>
class laplacian : laplacian_base<SPACE, TYPE> {
  M2_TYPEDEFS;

public:
  laplacian(asawa::surf<SPACE> *surf) : laplacian_base<SPACE, TYPE>(surf) {
    this->init();
  }

  ~laplacian() {}
  /*L = MinvC*/
  /*
    void initC() {
      _matC = this->build([](face_vertex_ptr fv, real &Km) {
        real cote = asawa::ci::cotan<SPACE>(fv->edge());
        assert(!isnan(cote));
        real K = cote;
        Km -= K;
        return K;
      });
      // std::cout << _matM << std::endl;
    }
  */
  void initC() {
    auto vertices = this->_surf->get_vertices();
    auto edges = this->_surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;

    for (auto v : this->_surf->get_vertices()) {

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            real K = asawa::ci::cotan<SPACE>(fv->edge());
            tripletList.push_back(triplet(i, j, K));
            Km -= K;
          });
      tripletList.push_back(triplet(i, i, Km));
      i++;
    }
    Eigen::SparseMatrix<real> mat(vertices.size(), vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    _matC = mat;
  }

  Eigen::SparseMatrix<real> initML() {
    auto vertices = this->_surf->get_vertices();
    auto edges = this->_surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;

    for (auto v : this->_surf->get_vertices()) {
      real A = 0.0;
      asawa::for_each_vertex<SPACE>(v, [&A, &tripletList](face_vertex_ptr fv) {
        real Aj = asawa::ci::area<SPACE>(fv->face());
        A += Aj / 3.0;
      });

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, A, &tripletList](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            real K = asawa::ci::cotan<SPACE>(fv->edge());
            K = K < 1e-16 ? 1.0 : K;
            K /= 2.0 * A;
            tripletList.push_back(triplet(i, j, K));
            Km -= K;
          });
      tripletList.push_back(triplet(i, i, Km));
      i++;
    }
    Eigen::SparseMatrix<real> mat(vertices.size(), vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
  }

  void initM() {
    auto vertices = this->_surf->get_vertices();
    auto edges = this->_surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;
    for (auto v : this->_surf->get_vertices()) {
      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList](face_vertex_ptr fv) {
            real aj = asawa::ci::area<SPACE>(fv->face());
            // real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
            real l = 0.33;
            Km += l * aj;
          });

      tripletList.push_back(triplet(i, i, Km));
      i++;
    }
    Eigen::SparseMatrix<real> mat(vertices.size(), vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    _matM = mat;
  }

  void init() {
    this->initM();
    this->initC();
  }

  std::vector<TYPE> multC(const std::vector<TYPE> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(U.data(), U.size());

    Eigen::VectorXd Le = _matC * Ue;
    return std::vector<real>(Le.data(), Le.data() + Le.rows() * Le.cols());
  }

  std::vector<TYPE> multM(const std::vector<TYPE> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(U.data(), U.size());

    Eigen::VectorXd Le = _matM * Ue;
    return std::vector<real>(Le.data(), Le.data() + Le.rows() * Le.cols());
  }

  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;

  vecX to(const std::vector<TYPE> &U) {
    vecX Ue = Eigen::Map<const vecX, Eigen::Unaligned>(U.data(), U.size());
    return Ue;
  }

  std::vector<TYPE> from(vecX U) {
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

  std::vector<TYPE> solve(const std::vector<TYPE> &f) {

    using namespace asawa;
    M2_TYPEDEFS;
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

  std::vector<TYPE> diffuse(const std::vector<TYPE> &f, const real &dt) {

    using namespace asawa;
    M2_TYPEDEFS;
    asawa::surf<SPACE> *surf = this->_surf;

    sparmat M = _matM;
    sparmat L = _matC;
    sparmat A = M - dt * L;
    vecX fe = to(f);

    return from(solve(A, fe));
  }

  std::vector<real> heatDist(const std::vector<real> &f, real dt) {

    asawa::surf<SPACE> *surf = this->_surf;
    std::vector<real> u = diffuse(f, dt);
    int i = 0;

    coordinate_array gradU = ci::gradient<SPACE>(u, surf);

    i = 0;

    std::for_each(gradU.begin(), gradU.end(), [](auto &v) {
      v.normalize();
      v *= -1;
    });

    std::vector<real> divu = ci::divergence<SPACE>(gradU, surf);

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
  Eigen::SparseMatrix<real> _matC;
  Eigen::SparseMatrix<real> _matM;
};

template <typename SPACE>
class laplacian3 : laplacian_base<SPACE, typename SPACE::coordinate_type> {
  M2_TYPEDEFS;

public:
  laplacian3(asawa::surf<SPACE> *surf)
      : laplacian_base<SPACE, coordinate_type>(surf) {
    this->init();
  }

  ~laplacian3() {}
  /*L = MinvC*/

  Eigen::SparseMatrix<real> build() {
    auto vertices = this->_surf->get_vertices();
    auto edges = this->_surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;

    for (auto v : this->_surf->get_vertices()) {

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            real K = asawa::ci::cotan<SPACE>(fv->edge());
            tripletList.push_back(triplet(3 * i + 0, 3 * j + 0, K));
            tripletList.push_back(triplet(3 * i + 1, 3 * j + 1, K));
            tripletList.push_back(triplet(3 * i + 2, 3 * j + 2, K));
            Km -= K;
          });

      tripletList.push_back(triplet(3 * i + 0, 3 * i + 0, Km));
      tripletList.push_back(triplet(3 * i + 1, 3 * i + 1, Km));
      tripletList.push_back(triplet(3 * i + 2, 3 * i + 2, Km));
      i++;
    }
    Eigen::SparseMatrix<real> mat(3 * vertices.size(), 3 * vertices.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
  }

  void init() { _mat = this->build(); }

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

  std::vector<coordinate_type> mult(const std::vector<coordinate_type> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ue = to(U);
    Eigen::VectorXd Le = _mat * Ue;
    std::vector<coordinate_type> Uc(U);
    from(Uc, Le);
    return Uc;
  }

  std::vector<coordinate_type> smooth(const std::vector<coordinate_type> &U,
                                      real c0, real c1) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd U0 = to(U);
    Eigen::SparseMatrix<real> I(_mat.rows(), _mat.cols());
    I.setIdentity();
    Eigen::VectorXd U1 = (I + c0 * _mat) * U0;
    Eigen::VectorXd U2 = (I - c1 * _mat) * U1;
    std::vector<coordinate_type> Uc(U);
    from(Uc, U2);
    return Uc;
  }

  bool inited = false;

private:
  Eigen::SparseMatrix<real> _mat;
};

template <typename SPACE, typename TYPE> class area_laplacian_0 {
  M2_TYPEDEFS;

public:
  area_laplacian_0(asawa::surf<SPACE> *surf) { _surf = surf; }

  ~area_laplacian_0() {}
  /*L = MinvC*/
  std::vector<TYPE> mult(const std::vector<TYPE> &U) {
    std::vector<TYPE> L(U.size());
    int i = 0;
    // real aAvg = 0;

    for (auto v : _surf->get_vertices()) {
      TYPE ui = U[i];
      TYPE u = z::zero<TYPE>();
      double area = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&u, &ui, &U, &area](face_vertex_ptr fv) {
            face_vertex_ptr fvp = fv->vprev()->next();
            face_vertex_ptr fvn = fv->vnext()->next();

            real aj = asawa::ci::area<SPACE>(fv->face()) + 1e-6;
            real l = fv->template get<real>(SPACE::face_vertex_index::BARY);

            int j = fv->next()->vertex()->position_in_set();

            TYPE uj = U[j];
            real K = l * aj;
            area += K;

            u += K * (uj - ui);
          });

      if (area < 1e-10)
        L[i] = z::zero<TYPE>();
      else
        L[i] = u / area;

      i++;
    }

    return L;
  }

private:
  asawa::surf<SPACE> *_surf;
};

template <typename SPACE, typename TYPE>
class area_laplacian : public laplacian_base<SPACE, TYPE> {
  M2_TYPEDEFS;

public:
  area_laplacian(asawa::surf<SPACE> *surf) : laplacian_base<SPACE, TYPE>(surf) {
    this->init();
  }

  ~area_laplacian() {}
  /*L = MinvC*/

  void initM() {
    _matM = this->build([](face_vertex_ptr fv, real &Km) {
      real aj = asawa::ci::area<SPACE>(fv->face());
      real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
      real K = l * aj + 1e-6;
      // K = 0.00001;
      Km -= K;
      return K;
    });
  }

  virtual void init() { this->initM(); }

  std::vector<coordinate_type> mult(const std::vector<coordinate_type> &U) {
    // Eigen::VectorXd test(U.data());
    Eigen::VectorXd Ux(U.size()), Uz(U.size()), Uy(U.size());

    int i = 0;

    for (auto v : U) {
      Ux[i] = v[0];
      Uy[i] = v[1];
      Uz[i] = v[2];
    }

    Eigen::VectorXd Lx = (_matM * Ux);
    Eigen::VectorXd Ly = (_matM * Uy);
    Eigen::VectorXd Lz = (_matM * Uz);

    std::vector<coordinate_type> Lxyz(Lx.size());
    auto diags = _matM.diagonal();

    for (int i = 0; i < Lx.size(); i++) {
      Lxyz[i] = coordinate_type(Lx[i], Ly[i], Lz[i]) / diags[i];
    }

    return Lxyz;
  }

private:
  Eigen::SparseMatrix<real> _matM;
};

template <typename SPACE> class cotan_curvature {
  M2_TYPEDEFS;

public:
  cotan_curvature(asawa::surf<SPACE> *surf) { _surf = surf; }

  ~cotan_curvature() {}
  /*L = MinvC*/
  std::vector<real> operator()() {
    int i = 0;

    std::vector<typename SPACE::real> K(_surf->get_vertices().size(), 0.0);
    for (auto vi : _surf->get_vertices()) {

      coordinate_type pi = ci::get_coordinate<SPACE>(vi);
      real k = 0.0;

      asawa::for_each_vertex<SPACE>(vi, [pi, &k, &K](face_vertex_ptr fv) {
        face_vertex_ptr fvp = fv->vprev()->next();
        face_vertex_ptr fvn = fv->vnext()->next();
        real cotp = asawa::ci::abs_cotan<SPACE>(fvp);
        real cotn = asawa::ci::abs_cotan<SPACE>(fvn);

        vertex_ptr vj = fv->next()->vertex();

        coordinate_type pj = ci::get_coordinate<SPACE>(vj);

        k += (cotp + cotn) * va::norm(coordinate_type(pj - pi));
      });
      K[i] = k;

      i++;
    }
    return K;
  }

private:
  asawa::surf<SPACE> *_surf;
};

} // namespace bontecou
#endif
