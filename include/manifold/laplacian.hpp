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

#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <math.h>

namespace m2 {
template <typename SPACE, typename TYPE> class laplacian_base {
  M2_TYPEDEFS;

public:
  laplacian_base(m2::surf<SPACE> *surf) : _surf(surf) { this->init(); }

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
      m2::for_each_vertex<SPACE>(
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
      m2::for_each_vertex<SPACE>(
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

  m2::surf<SPACE> *_surf;
};

template <typename SPACE, typename TYPE>
class laplacian : laplacian_base<SPACE, TYPE> {
  M2_TYPEDEFS;

public:
  laplacian(m2::surf<SPACE> *surf) : laplacian_base<SPACE, TYPE>(surf) {
    this->init();
  }

  ~laplacian() {}
  /*L = MinvC*/

  void initC() {
    _matC = this->build([](face_vertex_ptr fv, real &Km) {
      real cote = m2::ci::cotan<SPACE>(fv->edge());
      assert(!isnan(cote));
      real K = cote;
      Km -= K;
      return K;
    });
    // std::cout << _matM << std::endl;
  }

  void initM() {
    _matM = this->build([](face_vertex_ptr fv, real &Km) {
      real aj = m2::ci::area<SPACE>(fv->face());
      real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
      // real l = 0.33;
      Km += l * aj;
      return 0;
    });
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

  bool inited = false;

private:
  Eigen::SparseMatrix<real> _matC;
  Eigen::SparseMatrix<real> _matM;
};

template <typename SPACE>
class laplacian3 : laplacian_base<SPACE, typename SPACE::coordinate_type> {
  M2_TYPEDEFS;

public:
  laplacian3(m2::surf<SPACE> *surf)
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
      m2::for_each_vertex<SPACE>(v, [&Km, i, &tripletList](face_vertex_ptr fv) {
        int j = fv->next()->vertex()->position_in_set();
        real K = m2::ci::cotan<SPACE>(fv->edge());
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
  area_laplacian_0(m2::surf<SPACE> *surf) { _surf = surf; }

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
      m2::for_each_vertex<SPACE>(v, [&u, &ui, &U, &area](face_vertex_ptr fv) {
        face_vertex_ptr fvp = fv->vprev()->next();
        face_vertex_ptr fvn = fv->vnext()->next();

        real aj = m2::ci::area<SPACE>(fv->face()) + 1e-6;
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
  m2::surf<SPACE> *_surf;
};

template <typename SPACE, typename TYPE>
class area_laplacian : public laplacian_base<SPACE, TYPE> {
  M2_TYPEDEFS;

public:
  area_laplacian(m2::surf<SPACE> *surf) : laplacian_base<SPACE, TYPE>(surf) {
    this->init();
  }

  ~area_laplacian() {}
  /*L = MinvC*/

  void initM() {
    _matM = this->build([](face_vertex_ptr fv, real &Km) {
      real aj = m2::ci::area<SPACE>(fv->face());
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
  cotan_curvature(m2::surf<SPACE> *surf) { _surf = surf; }

  ~cotan_curvature() {}
  /*L = MinvC*/
  std::vector<real> operator()() {
    int i = 0;

    std::vector<typename SPACE::real> K(_surf->get_vertices().size(), 0.0);
    for (auto vi : _surf->get_vertices()) {

      coordinate_type pi = ci::get_coordinate<SPACE>(vi);
      real k = 0.0;

      m2::for_each_vertex<SPACE>(vi, [pi, &k, &K](face_vertex_ptr fv) {
        face_vertex_ptr fvp = fv->vprev()->next();
        face_vertex_ptr fvn = fv->vnext()->next();
        real cotp = m2::ci::abs_cotan<SPACE>(fvp);
        real cotn = m2::ci::abs_cotan<SPACE>(fvn);

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
  m2::surf<SPACE> *_surf;
};

} // namespace m2
#endif
