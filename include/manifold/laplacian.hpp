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

#include <cassert>
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>

#include "manifold/conj_grad.hpp"
namespace m2 {
template <typename SPACE, typename TYPE> class laplacian {
  M2_TYPEDEFS;

public:
  laplacian(m2::surf<SPACE> *surf) { _surf = surf; }

  ~laplacian() {}
  /*L = MinvC*/
  std::vector<TYPE> multC(const std::vector<TYPE> &U) {
    std::vector<TYPE> L(U.size());
    int i = 0;
    // real aAvg = 0;

    for (auto v : _surf->get_vertices()) {
      TYPE ui = U[i];
      TYPE u = z::zero<TYPE>();

      m2::for_each_vertex<SPACE>(v, [&u, &ui, &U](face_vertex_ptr fv) {
        face_vertex_ptr fvp = fv->vprev()->next();
        face_vertex_ptr fvn = fv->vnext()->next();

        int j = fv->next()->vertex()->position_in_set();
        TYPE uj = U[j];

        real cotp = m2::ci::abs_cotan<SPACE>(fvp);
        real cotn = m2::ci::abs_cotan<SPACE>(fvn);
        assert(!isnan(cotp));
        assert(!isnan(cotn));

        real K = (cotp + cotn);
        u += K * (uj - ui);
      });
      L[i] = u;

      i++;
    }

    return L;
  }

  std::vector<TYPE> multM(const std::vector<TYPE> &U) {
    std::vector<TYPE> L(U.size());
    int i = 0;
    // real aAvg = 0;

    for (auto v : _surf->get_vertices()) {
      TYPE ui = U[i];
      real area = 0.0;
      m2::for_each_vertex<SPACE>(v, [&area](face_vertex_ptr fv) {
        real aj = m2::ci::area<SPACE>(fv->face());
        real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
        //real l = 0.33;
        area += l * aj;
      });
      area = max(area, 1e-6);
      L[i] = area * ui;
      i++;
    }
    return L;
  }

  

private:
  m2::surf<SPACE> *_surf;
};

template <typename SPACE, typename TYPE> class area_laplacian {
  M2_TYPEDEFS;

public:
  area_laplacian(m2::surf<SPACE> *surf) { _surf = surf; }

  ~area_laplacian() {}
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

      L[i] = u / area;

      i++;
    }

    return L;
  }

private:
  m2::surf<SPACE> *_surf;
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
