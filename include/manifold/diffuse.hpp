/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __DIFFUSE_MAT__
#define __DIFFUSE_MAT__

#include <cassert>
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>

#include "bontecou/laplacian.hpp"
#include "conj_grad.hpp"
#include "manifold/asawa/m2.hpp"

namespace asawa {
template <typename SPACE> class diffuse {
  M2_TYPEDEFS;

public:
  diffuse(asawa::surf<SPACE> *surf) {
    _surf = surf;
    _M = new asawa::laplacian<SPACE, real>(surf);
  }

  ~diffuse() {}

  std::vector<typename SPACE::real>
  first_order(const std::vector<typename SPACE::real> &u,
              const std::vector<typename SPACE::real> &f,
              typename SPACE::real dt, typename SPACE::real C) {

    using namespace asawa;
    M2_TYPEDEFS;
    asawa::surf<SPACE> *surf = this->_surf;

    asawa::laplacian<SPACE, real> &M = *_M;

    auto diffMult = [&M, dt, C, surf](const std::vector<real> &X) {
      std::vector<real> MX = M.multM(X);
      std::vector<real> CX = M.multC(X);

      int i = 0;
      for (int i = 0; i < MX.size(); i++) {
        MX[i] = MX[i] - dt * C * CX[i];
      }
      return MX;
    };

    std::vector<typename SPACE::real> u_f(u);

    for (int i = 0; i < u_f.size(); i++) {
      u_f[i] = u[i] + dt * f[i];
    }

    std::vector<typename SPACE::real> au_f = M.multM(u_f);

    int its;
    asawa::gradient_descent<SPACE> solver;
    std::vector<real> x = solver.solveConjugateGradient(au_f, diffMult, its);

    return x;
  }

  std::vector<typename SPACE::real>
  second_order(const std::vector<typename SPACE::real> &u,
               const std::vector<typename SPACE::real> &f,
               typename SPACE::real dt, typename SPACE::real C) {

    using namespace asawa;
    M2_TYPEDEFS;
    asawa::surf<SPACE> *surf = this->_surf;

    asawa::laplacian<SPACE, real> &M = *_M;

    auto diffMult = [&M, dt, C, surf](const std::vector<real> &X) {
      std::vector<real> MX = M.multM(X);
      std::vector<real> CX = M.multC(X);

      int i = 0;
      for (int i = 0; i < MX.size(); i++) {
        MX[i] = MX[i] - 0.5 * dt * C * CX[i];
      }
      return MX;
    };

    std::vector<typename SPACE::real> u_f(u);

    for (int i = 0; i < u_f.size(); i++) {
      u_f[i] = u[i] + dt * f[i];
    }

    std::vector<typename SPACE::real> au_f = M.multM(u_f);
    std::vector<typename SPACE::real> lu = M.multC(u);

    for (int i = 0; i < u_f.size(); i++) {
      au_f[i] += 0.5 * dt * C * lu[i];
    }

    int its;
    asawa::gradient_descent<SPACE> solver;
    std::vector<real> x = solver.solveConjugateGradient(au_f, diffMult, its);

    return x;
  }
  asawa::surf<SPACE> *_surf;
  asawa::laplacian<SPACE, real> *_M;
};

} // namespace asawa
#endif