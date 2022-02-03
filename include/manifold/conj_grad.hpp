/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __CONJ_GRAD__
#define __CONJ_GRAD__

#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>

namespace m2 {
template <typename SPACE> class gradient_descent {
  M2_TYPEDEFS;

public:
  gradient_descent() {}
  ~gradient_descent() {}

  typename SPACE::real dot(const std::vector<typename SPACE::real> &A,
                           const std::vector<typename SPACE::real> &B) {

    real v = 0.0;
    for (int i = 0; i < A.size(); i++) {
      v += A[i] * B[i];
    }
    return v;
  }

  std::vector<typename SPACE::real>
  addScaledVec(const std::vector<typename SPACE::real> &X0,
               const std::vector<typename SPACE::real> &V,
               typename SPACE::real C) {

    std::vector<real> X(X0);
    for (int i = 0; i < X.size(); i++) {
      X[i] = X0[i] + C * V[i];
    }
    return X;
  }

  std::vector<typename SPACE::real>
  solveConjugateGradient(const std::vector<typename SPACE::real> &B,
                         function<std::vector<real>(const std::vector<real> &)> M,
                         int & its) {

    /*
     x =          initial guess for solution of Ax=b
     r = b - Ax = residual, to be made small
     p = r      = initial "search direction"

     do while ( new_r , new_r ) not small
        v = Ap                        ... matrix-vector multiply
        a = ( r , r ) / ( p , v )     ... dot product
        x = x + a*p                   ... updated approximate solution
        r_new = r - a*v               ... update the residual
        g = ( r_new , r_new ) / ( r , r )
        p = r_new + g*p               ... update search direction
        r = r_new
     enddo
    */

    // x =          initial guess for solution of Ax=b
    std::vector<real> X = B;

    // std::vector<real> MU = diffMult(surf, U0, dt, C);
    std::vector<real> MX = M(X);

    // r = b - Ax = residual, to be made small
    std::vector<real> R = addScaledVec(B, MX, -1);
    // p = r = initial "search direction"
    std::vector<real> P = R;
    real sig1 = dot(R, R);
    int ii = 0;
    while (ii < 512 && sig1 > 1e-16) {

      // v = Ap                        ... matrix-vector multiply
      std::vector<real> V = M(P);

      // a = ( r , r ) / ( p , v )     ... dot product
      real a = sig1 / dot(P, V);
      // std::cout << sig1 << std::endl;
      // x = x + a*p                   ... updated approximate solution
      X = addScaledVec(X, P, a);
      // r_new = r - a * v... update the residual
      std::vector<real> Rn = addScaledVec(R, V, -a);
      // g = ( r_new , r_new ) / ( r , r )
      real sig0 = sig1;
      sig1 = dot(Rn, Rn);
      real g = sig1 / sig0;

      // p = r_new + g *p... update search direction
      P = addScaledVec(Rn, P, g);

      // r = r_new
      R = Rn;
      ii++;
    }
    std::cout << "computed grad to: " << sig1 << " in " << ii << " iterations"
              << std::endl;
    its = ii;
    return X;
  }
}; // gradient_descent
} // namespace m2
#endif
