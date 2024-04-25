//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __CALDER_WEIGHT_FUNCTIONS__
#define __CALDER_WEIGHT_FUNCTIONS___

#include <gaudi/common.h>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi
{

  namespace calder
  {
    real calc_mollified(vec3 dp, real C, real p) {
      // mollified kernel
      real dist = dp.norm();
      real distp = pow(dist, p);
      real lp = pow(C, p);

      real kappa = (1.0 - exp(-distp / lp)) / distp;
      return kappa;
    };

    real calc_cauchy(vec3 dp, real C, real p) {
      // mollified kernel
      real dist = dp.norm();
      real kappa = pow(1.0 + pow(dist / C, 2.0), 2.0);
      return kappa;
    };

    real calc_inv_dist(vec3 dx, real eps, real p)
    {
      // std laplace kernel
      real dist = dx.norm();
      real distp = pow(dist, p);
      real kappa = 1.0 / (distp + eps);
      return kappa;
    };

    vec3 calc_d_inv_dist(vec3 dx, real eps, real p)
    {
      real dist = dx.norm();
      real distpm1 = pow(dist, p - 1);
      real distp = pow(dist, p);
      real distp_eps_2 = pow(distp + eps, 2.0);
      return -p * distpm1 / distp_eps_2 * dx / dist;
    };

    // calc w/dw but using gaussian kernel instead of std laplace
    real calc_gaussian(vec3 dx, real l)
    {
      real dist = dx.norm();
      real distp = pow(dist, 2.0);
      real lp = pow(l, 2.0);
      real C = 1.0 / std::sqrt(2.0 * M_PI);
      real expx2 = exp(-distp / lp);
      real kappa = C / l * expx2;
      // real kappa = exp(-distp / eps);
      return kappa;
    };

    vec3 calc_d_gaussian(vec3 dx, real l)
    {
      real dist = dx.norm();
      real distp = pow(dist, 2.0);
      real lp = pow(l, 2.0);
      real C = 1.0 / std::sqrt(2.0 * M_PI);
      real expx2 = exp(-distp / lp);

      return -2.0 * dx * C / l / lp * expx2;
    };

    real calc_tangent_point_radius(const vec3 &dp, const vec3 &N)
    {
      real ndp = dp.squaredNorm();
      real nPdp = (N * N.transpose() * dp).norm();
      return 0.5 * ndp / nPdp;
    };

    real calc_tangent_point_inverse_radius(const vec3 &dp, const vec3 &N,
                                              const real &l0, const double &p)
    {
      // computes the inverse radius of the tangent point to power p
      real ndp = dp.norm();
      real nPdp = (N * N.transpose() * dp).norm();
      real lp = pow(l0, p);
      real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);
      return k;
    };

    vec3 calc_tangent_point_radius_grad_0(const vec3 &dp0, const vec3 &N,
                                             const real &l0, const real &p)
    {
      real Ndp = dp0.dot(N);

      real ndp = dp0.norm();
      vec3 dp = ndp * dp0.normalized();

      mat3 P = N * N.transpose();
      vec3 Pdp = P * dp;

      real nPdp = (Pdp).norm();
      real lp = pow(l0, p);
      real l2 = pow(l0, 2.0);

      real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);
      vec3 dk0 = P * Pdp / (nPdp * nPdp + l2);
      vec3 dk1 = 2.0 * dp / (ndp * ndp + l2);
      vec3 dk = p * k * (dk0 - dk1);
      // vec3 dk = p * k * (dk0);

      return dk;
    };

    vec3 calc_tangent_point_radius_grad_1(const vec3 &dp0, const vec3 &N,
                                             const real &l0, const real &p)
    {
      real Ndp = dp0.dot(N);

      real fx = dp0.norm();
      vec3 dfx = fx * dp0.normalized();

      mat3 P = N * N.transpose();
      vec3 Pdp = P * dp0;

      real fPx = (Pdp).norm();
      vec3 dfPx = fPx * Pdp.normalized();

      real lp = pow(l0, p);

      real denom = (pow(fx, 2.0 * p) + lp);
      real k = pow(fPx, p) / denom;
      vec3 dk0 = P * pow(fPx, p - 1.0) * dfPx / denom;
      vec3 dk1 = 2.0 * k * pow(fx, 2.0 * p - 1.0) * dfx / denom;
      vec3 dk = p * (dk0 - dk1);
      // vec3 dk = p * k * (dk0);

      return -dk;
    };

    vec3 calc_tangent_point_radius_grad(const vec3 &dp, const vec3 &N,
                                           const real &l0, const real &p)
    {

      real fx = dp.norm();
      mat3 P = N * N.transpose();
      vec3 Px = P * dp;

      real fPx = Px.norm();
      // vec3 dfPx = fPx * Pdp.normalized();

      real lp = pow(l0, p);
      real l2 = pow(l0, 2.0);

      real k = pow(fPx, p) / (pow(fx, 2.0 * p) + lp);

      //  P = N*Nt => P * P = N*Nt*N*Nt = N*Nt = P
      vec3 dk = p * k * (Px / (Px.dot(Px) + l2) - 2.0 * dp / (dp.dot(dp) + l2));

      // vec3 dk = p * k * (dk0);

      return dk;
    };

  } // namespace calder
} // namespace gaudi
#endif