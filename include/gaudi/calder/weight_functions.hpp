//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __CALDER_WEIGHT_FUNCTIONS___
#define __CALDER_WEIGHT_FUNCTIONS___

#include <gaudi/common.h>
#include <gaudi/vec_addendum.h>
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

    real calc_inv_dist(vec3 dx, real l0, real p)
    {
      // Inverse-distance kernel with a length-scale regularizer. The
      // denominator has units length^p, so l0 must be raised to p as well.
      real dist = dx.norm();
      real distp = pow(dist, p);
      real lp = pow(l0, p);
      real kappa = 1.0 / (distp + lp);
      return kappa;
    };

    vec3 calc_d_inv_dist(vec3 dx, real l0, real p)
    {
      real dist = dx.norm();
      if (dist < 1e-16)
      {
        return vec3::Zero();
      }
      real distpm1 = pow(dist, p - 1);
      real distp = pow(dist, p);
      real lp = pow(l0, p);
      real denom2 = pow(distp + lp, 2.0);
      return -p * distpm1 / denom2 * dx / dist;
    };

    // PAT-style softmin / screened-Laplace kernel: exp(-β‖dp‖).
    // beta <= 0 → zero weight (disabled).
    real calc_exp_dist(vec3 dp, real beta)
    {
      if (!(beta > 0.0))
      {
        return 0.0;
      }
      const real dist = dp.norm();
      return std::exp(-beta * dist);
    }

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
      return va::tangent_point_radius(dp, N);
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

    /// Rosenhead / Cauchy filament TP density (N held fixed in ∇):
    ///   K = |N·dp|^p / (|dp|² + ε²)^p
    /// With rod N ∥ dp_⊥ this is ρ^p / (ρ² + s² + ε²)^p.
    inline real calc_tangent_point_inverse_radius_cauchy(const vec3 &dp,
                                                        const vec3 &N,
                                                        const real &eps,
                                                        const real &p) {
      const real f = std::abs(N.dot(dp));
      const real R2eps = dp.squaredNorm() + eps * eps;
      return std::pow(f, p) / std::pow(R2eps, p);
    }

    /// Exact ∇_dp of calc_tangent_point_inverse_radius_cauchy (N constant):
    ///   ∇K = p K ( Px / |Px|² - 2 dp / (|dp|² + ε²) ),  Px = (N·dp) N
    inline vec3 calc_tangent_point_radius_grad_cauchy(const vec3 &dp,
                                                     const vec3 &N,
                                                     const real &eps,
                                                     const real &p) {
      const real ndp = N.dot(dp);
      const vec3 Px = ndp * N;
      const real f2 = Px.squaredNorm(); // = ndp² for |N|=1
      if (f2 < real(1.0e-32))
        return vec3::Zero();

      const real R2eps = dp.squaredNorm() + eps * eps;
      const real k = std::pow(f2, real(0.5) * p) / std::pow(R2eps, p);
      return p * k * (Px / f2 - real(2.0) * dp / R2eps);
    }

    /// Stable log-sum-exp softmax: τ log(e^{a/τ}+e^{b/τ}).
    /// Stages: r = tan_point; rs = soft_floor(r); f = density(rs)
    /// (from scripts/sympy/tangent_point_soft.py)
    inline real tp_softmax(real a, real b, real tau) {
      const real M = std::max(a, b);
      return M + tau * std::log(std::exp((a - M) / tau) +
                                std::exp((b - M) / tau));
    }

    /// Geometric TP radius R = |dp|² / (2 |N·dp|).
    inline real tan_point_radius(const vec3 &dp, const vec3 &N) {
      const real f = std::abs(N.dot(dp));
      const real g2 = dp.squaredNorm();
      if (!(f > real(1.0e-16)) || !(g2 > real(1.0e-32)))
        return real(0.0);
      return real(0.5) * g2 / f;
    }

    /// Soft contact floor: rs = softmax_τ(r, R_min) ≈ max(r, R_min).
    inline real soft_floor_radius(real r, real R_min, real tau) {
      return tp_softmax(r, R_min, tau);
    }

    /// Classic density f = (1/rs)^p.
    inline real soft_tp_density(real rs, real p) {
      return std::pow(real(1.0) / rs, p);
    }

    /// ∂softmax(a,b)/∂a (stable).
    inline real soft_floor_dr(real a, real b, real tau) {
      const real inv_tau = real(1.0) / tau;
      const real M = std::max(a, b);
      const real ea = std::exp((a - M) * inv_tau);
      const real eb = std::exp((b - M) * inv_tau);
      return ea / (ea + eb);
    }

    /// ∂²softmax(a,b)/∂a² = w(1-w)/τ.
    inline real soft_floor_d2r(real a, real b, real tau) {
      const real w = soft_floor_dr(a, b, tau);
      return w * (real(1.0) - w) / tau;
    }

    /// ∇_dp R for R = |dp|²/(2|N·dp|) (N fixed).
    inline vec3 tan_point_radius_grad(const vec3 &dp, const vec3 &N) {
      const real f_s = N.dot(dp);
      const real f = std::abs(f_s);
      const real g2 = dp.squaredNorm();
      if (f < real(1.0e-16) || g2 < real(1.0e-32))
        return vec3::Zero();
      const real sgn = (f_s >= real(0.0)) ? real(1.0) : real(-1.0);
      return (real(2.0) * f * dp - g2 * (sgn * N)) / (real(2.0) * f * f);
    }

    /// ∇²_dp R for R = |dp|²/(2 N·dp) on the N·dp > 0 chart:
    ///   H = I/f_s − (dp⊗N + N⊗dp)/f_s² + (|dp|² N⊗N)/f_s³
    inline mat3 tan_point_radius_hess(const vec3 &dp, const vec3 &N) {
      const real f_s = N.dot(dp);
      const real g2 = dp.squaredNorm();
      if (!(std::abs(f_s) > real(1.0e-16)) || !(g2 > real(1.0e-32)))
        return mat3::Zero();
      const vec3 Ns = (f_s >= real(0.0)) ? N : vec3(-N);
      const real fs = std::abs(f_s);
      mat3 H = mat3::Identity() / fs;
      H.noalias() -= (dp * Ns.transpose() + Ns * dp.transpose()) / (fs * fs);
      H.noalias() += (g2 * (Ns * Ns.transpose())) / (fs * fs * fs);
      return H;
    }

    /// Softmax-floor TP density f = density(soft_floor(tan_point)).
    inline real calc_tangent_point_radius_soft(const vec3 &dp, const vec3 &N,
                                              real R_min, real tau, real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return real(0.0);
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return real(0.0);
      const real rs = soft_floor_radius(r, R_min, tau);
      return soft_tp_density(rs, p);
    }

    /// ∇_dp f (N fixed). ∇f = f'(r) ∇r.
    inline vec3 calc_tangent_point_radius_gradient_soft(const vec3 &dp,
                                                       const vec3 &N,
                                                       real R_min, real tau,
                                                       real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return vec3::Zero();
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return vec3::Zero();
      const real rs = soft_floor_radius(r, R_min, tau);
      const real f = soft_tp_density(rs, p);
      const real w = soft_floor_dr(r, R_min, tau);
      const real df_dr = (-p * f / rs) * w;
      const vec3 g = df_dr * tan_point_radius_grad(dp, N);
      if (!std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2]))
        return vec3::Zero();
      return g;
    }

    /// Full Hess_dp f (N fixed): H = f''(r) (∇r)(∇r)ᵀ + f'(r) ∇²r
    inline mat3 calc_tangent_point_radius_hessian_soft(const vec3 &dp,
                                                      const vec3 &N,
                                                      real R_min, real tau,
                                                      real p) {
      if (!(R_min > 0.0) || !(tau > 0.0))
        return mat3::Zero();
      const real r = tan_point_radius(dp, N);
      if (!(r > 0.0))
        return mat3::Zero();
      const real rs = soft_floor_radius(r, R_min, tau);
      const real f = soft_tp_density(rs, p);
      const real w = soft_floor_dr(r, R_min, tau);
      const real w2 = soft_floor_d2r(r, R_min, tau);
      const real df_drs = -p * f / rs;
      const real d2f_drs2 = p * (p + real(1.0)) * f / (rs * rs);
      const real df_dr = df_drs * w;
      const real d2f_dr2 = d2f_drs2 * w * w + df_drs * w2;
      const vec3 gr = tan_point_radius_grad(dp, N);
      const mat3 Hr = tan_point_radius_hess(dp, N);
      mat3 H = d2f_dr2 * (gr * gr.transpose()) + df_dr * Hr;
      if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) ||
          !std::isfinite(H(2, 2)))
        return mat3::Zero();
      return H;
    }

    // PAT-style Mahalanobis length from principal curvatures (mesh-free):
    //   M = (κ_min²⟨dp,e_min⟩² + κ_max²⟨dp,e_max⟩² + ε²⟨dp,n⟩²)^{1/2}
    inline real calc_mahalanobis_curvature(const vec3 &dp, real k_min,
                                           real k_max, const vec3 &e_min,
                                           const vec3 &e_max, const vec3 &n,
                                           real eps = 1e-3) {
      const real s_min = dp.dot(e_min);
      const real s_max = dp.dot(e_max);
      const real s_n = dp.dot(n);
      const real m2 = k_min * k_min * s_min * s_min +
                      k_max * k_max * s_max * s_max +
                      eps * eps * s_n * s_n;
      return std::sqrt(std::max(m2, real(0.0)));
    }

  } // namespace calder
} // namespace gaudi
#endif
