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

    /// κ = (C/l) exp(−(|dp|/l)^p). p = 2 matches calc_gaussian(dx, l).
    inline real calc_gaussian(const vec3 &dx, real l, real p) {
      if (!(l > real(1e-16)) || !(p > real(0)))
        return real(0);
      const real dist = dx.norm();
      const real u = std::pow(dist / l, p);
      const real C = real(1.0) / std::sqrt(real(2.0) * real(M_PI));
      return (C / l) * std::exp(-u);
    }

    /// ∂κ/∂dp for calc_gaussian(dx, l, p). At p = 2 this is calc_d_gaussian.
    inline vec3 calc_d_gaussian(const vec3 &dx, real l, real p) {
      if (!(l > real(1e-16)) || !(p > real(0)))
        return vec3::Zero();
      const real dist = dx.norm();
      if (!(dist > real(1e-16)))
        return vec3::Zero();
      const real kappa = calc_gaussian(dx, l, p);
      const real lp = std::pow(l, p);
      const real sp = std::pow(dist, p - real(2));
      return -kappa * (p / lp) * sp * dx;
    }

    /// Hess_dp of κ = (C/l) exp(−|dp|²/l²). Verified:
    ///   python3 scripts/sympy/tp_harmonic_compose.py
    inline mat3 calc_d2_gaussian(const vec3 &dp, real l) {
      const real lp = l * l;
      if (!(lp > real(1e-32)))
        return mat3::Zero();
      const real kappa = calc_gaussian(dp, l);
      const real c0 = real(-2.0) * kappa / lp;
      const real c1 = real(4.0) * kappa / (lp * lp);
      mat3 H = c0 * mat3::Identity();
      H += c1 * (dp * dp.transpose());
      return H;
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

    /// Exact ∇_dp of k = |P dp|^p / (|dp|^{2p}+l0^p), |N|=1:
    ///   ∇k = p k ( Px/f² - 2 |dp|^{2p-2} dp / (|dp|^{2p}+l0^p) )
    /// Commensurate with calc_tangent_point_inverse_radius / hessian of k.
    /// Near contact this ∇ is attractive (k→0 at coincidence).
    inline vec3 calc_tangent_point_radius_grad_exact(const vec3 &dp,
                                                    const vec3 &N,
                                                    const real &l0,
                                                    const real &p) {
      const real fx = dp.norm();
      const mat3 P = N * N.transpose();
      const vec3 Px = P * dp;
      const real f2 = Px.squaredNorm();
      if (!(f2 > real(1.0e-32)) || !(fx > real(1.0e-16)))
        return vec3::Zero();

      const real lp = std::pow(l0, p);
      const real g2p = std::pow(fx, real(2.0) * p);
      const real den = g2p + lp;
      if (!(den > real(1.0e-32)))
        return vec3::Zero();

      const real k = std::pow(f2, real(0.5) * p) / den;
      const real g2pm2 = std::pow(fx, real(2.0) * p - real(2.0));
      return p * k * (Px / f2 - real(2.0) * g2pm2 * dp / den);
    }

    /// Reg(∇E0): differentiate unregularized E0 = |Pdp|^p/|dp|^{2p}, then
    /// regularize denominators with l0²:
    ///   ∇_reg = p k0 ( Px/(|Px|²+l0²) - 2 dp/(|dp|²+l0²) ),  k0 = |Pdp|^p/|dp|^{2p}
    /// (Not ∇ of the l0^p-floor energy k.)
    vec3 calc_tangent_point_radius_grad(const vec3 &dp, const vec3 &N,
                                           const real &l0, const real &p)
    {
      real fx = dp.norm();
      mat3 P = N * N.transpose();
      vec3 Px = P * dp;

      real fPx = Px.norm();
      if (!(fPx > real(1.0e-32)) || !(fx > real(1.0e-16)))
        return vec3::Zero();

      real l2 = pow(l0, 2.0);
      real k0 = pow(fPx, p) / pow(fx, 2.0 * p);

      vec3 dk =
          p * k0 * (Px / (Px.dot(Px) + l2) - 2.0 * dp / (dp.dot(dp) + l2));

      return dk;
    };

    /// Exact Hess_dp of k = (N·dp)^p / (|dp|^{2p}+l0^p) on the N·dp>0, |N|=1
    /// chart. Generated from scripts/sympy/tangent_point_grad.py (--dump-hess).
    /// Caller should pass Ns = sign(N·dp)·N so N·dp > 0.
    /// Use for harmonic leaf (exact-of-E); regularized force uses *_hessian_reg.
    inline mat3 calc_tangent_point_radius_hessian_chart(
        const vec3 &dp, const vec3 &N, const real &l0, const real &p) {
      const real x = dp[0], y = dp[1], z = dp[2];
      const real nx = N[0], ny = N[1], nz = N[2];
      const real x0 = std::pow(nx, 2);
      const real x1 = std::pow(x, 2);
      const real x2 = std::pow(y, 2);
      const real x3 = std::pow(z, 2);
      const real x4 = x1 + x2 + x3;
      if (!(x4 > real(1.0e-32)))
        return mat3::Zero();
      const real x5 = std::pow(x4, p);
      const real x6 = std::pow(l0, p) + x5;
      if (!(x6 > real(1.0e-32)))
        return mat3::Zero();
      const real x7 = 1.0 / x6;
      const real x8 = nx * x;
      const real x9 = ny * y;
      const real x10 = nz * z;
      const real x11 = x10 + x8 + x9;
      if (!(x11 > real(1.0e-16)))
        return mat3::Zero();
      const real x12 = std::pow(x11, -2);
      const real x13 = x12 * x7;
      const real x14 = std::pow(x11, p);
      const real x15 = p * x14;
      const real x16 = std::pow(p, 2);
      const real x17 = 2 * p;
      const real x18 = std::pow(x6, -2);
      const real x19 = x18 * x5;
      const real x20 = x19 / x4;
      const real x21 = x14 * x17 * x20;
      const real x22 = std::pow(x4, -2);
      const real x23 = x14 * x16;
      const real x24 = x22 * x23;
      const real x25 = 4 * x19;
      const real x26 = std::pow(x4, x17);
      const real x27 = std::pow(x6, -3);
      const real x28 = x20 / x11;
      const real x29 = x23 * x28;
      const real x30 = 4 * x29;
      const real x31 = nx * x13 * x15;
      const real x32 = 2 * y;
      const real x33 = nx * x29;
      const real x34 = ny * x23;
      const real x35 = 2 * x * x28;
      const real x36 =
          nx * ny * x12 * x14 * x16 * x7 - ny * x31 +
          4 * p * x * x14 * x18 * x22 * x5 * y +
          8 * x * x14 * x16 * x22 * x26 * x27 * y - x * x24 * x25 * y -
          x32 * x33 - x34 * x35;
      const real x37 = 2 * z;
      const real x38 = nz * x23;
      const real x39 =
          nx * nz * x12 * x14 * x16 * x7 - nz * x31 +
          4 * p * x * x14 * x18 * x22 * x5 * z +
          8 * x * x14 * x16 * x22 * x26 * x27 * z - x * x24 * x25 * z -
          x33 * x37 - x35 * x38;
      const real x40 = std::pow(ny, 2);
      const real x41 =
          ny * nz * x12 * x14 * x16 * x7 - ny * nz * x13 * x15 +
          4 * p * x14 * x18 * x22 * x5 * y * z +
          8 * x14 * x16 * x22 * x26 * x27 * y * z - x24 * x25 * y * z -
          x28 * x32 * x38 - x28 * x34 * x37;
      const real x42 = std::pow(nz, 2);
      mat3 H;
      H(0, 0) = 4 * p * x1 * x14 * x18 * x22 * x5 + x0 * x12 * x14 * x16 * x7 -
                x0 * x13 * x15 + 8 * x1 * x14 * x16 * x22 * x26 * x27 -
                x1 * x24 * x25 - x21 - x30 * x8;
      H(0, 1) = x36;
      H(0, 2) = x39;
      H(1, 0) = x36;
      H(1, 1) = 4 * p * x14 * x18 * x2 * x22 * x5 +
                x12 * x14 * x16 * x40 * x7 - x13 * x15 * x40 +
                8 * x14 * x16 * x2 * x22 * x26 * x27 - x2 * x24 * x25 - x21 -
                x30 * x9;
      H(1, 2) = x41;
      H(2, 0) = x39;
      H(2, 1) = x41;
      H(2, 2) = 4 * p * x14 * x18 * x22 * x3 * x5 - x10 * x30 +
                x12 * x14 * x16 * x42 * x7 - x13 * x15 * x42 +
                8 * x14 * x16 * x22 * x26 * x27 * x3 - x21 - x24 * x25 * x3;
      if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) ||
          !std::isfinite(H(2, 2)))
        return mat3::Zero();
      return H;
    }

    /// Exact Hess_dp of regularized TP density k (N fixed). Flips to N·dp>0 chart.
    inline mat3 calc_tangent_point_radius_hessian(const vec3 &dp, const vec3 &N,
                                                  const real &l0,
                                                  const real &p) {
      const real fs = N.dot(dp);
      if (!(std::abs(fs) > real(1.0e-16)) || !(dp.squaredNorm() > real(1.0e-32)))
        return mat3::Zero();
      const vec3 Ns = (fs >= real(0.0)) ? N : vec3(-N);
      return calc_tangent_point_radius_hessian_chart(dp, Ns, l0, p);
    }

    /// Jacobian of Reg(∇E0) on N·dp>0 chart (commensurate with
    /// calc_tangent_point_radius_grad). Pass Ns = sign(N·dp)·N.
    inline mat3 calc_tangent_point_radius_hessian_reg_chart(
        const vec3 &dp, const vec3 &N, const real &l0, const real &p) {
      const real x = dp[0], y = dp[1], z = dp[2];
      const real nx = N[0], ny = N[1], nz = N[2];
      const real x0 = std::pow(l0, 2);
      const real x1 = std::pow(x, 2);
      const real x2 = std::pow(y, 2);
      const real x3 = std::pow(z, 2);
      const real x4 = x1 + x2 + x3;
      if (!(x4 > real(1.0e-32)))
        return mat3::Zero();
      const real x5 = x0 + x4;
      const real x6 = 2 / x5;
      const real x7 = nx * x + ny * y + nz * z;
      if (!(x7 > real(1.0e-16)))
        return mat3::Zero();
      const real x8 = std::pow(x7, 2);
      const real x9 = x0 + x8;
      const real x10 = 1.0 / x9;
      const real x11 = nx * x10;
      const real x12 = -x * x6 + x11 * x7;
      const real x13 = std::pow(x4, -p) * std::pow(x7, p);
      const real x14 = std::pow(p, 2) * x13;
      const real x15 = x14 / x7;
      const real x16 = x12 * x15;
      const real x17 = 2 * x14 / x4;
      const real x18 = x12 * x17;
      const real x19 = -x6;
      const real x20 = 4 / std::pow(x5, 2);
      const real x21 = std::pow(nx, 2);
      const real x22 = 2 * x8 / std::pow(x9, 2);
      const real x23 = p * x13;
      const real x24 = x * x20;
      const real x25 = nx * x22;
      const real x26 = x23 * (ny * x11 - ny * x25 + x24 * y);
      const real x27 = x23 * (nz * x11 - nz * x25 + x24 * z);
      const real x28 = ny * x10;
      const real x29 = x28 * x7 - x6 * y;
      const real x30 = x15 * x29;
      const real x31 = x17 * x29;
      const real x32 = std::pow(ny, 2);
      const real x33 = x23 * (-ny * nz * x22 + nz * x28 + x20 * y * z);
      const real x34 = nz * x10 * x7 - x6 * z;
      const real x35 = x15 * x34;
      const real x36 = x17 * x34;
      const real x37 = std::pow(nz, 2);
      mat3 H;
      H(0, 0) = nx * x16 - x * x18 +
                x23 * (x1 * x20 + x10 * x21 + x19 - x21 * x22);
      H(0, 1) = ny * x16 - x18 * y + x26;
      H(0, 2) = nz * x16 - x18 * z + x27;
      H(1, 0) = nx * x30 - x * x31 + x26;
      H(1, 1) = ny * x30 +
                x23 * (x10 * x32 + x19 + x2 * x20 - x22 * x32) - x31 * y;
      H(1, 2) = nz * x30 - x31 * z + x33;
      H(2, 0) = nx * x35 - x * x36 + x27;
      H(2, 1) = ny * x35 + x33 - x36 * y;
      H(2, 2) = nz * x35 +
                x23 * (x10 * x37 + x19 + x20 * x3 - x22 * x37) - x36 * z;
      if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) ||
          !std::isfinite(H(2, 2)))
        return mat3::Zero();
      return H;
    }

    /// Reg(H) for regularized TP force path (Jac of Reg(∇E0)).
    inline mat3 calc_tangent_point_radius_hessian_reg(const vec3 &dp,
                                                      const vec3 &N,
                                                      const real &l0,
                                                      const real &p) {
      const real fs = N.dot(dp);
      if (!(std::abs(fs) > real(1.0e-16)) ||
          !(dp.squaredNorm() > real(1.0e-32)))
        return mat3::Zero();
      const vec3 Ns = (fs >= real(0.0)) ? N : vec3(-N);
      return calc_tangent_point_radius_hessian_reg_chart(dp, Ns, l0, p);
    }

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

    /// Exact ∇_dp of Cauchy K = |N·dp|^p / (|dp|²+ε²)^p (N fixed).
    /// Closed form matches sympy on N·dp≠0; ε already regularizes the kernel.
    inline vec3 calc_tangent_point_radius_grad_cauchy(const vec3 &dp,
                                                     const vec3 &N,
                                                     const real &eps,
                                                     const real &p) {
      const real ndp = N.dot(dp);
      const vec3 Px = ndp * N;
      const real f2 = Px.squaredNorm();
      if (f2 < real(1.0e-32))
        return vec3::Zero();

      const real R2eps = dp.squaredNorm() + eps * eps;
      const real k = std::pow(f2, real(0.5) * p) / std::pow(R2eps, p);
      return p * k * (Px / f2 - real(2.0) * dp / R2eps);
    }

    /// Exact Hess_dp of Cauchy K=(N·dp)^p / (|dp|²+ε²)^p on N·dp>0, |N|=1
    /// chart. Generated from scripts/sympy/tangent_point_cauchy.py --dump-hess.
    /// Caller should pass Ns = sign(N·dp)·N so N·dp > 0.
    inline mat3 calc_tangent_point_radius_hessian_cauchy_chart(
        const vec3 &dp, const vec3 &N, const real &eps, const real &p) {
      const real x = dp[0], y = dp[1], z = dp[2];
      const real nx = N[0], ny = N[1], nz = N[2];
      const real x0 = nx * x;
      const real x1 = ny * y;
      const real x2 = nz * z;
      const real x3 = x0 + x1 + x2;
      if (!(x3 > real(1.0e-16)))
        return mat3::Zero();
      const real x4 = std::pow(x, 2);
      const real x5 = std::pow(y, 2);
      const real x6 = std::pow(z, 2);
      const real x7 = std::pow(eps, 2) + x4 + x5 + x6;
      if (!(x7 > real(1.0e-32)))
        return mat3::Zero();
      const real x8 = std::pow(x3, p) * std::pow(x7, -p);
      const real x9 = p * x8;
      const real x10 = 1.0 / x7;
      const real x11 = 2 * x10;
      const real x12 = -x11 * x9;
      const real x13 = std::pow(x3, -2);
      const real x14 = std::pow(nx, 2) * x13;
      const real x15 = std::pow(p, 2) * x8;
      const real x16 = 4 / std::pow(x7, 2);
      const real x17 = x16 * x4;
      const real x18 = 1.0 / x3;
      const real x19 = x15 * x18;
      const real x20 = 4 * x10 * x19;
      const real x21 = nx * x13;
      const real x22 = x21 * x9;
      const real x23 = ny * x15;
      const real x24 = x * x16;
      const real x25 = x24 * y;
      const real x26 = x11 * y;
      const real x27 = nx * x19;
      const real x28 = x * x11 * x18;
      const real x29 =
          -ny * x22 + x15 * x25 + x21 * x23 - x23 * x28 + x25 * x9 - x26 * x27;
      const real x30 = nz * x15;
      const real x31 = x24 * z;
      const real x32 = x11 * z;
      const real x33 =
          -nz * x22 + x15 * x31 + x21 * x30 - x27 * x32 - x28 * x30 + x31 * x9;
      const real x34 = std::pow(ny, 2) * x13;
      const real x35 = x16 * x5;
      const real x36 = nz * x13;
      const real x37 = x16 * y * z;
      const real x38 = -ny * x36 * x9 + x15 * x37 - x18 * x23 * x32 -
                      x18 * x26 * x30 + x23 * x36 + x37 * x9;
      const real x39 = std::pow(nz, 2) * x13;
      const real x40 = x16 * x6;
      mat3 H;
      H(0, 0) = -x0 * x20 + x12 + x14 * x15 - x14 * x9 + x15 * x17 + x17 * x9;
      H(0, 1) = x29;
      H(0, 2) = x33;
      H(1, 0) = x29;
      H(1, 1) = -x1 * x20 + x12 + x15 * x34 + x15 * x35 - x34 * x9 + x35 * x9;
      H(1, 2) = x38;
      H(2, 0) = x33;
      H(2, 1) = x38;
      H(2, 2) = x12 + x15 * x39 + x15 * x40 - x2 * x20 - x39 * x9 + x40 * x9;
      if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) ||
          !std::isfinite(H(2, 2)))
        return mat3::Zero();
      return H;
    }

    inline mat3 calc_tangent_point_radius_hessian_cauchy(const vec3 &dp,
                                                        const vec3 &N,
                                                        const real &eps,
                                                        const real &p) {
      const real fs = N.dot(dp);
      if (!(std::abs(fs) > real(1.0e-16)) ||
          !(dp.squaredNorm() > real(1.0e-32)))
        return mat3::Zero();
      const vec3 Ns = (fs >= real(0.0)) ? N : vec3(-N);
      return calc_tangent_point_radius_hessian_cauchy_chart(dp, Ns, eps, p);
    }

    /// Orient (Ni, Nj) onto the geometric circle chart:
    ///   Ni·dp > 0, Nj·(-dp) > 0. Returns false if either projection vanishes.
    inline bool circle_agree_orient(const vec3 &dp, const vec3 &Ni,
                                    const vec3 &Nj, vec3 &Nsi, vec3 &Nsj,
                                    real &ai, real &aj) {
      ai = Ni.dot(dp);
      aj = Nj.dot(-dp);
      if (!(std::abs(ai) > real(1.0e-16)) || !(std::abs(aj) > real(1.0e-16)) ||
          !(dp.squaredNorm() > real(1.0e-32)))
        return false;
      Nsi = (ai >= real(0.0)) ? Ni : vec3(-Ni);
      Nsj = (aj >= real(0.0)) ? Nj : vec3(-Nj);
      ai = Nsi.dot(dp);
      aj = Nsj.dot(-dp);
      return (ai > real(1.0e-16)) && (aj > real(1.0e-16));
    }

    /// Scale-free same-circle mismatch:
    ///   Ri=|dp|²/(2 Ni·dp), Rj=|dp|²/(2 Nj·(-dp)),
    ///   δc=-dp+Ri Ni-Rj Nj,  EΔ=½|δc|²/(|dp|²+ε²).
    inline real calc_circle_agree_energy(const vec3 &dp, const vec3 &Ni,
                                         const vec3 &Nj, const real &eps) {
      vec3 Nsi, Nsj;
      real ai, aj;
      if (!circle_agree_orient(dp, Ni, Nj, Nsi, Nsj, ai, aj))
        return real(0.0);
      const real r2 = dp.squaredNorm();
      const real Ri = real(0.5) * r2 / ai;
      const real Rj = real(0.5) * r2 / aj;
      const vec3 dc = -dp + Ri * Nsi - Rj * Nsj;
      const real D = r2 + eps * eps;
      return real(0.5) * dc.squaredNorm() / D;
    }

    /// Exact ∇_dp EΔ with Ni,Nj held fixed (oriented chart). From
    /// scripts/sympy/tangent_point_cauchy_agree.py.
    inline vec3 calc_circle_agree_grad(const vec3 &dp, const vec3 &Ni,
                                       const vec3 &Nj, const real &eps) {
      vec3 Nsi, Nsj;
      real ai, aj;
      if (!circle_agree_orient(dp, Ni, Nj, Nsi, Nsj, ai, aj))
        return vec3::Zero();
      const real r2 = dp.squaredNorm();
      const real Ri = real(0.5) * r2 / ai;
      const real Rj = real(0.5) * r2 / aj;
      const vec3 dc = -dp + Ri * Nsi - Rj * Nsj;
      const real D = r2 + eps * eps;
      const vec3 gRi = dp / ai - Ri * Nsi / ai;
      const vec3 gRj = dp / aj + Rj * Nsj / aj;
      const vec3 JTdc =
          -dc + Nsi.dot(dc) * gRi - Nsj.dot(dc) * gRj;
      return JTdc / D - dc.squaredNorm() * dp / (D * D);
    }

    // clang-format off
inline mat3 calc_circle_agree_hessian_chart(
    const vec3 &dp, const vec3 &Ni, const vec3 &Nj, const real &eps) {
  // generated: python3 scripts/sympy/tangent_point_cauchy_agree.py --dump-hess
  // Chart: Ni·dp > 0, Nj·(-dp) > 0.
  const real x = dp[0], y = dp[1], z = dp[2];
  const real nix = Ni[0], niy = Ni[1], niz = Ni[2];
  const real njx = Nj[0], njy = Nj[1], njz = Nj[2];
  const real x0 = std::pow(x, 2);
  const real x1 = std::pow(y, 2);
  const real x2 = std::pow(z, 2);
  const real x3 = x0 + x1 + x2;
  const real x4 = std::pow(eps, 2) + x3;
  const real x5 = std::pow(x4, -2);
  const real x6 = 2*x;
  const real x7 = njx*x6;
  const real x8 = 2*y;
  const real x9 = njy*x8;
  const real x10 = 2*z;
  const real x11 = njz*x10;
  const real x12 = -x11 - x7 - x9;
  const real x13 = 1.0/x12;
  const real x14 = njx*x13;
  const real x15 = nix*x6 + niy*x8 + niz*x10;
  const real x16 = 1.0/x15;
  const real x17 = nix*x16*x3 - x - x14*x3;
  const real x18 = njy*x13;
  const real x19 = niy*x16*x3 - x18*x3 - y;
  const real x20 = njz*x13;
  const real x21 = niz*x16*x3 - x20*x3 - z;
  const real x22 = (1.0/2.0)*std::pow(x17, 2) + (1.0/2.0)*std::pow(x19, 2) + (1.0/2.0)*std::pow(x21, 2);
  const real x23 = -2*x22*x5;
  const real x24 = 8*x22/std::pow(x4, 3);
  const real x25 = 4*x18;
  const real x26 = 4*x3;
  const real x27 = std::pow(x15, -2);
  const real x28 = nix*x27;
  const real x29 = niy*x28;
  const real x30 = std::pow(x12, -2);
  const real x31 = njx*x30;
  const real x32 = njy*x31;
  const real x33 = x26*x29 + x26*x32;
  const real x34 = 4*niy*x*x16 - x*x25 - x33;
  const real x35 = (1.0/2.0)*x19;
  const real x36 = 4*x20;
  const real x37 = niz*x28;
  const real x38 = njz*x31;
  const real x39 = x26*x37 + x26*x38;
  const real x40 = 4*niz*x*x16 - x*x36 - x39;
  const real x41 = (1.0/2.0)*x21;
  const real x42 = 4*x14;
  const real x43 = std::pow(nix, 2);
  const real x44 = x27*x43;
  const real x45 = std::pow(njx, 2);
  const real x46 = x30*x45;
  const real x47 = 4*nix*x*x16 - x*x42 - x26*x44 - x26*x46 - 2;
  const real x48 = (1.0/2.0)*x17;
  const real x49 = x5*(x34*x35 + x40*x41 + x47*x48);
  const real x50 = 1.0/x4;
  const real x51 = 2*x3;
  const real x52 = x29*x51 + x32*x51;
  const real x53 = (1.0/2.0)*x34;
  const real x54 = x37*x51 + x38*x51;
  const real x55 = (1.0/2.0)*x40;
  const real x56 = std::pow(x15, -3);
  const real x57 = 16*x;
  const real x58 = 16*x3/std::pow(x12, 3);
  const real x59 = -4*nix*x16 + x42;
  const real x60 = x28*x57;
  const real x61 = x31*x57;
  const real x62 = x45*x58;
  const real x63 = -16*niy*x3*x43*x56 + njy*x62;
  const real x64 = -4*niy*x16 + x25;
  const real x65 = -16*niz*x3*x43*x56 + njz*x62;
  const real x66 = -4*niz*x16 + x36;
  const real x67 = (1.0/2.0)*x47;
  const real x68 = x*x24;
  const real x69 = 4*nix*x16*y - x33 - x42*y;
  const real x70 = niy*niz*x27;
  const real x71 = njy*njz;
  const real x72 = x30*x71;
  const real x73 = x26*x70 + x26*x72;
  const real x74 = 4*niz*x16*y - x36*y - x73;
  const real x75 = std::pow(niy, 2);
  const real x76 = x27*x75;
  const real x77 = std::pow(njy, 2);
  const real x78 = x30*x77;
  const real x79 = 4*niy*x16*y - x25*y - x26*x76 - x26*x78 - 2;
  const real x80 = x35*x79 + x41*x74 + x48*x69;
  const real x81 = x5*x6;
  const real x82 = x51*x70 + x51*x72;
  const real x83 = 2*niz*x16*y - x20*x8 - x82;
  const real x84 = 2*niy*x16*y - x13*x9 - x51*x76 - x51*x78 - 1;
  const real x85 = 2*nix*x16*y - x14*x8 - x52;
  const real x86 = 8*x;
  const real x87 = -16*nix*niy*niz*x3*x56;
  const real x88 = njx*x58;
  const real x89 = x71*x88;
  const real x90 = x70*x86 + x72*x86 + x87 + x89;
  const real x91 = 8*y;
  const real x92 = x37*x91 + x38*x91;
  const real x93 = -16*nix*x3*x56*x75 + x77*x88;
  const real x94 = -x49*x8 + x50*(x35*(-x29*x91 - x32*x91 - x76*x86 - x78*x86 - x93) + x41*(-x90 - x92) + x48*(-x29*x86 - x32*x86 - x44*x91 - x46*x91 - x63) + x53*x84 + x55*x83 + x67*x85) + x68*y - x80*x81;
  const real x95 = 4*nix*x16*z - x39 - x42*z;
  const real x96 = 4*niy*x16*z - x25*z - x73;
  const real x97 = std::pow(niz, 2);
  const real x98 = x27*x97;
  const real x99 = std::pow(njz, 2);
  const real x100 = x30*x99;
  const real x101 = 4*niz*x16*z - x100*x26 - x26*x98 - x36*z - 2;
  const real x102 = x101*x41 + x35*x96 + x48*x95;
  const real x103 = 2*niy*x16*z - x10*x18 - x82;
  const real x104 = 2*niz*x16*z - x100*x51 - x11*x13 - x51*x98 - 1;
  const real x105 = 2*nix*x16*z - x10*x14 - x54;
  const real x106 = 8*z;
  const real x107 = x106*x29 + x106*x32;
  const real x108 = -16*nix*x3*x56*x97 + x88*x99;
  const real x109 = -x10*x49 - x102*x81 + x50*(x103*x53 + x104*x55 + x105*x67 + x35*(-x107 - x90) + x41*(-x100*x86 - x106*x37 - x106*x38 - x108 - x86*x98) + x48*(-x106*x44 - x106*x46 - x37*x86 - x38*x86 - x65)) + x68*z;
  const real x110 = x5*x80;
  const real x111 = (1.0/2.0)*x69;
  const real x112 = (1.0/2.0)*x74;
  const real x113 = 16*y;
  const real x114 = -16*niz*x3*x56*x75 + njz*x58*x77;
  const real x115 = (1.0/2.0)*x79;
  const real x116 = x102*x5;
  const real x117 = -16*niy*x3*x56*x97 + njy*x58*x99;
  const real x118 = -x10*x110 - x116*x8 + x24*y*z + x50*(x103*x115 + x104*x112 + x105*x111 + x35*(-x106*x76 - x106*x78 - x114 - x70*x91 - x72*x91) + x41*(-x100*x91 - x106*x70 - x106*x72 - x117 - x91*x98) + x48*(-x107 - x87 - x89 - x92));
  const real x119 = 16*z;
  mat3 H;
  H(0, 0) = -4*x*x49 + x0*x24 + x23 + x50*(x35*(-niy*x60 - njy*x61 - x63 - x64) + x41*(-niz*x60 - njz*x61 - x65 - x66) + x48*(16*std::pow(nix, 3)*x3*x56 - std::pow(njx, 3)*x58 - x44*x57 - x46*x57 - x59) + x53*(2*niy*x*x16 - x18*x6 - x52) + x55*(2*niz*x*x16 - x20*x6 - x54) + x67*(2*nix*x*x16 - x13*x7 - x44*x51 - x46*x51 - 1));
  H(0, 1) = x94;
  H(0, 2) = x109;
  H(1, 0) = x94;
  H(1, 1) = x1*x24 - 4*x110*y + x23 + x50*(x111*x85 + x112*x83 + x115*x84 + x35*(16*std::pow(niy, 3)*x3*x56 - std::pow(njy, 3)*x58 - x113*x76 - x113*x78 - x64) + x41*(-x113*x70 - x113*x72 - x114 - x66) + x48*(-x113*x29 - x113*x32 - x59 - x93));
  H(1, 2) = x118;
  H(2, 0) = x109;
  H(2, 1) = x118;
  H(2, 2) = -4*x116*z + x2*x24 + x23 + x50*((1.0/2.0)*x101*x104 + (1.0/2.0)*x103*x96 + (1.0/2.0)*x105*x95 + x35*(-x117 - x119*x70 - x119*x72 - x64) + x41*(16*std::pow(niz, 3)*x3*x56 - std::pow(njz, 3)*x58 - x100*x119 - x119*x98 - x66) + x48*(-x108 - x119*x37 - x119*x38 - x59));
  if (!std::isfinite(H(0, 0)) || !std::isfinite(H(1, 1)) || !std::isfinite(H(2, 2)))
    return mat3::Zero();
  return H;
}

    inline mat3 calc_circle_agree_hessian(const vec3 &dp, const vec3 &Ni,
                                          const vec3 &Nj, const real &eps) {
      vec3 Nsi, Nsj;
      real ai, aj;
      if (!circle_agree_orient(dp, Ni, Nj, Nsi, Nsj, ai, aj))
        return mat3::Zero();
      return calc_circle_agree_hessian_chart(dp, Nsi, Nsj, eps);
    }
    // clang-format on

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
