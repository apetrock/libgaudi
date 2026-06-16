
#ifndef ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#define ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "gaudi/common.h"
#include "ncls.hpp"
// stub: but least squares shape functions here
// sphere, cylinder, etc.

namespace gaudi
{
  namespace albers
  {

    TYPEDEF_VEC(14)
    TYPEDEF_MAT(14)
    TYPEDEF_MAT_NM(4, 14)
    // darboux cyclide is a type of implicit surface
    // X = (x*x + y*y + z*z)
    // L = m*x + n*y + k*z
    // Q = A*x*x + B*y*y + C*z*z + 2*D*x*y + 2*E*x*z + 2*F*y*z + 2*G*x + 2*H*y + 2*I*z + J
    // D = lambda * X * X + L * X + Q

    mat414 mk_darboux_A(vec3 dx)
    {
      mat414 A;
      real x = dx[0];
      real y = dx[1];
      real z = dx[2];

      real xx = x * x;
      real yy = y * y;
      real zz = z * z;
      real xy = 2.0 * x * y;
      real xz = 2.0 * x * z;
      real yz = 2.0 * y * z;

      real xyz = xx + yy + zz;
      real xyz2 = xyz * xyz;
      real xxyz = x * xyz;
      real yxyz = y * xyz;
      real zxyz = z * xyz;
      A.col(0) = vec4(xx, 2.0 * x, 0.0 * y, 0.0 * z); // A
      A.col(1) = vec4(yy, 0.0 * x, 2.0 * y, 0.0 * z); // B
      A.col(2) = vec4(zz, 0.0 * x, 0.0 * y, 2.0 * z); // C

      A.col(3) = vec4(xy, 2.0 * y, 2.0 * x, 0.0 * z); // D
      A.col(4) = vec4(xz, 2.0 * z, 0.0 * y, 2.0 * x); // E
      A.col(5) = vec4(yz, 0.0 * x, 2.0 * z, 2.0 * y); // F

      A.col(6) = vec4(2.0 * x, 2.0, 0.0, 0.0); // G
      A.col(7) = vec4(2.0 * y, 0.0, 2.0, 0.0); // H
      A.col(8) = vec4(2.0 * z, 0.0, 0.0, 2.0); // I

      A.col(9) = vec4(1.0, 0.0, 0.0, 0.0); // J

      A.col(10) = vec4(xyz2, 4.0 * xxyz, 4.0 * yxyz, 4.0 * zxyz); // lambda
      A.col(11) = vec4(x * xyz, 2.0 * xx + xyz, xy, xz);          // mu
      A.col(12) = vec4(y * xyz, xy, 2.0 * yy + xyz, yz);          // nu
      A.col(13) = vec4(z * xyz, xz, yz, 2.0 * zz + xyz);          // kappa

      return A;
    }

    void unpack_darboux(const vec14 &Q, real &A, real &B, real &C, real &D, real &E, real &F, real &G, real &H, real &I, real &J, real &lambda, real &mu, real &nu, real &kappa)
    {
      A = Q[0];
      B = Q[1];
      C = Q[2];
      D = Q[3];
      E = Q[4];
      F = Q[5];
      G = Q[6];
      H = Q[7];
      I = Q[8];
      J = Q[9];
      lambda = Q[10];
      mu = Q[11];
      nu = Q[12];
      kappa = Q[13];
    }

    real eval_darboux(const vec14 &Q, const vec3 x_v)
    {
      vec4 p = {x_v[0], x_v[1], x_v[2], 1.0};

      real x = x_v[0], y = x_v[1], z = x_v[2];
      real A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa;

      unpack_darboux(Q, A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa);
      mat4 Qm = mat4::Zero();
      Qm << A, D, E, G,
          D, B, F, H,
          E, F, C, I,
          G, H, I, J;

      real X = x_v.dot(x_v);
      real L = vec3(mu, nu, kappa).dot(x_v);
      real Q_ = p.transpose() * Qm * p;

      real D_ = lambda * X * X + L * X + Q_;
      return D_;
    }

    vec3 darboux_grad(const vec14 &Q, const vec3 x_v)
    {
      real x = x_v[0], y = x_v[1], z = x_v[2];

      real A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa;
      unpack_darboux(Q, A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa);
      real X = x_v.dot(x_v);
      real L = vec3(mu, nu, kappa).dot(x_v);
      real dx = 2.0 * A * x + 2.0 * D * y + 2.0 * E * z + 2.0 * G + 4.0 * lambda * x * X + mu * X + 2.0 * x * L;
      real dy = 2.0 * B * y + 2.0 * D * x + 2.0 * F * z + 2.0 * H + 4.0 * lambda * y * X + nu * X + 2.0 * y * L;
      real dz = 2.0 * C * z + 2.0 * E * x + 2.0 * F * y + 2.0 * I + 4.0 * lambda * z * X + kappa * X + 2.0 * z * L;
      return vec3(dx, dy, dz);
    }

    vec3 darboux_center(const vec14 &Q)
    {
      // the machine spit this out, this is probably wrongo

      real A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa;
      unpack_darboux(Q, A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa);
      vec3 center = vec3(-G / A, -H / B, -I / C);
      return center;
    }

    mat3 darboux_hessian(const vec14 &Q, const vec3 x_v)
    {
      real x = x_v[0], y = x_v[1], z = x_v[2];
      real A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa;
      unpack_darboux(Q, A, B, C, D, E, F, G, H, I, J, lambda, mu, nu, kappa);
      real X = x_v.dot(x_v);
      real L = vec3(mu, nu, kappa).dot(x_v);
      mat3 M = mat3::Zero();
      M(0, 0) = 2.0 * A + 8.0 * lambda * x * x + 4.0 * lambda * X + 6.0 * mu * x + 2.0 * nu * y + 2.0 * kappa * z;
      M(1, 1) = 2.0 * B + 8.0 * lambda * y * y + 4.0 * lambda * X + 2.0 * mu * x + 6.0 * nu * y + 2.0 * kappa * z;
      M(2, 2) = 2.0 * C + 8.0 * lambda * z * z + 4.0 * lambda * X + 2.0 * mu * x + 2.0 * nu * y + 6.0 * kappa * z;

      M(0, 1) = M(1, 0) = 2.0 * D + 8.0 * lambda * x * y + 2.0 * mu * y + 2.0 * nu * x;
      M(0, 2) = M(2, 0) = 2.0 * E + 8.0 * lambda * x * z + 2.0 * mu * z + 2.0 * kappa * x;
      M(1, 2) = M(2, 1) = 2.0 * F + 8.0 * lambda * y * z + 2.0 * nu * z + 2.0 * kappa * y;

      return M;
    }

    enum class darboux_ridge_failure {
      none,
      invalid_direction,
      surface_projection_failed,
      small_denominator,
      nonfinite,
      max_travel,
      not_converged
    };

    inline const char *darboux_ridge_failure_name(darboux_ridge_failure failure) {
      switch (failure) {
      case darboux_ridge_failure::none:
        return "none";
      case darboux_ridge_failure::invalid_direction:
        return "invalid_direction";
      case darboux_ridge_failure::surface_projection_failed:
        return "surface_projection_failed";
      case darboux_ridge_failure::small_denominator:
        return "small_denominator";
      case darboux_ridge_failure::nonfinite:
        return "nonfinite";
      case darboux_ridge_failure::max_travel:
        return "max_travel";
      case darboux_ridge_failure::not_converged:
        return "not_converged";
      }
      return "unknown";
    }

    struct darboux_ridge_estimate {
      vec3 center = vec3::Zero(); // local medial-axis point
      real travel = 0.0;          // |center - start|
      real residual = 0.0;        // |grad D(center) · dir|
      real surface_residual = 0.0;
      int iterations = 0;
      int clamped_steps = 0;
      darboux_ridge_failure failure = darboux_ridge_failure::not_converged;
      darboux_ridge_failure projection_failure = darboux_ridge_failure::none;
      bool projection_attempted = false;
      bool projection_converged = false;
      bool accepted = false;
    };

    inline real darboux_directional_deriv(const vec14 &Q, const vec3 &x,
                                          const vec3 &dir) {
      return darboux_grad(Q, x).dot(dir);
    }

    inline real darboux_directional_second_deriv(const vec14 &Q, const vec3 &x,
                                                 const vec3 &dir) {
      return dir.transpose() * darboux_hessian(Q, x) * dir;
    }

    // Unconstrained Newton-Raphson onto the fitted zero level set D(x) = 0.
    //
    // We solve D(x) = 0 by minimizing f(x) = D(x)^2 with full Newton steps
    // (no projection onto an input direction, no KKT/closest-point objective):
    //   grad f = 2 D grad D
    //   hess f = 2 (grad D grad D^T + D H_D)
    //   x <- x - (grad D grad D^T + D H_D)^{-1} (D grad D)
    // The grad D grad D^T term alone is the Gauss-Newton approximation; adding
    // the D H_D curvature term makes this the true Newton solve. Near D=0 this
    // Hessian can be singular, so we fall back to the minimum-norm scalar root
    // correction in the gradient direction.
    inline bool darboux_surface_point_newton(
        const vec14 &Q, const vec3 &x0, int max_iters, real tol, vec3 *x_out,
        int *iters_out = nullptr, std::vector<vec3> *trace = nullptr,
        real max_travel = std::numeric_limits<real>::infinity(),
        const vec3 *travel_origin = nullptr,
        darboux_ridge_failure *failure_out = nullptr) {
      vec3 x = x0;
      const vec3 origin = travel_origin != nullptr ? *travel_origin : x0;
      if (failure_out != nullptr) {
        *failure_out = darboux_ridge_failure::not_converged;
      }

      for (int iter = 0; iter < std::max(1, max_iters); ++iter) {
        const real phi = eval_darboux(Q, x);
        const real surface_tol = std::max(tol, std::sqrt(tol));
        if (!std::isfinite(phi) || !x.allFinite()) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::nonfinite;
          }
          return false;
        }
        if (std::abs(phi) <= surface_tol) {
          if (x_out != nullptr) {
            *x_out = x;
          }
          if (iters_out != nullptr) {
            *iters_out = iter + 1;
          }
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::none;
          }
          return true;
        }

        const vec3 g = darboux_grad(Q, x);
        const mat3 Hd = darboux_hessian(Q, x);
        if (!g.allFinite() || !Hd.allFinite()) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::nonfinite;
          }
          return false;
        }

        const mat3 Hf = g * g.transpose() + phi * Hd;
        const vec3 gf = phi * g;
        Eigen::SelfAdjointEigenSolver<mat3> eig(Hf);
        if (eig.info() != Eigen::Success ||
            !eig.eigenvalues().allFinite()) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::small_denominator;
          }
          return false;
        }

        const vec3 evals = eig.eigenvalues().cwiseAbs();
        const real min_eval = evals.minCoeff();
        const real max_eval = evals.maxCoeff();
        const bool ill_conditioned =
            max_eval <= 0.0 || min_eval <= real(1e-6) * max_eval;

        auto gradient_root_step = [&]() -> vec3 {
          const real g2 = g.squaredNorm();
          if (g2 < 1e-24) {
            return vec3::Constant(std::numeric_limits<real>::quiet_NaN());
          }
          return -(phi / g2) * g;
        };

        vec3 dx = vec3::Zero();
        if (!ill_conditioned) {
          Eigen::ColPivHouseholderQR<mat3> qr(Hf);
          dx = qr.solve(-gf);
        } else {
          dx = gradient_root_step();
        }

        const vec3 x_trial = x + dx;
        if (dx.allFinite() && x_trial.allFinite() &&
            (x_trial - origin).norm() <= max_travel) {
          const real trial_phi = eval_darboux(Q, x_trial);
          if (!std::isfinite(trial_phi) || std::abs(trial_phi) > std::abs(phi)) {
            dx = gradient_root_step();
          }
        } else {
          dx = gradient_root_step();
        }
        if (!dx.allFinite()) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::small_denominator;
          }
          return false;
        }

        x += dx;
        if (!x.allFinite()) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::nonfinite;
          }
          return false;
        }
        if ((x - origin).norm() > max_travel) {
          if (failure_out != nullptr) {
            *failure_out = darboux_ridge_failure::max_travel;
          }
          return false;
        }
        if (trace != nullptr) {
          trace->push_back(x);
        }
      }

      const real phi = eval_darboux(Q, x);
      if (x.allFinite() && std::abs(phi) <= std::sqrt(tol) * 10.0) {
        if (x_out != nullptr) {
          *x_out = x;
        }
        if (iters_out != nullptr) {
          *iters_out = max_iters;
        }
        if (failure_out != nullptr) {
          *failure_out = darboux_ridge_failure::none;
        }
        return true;
      }
      return false;
    }

    // Newton ridge search on a fixed inward ray through a Darboux signed-value
    // field. Parameterize x(t) = x_base + t * dir and solve grad D(x) · dir = 0.
    //
    // If the start point is outside the fitted field (D > 0), first solve down
    // the gradient to the zero level set, then continue the ridge search
    // inward from that surface point. Inside points (D <= 0) march inward
    // directly from the POV.
    inline darboux_ridge_estimate estimate_center_ridge(
        const vec14 &Q, const vec3 &start, const vec3 &inward_dir,
        int max_iters = 12, real tol = 1e-8,
        std::vector<vec3> *trace = nullptr,
        real max_travel = std::numeric_limits<real>::infinity()) {
      darboux_ridge_estimate out;
      if (inward_dir.norm() < 1e-12) {
        out.failure = darboux_ridge_failure::invalid_direction;
        return out;
      }
      const vec3 dir = inward_dir.normalized();
      vec3 x_base = start;
      real t = 0.0;
      if (trace != nullptr) {
        trace->clear();
        trace->push_back(start);
      }

      const real d0 = eval_darboux(Q, start);
      out.surface_residual = std::abs(d0);
      const vec3 g0 = darboux_grad(Q, start);
      const real g0_norm = g0.norm();
      const real surface_dist =
          g0_norm > 1e-12 ? std::abs(d0) / g0_norm
                           : std::numeric_limits<real>::infinity();
      const real surface_dist_tol =
          std::max(std::sqrt(tol), real(0.01) * max_travel);
      if (d0 > 0.0 && surface_dist > surface_dist_tol) {
        out.projection_attempted = true;
        vec3 x_surface = start;
        int root_iters = 0;
        darboux_ridge_failure projection_failure =
            darboux_ridge_failure::not_converged;
        if (darboux_surface_point_newton(Q, start, max_iters, tol, &x_surface,
                                         &root_iters, trace, max_travel,
                                         &start, &projection_failure)) {
          x_base = x_surface;
          out.projection_converged = true;
          out.surface_residual = std::abs(eval_darboux(Q, x_base));
          out.iterations = root_iters;
        } else {
          out.projection_failure = projection_failure;
        }
      }

      auto over_travel = [&](const vec3 &x) {
        return (x - start).norm() > max_travel;
      };

      for (int iter = 0; iter < std::max(1, max_iters); ++iter) {
        t = std::max(t, real(0.0));
        const vec3 x = x_base + t * dir;
        const real ridge = darboux_directional_deriv(Q, x, dir);
        out.residual = std::abs(ridge);
        out.iterations += 1;

        if (!std::isfinite(ridge) || !x.allFinite()) {
          out.failure = darboux_ridge_failure::nonfinite;
          return out;
        }
        if (out.residual <= tol) {
          out.center = x;
          out.travel = (x - start).norm();
          out.failure = darboux_ridge_failure::none;
          out.accepted = out.center.allFinite() && out.travel > 1e-12;
          return out;
        }

        const real denom = darboux_directional_second_deriv(Q, x, dir);
        if (!std::isfinite(denom) || std::abs(denom) < 1e-12) {
          out.failure = darboux_ridge_failure::small_denominator;
          return out;
        }

        const real dt = -ridge / denom;
        if (!std::isfinite(dt)) {
          out.failure = darboux_ridge_failure::nonfinite;
          return out;
        }

        const real t_next = t + dt;
        if (t_next < 0.0) {
          ++out.clamped_steps;
        }
        t = std::max(t_next, real(0.0));
        if (over_travel(x_base + t * dir)) {
          const vec3 x_reject = x_base + t * dir;
          out.center = x_reject;
          out.travel = (x_reject - start).norm();
          out.failure = darboux_ridge_failure::max_travel;
          if (trace != nullptr && x_reject.allFinite()) {
            trace->push_back(x_reject);
          }
          return out;
        }
        if (trace != nullptr) {
          trace->push_back(x_base + t * dir);
        }
      }

      t = std::max(t, real(0.0));
      const vec3 x = x_base + t * dir;
      out.center = x;
      out.travel = (x - start).norm();
      out.residual = std::abs(darboux_directional_deriv(Q, x, dir));
      out.failure = over_travel(x) ? darboux_ridge_failure::max_travel
                                   : darboux_ridge_failure::not_converged;
      out.accepted = false;
      return out;
    }

    class darboux_cyclide
    {
    public:
      using coefficients = vec14;

      darboux_cyclide()
      {
        A = mat14::Zero();
        b = vec14::Zero();
      }

      void accumulate(real w, const vec3 &x, const vec3 &N)
      {
        mat414 Ab = mk_darboux_A(x);
        vec4 Nb = mk_N(N);
        b += w * Ab.transpose() * Nb;
        A += w * Ab.transpose() * Ab;
      }

      vec14 solve()
      {
        vec14 x = A.colPivHouseholderQr().solve(b);
        return x;
      }

      mat14 A;
      vec14 b;
    };

  }
}
#endif