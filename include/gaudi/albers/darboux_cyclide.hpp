
#ifndef ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#define ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
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
      int negative_steps = 0;
      int backward_clamps = 0;
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

    struct darboux_line_newton_result {
      vec3 x = vec3::Zero();
      real t = 0.0;
      real residual = 0.0;
      int iterations = 0;
      int clamped_steps = 0;
      int negative_steps = 0;
      int backward_clamps = 0;
      darboux_ridge_failure failure = darboux_ridge_failure::not_converged;
      bool converged = false;
    };

    // Shared 1D Newton machinery for scalar equations along a fixed ray:
    //   x(t) = start + t * dir, F(t) = 0, t <- t - F/F'.
    // The caller supplies F and F' at the current point. This keeps the
    // surface projection and ridge solve mechanically identical while leaving
    // their scalar equations explicit at the call site.
    template <typename F_TYPE>
    inline darboux_line_newton_result darboux_line_newton_on_ray(
        const vec3 &start, const vec3 &dir_in, int max_iters, real tol,
        F_TYPE &&eval_value_deriv, std::vector<vec3> *trace = nullptr,
        real max_travel = std::numeric_limits<real>::infinity(),
        const vec3 *travel_origin = nullptr, bool clamp_nonnegative = false,
        real final_tol = -1.0,
        real max_step = std::numeric_limits<real>::infinity()) {
      darboux_line_newton_result out;
      out.x = start;
      const vec3 origin = travel_origin != nullptr ? *travel_origin : start;
      if (dir_in.norm() < 1e-12) {
        out.failure = darboux_ridge_failure::invalid_direction;
        return out;
      }

      const vec3 dir = dir_in.normalized();
      real t = 0.0;
      for (int iter = 0; iter < std::max(1, max_iters); ++iter) {
        if (clamp_nonnegative) {
          t = std::max(t, real(0.0));
        }
        const vec3 x = start + t * dir;
        const std::pair<real, real> value_deriv = eval_value_deriv(x, dir);
        const real value = value_deriv.first;
        const real deriv = value_deriv.second;
        out.x = x;
        out.t = t;
        out.residual = std::abs(value);
        out.iterations = iter + 1;

        if (!x.allFinite() || !std::isfinite(value)) {
          out.failure = darboux_ridge_failure::nonfinite;
          return out;
        }
        if (out.residual <= tol) {
          out.failure = darboux_ridge_failure::none;
          out.converged = true;
          return out;
        }
        if (!std::isfinite(deriv) || std::abs(deriv) < 1e-12) {
          out.failure = darboux_ridge_failure::small_denominator;
          return out;
        }

        real dt = -value / deriv;
        if (!std::isfinite(dt)) {
          out.failure = darboux_ridge_failure::nonfinite;
          return out;
        }
        if (std::isfinite(max_step) && max_step > 0.0 &&
            std::abs(dt) > max_step) {
          dt = std::copysign(max_step, dt);
          ++out.clamped_steps;
        }
        if (dt < 0.0) {
          ++out.negative_steps;
        }

        const real t_next_raw = t + dt;
        if (clamp_nonnegative && t_next_raw < 0.0) {
          ++out.clamped_steps;
          ++out.backward_clamps;
        }
        t = clamp_nonnegative ? std::max(t_next_raw, real(0.0)) : t_next_raw;
        const vec3 x_next = start + t * dir;
        if (!x_next.allFinite()) {
          out.failure = darboux_ridge_failure::nonfinite;
          return out;
        }
        if ((x_next - origin).norm() > max_travel) {
          out.x = x_next;
          out.t = t;
          out.failure = darboux_ridge_failure::max_travel;
          if (trace != nullptr) {
            trace->push_back(x_next);
          }
          return out;
        }
        if (trace != nullptr) {
          trace->push_back(x_next);
        }
      }

      out.x = start + t * dir;
      out.t = t;
      if (out.x.allFinite()) {
        const std::pair<real, real> value_deriv = eval_value_deriv(out.x, dir);
        out.residual = std::abs(value_deriv.first);
        if (final_tol >= 0.0 && out.residual <= final_tol) {
          out.failure = darboux_ridge_failure::none;
          out.converged = true;
          return out;
        }
      }
      out.failure = (out.x - origin).norm() > max_travel
                        ? darboux_ridge_failure::max_travel
                        : darboux_ridge_failure::not_converged;
      return out;
    }

    // 1D Newton-Raphson to find the zero level set along a fixed ray.
    //
    // Parameterize x(t) = start + t * dir and solve D(x(t)) = 0 with:
    //   t <- t - D(x) / (grad D(x) dot dir)
    // This is intentionally constrained to the supplied direction, so it avoids
    // tangential motion from an underdetermined 3D scalar solve near D = 0.
    inline bool darboux_surface_point_line_newton(
        const vec14 &Q, const vec3 &start, const vec3 &dir_in, int max_iters,
        real tol, vec3 *x_out, int *iters_out = nullptr,
        std::vector<vec3> *trace = nullptr,
        real max_travel = std::numeric_limits<real>::infinity(),
        darboux_ridge_failure *failure_out = nullptr) {
      const real surface_tol = std::max(tol, std::sqrt(tol));
      const darboux_line_newton_result result = darboux_line_newton_on_ray(
          start, dir_in, max_iters, surface_tol,
          [&](const vec3 &x, const vec3 &dir) {
            return std::make_pair(eval_darboux(Q, x), darboux_grad(Q, x).dot(dir));
          },
          trace, max_travel, nullptr, false, std::sqrt(tol) * 10.0);
      if (x_out != nullptr) {
        *x_out = result.x;
      }
      if (iters_out != nullptr) {
        *iters_out = result.iterations;
      }
      if (failure_out != nullptr) {
        *failure_out = result.failure;
      }
      return result.converged;
    }

    // Newton ridge search on a fixed ray through a Darboux signed-value field.
    // Parameterize x(t) = start + t * dir and solve grad D(x) · dir = 0.
    // Surface projection, if needed, should happen before calling this helper.
    inline darboux_ridge_estimate estimate_center_ridge(
        const vec14 &Q, const vec3 &start, const vec3 &inward_dir,
        int max_iters = 12, real tol = 1e-8,
        std::vector<vec3> *trace = nullptr,
        real max_travel = std::numeric_limits<real>::infinity(),
        real max_step = std::numeric_limits<real>::infinity()) {
      darboux_ridge_estimate out;
      if (inward_dir.norm() < 1e-12) {
        out.failure = darboux_ridge_failure::invalid_direction;
        return out;
      }
      const vec3 dir = inward_dir.normalized();
      vec3 x_base = start;
      if (trace != nullptr) {
        trace->clear();
        trace->push_back(start);
      }

      const real d0 = eval_darboux(Q, start);
      out.surface_residual = std::abs(d0);

      const darboux_line_newton_result ridge_result = darboux_line_newton_on_ray(
          x_base, dir, max_iters, tol,
          [&](const vec3 &x, const vec3 &ray_dir) {
            return std::make_pair(
                darboux_directional_deriv(Q, x, ray_dir),
                darboux_directional_second_deriv(Q, x, ray_dir));
          },
          trace, max_travel, &start, true, -1.0, max_step);
      out.center = ridge_result.x;
      out.travel = (out.center - start).norm();
      out.residual = ridge_result.residual;
      out.iterations += ridge_result.iterations;
      out.clamped_steps += ridge_result.clamped_steps;
      out.negative_steps += ridge_result.negative_steps;
      out.backward_clamps += ridge_result.backward_clamps;
      out.failure = ridge_result.failure;
      out.accepted = ridge_result.converged && out.center.allFinite() &&
                     out.travel > 1e-12;
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

        const vec3 n = N.normalized();
        const vec3 seed =
            std::abs(n.dot(vec3::UnitZ())) < 0.9 ? vec3::UnitZ() : vec3::UnitX();
        const vec3 t1 = seed.cross(n).normalized();
        const vec3 t2 = n.cross(t1).normalized();

        const auto value_row = Ab.row(0);
        const auto tangent_row_1 =
            t1[0] * Ab.row(1) + t1[1] * Ab.row(2) + t1[2] * Ab.row(3);
        const auto tangent_row_2 =
            t2[0] * Ab.row(1) + t2[1] * Ab.row(2) + t2[2] * Ab.row(3);

        A += w * value_row.transpose() * value_row;
        A += w * tangent_row_1.transpose() * tangent_row_1;
        A += w * tangent_row_2.transpose() * tangent_row_2;

        // Orientation only for the homogeneous solve.
        b += w * Ab.transpose() * Nb;
      }

      vec14 solve()
      {
        Eigen::SelfAdjointEigenSolver<mat14> es(A);
        vec14 x = es.eigenvectors().col(0);
        if (x.dot(b) < 0.0)
        {
          x *= -1.0;
        }
        return x;
      }

      mat14 A;
      vec14 b;
    };

    class normal_constrained_darboux_cyclide
    {
    public:
      using coefficients = vec14;

      normal_constrained_darboux_cyclide()
      {
        A = mat14::Zero();
        b = vec14::Zero();
      }

      void accumulate(real w, const vec3 &x, const vec3 &N)
      {
        mat414 Ab = mk_darboux_A(x);
        vec4 Nb = mk_N(N);

        A += w * Ab.transpose() * Ab;
        b += w * Ab.transpose() * Nb;
      }

      vec14 solve()
      {
        Eigen::ColPivHouseholderQR<mat14> qr(A);
        return qr.solve(b);
      }

      mat14 A;
      vec14 b;
    };

    using tangent_plane_darboux_cyclide = darboux_cyclide;

  }
}
#endif