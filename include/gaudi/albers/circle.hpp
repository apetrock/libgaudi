

#ifndef ALBERS_CIRCLES_H
#define ALBERS_CIRCLES_H
#include <Eigen/Dense>
#include "gaudi/common.h"
// stub: but least squares shape functions here
// sphere, cylinder, etc.

namespace gaudi
{
  namespace albers
  {

    mat3 mk_circle_A(vec2 dx, const vec2 &N)
    {
      mat3 A;
      real x = dx[0];
      real y = dx[1];
      real x2 = 2.0 * x;
      real y2 = 2.0 * y;
      real Nx = N[0];
      real Ny = N[1];

      real Nx1 = sqrt(1.0 - Nx * Nx);
      real Ny1 = sqrt(1.0 - Ny * Ny);

      A.row(0) = vec3(x2, y2, 1.0);
      // A = X^2 * I + sort of Skew(X) + 1
      A.row(1) = vec3(Nx1, -Nx, 0.0);
      A.row(2) = vec3(-Ny, Ny1, 0.0);

      return A;
    }

    vec3 mk_circle_b(const vec2 dx, const vec2 &N)
    {
      real x = dx[0];
      real y = dx[1];
      real xx = x * x;
      real yy = y * y;
      real Nx = N[0];
      real Ny = N[1];

      real Nx1 = sqrt(1.0 - Nx * Nx);
      real Ny1 = sqrt(1.0 - Ny * Ny);
      real b0 = xx + yy;
      real b1 = Nx1 * x - Nx * y;
      real b2 = -Ny * x + Ny1 * y;
      return vec3(b0, b1, b2);
    }

    vec3 mk_circle_A(vec2 dx)
    {
      mat3 A;
      real x = dx[0];
      real y = dx[1];
      real x2 = 2.0 * x;
      real y2 = 2.0 * y;
      return vec3(x2, y2, 1.0);
    }

    real mk_circle_b(const vec2 dx)
    {
      real x = dx[0];
      real y = dx[1];
      real xx = x * x;
      real yy = y * y;
      return xx + yy;
    }

    class circle
    {
    public:
      using coefficients = vec3;

      circle()
      {
        A = mat3::Zero();
        b = vec3::Zero();
      }

      void accumulate(real w, const vec2 &x, const vec2 &N)
      {
        vec3 Ai = mk_circle_A(x);
        real bi = mk_circle_b(x);
        b += w * Ai.transpose() * bi;
        A += w * Ai * Ai.transpose();
      }

      vec3 solve()
      {
        return A.colPivHouseholderQr().solve(b);
      }

      mat3 A;
      vec3 b;
    };

    class constrained_circle
    {
    public:
      using coefficients = vec3;

      constrained_circle()
      {
        A = mat3::Zero();
        b = vec3::Zero();
      }

      void accumulate(real w, const vec2 &x, const vec2 &N)
      {
        mat3 Ab = mk_circle_A(x, N);
        vec3 Nb = mk_circle_b(x, N);
        b += w * Ab.transpose() * Nb;
        A += w * Ab.transpose() * Ab;
      }

      vec3 solve()
      {
        return A.colPivHouseholderQr().solve(b);
      }

      mat3 A;
      vec3 b;
    };
  }
}

#endif