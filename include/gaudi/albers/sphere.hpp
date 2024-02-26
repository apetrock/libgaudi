#ifndef ALBERS_SPHERES_H
#define ALBERS_SPHERES_H
#include <Eigen/Dense>
#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "ncls.hpp"

// https://arxiv.org/pdf/1506.02776.pdf
namespace gaudi
{
  namespace albers
  {
    TYPEDEF_VEC(18)
    vec18 mk_sphere_v(const vec3 &xi)
    {
      vec18 S = vec18::Zero();
      real x = xi[0];
      real y = xi[1];
      real z = xi[2];
      S[0] = x;
      S[1] = y;
      S[2] = z;
      S[3] = x * x, S[4] = y * y, S[5] = z * z;
      S[6] = x * y, S[7] = x * z, S[8] = y * z;
      S[9] = x * x * x, S[10] = x * x * y, S[11] = x * x * z;
      S[12] = y * y * x, S[13] = y * y * y, S[14] = y * y * z;
      S[15] = z * z * x, S[16] = z * z * y, S[17] = z * z * z;
      return S;
    }

    vec4 calc_sphere(real W, const vec18 &S)
    {
      real Sx = S[0], Sy = S[1], Sz = S[2];
      real Sxx = S[3], Syy = S[4], Szz = S[5];
      real Sxy = S[6], Sxz = S[7], Syz = S[8];
      real Sxxx = S[9], Sxxy = S[10], Sxxz = S[11];
      real Syyx = S[12], Syyy = S[13], Syyz = S[14];
      real Szzx = S[15], Szzy = S[16], Szzz = S[17];

      real A1 = Sxx + Syy + Szz;
      real a = 2.0 * Sx * Sx - 2.0 * W * Sxx;
      real b = 2.0 * Sx * Sy - 2.0 * W * Sxy;
      real c = 2.0 * Sx * Sz - 2.0 * W * Sxz;
      real d = -W * (Sxxx + Syyx + Szzx) + A1 * Sx;

      real e = 2.0 * Sx * Sy - 2.0 * W * Sxy;
      real f = 2.0 * Sy * Sy - 2.0 * W * Syy;
      real g = 2.0 * Sy * Sz - 2.0 * W * Syz;
      real h = -W * (Sxxy + Syyy + Szzy) + A1 * Sy;
      real j = 2.0 * Sx * Sz - 2.0 * W * Sxz;
      real k = 2.0 * Sy * Sz - 2.0 * W * Syz;
      real l = 2.0 * Sz * Sz - 2.0 * W * Szz;
      real m = -W * (Sxxz + Syyz + Szzz) + A1 * Sz;

      real delta = a * (f * l - g * k) - e * (b * l - c * k) + j * (b * g - c * f);
      vec3 X = vec3::Zero();
      X[0] =
          (d * (f * l - g * k) - h * (b * l - c * k) + m * (b * g - c * f)) / delta;
      X[1] =
          (a * (h * l - m * g) - e * (d * l - m * c) + j * (d * g - h * c)) / delta;
      X[2] =
          (a * (f * m - h * k) - e * (b * m - d * k) + j * (b * h - d * f)) / delta;

      real R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2] +
                    (A1 - 2.0 * (X[0] * Sx + X[1] * Sy + X[2] * Sz)) / W);
      return vec4(X[0], X[1], X[2], R);
    }

    class sphere
    {
    public:
      using coefficients = vec4;

      sphere()
      {
        A = vec18::Zero();
        w = 0.0;
      }

      void accumulate(real w, const vec3 &x, const vec3 &N)
      {
        A += w * mk_sphere_v(x);
        w += w;
      }

      vec4 solve()
      {
        return calc_sphere(w, A);
      }
      real w;
      vec18 A;
    };

    TYPEDEF_VEC(5)
    TYPEDEF_MAT(5)
    TYPEDEF_MAT_NM(4, 5)

    mat45 mk_sphere_A(vec3 dx, const vec3 &N)
    {
      mat45 A;
      real x = dx[0];
      real y = dx[1];
      real z = dx[2];

      real xx = dx.dot(dx);

      A.row(0) = vec5(xx, x, y, z, 1.0);
      // A = X^2 * I + sort of Skew(X) + 1
      A.row(1) = vec5(2.0 * x, 1.0, 0.0, 0.0, 0.0);
      A.row(2) = vec5(2.0 * y, 0.0, 1.0, 0.0, 0.0);
      A.row(3) = vec5(2.0 * z, 0.0, 0.0, 1.0, 0.0);

      return A;
    }

    class constrained_sphere
    {
    public:

      /*
       * This code implements the sphere fitting method described in:
       * Pratt, V. "Direct least-squares fitting of algebraic surfaces."
       * Computer Graphics, Vol. 21, pages 145-152 (1987).
       * 
       * This method fits a sphere to a set of points by minimizing the algebraic distance.
       * I've augmented the method to include a normal constraint.
       * 
       * the pratt method assumes a normalization, which I don't do here.
       * also an iterative update of the center is not implemented, but I could add 
       * a constructor that would allow you to use this accumulator iteratively
       */
      using coefficients = vec4;

      constrained_sphere()
      {
        A = mat5::Zero();
        b = vec5::Zero();
      }

      void accumulate(real w, const vec3 &x, const vec3 &N)
      {
        vec3 dx = x - vec3::Zero();
        // placeholder for center, somehow one should be able to pass
        // center for iterative conditioning

        mat45 Ab = mk_sphere_A(dx, N);
        vec4 Nb = mk_N(N);
        b += w * Ab.transpose() * Nb;
        A += w * Ab.transpose() * Ab;
      }

      vec4 solve()
      {
        real det = A.determinant();

        if (det < 1e-10)
        {
          return vec4::Zero();
        }

        vec5 x = A.colPivHouseholderQr().solve(b);
        real A2 = 2.0 * x[0];
        vec3 BCD = x.segment(1, 3);
        vec3 cen = -BCD / A2;
        real r = sqrt(2.0 * A2 + BCD.dot(BCD)) / A2;
        return vec4(cen[0], cen[1], cen[2], r);
      }

      mat5 A;
      vec5 b;
    };

  }
}

#endif