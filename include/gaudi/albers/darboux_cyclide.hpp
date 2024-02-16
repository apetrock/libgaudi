
#ifndef ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#define ALBERS_NORMAL_CONSTRAINED_CYCLIDE_H
#include <Eigen/Dense>
#include "gaudi/common.h"
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

      A.col(10) = vec4(xyz2, 4.0 * xxyz, 4.0 * yxyz, 4.0 * zxyz);    // lambda
      A.col(11) = vec4(x * xyz, 2.0 * xx + xyz, 2.0 * xy, 2.0 * xz); // mu
      A.col(12) = vec4(y * xyz, 2.0 * xy, 2.0 * yy + xyz, 2.0 * yz); // nu
      A.col(13) = vec4(z * xyz, 2.0 * xz, 2.0 * yz, 2.0 * zz + xyz); // kappa

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