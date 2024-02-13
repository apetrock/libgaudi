
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