#ifndef ALBERS_NORMAL_CONSTRAINED_QUADRICS_H
#define ALBERS_NORMAL_CONSTRAINED_QUADRICS_H
#include <Eigen/Dense>
#include "gaudi/common.h"
#include "gaudi/vec_addendum.h"
#include "ncls.hpp"
// stub: but least squares shape functions here
// sphere, cylinder, etc.

namespace gaudi
{
    namespace albers
    {


        TYPEDEF_VEC(7)
        TYPEDEF_VEC(10)
        TYPEDEF_MAT(10)
        TYPEDEF_MAT_NM(4, 10)

        mat410 mk_quad_A(vec3 dx)
        {
            mat410 A;
            real x = dx[0];
            real y = dx[1];
            real z = dx[2];
            real x2 = 2.0 * x;
            real y2 = 2.0 * y;
            real z2 = 2.0 * z;

            real xx = x * x;
            real yy = y * y;
            real zz = z * z;
            real xy = 2.0 * x * y;
            real xz = 2.0 * x * z;
            real yz = 2.0 * y * z;

            A.row(0) = vec10(xx, yy, zz, xy, xz, yz, x2, y2, z2, 1.0);
            // A = X^2 * I + sort of Skew(X) + 1
            A.row(1) = vec10(x2, 0.0, 0.0, y2, z2, 0.0, 2.0, 0.0, 0.0, 0.0);
            A.row(2) = vec10(0.0, y2, 0.0, x2, 0.0, z2, 0.0, 2.0, 0.0, 0.0);
            A.row(3) = vec10(0.0, 0.0, z2, 0.0, x2, y2, 0.0, 0.0, 2.0, 0.0);

            return A;
        }


        mat3 quadric_hessian(const vec10 &Q)
        {
            mat3 A = mat3::Zero();
            A.row(0) = vec3(Q[0], Q[3], Q[4]);
            A.row(1) = vec3(Q[3], Q[1], Q[5]);
            A.row(2) = vec3(Q[4], Q[5], Q[2]);
            return A;
        }
        
        vec3 quadric_grad(const vec10 &Q, const vec3 x)
        {

            vec3 g = {
                2.0 * Q[0] * x[0] + 2.0 * Q[3] * x[1] + 2.0 * Q[4] * x[2] + 2.0 * Q[6], //
                2.0 * Q[3] * x[0] + 2.0 * Q[1] * x[1] + 2.0 * Q[5] * x[2] + 2.0 * Q[7], //
                2.0 * Q[4] * x[0] + 2.0 * Q[5] * x[1] + 2.0 * Q[2] * x[2] + 2.0 * Q[8]  //
            };
            return g;
        }

        vec3 quadric_center(const vec10 &Q)
        {
            mat3 A = quadric_hessian(Q);
            vec3 b = -vec3(Q[6], Q[7], Q[8]);
            //hmmm... what if we solve in least squares sense?
            A =  A + 1e-8 * mat3::Identity();
            Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
            int rank = lu_decomp.rank();
            if (rank < 3)
            {
                return vec3::Zero();
            }

            vec3 c = lu_decomp.solve(b);
            return c;
            
        }



        class quadric
        {
        public:

        using coefficients = vec10;
        
            quadric()
            {
                A = mat10::Zero();
                b = vec10::Zero();
            }

            void accumulate(real w, const vec3 &x, const vec3 &N)
            {
                mat410 Ab = mk_quad_A(x);
                vec4 Nb = mk_N(N);
                b += w * Ab.transpose() * Nb;
                A += w * Ab.transpose() * Ab;
            }
            
            vec10 solve()
            {
                vec10 x = A.colPivHouseholderQr().solve(b);
                return x;
            }

            mat10 A;
            vec10 b;
        };


    }
}

#endif