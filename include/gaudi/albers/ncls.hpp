#ifndef NORMAL_CONSTRAINED_LEAST_SQUARES_H
#define NORMAL_CONSTRAINED_LEAST_SQUARES_H
#include <Eigen/Dense>
#include "gaudi/common.h"

namespace gaudi
{
    namespace albers
    {
        vec4 mk_N(const vec3 &N) { return vec4(0.0, N[0], N[1], N[2]); }

        inline vec6 vech3(const mat3 &M)
        {
            return vec6(M(0, 0), M(1, 1), M(2, 2), M(0, 1), M(0, 2), M(1, 2));
        }

        inline vec3 vech2(const mat2 &M)
        {
            return vec3(M(0, 0), M(1, 1), M(0, 1));
        }

        inline mat3 unvech3(const vec6 &v)
        {
            mat3 M;
            M << v[0], v[3], v[4], v[3], v[1], v[5], v[4], v[5], v[2];
            return M;
        }

        inline const vec6 &hessian_frobenius_weight3()
        {
            static const vec6 W(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
            return W;
        }

        inline const vec3 &hessian_frobenius_weight2()
        {
            static const vec3 W(1.0, 1.0, 2.0);
            return W;
        }

        template <typename MatH, typename MatA, typename VecB>
        inline void accumulate_hessian_block3(real hessian_weight, real w,
                                            const MatH &H, const mat3 &S,
                                            MatA &A, VecB &b)
        {
            if (hessian_weight <= 0.0)
            {
                return;
            }
            const vec6 s = vech3(S);
            const vec6 &W = hessian_frobenius_weight3();
            const vec6 Ws = W.cwiseProduct(s);
            A += hessian_weight * w * (H.transpose() * W.asDiagonal() * H);
            b += hessian_weight * w * (H.transpose() * Ws);
        }

        template <typename MatH, typename MatA, typename VecB>
        inline void accumulate_hessian_block2(real hessian_weight, real w,
                                              const MatH &H, const mat2 &S,
                                              MatA &A, VecB &b)
        {
            if (hessian_weight <= 0.0)
            {
                return;
            }
            const vec3 s = vech2(S);
            const vec3 &W = hessian_frobenius_weight2();
            const vec3 Ws = W.cwiseProduct(s);
            A += hessian_weight * w * (H.transpose() * W.asDiagonal() * H);
            b += hessian_weight * w * (H.transpose() * Ws);
        }
    }
}

#endif
