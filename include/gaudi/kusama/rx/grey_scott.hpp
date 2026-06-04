#ifndef GAUDI_KUSAMA_RX_GREY_SCOTT_HPP
#define GAUDI_KUSAMA_RX_GREY_SCOTT_HPP

/// Grey–Scott reaction kinetics (mesh-agnostic).

#include "gaudi/common.h"
#include <array>
#include <cmath>

#include <Eigen/LU>

namespace gaudi {
namespace kusama {
namespace rx {
namespace grey_scott {

/// Backward–Euler reaction substep: solve \((u,v) - (u_0,v_0) = h G(u,v)\) with Newton.
inline std::array<real, 2> backward_euler(std::array<real, 2> uv0a, const real &f,
                                            const real &k, const real &h,
                                            const int N = 40) {
  vec2 uv0 = vec2(uv0a[0], uv0a[1]);
  vec2 uvi = uv0;
  for (int j = 0; j < N; j++) {
    real u = uvi[0];
    real v = uvi[1];
    real uv2 = u * v * v;
    real v2 = v * v;
    real uv = u * v;
    vec2 G;
    G(0) = (-uv2 + f * (1.0 - u));
    G(1) = (uv2 - (f + k) * v);

    vec2 F = uvi - uv0 - h * G;
    mat2 dG;
    dG(0, 0) = -(v2 + f);
    dG(0, 1) = -2.0 * uv;
    dG(1, 0) = v2;
    dG(1, 1) = 2.0 * uv - (f + k);

    mat2 dF = mat2::Identity() - h * dG;
    Eigen::PartialPivLU<mat2> lu = dF.partialPivLu();
    uvi += lu.solve(-F);

    if (F.norm() < 1.0e-12)
      break;
  }
  return {uvi[0], uvi[1]};
}

/// Trapezoidal reaction substep: uses average of \((G+G_0)\) in the implicit map.
inline std::array<real, 2> trapezoid(std::array<real, 2> uv0a, const real &f, const real &k,
                                    const real &h, const int N = 40) {
  vec2 uv0 = vec2(uv0a[0], uv0a[1]);
  real u0 = uv0[0];
  real v0 = uv0[1];
  real uv20 = u0 * v0 * v0;
  vec2 G0;
  G0(0) = (-uv20 + f * (1.0 - u0));
  G0(1) = (uv20 - (f + k) * v0);

  vec2 uvi = uv0;

  for (int j = 0; j < N; j++) {
    real u = uvi[0];
    real v = uvi[1];
    real uv2 = u * v * v;

    real v2 = v * v;
    real uv = u * v;
    vec2 G;
    G(0) = (-uv2 + f * (1.0 - u));
    G(1) = (uv2 - (f + k) * v);

    vec2 F = uvi - uv0 - 0.5 * h * (G + G0);
    mat2 dG;
    dG(0, 0) = -(v2 + f);
    dG(0, 1) = -2.0 * uv;
    dG(1, 0) = v2;
    dG(1, 1) = 2.0 * uv - (f + k);

    mat2 dF = mat2::Identity() - 0.5 * h * dG;
    Eigen::PartialPivLU<mat2> lu = dF.partialPivLu();
    uvi += lu.solve(-F);

    if (F.norm() < 1.0e-12)
      break;
  }
  return {uvi[0], uvi[1]};
}

} // namespace grey_scott
} // namespace rx
} // namespace kusama
} // namespace gaudi

#endif
