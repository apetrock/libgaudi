#ifndef GAUDI_KUSAMA_RX_GINZBURG_LANDAU_HPP
#define GAUDI_KUSAMA_RX_GINZBURG_LANDAU_HPP

/// Ginzburg–Landau cubic reaction in real (u, v) form, mesh-agnostic.

#include "gaudi/common.h"
#include <cmath>

namespace gaudi {
namespace kusama {
namespace rx {
namespace ginzburg_landau {

/// Newton step for implicit reaction: state minus previous minus `h·c·G(u,v)`.
inline void reaction_newton_2d(real &u, real &v, real u0, real v0, real c, real b,
                               real h) {
  for (int it = 0; it < 25; ++it) {
    real r2 = u * u + v * v;
    if (!std::isfinite(r2))
      break;
    real gu = u - r2 * (u - b * v);
    real gv = v - r2 * (v + b * u);
    real fu = u - u0 - h * c * gu;
    real fv = v - v0 - h * c * gv;
    if (fu * fu + fv * fv < 1e-20 * (1.0 + u0 * u0 + v0 * v0)) {
      return;
    }
    real dgu_du = 1.0 - (2.0 * u) * (u - b * v) - r2;
    real dgu_dv = -(2.0 * v) * (u - b * v) + r2 * b;
    real dgv_du = -(2.0 * u) * (v + b * u) - r2 * b;
    real dgv_dv = 1.0 - (2.0 * v) * (v + b * u) - r2;
    const real J00 = 1.0 - h * c * dgu_du;
    const real J01 = -h * c * dgu_dv;
    const real J10 = -h * c * dgv_du;
    const real J11 = 1.0 - h * c * dgv_dv;
    const real det = J00 * J11 - J01 * J10;
    if (std::abs(det) < 1e-20)
      break;
    u -= (J11 * fu - J01 * fv) / det;
    v -= (-J10 * fu + J00 * fv) / det;
  }
}

} // namespace ginzburg_landau
} // namespace rx
} // namespace kusama
} // namespace gaudi

#endif
