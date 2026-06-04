#ifndef GAUDI_KUSAMA_RX_SWIFT_HOHENBERG_HPP
#define GAUDI_KUSAMA_RX_SWIFT_HOHENBERG_HPP

/// Swift–Hohenberg local reaction (scalar cubic), mesh-agnostic.

#include "gaudi/common.h"
#include <cmath>

namespace gaudi {
namespace kusama {
namespace rx {
namespace swift_hohenberg {

enum class reaction_step_mode { forward_euler, newton };

/// Newton solve for `u − u0 − h·c·(ε·u − g·u³) = 0`.
inline real reaction_newton_solve(real u0, real c, real eps, real g, real h) {
  auto F = [u0, c, h, eps, g](real u) {
    return u - u0 - h * c * (eps * u - g * u * u * u);
  };
  auto dF = [c, h, eps, g](real u) {
    return 1.0 - h * c * (eps - 3.0 * g * u * u);
  };
  real u = u0;
  for (int it = 0; it < 20; ++it) {
    const real Fu = F(u);
    if (std::abs(Fu) < 1e-12 * (1.0 + std::abs(u0)))
      return u;
    const real J = dF(u);
    if (std::abs(J) < 1e-20)
      break;
    u -= Fu / J;
  }
  return u;
}

} // namespace swift_hohenberg
} // namespace rx
} // namespace kusama
} // namespace gaudi

#endif
