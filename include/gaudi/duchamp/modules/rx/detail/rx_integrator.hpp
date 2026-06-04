#ifndef GAUDI_DUCHAMP_RX_DETAIL_RX_INTEGRATOR_HPP
#define GAUDI_DUCHAMP_RX_DETAIL_RX_INTEGRATOR_HPP

/// Generic local timestep helpers (not shell-specific).

#include "gaudi/common.h"
#include <functional>

namespace gaudi {
namespace duchamp {
namespace rx {
namespace detail {

/// Forward Euler: `y + h * f(y, p)`.
template <class State, class Params, class RHS>
inline State step_forward_euler(State y, const Params &p, real h, RHS f) {
  return y + h * f(y, p);
}

/// Midpoint (RK2): two stages, O(h^2) per step for smooth RHS.
template <class State, class Params, class RHS>
inline State step_rk2(State y, const Params &p, real h, RHS f) {
  State k1 = f(y, p);
  State y_mid = y + (0.5 * h) * k1;
  State k2 = f(y_mid, p);
  return y + h * k2;
}

/// Classical RK4 per vertex / per degree of freedom when `State` supports scalar ops.
template <class State, class Params, class RHS>
inline State step_rk4(State y, const Params &p, real h, RHS f) {
  State k1 = f(y, p);
  State k2 = f(y + (0.5 * h) * k1, p);
  State k3 = f(y + (0.5 * h) * k2, p);
  State k4 = f(y + h * k3, p);
  return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

} // namespace detail
} // namespace rx
} // namespace duchamp
} // namespace gaudi

#endif
