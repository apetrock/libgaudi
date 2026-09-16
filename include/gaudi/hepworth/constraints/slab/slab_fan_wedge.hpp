#ifndef GAUDI_HEPWORTH_SLAB_FAN_WEDGE_HPP
#define GAUDI_HEPWORTH_SLAB_FAN_WEDGE_HPP

#include "gaudi/common.h"
#include "trivector3.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// One hub fan triangle extruded along spine d: τ = (e₀∧e₁)∧d.
struct slab_fan_wedge {
  index_t hub = -1;
  index_t i0 = -1;
  index_t i1 = -1;
  real weight = 1.0;
  trivector3 tau_rest;
  real omega_rest = 0.0;
  bool bottom_cap = false;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_FAN_WEDGE_HPP
