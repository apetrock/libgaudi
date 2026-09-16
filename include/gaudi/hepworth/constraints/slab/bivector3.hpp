#ifndef GAUDI_HEPWORTH_BIVECTOR3_HPP
#define GAUDI_HEPWORTH_BIVECTOR3_HPP

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {
namespace block {

/// Grade-2 element of Cl(3,0): e₀∧e₁ ↔ Hodge dual vector (e₀×e₁) in R³.
struct bivector3 {
  vec3 dual;

  static bivector3 wedge(const vec3 &a, const vec3 &b) { return {a.cross(b)}; }

  bivector3 operator*(real s) const { return {dual * s}; }
};

inline bivector3 operator*(real s, const bivector3 &b) { return b * s; }

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_BIVECTOR3_HPP
