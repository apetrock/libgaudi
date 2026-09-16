#ifndef GAUDI_HEPWORTH_TRIVECTOR3_HPP
#define GAUDI_HEPWORTH_TRIVECTOR3_HPP

#include "gaudi/common.h"

#include "bivector3.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Oriented volume 3-form in R³ (grade-3 of Cl(3,0)); Λ³≅R via pseudoscalar I.
struct trivector3 {
  real coeff = 0.0;

  /// e₀∧e₁∧e₂ coefficient (same as scalar triple product).
  static trivector3 wedge(const vec3 &a, const vec3 &b, const vec3 &c) {
    return {a.dot(b.cross(c))};
  }

  /// (e₀∧e₁)∧d — bivector–vector wedge in 3D.
  static trivector3 wedge(const bivector3 &B, const vec3 &d) {
    return {B.dual.dot(d)};
  }

  trivector3 &operator+=(const trivector3 &o) {
    coeff += o.coeff;
    return *this;
  }

  trivector3 operator+(const trivector3 &o) const { return {coeff + o.coeff}; }

  trivector3 operator-(const trivector3 &o) const { return {coeff - o.coeff}; }

  trivector3 operator*(real s) const { return {coeff * s}; }

  real squared_norm() const { return coeff * coeff; }
};

inline trivector3 operator*(real s, const trivector3 &t) { return t * s; }

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_TRIVECTOR3_HPP
