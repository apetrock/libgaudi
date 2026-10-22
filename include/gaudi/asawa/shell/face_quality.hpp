#ifndef GAUDI_ASAWA_SHELL_FACE_QUALITY_HPP
#define GAUDI_ASAWA_SHELL_FACE_QUALITY_HPP

#include "gaudi/common.h"

#include <cmath>

namespace gaudi {
namespace asawa {
namespace shell {

/// `4 sqrt(3) A / (l0^2 + l1^2 + l2^2)`: 1 for an equilateral triangle, 0 for
/// a degenerate one (needle or cap).
inline real face_quality(const vec3 &a, const vec3 &b, const vec3 &c) {
  const real A = real(0.5) * (b - a).cross(c - a).norm();
  const real s = (b - a).squaredNorm() + (c - b).squaredNorm() + (a - c).squaredNorm();
  if (!std::isfinite(A) || !std::isfinite(s) || !(s > real(0)))
    return real(0);
  return real(4) * std::sqrt(real(3)) * A / s;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi

#endif
