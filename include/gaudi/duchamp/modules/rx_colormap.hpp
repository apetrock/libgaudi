#ifndef __GAUDI_DUCHAMP_RX_COLORMAP__
#define __GAUDI_DUCHAMP_RX_COLORMAP__

#include "gaudi/common.h" // vec4, real
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace gaudi {
namespace duchamp {

/// Min–max normalize \p field to \f$[0,1]\f$ and map through \p gradient(t).
/// When max ≈ min, uses \p t = 0.5.
inline std::vector<vec4>
field_colors_minmax(const std::vector<real> &field,
                    const std::function<vec4(real t01)> &gradient) {
  std::vector<vec4> out(field.size(), vec4(0.5, 0.5, 0.5, 1.0));
  if (field.empty())
    return out;

  real mn = std::numeric_limits<real>::infinity();
  real mx = -std::numeric_limits<real>::infinity();
  for (real v : field) {
    if (std::isfinite(v)) {
      mn = std::min(mn, v);
      mx = std::max(mx, v);
    }
  }
  real denom = mx - mn;
  const real eps = 1.0e-30;
  for (std::size_t i = 0; i < field.size(); ++i) {
    real v = field[i];
    real t = 0.5;
    if (std::isfinite(v) && denom > eps)
      t = (v - mn) / denom;
    else if (!std::isfinite(v))
      t = 0.0;
    t = std::max(real(0), std::min(real(1), t));
    out[i] = gradient(t);
  }
  return out;
}

/// Map phase in radians to \f$[0,1]\f$ for a cyclic gradient (wrap to one period).
inline std::vector<vec4>
field_colors_phase(const std::vector<real> &phase_rad,
                   const std::function<vec4(real t01)> &gradient) {
  const real two_pi = 8.0 * std::atan(1.0);
  std::vector<vec4> out(phase_rad.size(), vec4(0.5, 0.5, 0.5, 1.0));
  for (std::size_t i = 0; i < phase_rad.size(); ++i) {
    real p = phase_rad[i];
    real t = 0.5;
    if (std::isfinite(p)) {
      real a = std::atan2(std::sin(p), std::cos(p));
      t = (a + two_pi / 2.0) / two_pi;
      if (t < 0)
        t += 1.0;
      if (t >= 1)
        t -= 1.0;
    }
    out[i] = gradient(t);
  }
  return out;
}

} // namespace duchamp
} // namespace gaudi

#endif
