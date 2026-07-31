#ifndef __GAUDI_SYMMETRIC_BEZIER_TESTS_HPP__
#define __GAUDI_SYMMETRIC_BEZIER_TESTS_HPP__

#include <cmath>
#include <vector>

#include "gaudi/albers/symmetric_bezier.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

// Screw-symmetric S-curve: central inversion through origin
// (s, v) → (-s, -v) — coil/flip symmetry, not planar mirror.
inline std::vector<vec3> make_screw_symmetric_arc() {
  return {
      vec3(-0.40, 0.00, -0.10),
      vec3(-0.20, 0.05, -0.15),
      vec3(0.00, 0.00, 0.00),
      vec3(0.20, -0.05, 0.15),
      vec3(0.40, 0.00, 0.10),
  };
}

GAUDI_TEST(symmetric_bezier_exact_screw_fit) {
  const std::vector<vec3> pts = make_screw_symmetric_arc();
  const vec3 T_curve = vec3(1.0, 0.0, 0.0);

  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_cubic_bezier(pts, T_curve);
  GAUDI_ASSERT(fit.ok);

  // Curve should obey B(1-t) ≈ screw_flip(B(t)).
  real max_sym = 0.0;
  for (int i = 1; i < 10; ++i) {
    const real t = real(i) / 10.0;
    const vec3 a = albers::eval_cubic_bezier(fit.P, t);
    const vec3 b = albers::eval_cubic_bezier(fit.P, 1.0 - t);
    const vec3 bf = albers::screw_flip(b, fit.frame.origin);
    max_sym = std::max(max_sym, (a - bf).norm());
  }
  GAUDI_EXPECT(max_sym < 5e-2);

  real max_dist = 0.0;
  for (const vec3 &p : pts)
    max_dist = std::max(max_dist, (fit.eval_at_point(p) - p).norm());
  GAUDI_EXPECT(max_dist < 5e-2);
}

GAUDI_TEST(symmetric_bezier_noisy_restores_screw) {
  std::vector<vec3> pts = make_screw_symmetric_arc();
  pts[1] += vec3(0.0, 0.03, 0.0);
  pts[3] += vec3(0.0, 0.01, -0.02);

  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_cubic_bezier(pts, vec3(1.0, 0.0, 0.0));
  GAUDI_ASSERT(fit.ok);

  real max_sym = 0.0;
  for (int i = 1; i < 10; ++i) {
    const real t = real(i) / 10.0;
    const vec3 a = albers::eval_cubic_bezier(fit.P, t);
    const vec3 b = albers::eval_cubic_bezier(fit.P, 1.0 - t);
    const vec3 bf = albers::screw_flip(b, fit.frame.origin);
    max_sym = std::max(max_sym, (a - bf).norm());
  }
  GAUDI_EXPECT(max_sym < 0.12);
}

GAUDI_TEST(symmetric_bezier_degenerate_spectrum_guard_rejects) {
  std::vector<vec3> pts;
  for (int i = 0; i < 5; ++i) {
    const real th = 2.0 * M_PI * real(i) / 5.0;
    pts.emplace_back(std::cos(th), std::sin(th), 0.0);
  }
  // Isotropic in-plane: both ratios near 1 → guard rejects.
  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_cubic_bezier(pts, vec3(1.0, 0.0, 0.0), 1.15);
  GAUDI_EXPECT(!fit.ok);
}

} // namespace test
} // namespace gaudi

#endif
