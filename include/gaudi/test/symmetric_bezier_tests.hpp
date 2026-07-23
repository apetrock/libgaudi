#ifndef __GAUDI_SYMMETRIC_BEZIER_TESTS_HPP__
#define __GAUDI_SYMMETRIC_BEZIER_TESTS_HPP__

#include <cmath>
#include <vector>

#include "gaudi/albers/symmetric_bezier.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

inline std::vector<vec3> make_symmetric_c_arc() {
  // Mild C in the xy-plane, mirror across yz (normal = x).
  // Points are even in x about the apex at x=0.
  return {
      vec3(-0.40, 0.00, 0.0),
      vec3(-0.20, 0.15, 0.0),
      vec3(0.00, 0.22, 0.0),
      vec3(0.20, 0.15, 0.0),
      vec3(0.40, 0.00, 0.0),
  };
}

GAUDI_TEST(symmetric_bezier_exact_mirror_fit) {
  const std::vector<vec3> pts = make_symmetric_c_arc();
  // Along-arc / left-right for this C is +x.
  const vec3 T_curve = vec3(1.0, 0.0, 0.0);

  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_quartic_bezier(pts, T_curve);
  GAUDI_ASSERT(fit.ok);

  const vec3 &n = fit.frame.n_mirror;
  // Mirror normal should align with ±x for this cloud.
  GAUDI_EXPECT(std::abs(std::abs(n.dot(vec3::UnitX())) - 1.0) < 0.15);

  // Curve should be nearly reflection-symmetric.
  real max_sym = 0.0;
  for (int i = 1; i < 10; ++i) {
    const real t = real(i) / 10.0;
    const vec3 a = albers::eval_quartic_bezier(fit.P, t);
    const vec3 b = albers::eval_quartic_bezier(fit.P, 1.0 - t);
    const vec3 br =
        albers::reflect_plane(b, fit.frame.origin, fit.frame.n_mirror);
    max_sym = std::max(max_sym, (a - br).norm());
  }
  GAUDI_EXPECT(max_sym < 5e-2);

  // Samples near the curve.
  real max_dist = 0.0;
  for (const vec3 &p : pts)
    max_dist = std::max(max_dist, (fit.eval_at_point(p) - p).norm());
  GAUDI_EXPECT(max_dist < 5e-2);
}

GAUDI_TEST(symmetric_bezier_noisy_restores_evenness) {
  std::vector<vec3> pts = make_symmetric_c_arc();
  // Antisymmetric noise (pushes one side only).
  pts[1] += vec3(0.0, 0.0, 0.04);
  pts[3] += vec3(0.0, 0.0, -0.01);

  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_quartic_bezier(pts, vec3(1.0, 0.0, 0.0));
  GAUDI_ASSERT(fit.ok);

  real max_sym = 0.0;
  for (int i = 1; i < 10; ++i) {
    const real t = real(i) / 10.0;
    const vec3 a = albers::eval_quartic_bezier(fit.P, t);
    const vec3 b = albers::eval_quartic_bezier(fit.P, 1.0 - t);
    const vec3 br =
        albers::reflect_plane(b, fit.frame.origin, fit.frame.n_mirror);
    max_sym = std::max(max_sym, (a - br).norm());
  }
  // Looser than exact, but should still be fairly even.
  GAUDI_EXPECT(max_sym < 0.12);
}

GAUDI_TEST(symmetric_bezier_degenerate_spectrum_guard_rejects) {
  // Near-isotropic in-plane pentagon: two large singular values nearly equal
  // → separation guard should reject.
  std::vector<vec3> pts;
  for (int i = 0; i < 5; ++i) {
    const real th = 2.0 * M_PI * real(i) / 5.0;
    pts.emplace_back(std::cos(th), std::sin(th), 0.0);
  }
  const albers::symmetric_bezier_fit fit =
      albers::fit_symmetric_quartic_bezier(pts, vec3(1.0, 0.0, 0.0), 1.15);
  GAUDI_EXPECT(!fit.ok);
}

} // namespace test
} // namespace gaudi

#endif
