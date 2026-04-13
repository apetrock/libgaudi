#ifndef __GAUDI_CALDER_WINDING_TESTS_HPP__
#define __GAUDI_CALDER_WINDING_TESTS_HPP__

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/calder/integrators.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/test/bvh_tests.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace gaudi {
namespace test {

namespace {

real winding_brute(asawa::shell::shell &M, const std::vector<vec3> &x,
                   const vec3 &pi) {
  std::vector<index_t> fv = M.get_face_vert_ids();
  const size_t ntri = fv.size() / 3;
  real s = 0;
  for (size_t t = 0; t < ntri; ++t) {
    vec3 p0 = x[fv[3 * t + 0]];
    vec3 p1 = x[fv[3 * t + 1]];
    vec3 p2 = x[fv[3 * t + 2]];
    s += 0.25 / M_PI * va::solidAngle(pi, p0, p1, p2);
  }
  return s;
}

real max_vertex_radius(const std::vector<vec3> &x, const vec3 &c) {
  real r = 0;
  for (const auto &v : x)
    r = std::max(r, (v - c).norm());
  return r;
}

} // namespace

GAUDI_TEST(calder_fast_winding_sphere_matches_brute) {
  MeshData mesh = load_sphere_mesh();
  asawa::shell::shell &M = *mesh.shell;
  std::vector<vec3> &x = mesh.vertices;

  vec3 c(0, 0, 0);
  for (const auto &v : x)
    c += v;
  c /= static_cast<real>(x.size());

  const real r = max_vertex_radius(x, c);
  GAUDI_ASSERT(r > 1e-12);

  // Ray from centroid along +X: t=0 at center, t=1 near the surface, t>1 outside.
  // Sample every (1/16) r from 0 to 2r inclusive.
  const vec3 axis = vec3(1, 0, 0);
  constexpr int k_max = 32;
  constexpr real dt_r = 1.0 / 16.0;

  std::vector<vec3> pov;
  for (int k = 0; k <= k_max; ++k) {
    const real t = dt_r * static_cast<real>(k);
    pov.push_back(c + t * r * axis);
  }

  const real l0 = asawa::shell::avg_length(M, x);
  std::vector<real> w_fast = calder::fast_winding(M, x, pov, l0);
  GAUDI_ASSERT(w_fast.size() == pov.size());

  // fast_winding uses fast_summation with traversal split 0.5; single loose bound
  // vs brute solid-angle sum (BH is first-order).
  const real eps = 0.5;
 
  std::cerr << "\n[calder_fast_winding] radial sweep along +X from centroid; t = distance / r\n"
            << "  r (max vertex radius) = " << std::fixed << std::setprecision(6) << r
            << "  l0 = " << l0 << "\n"
            << "  t/r       w_fast      w_brute     |delta|\n";

  for (size_t i = 0; i < pov.size(); ++i) {
    const real t = dt_r * static_cast<real>(static_cast<int>(i));
    const real w_brute = winding_brute(M, x, pov[i]);
    const real delta = std::abs(w_fast[i] - w_brute);
    std::cerr << "  " << std::setw(6) << t << "  " << std::setw(10) << w_fast[i]
              << "  " << std::setw(10) << w_brute << "  " << std::setw(10) << delta
              << "\n";
    GAUDI_EXPECT(delta < eps);
  }
}

} // namespace test
} // namespace gaudi

#endif
