#ifndef GAUDI_TEST_STRIPES_TESTS_HPP
#define GAUDI_TEST_STRIPES_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/kusama/stripes.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <random>
#include <vector>

namespace gaudi {
namespace test {

namespace stripes_detail {

/// Unit tangent of the (1, q) torus-knot family about the z axis through
/// each vertex of a torus with core radius `R`: `rho phi_hat + q r theta_hat`.
inline std::vector<vec3> helix(const std::vector<vec3> &x, real R, int q) {
  std::vector<vec3> u(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    const vec3 &p = x[i];
    const real rho = std::hypot(p.x(), p.y());
    const vec3 phi(-p.y() / rho, p.x() / rho, 0);
    const vec3 c(p.x() / rho * R, p.y() / rho * R, 0);
    const vec3 n = (p - c).normalized();
    const vec3 theta = phi.cross(n);
    const real r = (p - c).norm();
    u[i] = (rho * phi + real(q) * r * theta).normalized();
  }
  return u;
}

inline std::vector<vec3> across(const asawa::shell::shell &M, const std::vector<vec3> &x,
                                const std::vector<vec3> &u) {
  std::vector<vec3> X(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    vec3 n = asawa::shell::vert_normal(M, asawa::shell::vert_id(static_cast<int>(i)), x);
    n.normalize();
    X[i] = n.cross(u[i]).normalized();
  }
  return X;
}

inline real mean_angle(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real s = 0;
  int n = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!(a[i].squaredNorm() > 0) || !(b[i].squaredNorm() > 0))
      continue;
    s += std::acos(std::clamp(std::abs(a[i].dot(b[i])), real(0), real(1)));
    ++n;
  }
  return n ? s / real(n) : real(M_PI);
}

} // namespace stripes_detail

GAUDI_TEST(stripes_level_lines_follow_a_torus_knot_field) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 48, 24);
  const std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const std::vector<vec3> u = stripes_detail::helix(x, real(1), 5);
  const real l = asawa::shell::avg_length(*M, x);
  const real nu = real(2) * real(M_PI) / (real(8) * l);
  const kusama::stripes::system S =
      kusama::stripes::build(*M, x, stripes_detail::across(*M, x, u), nu);
  vecX psi;
  const real lambda = kusama::stripes::solve(S, psi, 40);
  int skipped = -1;
  const std::vector<vec3> d = kusama::stripes::level_directions(*M, x, S, psi, &skipped);
  const real err = stripes_detail::mean_angle(d, u);
  std::cout << "  stripes lambda " << lambda << " / nu^2 " << lambda / (nu * nu) << "  skipped "
            << skipped << "  mean angle to the field " << err << " rad" << std::endl;
  GAUDI_EXPECT(std::isfinite(lambda) && lambda >= real(0));
  GAUDI_EXPECT(skipped == 0);
  GAUDI_EXPECT(err < real(0.1));
}

/// After a small change of shape a warm start settles in one or two
/// backsolves and keeps its pattern; a cold start may find a neighbor of
/// nearly the same energy (another winding), within a small angle of it.
GAUDI_TEST(stripes_warm_start_keeps_its_pattern) {
  asawa::shell::shell::ptr M = asawa::shell::load_torus(1.0, 0.33, 48, 24);
  std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
  const std::vector<vec3> u = stripes_detail::helix(x, real(1), 5);
  const real l = asawa::shell::avg_length(*M, x);
  const real nu = real(2) * real(M_PI) / (real(8) * l);
  vecX psi;
  kusama::stripes::solve(kusama::stripes::build(*M, x, stripes_detail::across(*M, x, u), nu),
                         psi, 40);
  // A small smooth change of shape, as between two frames.
  for (vec3 &p : x)
    p.z() *= real(1.02);
  const kusama::stripes::system S =
      kusama::stripes::build(*M, x, stripes_detail::across(*M, x, u), nu);
  vecX cold;
  const real l_cold = kusama::stripes::solve(S, cold, 200, 9);
  const std::vector<vec3> d_cold = kusama::stripes::level_directions(*M, x, S, cold);
  real l_first = 0, l_warm = 0, gap = 0;
  for (const int its : {1, 2, 8}) {
    vecX warm = psi;
    l_warm = kusama::stripes::solve(S, warm, its);
    if (its == 1)
      l_first = l_warm;
    gap = stripes_detail::mean_angle(kusama::stripes::level_directions(*M, x, S, warm), d_cold);
    std::cout << "  warm " << its << " backsolves: lambda " << l_warm << " (cold " << l_cold
              << ")  direction gap " << gap << " rad" << std::endl;
  }
  GAUDI_EXPECT(std::abs(l_warm - l_first) < real(2e-3) * l_warm);
  GAUDI_EXPECT(l_cold <= l_warm * (real(1) + real(1e-9)));
  GAUDI_EXPECT(gap < real(0.03));
}

} // namespace test
} // namespace gaudi

#endif
