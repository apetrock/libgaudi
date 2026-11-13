#ifndef GAUDI_TEST_GROWTH_TESTS_HPP
#define GAUDI_TEST_GROWTH_TESTS_HPP

/// Tests for the Moebius growth (`mobius_growth.hpp`) and the flat-direction
/// growth (`flat_growth.hpp`) on the analytic torus charts.

#include "gaudi/hepworth/block/flat_growth.hpp"
#include "gaudi/hepworth/block/mobius_growth.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/dupin_circle_tests.hpp"
#include "gaudi/test/test.hpp"

#include <cstdio>
#include <vector>

namespace gaudi {
namespace test {
namespace growth_t {

namespace hb = hepworth::block;

struct rig {
  std::vector<vec3> x;
  dupin_circle_t::tri_list tris;
  std::vector<hb::vertex_chart> charts;
  kusama::stencil_table t;
  rig(int nu = 32, int nv = 16) {
    dupin_circle_t::torus(2.0, 0.7, nu, nv, x, tris, charts);
    t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  }
  real volume(const std::vector<vec3> &p) const {
    real V = 0;
    for (const auto &tr : tris)
      V += p[static_cast<size_t>(tr[0])].dot(p[static_cast<size_t>(tr[1])].cross(p[static_cast<size_t>(tr[2])])) / 6;
    return V;
  }
};

} // namespace growth_t

/// Mean-curvature source: convex grows (`-H` 0.90 on the outer equator,
/// 0.33 on the inner), `balanced` removes the mean; the flat source is 1 at
/// the top of the tube (ring curvature 0) and lower on the equators.
GAUDI_TEST(mobius_sources_on_the_torus) {
  growth_t::rig r;
  hepworth::block::mobius_target spec;
  spec.alpha = 0.05;
  spec.max_u = 0;
  spec.source = hepworth::block::mobius_source::mean_curvature;
  const std::vector<real> u = hepworth::block::mobius_log_scale(r.t, r.x, r.charts, spec);
  real u_out = 0, u_in = 0;
  int n_out = 0, n_in = 0;
  for (size_t i = 0; i < r.x.size(); ++i) {
    const real rho = std::hypot(r.x[i][0], r.x[i][1]);
    if (rho > 2.65) {
      u_out += u[i];
      ++n_out;
    } else if (rho < 1.35) {
      u_in += u[i];
      ++n_in;
    }
  }
  u_out /= std::max(n_out, 1);
  u_in /= std::max(n_in, 1);
  spec.balanced = true;
  const std::vector<real> ub = hepworth::block::mobius_log_scale(r.t, r.x, r.charts, spec);
  const std::vector<real> a = hepworth::block::growth_detail::vertex_areas(r.t, r.x);
  real mb = 0, sa = 0;
  for (size_t i = 0; i < ub.size(); ++i) {
    mb += a[i] * ub[i];
    sa += a[i];
  }
  mb /= sa;
  spec.balanced = false;
  spec.alpha = 1;
  spec.source = hepworth::block::mobius_source::flat;
  const std::vector<real> uf = hepworth::block::mobius_log_scale(r.t, r.x, r.charts, spec);
  real top = 0, eq = 0;
  for (size_t i = 0; i < r.x.size(); ++i) {
    const real rho = std::hypot(r.x[i][0], r.x[i][1]);
    if (std::abs(r.x[i][2]) > 0.69)
      top = std::max(top, uf[i]);
    if (rho > 2.65)
      eq = std::max(eq, uf[i]);
  }
  std::printf("  mobius mean curvature: u outer %.4f inner %.4f, balanced mean %.1e; flat source top %.3f outer %.3f\n",
              u_out, u_in, mb, top, eq);
  GAUDI_EXPECT(u_out > 2 * u_in);
  GAUDI_EXPECT(std::abs(u_out - 0.05 * 0.90) < 0.005);
  GAUDI_EXPECT(std::abs(mb) < 1e-12);
  GAUDI_EXPECT(top > 0.99);
  GAUDI_EXPECT(eq < 0.8);
}

/// Driven by the mean curvature, a convex body grows: the edge rows alone
/// raise the torus's volume.
GAUDI_TEST(mobius_mean_curvature_inflates_volume) {
  growth_t::rig r;
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  const real V0 = r.volume(x);
  hepworth::block::mobius_target spec;
  spec.alpha = 0.02;
  spec.source = hepworth::block::mobius_source::mean_curvature;
  hepworth::block::projection_solver solver;
  for (int frame = 0; frame < 3; ++frame) {
    std::vector<hepworth::projection_constraint::ptr> cs;
    hepworth::block::init_mobius_growth(r.t, x, r.charts, spec, real(1), cs, blocks);
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
  }
  const real V = r.volume(x);
  std::printf("  mobius mean curvature, 3 frames at alpha 0.02: volume %+.2f%%\n", 100 * (V / V0 - 1));
  GAUDI_EXPECT(V > V0 * 1.01);
}

/// Directional growth on flatness (power 8): on the torus the loose
/// direction is the hoop; the top of the tube (hoop curvature 0) gets far
/// more weight than the outer equator, and the realized hoop edges lengthen
/// accordingly.
GAUDI_TEST(flat_growth_lengthens_the_top_hoop) {
  growth_t::rig r;
  hepworth::block::flat_growth_target spec;
  spec.g = 0.1;
  spec.power = 8;
  std::vector<real> w;
  std::vector<vec3> dir;
  hepworth::block::flat_growth_field(r.charts, spec, w, dir);
  real w_top = 0, w_eq = 1, align = 1;
  for (size_t i = 0; i < r.x.size(); ++i) {
    const real rho = std::hypot(r.x[i][0], r.x[i][1]);
    const vec3 hoop(-r.x[i][1] / rho, r.x[i][0] / rho, 0);
    align = std::min(align, std::abs(dir[i].dot(hoop)));
    if (std::abs(r.x[i][2]) > 0.69)
      w_top = std::max(w_top, w[i]);
    if (rho > 2.65)
      w_eq = std::min(w_eq, w[i]);
  }
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  std::vector<hepworth::projection_constraint::ptr> cs;
  hepworth::block::init_flat_growth(r.t, x, r.charts, spec, real(1), cs, blocks);
  hepworth::block::projection_solver solver;
  solver.set_constraints(cs);
  solver.step(blocks, real(0.1), real(1), 30);
  // Hoop edge length at the top row and on the outer equator, before and after.
  auto hoop_len = [&](const std::vector<vec3> &p, bool top) {
    real s = 0;
    int c = 0;
    for (int e = 0; e < r.t.n_edges; ++e) {
      const auto ab = hepworth::block::growth_detail::edge_ends(r.t, e);
      const vec3 a = r.x[static_cast<size_t>(ab[0])], b = r.x[static_cast<size_t>(ab[1])];
      if (std::abs(a[2] - b[2]) > 1e-9)
        continue;
      const bool is_top = std::abs(a[2]) > 0.69, is_eq = std::hypot(a[0], a[1]) > 2.65;
      if ((top && !is_top) || (!top && !is_eq))
        continue;
      s += (p[static_cast<size_t>(ab[1])] - p[static_cast<size_t>(ab[0])]).norm();
      ++c;
    }
    return s / std::max(c, 1);
  };
  const real gt = hoop_len(x, true) / hoop_len(r.x, true) - 1, ge = hoop_len(x, false) / hoop_len(r.x, false) - 1;
  std::printf("  flat growth: weight top %.3f outer %.3f, hoop alignment %.3f; hoop growth top %+.3f outer %+.3f\n",
              w_top, w_eq, align, gt, ge);
  GAUDI_EXPECT(w_top > 5 * w_eq);
  GAUDI_EXPECT(align > 0.99);
  // One short solve with edge rows only: the top's extra length spreads into
  // the connected circles, so the realized differential is smaller.
  GAUDI_EXPECT(gt > real(1.3) * ge);
}

/// The growth axis turned a quarter (`angle` pi / 2) is the tight direction:
/// on the torus the tube circle, everywhere.
GAUDI_TEST(flat_growth_angle_turns_the_axis) {
  growth_t::rig r;
  hepworth::block::flat_growth_target spec;
  spec.power = 0;
  spec.angle = real(0.5 * M_PI);
  std::vector<real> w;
  std::vector<vec3> dir;
  hepworth::block::flat_growth_field(r.charts, spec, w, dir);
  real align = 1;
  for (size_t i = 0; i < r.x.size(); ++i) {
    const real rho = std::hypot(r.x[i][0], r.x[i][1]);
    const vec3 hoop(-r.x[i][1] / rho, r.x[i][0] / rho, 0);
    const vec3 tube = r.charts[i].N.cross(hoop);
    if (w[i] > 0)
      align = std::min(align, std::abs(dir[i].dot(tube)));
  }
  std::printf("  flat growth at angle pi/2: min alignment with the tube circle %.3f\n", align);
  GAUDI_EXPECT(align > 0.99);
}

/// `iso = 1, power = 0`: every edge grows by the same fraction (a uniform
/// scale growth).
GAUDI_TEST(flat_growth_iso_is_uniform) {
  growth_t::rig r;
  hepworth::block::flat_growth_target spec;
  spec.g = 0.1;
  spec.power = 0;
  spec.iso = 1;
  std::vector<hepworth::projection_constraint::ptr> cs;
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-8));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  hepworth::block::init_flat_growth(r.t, x, r.charts, spec, real(1), cs, blocks);
  hepworth::block::projection_solver solver;
  solver.set_constraints(cs);
  solver.step(blocks, real(0.1), real(1), 30);
  real lo = 1e30, hi = -1e30;
  for (int e = 0; e < r.t.n_edges; ++e) {
    const auto ab = hepworth::block::growth_detail::edge_ends(r.t, e);
    const real g = (x[static_cast<size_t>(ab[1])] - x[static_cast<size_t>(ab[0])]).norm() /
                       (r.x[static_cast<size_t>(ab[1])] - r.x[static_cast<size_t>(ab[0])]).norm() -
                   1;
    lo = std::min(lo, g);
    hi = std::max(hi, g);
  }
  std::printf("  iso growth: edge growth %.4f .. %.4f (asked 0.1)\n", lo, hi);
  GAUDI_EXPECT(std::abs(lo - 0.1) < 0.01);
  GAUDI_EXPECT(std::abs(hi - 0.1) < 0.01);
}

} // namespace test
} // namespace gaudi

#endif
