#ifndef GAUDI_TEST_KITE_GROWTH_TESTS_HPP
#define GAUDI_TEST_KITE_GROWTH_TESTS_HPP

/// Tests for the Moebius kite growth (`kite_growth.hpp`).

#include "gaudi/hepworth/block/kite_growth.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/dupin_circle_tests.hpp"
#include "gaudi/test/test.hpp"

#include <cstdio>
#include <vector>

namespace gaudi {
namespace test {
namespace kite_t {

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
  std::array<int, 2> ends(int e) const {
    const int f = t.edge_faces[static_cast<size_t>(e)][0], i = t.edge_local[static_cast<size_t>(e)][0];
    const auto &id = t.ids[static_cast<size_t>(f)];
    return {id[static_cast<size_t>(i)], id[static_cast<size_t>((i + 1) % 3)]};
  }
  vec3 mid(int e) const {
    const auto v = ends(e);
    return real(0.5) * (x[static_cast<size_t>(v[0])] + x[static_cast<size_t>(v[1])]);
  }
};

} // namespace kite_t

/// A constant vector is the gradient of a linear function, whose piecewise
/// linear interpolant reproduces it on every face: the Poisson solve must
/// return that function.
GAUDI_TEST(kite_poisson_recovers_a_gradient) {
  kite_t::rig r;
  const vec3 c(0.3, -0.2, 0.9);
  std::vector<vec3> g(static_cast<size_t>(r.t.n_edges), c);
  const std::vector<real> u = hepworth::block::kite_poisson(r.t, r.x, g);
  real mean = 0;
  for (const vec3 &p : r.x)
    mean += c.dot(p);
  mean /= real(r.x.size());
  real ref = 0;
  for (size_t i = 0; i < r.x.size(); ++i)
    ref = std::max(ref, std::abs(c.dot(r.x[i]) - mean));
  // Up to the constant (area-weighted mean zero against the plain mean).
  real err = 0, shift = 0;
  for (size_t i = 0; i < r.x.size(); ++i)
    shift += u[i] - (c.dot(r.x[i]) - mean);
  shift /= real(r.x.size());
  for (size_t i = 0; i < r.x.size(); ++i)
    err = std::max(err, std::abs(u[i] - (c.dot(r.x[i]) - mean) - shift));
  std::printf("  kite poisson: linear function recovered to %.1e (range %.2f)\n", err, ref);
  GAUDI_EXPECT(err < 1e-6 * ref);
}

/// The ring direction on a torus is divergence free (it follows the
/// rotation): no integrable part, so no growth.
GAUDI_TEST(kite_poisson_drops_the_rotational_part) {
  kite_t::rig r;
  std::vector<vec3> g(static_cast<size_t>(r.t.n_edges));
  for (int e = 0; e < r.t.n_edges; ++e) {
    const vec3 p = r.mid(e);
    g[static_cast<size_t>(e)] = vec3(-p[1], p[0], 0).normalized();
  }
  const std::vector<real> u = hepworth::block::kite_poisson(r.t, r.x, g);
  real mx = 0;
  for (real v : u)
    mx = std::max(mx, std::abs(v));
  std::printf("  kite poisson: ring field (unit, divergence free) gives max |u| %.1e\n", mx);
  GAUDI_EXPECT(mx < 0.05);
}

/// Zero rate: every edge row asks for its frame-start length, the frame
/// start is fixed.
GAUDI_TEST(kite_growth_frame_start_is_fixed_at_zero_rate) {
  kite_t::rig r;
  std::vector<vec3> dirs(static_cast<size_t>(r.t.n_edges), vec3(0.3, -0.2, 0.9)),
      store(static_cast<size_t>(r.t.n_edges), vec3(1, 0, -1));
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::projection_constraint::ptr> cs;
  hepworth::block::kite_target spec;
  spec.alpha = 0;
  hepworth::block::init_kite_growth(r.t, x, dirs, {}, store, spec, real(1), cs, {pos});
  hepworth::block::projection_solver solver;
  solver.set_constraints(cs);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  solver.step(blocks, real(0.1), real(1), 10);
  real moved = 0;
  for (size_t i = 0; i < x.size(); ++i)
    moved = std::max(moved, (x[i] - r.x[i]).norm());
  std::printf("  kite growth at zero rate: moved %.1e\n", moved);
  GAUDI_EXPECT(moved < 1e-10);
}

/// The edge rows realize the discrete conformal scaling: per edge, the log
/// length ratio after a few frames follows `(u_i + u_j) / 2` frame by
/// frame, triangles change shape where their vertices' `u` differ, and the
/// total area holds (mean-zero `u`).
GAUDI_TEST(kite_growth_realizes_the_edge_scaling) {
  kite_t::rig r;
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  hepworth::block::projection_solver solver;
  const std::vector<vec3> store(static_cast<size_t>(r.t.n_edges), vec3(1, 0, -1)); // |z| = 1
  hepworth::block::kite_target spec;
  spec.alpha = 0.05;
  spec.angle = 0;
  spec.source = hepworth::block::kite_source::tangent;
  const int nE = r.t.n_edges;
  auto len = [&](const std::vector<vec3> &p, int e) {
    const auto ab = r.ends(e);
    return (p[static_cast<size_t>(ab[1])] - p[static_cast<size_t>(ab[0])]).norm();
  };
  auto total_area = [&](const std::vector<vec3> &p) {
    real A = 0;
    for (int fi = 0; fi < r.t.n_faces; ++fi) {
      const auto &id = r.t.ids[static_cast<size_t>(fi)];
      A += real(0.5) * (p[static_cast<size_t>(id[1])] - p[static_cast<size_t>(id[0])])
                           .cross(p[static_cast<size_t>(id[2])] - p[static_cast<size_t>(id[0])])
                           .norm();
    }
    return A;
  };
  std::vector<real> L0(static_cast<size_t>(nE)), target(static_cast<size_t>(nE), real(0));
  for (int e = 0; e < nE; ++e)
    L0[static_cast<size_t>(e)] = len(x, e);
  const real total0 = total_area(x);
  for (int frame = 0; frame < 3; ++frame) {
    std::vector<vec3> dirs(static_cast<size_t>(nE), vec3(0.3, -0.2, 0.9)); // has a gradient part
    std::vector<hepworth::projection_constraint::ptr> cs;
    const std::vector<real> u = hepworth::block::init_kite_growth(r.t, x, dirs, {}, store, spec, real(1), cs, blocks);
    for (int e = 0; e < nE; ++e) {
      const auto ab = r.ends(e);
      target[static_cast<size_t>(e)] += real(0.5) * (u[static_cast<size_t>(ab[0])] + u[static_cast<size_t>(ab[1])]);
    }
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
  }
  real sxy = 0, sxx = 0, syy = 0, mx = 0, my = 0;
  std::vector<real> got(static_cast<size_t>(nE));
  for (int e = 0; e < nE; ++e) {
    got[static_cast<size_t>(e)] = std::log(len(x, e) / L0[static_cast<size_t>(e)]);
    mx += got[static_cast<size_t>(e)];
    my += target[static_cast<size_t>(e)];
  }
  mx /= nE;
  my /= nE;
  for (int e = 0; e < nE; ++e) {
    const real a = got[static_cast<size_t>(e)] - mx, b = target[static_cast<size_t>(e)] - my;
    sxy += a * b;
    sxx += a * a;
    syy += b * b;
  }
  const real corr = sxy / std::sqrt(std::max(sxx * syy, real(1e-300)));
  const real gain = sxy / std::max(syy, real(1e-300));
  const real total = total_area(x);
  std::printf("  kite growth: log edge ratio vs (u_i + u_j) / 2 over 3 frames: corr %.3f, gain %.2f; total area %+.2f%%\n",
              corr, gain, 100 * (total / total0 - 1));
  GAUDI_EXPECT(corr > 0.9);
  GAUDI_EXPECT(std::abs(total / total0 - 1) < 0.02);
}

namespace kite_t {

/// The exact torus curvature per edge of the rig, in the edge frame
/// `(t, n x t)`, outward normals (`k > 0` toward the normal), and the ring
/// direction per edge.
inline void torus_store(const rig &r, std::vector<vec3> &store, std::vector<vec3> &ring) {
  const int nE = r.t.n_edges;
  store.resize(static_cast<size_t>(nE));
  ring.resize(static_cast<size_t>(nE));
  for (int e = 0; e < nE; ++e) {
    const vec3 p = r.mid(e);
    const real rho = std::hypot(p[0], p[1]);
    const vec3 er(-p[1] / rho, p[0] / rho, 0);
    const vec3 c(p[0] / rho * 2.0, p[1] / rho * 2.0, 0);
    const vec3 n = (p - c).normalized();
    const vec3 et = n.cross(er);
    const real cv = n.dot(vec3(p[0] / rho, p[1] / rho, 0));
    const real k_ring = -cv / rho, k_tube = -1 / 0.7;
    const auto ab = r.ends(e);
    vec3 te = r.x[static_cast<size_t>(ab[1])] - r.x[static_cast<size_t>(ab[0])];
    te = (te - n * n.dot(te)).normalized();
    const vec3 be = n.cross(te);
    auto kk = [&](const vec3 &a, const vec3 &b) {
      return k_ring * a.dot(er) * b.dot(er) + k_tube * a.dot(et) * b.dot(et);
    };
    store[static_cast<size_t>(e)] = vec3(kk(te, te), kk(te, be), kk(be, be));
    ring[static_cast<size_t>(e)] = er;
  }
}

} // namespace kite_t

/// The mean-curvature source: each element's Moebius map anchored on its
/// mean-curvature sphere, so convex grows. On the torus both sides are
/// convex on average (`-H` 0.90 outside, 0.33 inside): `u` is positive
/// everywhere and largest on the outer equator; `balanced` removes the mean.
/// The ring tangents circulate (no potential) and give almost no `u`.
GAUDI_TEST(kite_mean_curvature_grows_convex) {
  kite_t::rig r;
  std::vector<vec3> store, ring;
  kite_t::torus_store(r, store, ring);
  hepworth::block::kite_target spec;
  spec.alpha = 0.05;
  spec.max_u = 0;
  spec.source = hepworth::block::kite_source::mean_curvature;
  const std::vector<real> u = hepworth::block::kite_log_scale(r.t, r.x, {}, {}, store, spec);
  real u_out = 0, u_in = 0, u_min = 1e30;
  int n_out = 0, n_in = 0;
  for (size_t i = 0; i < r.x.size(); ++i) {
    const real rho = std::hypot(r.x[i][0], r.x[i][1]);
    u_min = std::min(u_min, u[i]);
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
  const std::vector<real> ub = hepworth::block::kite_log_scale(r.t, r.x, {}, {}, store, spec);
  const std::vector<real> a = hepworth::block::kite_detail::vertex_areas(r.t, r.x);
  real mb = 0, sa = 0;
  for (size_t i = 0; i < ub.size(); ++i) {
    mb += a[i] * ub[i];
    sa += a[i];
  }
  mb /= sa;
  spec.balanced = false;
  spec.source = hepworth::block::kite_source::tangent;
  spec.angle = 0;
  const std::vector<real> ut = hepworth::block::kite_log_scale(r.t, r.x, ring, {}, store, spec);
  const auto mm = std::minmax_element(ut.begin(), ut.end());
  std::printf("  kite mean curvature: u outer %.4f, inner %.4f (expect 0.05 x 0.90 / 0.33), min %.4f; balanced mean "
              "%.1e; ring tangents u range %.1e\n",
              u_out, u_in, u_min, mb, *mm.second - *mm.first);
  GAUDI_EXPECT(u_min > 0);
  GAUDI_EXPECT(u_out > 2 * u_in);
  GAUDI_EXPECT(std::abs(u_out - 0.05 * 0.90) < 0.005);
  GAUDI_EXPECT(std::abs(mb) < 1e-12);
  GAUDI_EXPECT(*mm.second - *mm.first < 0.01);
}

/// Driven by the mean curvature, a convex body grows: the edge rows alone
/// raise the torus's volume and area.
GAUDI_TEST(kite_mean_curvature_inflates_volume) {
  kite_t::rig r;
  std::vector<vec3> store, ring;
  kite_t::torus_store(r, store, ring);
  std::vector<vec3> x = r.x, v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  auto volume = [&](const std::vector<vec3> &p) {
    real V = 0;
    for (const auto &t : r.tris)
      V += p[static_cast<size_t>(t[0])].dot(p[static_cast<size_t>(t[1])].cross(p[static_cast<size_t>(t[2])])) / 6;
    return V;
  };
  const real V0 = volume(x);
  hepworth::block::kite_target spec;
  spec.alpha = 0.02;
  spec.source = hepworth::block::kite_source::mean_curvature;
  hepworth::block::projection_solver solver;
  for (int frame = 0; frame < 3; ++frame) {
    std::vector<hepworth::projection_constraint::ptr> cs;
    hepworth::block::init_kite_growth(r.t, x, {}, {}, store, spec, real(1), cs, blocks);
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
  }
  const real V = volume(x);
  std::printf("  kite mean curvature, 3 frames at alpha 0.02: volume %+.2f%%\n", 100 * (V / V0 - 1));
  GAUDI_EXPECT(V > V0 * 1.01);
}

} // namespace test
} // namespace gaudi

#endif
