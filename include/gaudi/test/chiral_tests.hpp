#ifndef GAUDI_TEST_CHIRAL_TESTS_HPP
#define GAUDI_TEST_CHIRAL_TESTS_HPP

/// Tests for the chiral curvature drive (`chiral_drive.hpp`, MATH.md 19.22):
/// the spoke targets carry the target tensor exactly, and on a strip with a
/// held metric a narrow strip twists while a wide one twists less and bends.

#include "gaudi/hepworth/block/chiral_drive.hpp"
#include "gaudi/hepworth/block/contact_net_gn.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/growth_tests.hpp"
#include "gaudi/test/test.hpp"

#include <cstdio>
#include <vector>

namespace gaudi {
namespace test {
namespace chiral_t {

namespace hb = hepworth::block;
namespace cnet = kusama::contact_net;

struct strip_result {
  real twist = 0, H = 0, chord = 0;
};

/// A flat strip `length` by `width` at spacing `h` (x along the length),
/// driven by the held chiral drive (`k0`, frame axis along x) against its
/// frame-0 edge lengths for `frames` frames.
inline strip_result run_strip(real width, real k0, int frames = 8, real length = 2, real h = 0.1) {
  const int ns = static_cast<int>(std::lround(length / h)) + 1;
  const int nz = std::max(3, static_cast<int>(std::lround(width / h)) + 1) | 1;
  std::vector<vec3> x;
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < nz; ++j)
      x.emplace_back(-length / 2 + length * i / (ns - 1), -width / 2 + width * j / (nz - 1), real(0));
  auto id = [&](int i, int j) { return i * nz + j; };
  dupin_circle_t::tri_list tris;
  for (int i = 0; i + 1 < ns; ++i)
    for (int j = 0; j + 1 < nz; ++j) {
      tris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      tris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
  const int n = static_cast<int>(x.size());
  // Break the flat symmetry.
  for (int i = 0; i < n; ++i)
    x[static_cast<size_t>(i)][2] += real(1e-4) * h * std::sin(real(12.9898) * i);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  const cnet::topology tp = cnet::build_topology(t);
  std::vector<real> l0(static_cast<size_t>(t.n_edges), real(0));
  for (int e = 0; e < t.n_edges; ++e) {
    const std::array<int, 2> v = hb::growth_detail::edge_ends(t, e);
    if (v[0] >= 0)
      l0[static_cast<size_t>(e)] = (x[static_cast<size_t>(v[1])] - x[static_cast<size_t>(v[0])]).norm();
  }
  std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Constant(real(1e-4)));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<vec3> gauge(static_cast<size_t>(n), vec3::UnitX());
  hb::projection_solver solver;
  solver.gn_refresh_every = 5;
  solver.gn_backtracks = 8;
  hb::chiral_target spec;
  spec.k0 = k0;
  spec.axis = hb::shear_axis::frame;
  spec.mode = hb::chiral_mode::held;
  std::vector<vec3> href(static_cast<size_t>(n), vec3::Zero());
  std::vector<hb::frame_detail::frame3> F(static_cast<size_t>(n));
  std::vector<hb::vertex_chart> charts(static_cast<size_t>(n));
  for (int k = 0; k < frames; ++k) {
    const std::vector<vec3> N = cnet::vertex_normals(tp, x);
    std::vector<vec4> Y(static_cast<size_t>(n)), ym(static_cast<size_t>(n), vec4::Constant(real(1e-8)));
    for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
      vec3 tt = gauge[i] - N[i] * N[i].dot(gauge[i]);
      tt = tt.norm() > real(1e-9) ? vec3(tt.normalized()) : vec3::UnitX();
      mat3 R;
      R.col(0) = tt;
      R.col(1) = N[i].cross(tt);
      R.col(2) = N[i];
      Y[i] = hb::quat_ops::quat_from_matrix(R);
      F[i] = hb::frame_detail::frame_of(Y[i]);
    }
    for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
      charts[i] = {F[i].T, F[i].B, F[i].N, hb::frame_detail::fit_vertex_tensor(tp, static_cast<int>(i), x, F), true};
      charts[i].ok = charts[i].h.allFinite();
    }
    auto fb = hepworth::vec4_block::create(ym, Y);
    std::vector<hepworth::sim_block::ptr> blocks{pos, fb};
    auto net = hb::contact_net_gn::create(&t, x, N, std::vector<vec3>{}, real(1), blocks);
    hb::frame_rows_config fc;
    fc.w_target = 1;
    auto rows = hb::frame_rows::create(&t, x, Y, fc, blocks);
    std::vector<real> kt, taut;
    hb::chiral_spoke_targets(tp, x, charts, &href, spec, kt, taut);
    rows->set_spoke_targets(kt, taut);
    std::vector<hepworth::projection_constraint::ptr> cs{net, rows};
    for (int e = 0; e < t.n_edges; ++e) {
      const std::array<int, 2> ev = hb::growth_detail::edge_ends(t, e);
      if (ev[0] < 0)
        continue;
      auto row = hb::edge_strain::create(std::vector<index_t>{ev[0], ev[1]}, x, real(1),
                                         std::vector<hepworth::sim_block::ptr>{pos});
      row->set_goal_length(l0[static_cast<size_t>(e)]);
      cs.push_back(row);
    }
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 30);
    for (size_t i = 0; i < static_cast<size_t>(n); ++i)
      gauge[i] = hb::quat_ops::qrot(Y[i], vec3::UnitX());
  }
  // Read-outs along the centerline row and over the interior.
  strip_result r;
  const std::vector<vec3> N = cnet::vertex_normals(tp, x);
  real arc = 0, tw = 0;
  for (int i = 0; i + 1 < ns; ++i) {
    const size_t a = static_cast<size_t>(id(i, nz / 2)), b = static_cast<size_t>(id(i + 1, nz / 2));
    const vec3 d = x[b] - x[a];
    arc += d.norm();
    tw += N[a].cross(N[b]).dot(d.normalized());
  }
  r.twist = tw / arc;
  r.chord = (x[static_cast<size_t>(id(ns - 1, nz / 2))] - x[static_cast<size_t>(id(0, nz / 2))]).norm() / arc;
  for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
    vec3 tt = gauge[i] - N[i] * N[i].dot(gauge[i]);
    tt.normalize();
    F[i] = {tt, N[i].cross(tt), N[i]};
  }
  int nh = 0;
  for (int i = 2; i + 2 < ns; ++i)
    for (int j = 1; j + 1 < nz; ++j) {
      const vec3 s = hb::frame_detail::fit_vertex_tensor(tp, id(i, j), x, F);
      r.H += std::abs(real(0.5) * (s[0] + s[2]));
      ++nh;
    }
  r.H /= std::max(nh, 1);
  return r;
}

} // namespace chiral_t

/// The spoke targets are the target tensor's readings: the per-vertex
/// least-squares fit of a ring's target `(k, tau)` gives back
/// `h + k0 (e+ e+^T - e- e-^T)` exactly (analytic torus charts).
GAUDI_TEST(chiral_targets_fit_back_to_the_tensor) {
  namespace hb = hepworth::block;
  growth_t::rig r;
  const kusama::contact_net::topology tp = kusama::contact_net::build_topology(r.t);
  hb::chiral_target spec;
  spec.k0 = 0.7;
  spec.angle = 0.3;
  spec.axis = hb::shear_axis::principal;
  const std::vector<vec3> b = hb::chiral_tensors(r.charts, nullptr, spec);
  std::vector<real> k, tau;
  hb::chiral_spoke_targets(tp, r.x, r.charts, nullptr, spec, k, tau);
  real err = 0, mag = 0;
  for (size_t v = 0; v < r.x.size(); ++v) {
    const hb::vertex_chart &c = r.charts[v];
    if (!c.ok)
      continue;
    Eigen::Matrix<real, 3, 3> M = Eigen::Matrix<real, 3, 3>::Zero();
    vec3 rhs = vec3::Zero();
    for (int s = tp.spoke_begin[v]; s < tp.spoke_begin[v + 1]; ++s) {
      const size_t j = static_cast<size_t>(tp.spokes[static_cast<size_t>(s)][1]);
      const vec3 dx = r.x[j] - r.x[v], nb = (c.N + r.charts[j].N).normalized();
      const vec3 dxt = dx - nb * dx.dot(nb);
      real a = dxt.dot(c.T), cc = dxt.dot(c.B);
      const real rr = std::hypot(a, cc);
      a /= rr;
      cc /= rr;
      const vec3 r1(a * a, 2 * a * cc, cc * cc), r2(-a * cc, a * a - cc * cc, a * cc);
      M += r1 * r1.transpose() + r2 * r2.transpose();
      rhs += r1 * k[static_cast<size_t>(s)] + r2 * tau[static_cast<size_t>(s)];
    }
    const vec3 fit = M.ldlt().solve(rhs);
    err = std::max(err, (fit - b[v]).norm());
    mag = std::max(mag, b[v].norm());
  }
  std::printf("  chiral targets: fit back to the tensor, max error %.2e (|b| up to %.2f)\n", err, mag);
  GAUDI_EXPECT(err < real(1e-9) * mag);
}

/// Armon's incompatible strip on the E3 net: under the held drive with a
/// frame-0 metric a narrow strip twists at nearly `k0` and stays nearly
/// straight with `H ~ 0`; a wide one twists less and bends more (20 frames:
/// twist 0.82 / 0.60 of `k0`, `|H|` 0.04 / 0.14 of `k0`; env `CHIRAL_K0`,
/// `CHIRAL_FRAMES`, `CHIRAL_WIDE`).
GAUDI_TEST(chiral_strip_twists_when_narrow) {
  const char *ek = std::getenv("CHIRAL_K0"), *ef = std::getenv("CHIRAL_FRAMES"), *ew = std::getenv("CHIRAL_WIDE");
  const real k0 = ek ? real(std::atof(ek)) : real(2);
  const int frames = ef ? std::atoi(ef) : 20;
  const real wide = ew ? real(std::atof(ew)) : real(1.0);
  const chiral_t::strip_result a = chiral_t::run_strip(0.2, k0, frames), b = chiral_t::run_strip(wide, k0, frames);
  std::printf("  chiral strip k0 %.1f: narrow twist/k0 %.3f chord %.3f |H|/k0 %.3f; wide twist/k0 %.3f chord %.3f "
              "|H|/k0 %.3f\n",
              k0, a.twist / k0, a.chord, a.H / k0, b.twist / k0, b.chord, b.H / k0);
  GAUDI_EXPECT(std::abs(a.twist) > real(0.7) * k0);
  GAUDI_EXPECT(a.chord > real(0.9));
  GAUDI_EXPECT(a.H < real(0.1) * k0);
  GAUDI_EXPECT(std::abs(b.twist) < std::abs(a.twist));
  GAUDI_EXPECT(b.H > real(2) * a.H);
}

} // namespace test
} // namespace gaudi

#endif
