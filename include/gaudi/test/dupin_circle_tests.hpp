#ifndef GAUDI_TEST_DUPIN_CIRCLE_TESTS_HPP
#define GAUDI_TEST_DUPIN_CIRCLE_TESTS_HPP

/// Tests for the position-space Dupin circle rows (`dupin_circle.hpp`).

#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/dupin_circle.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/test.hpp"

#include <complex>
#include <cstdio>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace dupin_circle_t {

namespace hb = hepworth::block;
using tri_list = std::vector<std::array<int, 3>>;

/// Four points on a circle (centre `c`, radius `r`) in the plane spanned by
/// `u`, `v`, at angles `th`.
inline std::array<vec3, 4> on_circle(const vec3 &c, real r, const vec3 &u, const vec3 &v,
                                     const std::array<real, 4> &th) {
  std::array<vec3, 4> p;
  for (int k = 0; k < 4; ++k)
    p[k] = c + r * (std::cos(th[k]) * u + std::sin(th[k]) * v);
  return p;
}

/// The four points' complex cross ratio in their plane.
inline std::complex<real> cross_ratio(const std::array<vec3, 4> &p) {
  const hb::dupin_circle_detail::fit4 f = hb::dupin_circle_detail::fit(p);
  const vec3 e1 = f.dir, e2 = f.nu.cross(e1);
  std::array<std::complex<real>, 4> z;
  for (int k = 0; k < 4; ++k)
    z[k] = {(p[k] - f.m).dot(e1), (p[k] - f.m).dot(e2)};
  return (z[0] - z[1]) / (z[1] - z[2]) * (z[2] - z[3]) / (z[3] - z[0]);
}

/// A rig of one `circle4` row over four given points.
struct one {
  std::vector<vec3> x, v, f, m;
  hepworth::vec3_block::ptr pos;
  hb::circle4::ptr row;
  vecX q;
  one(const std::array<vec3, 4> &p, const hb::dupin_circle_spec &spec, int sigma, const vec3 &sc = vec3::Zero(),
      real sr = 0, const vec3 &n = vec3::Zero())
      : x(p.begin(), p.end()), v(4, vec3::Zero()), f(v), m(4, vec3::Constant(1)) {
    pos = hepworth::vec3_block::create(m, x, v, f);
    row = hb::circle4::create(std::array<index_t, 4>{0, 1, 2, 3}, x, spec, real(1), sigma, sc, sr, n,
                              std::vector<hepworth::sim_block::ptr>{pos});
    pos->map_to_x(q);
  }
  /// The row's target points (relative to their centroid).
  std::array<vec3, 4> targets() {
    std::vector<trip> T;
    index_t id0 = 0;
    row->fill_A(id0, T);
    vecX p = vecX::Zero(12);
    row->project(q, p);
    std::array<vec3, 4> t;
    for (int k = 0; k < 4; ++k)
      t[k] = p.segment(3 * k, 3);
    return t;
  }
};

inline real angle_on(const hb::dupin_circle_detail::fit4 &f, const vec3 &u, const vec3 &v, const vec3 &q) {
  return std::atan2((q - f.c).dot(v), (q - f.c).dot(u));
}

inline void torus(real R, real r, int nu, int nv, std::vector<vec3> &x, tri_list &tris,
                  std::vector<hb::vertex_chart> &charts) {
  x.clear();
  tris.clear();
  charts.clear();
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      const real u = real(2) * real(M_PI) * i / nu, v = real(2) * real(M_PI) * j / nv;
      const real rho = R + r * std::cos(v);
      x.push_back(vec3(rho * std::cos(u), rho * std::sin(u), r * std::sin(v)));
      hb::vertex_chart c;
      c.N = vec3(std::cos(v) * std::cos(u), std::cos(v) * std::sin(u), std::sin(v));
      c.T = vec3(-std::sin(u), std::cos(u), 0);
      c.B = c.N.cross(c.T);
      // Outward normal, `k > 0` toward the normal: ring `-cos v / rho`, tube `-1 / r`.
      c.h = vec3(-std::cos(v) / rho, 0, -real(1) / r);
      c.ok = true;
      charts.push_back(c);
    }
  auto id = [&](int i, int j) { return ((i % nu + nu) % nu) * nv + ((j % nv + nv) % nv); };
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      tris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      tris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
}

} // namespace dupin_circle_t

GAUDI_TEST(dupin_circle_fit_and_closest_projection) {
  using namespace dupin_circle_t;
  const vec3 u = vec3(1, 2, 0.5).normalized(), v = u.cross(vec3(0.3, -1, 2)).normalized();
  const std::array<vec3, 4> p = on_circle(vec3(0.4, -1, 2), 1.7, u, v, {-0.6, 0.0, 0.5, 1.2});
  one rig(p, hb::dupin_circle_spec{}, 0);
  const std::array<vec3, 4> t = rig.targets();
  const vec3 m = real(0.25) * (p[0] + p[1] + p[2] + p[3]);
  real fixed = 0;
  for (int k = 0; k < 4; ++k)
    fixed = std::max(fixed, (t[k] - (p[k] - m)).norm());
  std::array<vec3, 4> pp = p;
  pp[2] += vec3(0.05, -0.03, 0.04);
  one rig2(pp, hb::dupin_circle_spec{}, 0);
  const std::array<vec3, 4> t2 = rig2.targets();
  const real after = hb::dupin_circle_detail::residual(t2);
  const std::array<vec3, 4> line{vec3(0, 0, 0), vec3(1, 1e-7, 0), vec3(2, 1e-7, 0), vec3(3, 0, 0)};
  const hb::dupin_circle_detail::fit4 fl = hb::dupin_circle_detail::fit(line);
  std::printf("  circle4: points on a circle fixed %.1e, perturbed onto a circle %.1e (from %.1e), line %d\n", fixed,
              after, hb::dupin_circle_detail::residual(pp), int(fl.line));
  GAUDI_EXPECT(fixed < 1e-12);
  GAUDI_EXPECT(after < 1e-12);
  GAUDI_EXPECT(fl.line);
}

GAUDI_TEST(dupin_circle_modes_local_step) {
  using namespace dupin_circle_t;
  using mode = hb::dupin_circle_spec::mode;
  const vec3 c(0.4, -1, 2), u = vec3(1, 2, 0.5).normalized(), v = u.cross(vec3(0.3, -1, 2)).normalized();
  const real r = 1.7;
  const std::array<real, 4> th{-0.6, 0.0, 0.5, 1.2};
  const std::array<vec3, 4> p = on_circle(c, r, u, v, th);
  const real l = ((p[1] - p[0]).norm() + (p[2] - p[1]).norm() + (p[3] - p[2]).norm()) / 3;
  auto seeds_of = [&](hb::dupin_circle_spec s, int sigma) { return one(p, s, sigma).row->seeds(); };
  const hb::dupin_circle_detail::fit4 f = hb::dupin_circle_detail::fit(p);
  vec3 fu = (p[1] - f.c).normalized(), fv = f.nu.cross(fu);
  if (angle_on(f, fu, fv, p[2]) < 0)
    fv = -fv;
  real worst = 0;
  {
    hb::dupin_circle_spec s;
    s.kind = mode::slide;
    s.delta = 0.2;
    const auto sd = seeds_of(s, 1);
    for (int k = 0; k < 4; ++k)
      worst = std::max(worst, std::abs(angle_on(f, fu, fv, sd[k]) - (angle_on(f, fu, fv, p[k]) + 0.2 * l / r)));
  }
  {
    hb::dupin_circle_spec s;
    s.kind = mode::even;
    const auto sd = seeds_of(s, 0);
    const real a0 = angle_on(f, fu, fv, sd[0]), a1 = angle_on(f, fu, fv, sd[1]), a2 = angle_on(f, fu, fv, sd[2]),
               a3 = angle_on(f, fu, fv, sd[3]);
    worst = std::max({worst, std::abs((a1 - a0) - (a2 - a1)), std::abs((a2 - a1) - (a3 - a2)),
                      std::abs(a0 - th[0]), std::abs(a3 - th[3])});
  }
  {
    hb::dupin_circle_spec s;
    s.kind = mode::graded;
    s.g = 2;
    const auto sd = seeds_of(s, 1);
    const real a0 = angle_on(f, fu, fv, sd[0]), a1 = angle_on(f, fu, fv, sd[1]), a2 = angle_on(f, fu, fv, sd[2]),
               a3 = angle_on(f, fu, fv, sd[3]);
    worst = std::max({worst, std::abs((a2 - a1) - 2 * (a1 - a0)), std::abs((a3 - a2) - 4 * (a1 - a0))});
  }
  real cr = 0, moved = 1e9;
  for (mode md : {mode::hyperbolic, mode::parabolic}) {
    hb::dupin_circle_spec s;
    s.kind = md;
    s.lambda = 1.5;
    s.beta = 0.2;
    const auto sd = seeds_of(s, 1);
    cr = std::max(cr, std::abs(cross_ratio(sd) - cross_ratio(p)));
    real mv = 0;
    for (int k = 0; k < 4; ++k)
      mv = std::max(mv, (sd[k] - p[k]).norm());
    moved = std::min(moved, mv);
  }
  // Loxodromic: points on a circle of a sphere through x_i.
  const vec3 C(0.2, 0.1, -0.3);
  const real R = 1.4;
  const vec3 xi = C + R * vec3(0, 0, 1);
  const vec3 nrm(0, 0, 1);
  // A tilted circle on the sphere through xi.
  const vec3 axis = vec3(0.3, 0.8, 0.5).normalized();
  const real h = (xi - C).dot(axis);
  const vec3 cc = C + h * axis;
  const real rc = std::sqrt(R * R - h * h);
  const vec3 cu = (xi - cc).normalized(), cv = axis.cross(cu);
  const std::array<vec3, 4> ps = on_circle(cc, rc, cu, cv, {-0.4, 0.0, 0.35, 0.8});
  hb::dupin_circle_spec sl;
  sl.kind = mode::loxodromic;
  sl.lambda = 1.0;
  sl.psi = 0.3;
  const auto sd = one(ps, sl, 0, C, R, nrm).row->seeds();
  real on_sphere = 0;
  for (const vec3 &q : sd)
    on_sphere = std::max(on_sphere, std::abs((q - C).norm() - R));
  const real fixed_i = (sd[1] - xi).norm();
  const real conc = hb::dupin_circle_detail::residual(sd);
  const hb::dupin_circle_detail::fit4 f0 = hb::dupin_circle_detail::fit(ps), f1 = hb::dupin_circle_detail::fit(sd);
  const vec3 t0 = f0.nu.cross(xi - f0.c).normalized(), t1 = f1.nu.cross(xi - f1.c).normalized();
  const real turn = std::acos(std::clamp(std::abs(t0.dot(t1)), real(0), real(1)));
  // No drift: two local steps at the same configuration agree.
  hb::dupin_circle_spec sd2;
  sd2.kind = mode::slide;
  one rig(p, sd2, 1);
  const auto ta = rig.targets(), tb = rig.targets();
  real drift = 0;
  for (int k = 0; k < 4; ++k)
    drift = std::max(drift, (ta[k] - tb[k]).norm());
  std::printf("  modes: slide / even / graded %.1e, cross ratio kept %.1e (moved %.2f), loxodromic on sphere %.1e "
              "x_i %.1e concyclic %.1e turn %.4f (psi 0.3), drift %.1e\n",
              worst, cr, moved, on_sphere, fixed_i, conc, turn, drift);
  GAUDI_EXPECT(worst < 1e-10);
  GAUDI_EXPECT(cr < 1e-10);
  GAUDI_EXPECT(moved > 1e-3);
  GAUDI_EXPECT(on_sphere < 1e-10);
  GAUDI_EXPECT(fixed_i < 1e-10);
  GAUDI_EXPECT(conc < 1e-10);
  GAUDI_EXPECT(std::abs(turn - 0.3) < 1e-6);
  GAUDI_EXPECT(drift == 0);
}

GAUDI_TEST(dupin_circle_torus_is_dupin_and_noise_falls) {
  using namespace dupin_circle_t;
  std::vector<vec3> x;
  tri_list tris;
  std::vector<hb::vertex_chart> charts;
  torus(2.0, 0.7, 32, 16, x, tris, charts);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  const real clean = hb::dupin_circle_residual(t, x, charts);
  size_t ngroups = 0;
  for (int fam = 1; fam <= 2; ++fam)
    ngroups += hb::select_dupin_groups(t, x, charts, fam, real(25.0 * M_PI / 180.0), real(0.05)).size();
  std::mt19937 rng(7);
  std::normal_distribution<real> nd(0, 0.01);
  std::vector<vec3> xn = x;
  for (vec3 &p : xn)
    p += vec3(nd(rng), nd(rng), nd(rng));
  const real noisy = hb::dupin_circle_residual(t, xn, charts);
  std::vector<vec3> v(xn.size(), vec3::Zero()), f = v, m(xn.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, xn, v, f);
  std::vector<hepworth::projection_constraint::ptr> cs;
  hb::dupin_circle_spec spec;
  const int n = hb::init_dupin_circles(t, xn, charts, {spec}, real(1), nullptr, cs,
                                       std::vector<hepworth::sim_block::ptr>{pos});
  hb::projection_solver solver;
  solver.set_constraints(cs);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  for (int k = 0; k < 3; ++k) {
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
  }
  const real after = hb::dupin_circle_residual(t, xn, charts);
  std::printf("  torus: %zu groups, clean residual %.1e; noisy %.2e -> %.2e after the rows (%d rows)\n", ngroups,
              clean, noisy, after, n);
  GAUDI_EXPECT(ngroups > x.size());
  GAUDI_EXPECT(clean < 1e-12);
  GAUDI_EXPECT(after < 0.5 * noisy);
}

/// The circle condition is scale-free, so the rows alone should not change
/// the enclosed volume (beyond what removing noise does).
GAUDI_TEST(dupin_circle_rows_keep_volume) {
  using namespace dupin_circle_t;
  std::vector<vec3> x;
  tri_list tris;
  std::vector<hb::vertex_chart> charts;
  torus(2.0, 0.7, 32, 16, x, tris, charts);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  auto volume = [&](const std::vector<vec3> &p) {
    real V = 0;
    for (const auto &f : tris)
      V += p[static_cast<size_t>(f[0])].dot(p[static_cast<size_t>(f[1])].cross(p[static_cast<size_t>(f[2])])) / 6;
    return V;
  };
  std::mt19937 rng(11);
  std::normal_distribution<real> nd(0, 0.01);
  std::vector<vec3> xn = x;
  for (vec3 &p : xn)
    p += vec3(nd(rng), nd(rng), nd(rng));
  const real V0 = volume(xn);
  std::vector<vec3> v(xn.size(), vec3::Zero()), f = v, m(xn.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, xn, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  hb::projection_solver solver;
  std::vector<real> ratio;
  for (int k = 0; k < 20; ++k) {
    std::vector<hepworth::projection_constraint::ptr> cs;
    hb::init_dupin_circles(t, xn, charts, {hb::dupin_circle_spec{}}, real(1), nullptr, cs, blocks);
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
    if (k == 0 || k == 4 || k == 19)
      ratio.push_back(volume(xn) / V0);
  }
  std::printf("  circle rows alone, volume / start after 1 / 5 / 20 frames: %.4f / %.4f / %.4f (clean %.4f)\n",
              ratio[0], ratio[1], ratio[2], volume(x) / V0);
  GAUDI_EXPECT(std::abs(ratio[2] - volume(x) / V0) < 0.01);
}

/// Groups that zig-zag sideways about their line (as on meshes whose edges
/// do not follow the curvature lines) must not flatten the line: on a fine
/// torus the ring arcs' sagitta is below the zig-zag, where a least-squares
/// plane of the four points can lie in the tangent plane. In the Python
/// replica (`dupin_circle_study.py`, every grid group at equal weight) that
/// crumples the torus; here the 25 degree alignment test drops the groups
/// where the zig-zag dwarfs the sagitta (the inner rings), and this guards it.
GAUDI_TEST(dupin_circle_rows_keep_volume_under_zigzag) {
  using namespace dupin_circle_t;
  const int nu = 64, nv = 32;
  std::vector<vec3> x;
  tri_list tris;
  std::vector<hb::vertex_chart> charts;
  torus(2.0, 0.7, nu, nv, x, tris, charts);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  auto volume = [&](const std::vector<vec3> &p) {
    real V = 0;
    for (const auto &f : tris)
      V += p[static_cast<size_t>(f[0])].dot(p[static_cast<size_t>(f[1])].cross(p[static_cast<size_t>(f[2])])) / 6;
    return V;
  };
  const real V0 = volume(x);
  // Ring lines offset along the tube direction, alternating by vertex.
  const real h = real(2) * real(M_PI) * real(0.7) / nv, amp = real(0.25) * h;
  std::vector<vec3> xz = x;
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      const size_t k = static_cast<size_t>(i * nv + j);
      xz[k] += (i % 2 == 0 ? amp : -amp) * charts[k].B;
    }
  std::vector<vec3> v(xz.size(), vec3::Zero()), f = v, m(xz.size(), vec3::Constant(1e-3));
  auto pos = hepworth::vec3_block::create(m, xz, v, f);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  hb::projection_solver solver;
  std::vector<hepworth::projection_constraint::ptr> c0;
  const int n_clean = hb::init_dupin_circles(t, x, charts, {hb::dupin_circle_spec{}}, real(1), nullptr, c0, blocks);
  const int n_zig = hb::init_dupin_circles(t, xz, charts, {hb::dupin_circle_spec{}}, real(1), nullptr, c0, blocks);
  std::printf("  groups: clean %d, zig-zag %d\n", n_clean, n_zig);
  for (int k = 0; k < 6; ++k) {
    std::vector<hepworth::projection_constraint::ptr> cs;
    // Rows at 10: well above the mass (`1e-3 / h^2`), as in the demo where
    // `w_dupin_circle` competes with the coupling rather than inertia.
    hb::init_dupin_circles(t, xz, charts, {hb::dupin_circle_spec{}}, real(10), nullptr, cs, blocks);
    solver.set_constraints(cs);
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 15);
  }
  const real ratio = volume(xz) / V0;
  real A0 = 0, A = 0;
  for (const auto &f : tris) {
    const auto area = [&](const std::vector<vec3> &p) {
      return real(0.5) * (p[static_cast<size_t>(f[1])] - p[static_cast<size_t>(f[0])])
                             .cross(p[static_cast<size_t>(f[2])] - p[static_cast<size_t>(f[0])])
                             .norm();
    };
    A0 += area(x);
    A += area(xz);
  }
  std::printf("  zig-zag 0.25 h on a 64 x 32 torus, circle rows at 10, volume / area against clean after 6 frames: "
              "%.4f / %.4f\n",
              ratio, A / A0);
  GAUDI_EXPECT(std::abs(ratio - real(1)) < 0.01);
  GAUDI_EXPECT(std::abs(A / A0 - real(1)) < 0.02);
}

} // namespace test
} // namespace gaudi

#endif
