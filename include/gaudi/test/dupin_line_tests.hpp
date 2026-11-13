#ifndef GAUDI_TEST_DUPIN_LINE_TESTS_HPP
#define GAUDI_TEST_DUPIN_LINE_TESTS_HPP

/// Tests for the straightest walk (`surface_walk.hpp`) and the long-stencil
/// Dupin store rows (`edge_fair_rows::dupin_line`).

#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/surface_walk.hpp"
#include "gaudi/hepworth/block/variational_audit.hpp"
#include "gaudi/kusama/weingarten.hpp"
#include "gaudi/test/dupin_circle_tests.hpp"
#include "gaudi/test/test.hpp"

#include <cstdio>
#include <vector>

namespace gaudi {
namespace test {
namespace dupin_line_t {

namespace hb = hepworth::block;
namespace wg = kusama::weingarten;

struct torus_mesh {
  std::vector<vec3> x;
  dupin_circle_t::tri_list tris;
  kusama::stencil_table t;
  torus_mesh(real R, real r, int nu, int nv, real bump = 0) {
    std::vector<hb::vertex_chart> charts;
    dupin_circle_t::torus(R, r, nu, nv, x, tris, charts);
    if (bump != real(0))
      for (vec3 &p : x) {
        const real phi = std::atan2(p[1], p[0]);
        const real rho = std::hypot(p[0], p[1]) - R;
        const real theta = std::atan2(p[2], rho);
        const real rt = r * (real(1) + bump * std::sin(real(3) * phi) * std::cos(theta));
        p = vec3((R + rt * std::cos(theta)) * std::cos(phi), (R + rt * std::cos(theta)) * std::sin(phi),
                 rt * std::sin(theta));
      }
    t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  }
};

/// `|A h| / |h|` of a row set on the store `h`.
inline real rows_ratio(const kusama::stencil_table &t, const std::vector<vec3> &x, std::vector<vec3> h,
                       hb::edge_fair_rows rows, hb::dupin_line_spec line = {}) {
  const std::vector<vec3> h0 = h;
  std::vector<vec3> hm(h.size(), vec3::Ones()), hv(h.size(), vec3::Zero()), hf = hv;
  auto blk = hepworth::vec3_block::create(hm, h, hv, hf);
  auto c = hb::edge_shape_fair::create(&t, x, rows, real(1), std::vector<hepworth::sim_block::ptr>{blk},
                                       hb::goal_zero(), real(1), &h0, false, real(0), line);
  vecX q;
  blk->map_to_x(q);
  index_t id0 = 0;
  std::vector<trip> trips;
  c->fill_A(id0, trips);
  vecX Aq = vecX::Zero(id0);
  for (const trip &tr : trips)
    Aq[tr.row()] += tr.value() * q[tr.col()];
  return Aq.norm() / std::max(q.norm(), real(1e-30));
}

} // namespace dupin_line_t

/// A meridian of a torus is a geodesic: the straightest walk along it from a
/// face centroid comes back near its start after one turn, and every sample
/// stays on the (faceted) surface.
GAUDI_TEST(surface_walk_closes_a_meridian) {
  using namespace dupin_line_t;
  const real R = 2.0, r = 0.7;
  torus_mesh m(R, r, 48, 24);
  const int f = 0;
  const auto &id = m.t.ids[static_cast<size_t>(f)];
  const vec3 c = (m.x[static_cast<size_t>(id[0])] + m.x[static_cast<size_t>(id[1])] + m.x[static_cast<size_t>(id[2])]) / 3;
  const real phi = std::atan2(c[1], c[0]), theta = std::atan2(c[2], std::hypot(c[0], c[1]) - R);
  const vec3 e_tube(-std::sin(theta) * std::cos(phi), -std::sin(theta) * std::sin(phi), std::cos(theta));
  const real L = real(2) * real(M_PI) * r;
  const std::vector<hb::walk_sample> w = hb::trace_on_surface(m.t, m.x, f, c, e_tube, L, real(0.05));
  real off = 0;
  for (const hb::walk_sample &s : w) {
    const real d = std::hypot(std::hypot(s.p[0], s.p[1]) - R, s.p[2]) - r;
    off = std::max(off, std::abs(d));
  }
  const real gap = w.empty() ? real(1e9) : (w.back().p - c).norm();
  std::printf("  surface walk: meridian of %.3f in %zu samples, end gap %.4f, max off-surface %.4f\n", L, w.size(),
              gap, off);
  GAUDI_EXPECT(w.size() > 80);
  GAUDI_EXPECT(gap < real(0.05) * r);
  GAUDI_EXPECT(off < real(0.02));
}

/// The long-stencil Dupin rows vanish (to discretization) on a torus and not
/// on a bumpy one; with a longer reach the bumpy torus stands out more.
GAUDI_TEST(dupin_line_rows_hold_a_torus) {
  using namespace dupin_line_t;
  torus_mesh clean(2.0, 0.7, 48, 24), bumpy(2.0, 0.7, 48, 24, real(0.2));
  const std::vector<vec3> hc = wg::edge_readings(clean.t, clean.x), hb_ = wg::edge_readings(bumpy.t, bumpy.x);
  const real face_c = rows_ratio(clean.t, clean.x, hc, hb::edge_fair_rows::dupin);
  const real face_b = rows_ratio(bumpy.t, bumpy.x, hb_, hb::edge_fair_rows::dupin);
  std::printf("  dupin face rows: torus %.2e, bumpy %.2e (ratio %.3f)\n", face_c, face_b, face_c / face_b);
  for (const real reach : {real(1), real(3), real(6)}) {
    const hb::dupin_line_spec spec{reach, 0};
    const real rc = rows_ratio(clean.t, clean.x, hc, hb::edge_fair_rows::dupin_line, spec);
    const real rb = rows_ratio(bumpy.t, bumpy.x, hb_, hb::edge_fair_rows::dupin_line, spec);
    std::printf("  dupin line rows, reach %.0f: torus %.2e, bumpy %.2e (ratio %.3f)\n", reach, rc, rb, rc / rb);
    GAUDI_EXPECT(rc < real(0.25) * rb);
  }
}

/// A noisy torus store relaxed under the rows (pin to the noisy store at a
/// low weight plus the Dupin rows): the long stencil recovers the clean
/// store at least as well as the face rows, and the rows are variational.
GAUDI_TEST(dupin_line_rows_recover_a_noisy_store) {
  using namespace dupin_line_t;
  torus_mesh m(2.0, 0.7, 48, 24);
  const std::vector<vec3> clean = wg::edge_readings(m.t, m.x);
  std::mt19937 rng(5);
  real rms = 0;
  for (const vec3 &h : clean)
    rms += h.squaredNorm();
  rms = std::sqrt(rms / clean.size());
  std::normal_distribution<real> nd(0, 0.1 * rms);
  std::vector<vec3> noisy = clean;
  for (vec3 &h : noisy)
    h += vec3(nd(rng), nd(rng), nd(rng));
  auto relax = [&](hb::edge_fair_rows rows, hb::dupin_line_spec spec, real *audit) {
    std::vector<vec3> h = noisy, hm(h.size(), vec3::Constant(1e-6)), hv(h.size(), vec3::Zero()), hf = hv;
    auto blk = hepworth::vec3_block::create(hm, h, hv, hf);
    std::vector<hepworth::sim_block::ptr> blocks{blk};
    std::vector<hepworth::projection_constraint::ptr> cs;
    cs.push_back(hb::edge_shape_fair::create(&m.t, m.x, hb::edge_fair_rows::pin, real(0.1), blocks,
                                             hb::goal_fixed(noisy)));
    cs.push_back(hb::edge_shape_fair::create(&m.t, m.x, rows, real(1), blocks, hb::goal_zero(), real(1), &noisy,
                                             false, real(0), spec));
    if (audit) {
      vecX q;
      blk->map_to_x(q);
      const auto table = hb::audit_variational(std::vector<hepworth::projection_constraint::ptr>{cs[1]}, q);
      *audit = table.empty() ? real(0) : table.front().max;
    }
    hb::projection_solver solver;
    solver.set_constraints(cs);
    solver.step(blocks, real(0.1), real(1), 5);
    real err = 0, ref = 0;
    for (size_t e = 0; e < h.size(); ++e) {
      err += (h[e] - clean[e]).squaredNorm();
      ref += (noisy[e] - clean[e]).squaredNorm();
    }
    return std::sqrt(err / ref);
  };
  real audit = 0;
  const real face = relax(hb::edge_fair_rows::dupin, {}, nullptr);
  const real line1 = relax(hb::edge_fair_rows::dupin_line, {1, 0}, nullptr);
  const real line3 = relax(hb::edge_fair_rows::dupin_line, {3, 0}, &audit);
  std::printf("  noisy torus store, error / noise after the rows: face %.3f, line reach 1 %.3f, reach 3 %.3f; audit %.1e\n",
              face, line1, line3, audit);
  GAUDI_EXPECT(line3 < real(1));
  GAUDI_EXPECT(line3 <= face * real(1.05));
  GAUDI_EXPECT(audit < real(1e-4));
}

} // namespace test
} // namespace gaudi

#endif
