#ifndef GAUDI_TEST_AMAP_TESTS_HPP
#define GAUDI_TEST_AMAP_TESTS_HPP

/// Tests for the AMAP rows (`amap_gn.hpp`, Vaxman, Mueller, Weber 2015) and
/// the surface-tension forms (`goal_mean_tension`, zero-area rows).

#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/amap_gn.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"
#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/test/test.hpp"

#include <Eigen/SparseCholesky>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace gaudi {
namespace test {
namespace amap_detail_t {

using tri_list = std::vector<std::array<int, 3>>;
namespace hb = hepworth::block;

inline void torus(real R, real r, int nu, int nv, std::vector<vec3> &x, tri_list &tris) {
  x.clear();
  tris.clear();
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      const real u = real(2) * real(M_PI) * i / nu, v = real(2) * real(M_PI) * j / nv;
      const real rho = R + r * std::cos(v);
      x.push_back(vec3(rho * std::cos(u), rho * std::sin(u), r * std::sin(v)));
    }
  auto id = [&](int i, int j) { return ((i % nu + nu) % nu) * nv + ((j % nv + nv) % nv); };
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      tris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
      tris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
    }
}

inline vec3 invert(const vec3 &p, const vec3 &c, real rho) {
  const vec3 d = p - c;
  return c + (rho * rho / d.squaredNorm()) * d;
}

/// The AMAP rows on rest `x0`, positions `x`, quaternions `Y`, and the stacked
/// iterate `q` (positions then quaternions).
struct rig {
  std::vector<vec3> x, v, f, m;
  std::vector<vec4> Y, mY;
  hepworth::vec3_block::ptr pos;
  hepworth::vec4_block::ptr yb;
  hb::amap_gn::ptr row;
  vecX q;
  rig(const kusama::stencil_table &t, const std::vector<vec3> &x0, const std::vector<vec3> &xl, real alpha)
      : x(xl), v(xl.size(), vec3::Zero()), f(v), m(xl.size(), vec3::Constant(1)), Y(xl.size(), vec4(1, 0, 0, 0)),
        mY(xl.size(), vec4::Constant(1)) {
    pos = hepworth::vec3_block::create(m, x, v, f);
    yb = hepworth::vec4_block::create(mY, Y);
    row = hb::amap_gn::create(&t, x0, real(1), alpha, std::vector<hepworth::sim_block::ptr>{pos, yb});
    restack();
  }
  void restack() {
    q = vecX();
    pos->map_to_x(q);
    yb->map_to_x(q);
  }
  /// Gauss-Newton on the quaternions only (positions held): the least AMAP
  /// energy the deformed positions allow.
  real fit_quaternions(int iterations = 20) {
    const index_t ny = static_cast<index_t>(4 * Y.size()), off = yb->_offset;
    for (int it = 0; it < iterations; ++it) {
      std::vector<trip> T;
      index_t id0 = 0;
      row->fill_J(q, id0, T);
      std::vector<trip> TY;
      for (const trip &t : T)
        if (t.col() >= off)
          TY.emplace_back(t.row(), t.col() - off, t.value());
      Eigen::SparseMatrix<real> J(id0, ny);
      J.setFromTriplets(TY.begin(), TY.end());
      vecX g = vecX::Zero(q.size());
      row->add_gradient(q, g);
      Eigen::SparseMatrix<real> H = J.transpose() * J;
      for (index_t k = 0; k < ny; ++k)
        H.coeffRef(k, k) += real(1e-12);
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>> ldlt(H);
      const vecX dy = ldlt.solve(-g.segment(off, ny));
      q.segment(off, ny) += dy;
      if (dy.lpNorm<Eigen::Infinity>() < real(1e-14))
        break;
    }
    return row->residual_sq(q);
  }
};

} // namespace amap_detail_t

GAUDI_TEST(amap_frame_start_is_a_fixed_point) {
  using namespace amap_detail_t;
  std::vector<vec3> x;
  tri_list tris;
  torus(2.0, 0.7, 24, 12, x, tris);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  rig r(t, x, x, real(0.1));
  const real e = r.row->residual_sq(r.q);
  std::printf("  amap at the frame start: %.1e\n", e);
  GAUDI_EXPECT(e < 1e-28);
}

GAUDI_TEST(amap_gradient_matches_differences) {
  using namespace amap_detail_t;
  std::vector<vec3> x;
  tri_list tris;
  torus(2.0, 0.7, 12, 8, x, tris);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  std::mt19937 rng(4);
  std::normal_distribution<real> nd(0, 0.05);
  std::vector<vec3> xl = x;
  for (vec3 &p : xl)
    p += vec3(nd(rng), nd(rng), nd(rng));
  rig r(t, x, xl, real(0.3));
  for (index_t k = r.yb->_offset; k < r.q.size(); ++k)
    r.q[k] += nd(rng);
  vecX g = vecX::Zero(r.q.size());
  r.row->add_gradient(r.q, g);
  real worst = 0;
  const real eps = 1e-6;
  for (index_t k = 0; k < r.q.size(); k += 7) {
    vecX qp = r.q, qm = r.q;
    qp[k] += eps;
    qm[k] -= eps;
    const real fd = (r.row->residual_sq(qp) - r.row->residual_sq(qm)) / (4 * eps);
    worst = std::max(worst, std::abs(fd - g[k]) / std::max(real(1), std::abs(g[k])));
  }
  std::printf("  amap gradient vs differences: %.1e\n", worst);
  GAUDI_EXPECT(worst < 1e-6);
}

/// The defining property: a Moebius image of the mesh (two inversions, so
/// orientation is kept) admits quaternions with zero AMAP energy; an
/// anisotropic stretch does not.
GAUDI_TEST(amap_admits_moebius_maps_exactly) {
  using namespace amap_detail_t;
  std::vector<vec3> x;
  tri_list tris;
  torus(2.0, 0.7, 24, 12, x, tris);
  const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::triangle_walker(x.size(), tris));
  std::vector<vec3> xm(x.size()), xs(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    xm[i] = invert(invert(x[i], vec3(4.0, 1.0, 2.0), 3.0), vec3(-3.0, 0.5, -2.5), 4.0);
    xs[i] = vec3(1.3 * x[i][0], x[i][1], x[i][2]);
  }
  rig rm(t, x, xm, real(0));
  rig rs(t, x, xs, real(0));
  real scale = 0;
  for (size_t i = 0; i < xm.size(); ++i)
    scale += (xm[i] - xm[(i + 1) % xm.size()]).squaredNorm();
  const real em = rm.fit_quaternions(), es = rs.fit_quaternions();
  std::printf("  amap least energy: moebius image %.1e, stretch 1.3 %.1e\n", em, es);
  GAUDI_EXPECT(em < 1e-14);
  GAUDI_EXPECT(es > 1e-4);
}

GAUDI_TEST(tension_goal_moves_the_mean_and_keeps_the_trace_free_part) {
  namespace hb = hepworth::block;
  std::vector<vec3> h{vec3(1.0, 0.2, -0.4), vec3(-0.5, 0.0, 0.7), vec3(2.0, -0.3, 1.0)};
  std::vector<vec3> hv(h.size(), vec3::Zero()), hf = hv, hm(h.size(), vec3::Constant(1));
  auto store = hepworth::vec3_block::create(hm, h, hv, hf);
  vecX q;
  store->map_to_x(q);
  std::vector<vec3> g(h.size());
  hb::goal_mean_tension(store, {hb::tension_target::mode::constant, real(-0.25)})(q, g);
  real worst = 0;
  for (size_t e = 0; e < h.size(); ++e) {
    worst = std::max(worst, std::abs(real(0.5) * (g[e][0] + g[e][2]) + real(0.25)));
    worst = std::max(worst, std::abs((g[e][0] - g[e][2]) - (h[e][0] - h[e][2])) + std::abs(g[e][1] - h[e][1]));
  }
  hb::goal_mean_tension(store, {hb::tension_target::mode::mean, real(0)})(q, g);
  real mean = 0;
  for (const vec3 &v : h)
    mean += real(0.5) * (v[0] + v[2]) / real(h.size());
  for (const vec3 &v : g)
    worst = std::max(worst, std::abs(real(0.5) * (v[0] + v[2]) - mean));
  std::printf("  tension goal: %.1e\n", worst);
  GAUDI_EXPECT(worst < 1e-14);
}

GAUDI_TEST(zero_area_rows_shrink_area) {
  using namespace amap_detail_t;
  std::vector<vec3> x;
  tri_list tris;
  torus(2.0, 0.7, 24, 12, x, tris);
  auto area_of = [&](const std::vector<vec3> &p) {
    real a = 0;
    for (const auto &f : tris)
      a += real(0.5) * (p[static_cast<size_t>(f[1])] - p[static_cast<size_t>(f[0])])
                           .cross(p[static_cast<size_t>(f[2])] - p[static_cast<size_t>(f[0])])
                           .norm();
    return a;
  };
  const real a0 = area_of(x);
  std::vector<vec3> v(x.size(), vec3::Zero()), f = v, m(x.size(), vec3::Constant(1e-2));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  std::vector<hepworth::projection_constraint::ptr> cs;
  for (const auto &tr : tris)
    cs.push_back(hb::area::create(std::vector<index_t>{tr[0], tr[1], tr[2]}, x, real(0), real(0), real(1),
                                  std::vector<hepworth::sim_block::ptr>{pos}, true));
  hb::projection_solver solver;
  solver.set_constraints(cs);
  std::vector<hepworth::sim_block::ptr> blocks{pos};
  for (int k = 0; k < 3; ++k) {
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 10);
  }
  const real a1 = area_of(x);
  std::printf("  zero-area rows: area %.3f -> %.3f\n", a0, a1);
  GAUDI_EXPECT(a1 < 0.95 * a0);
}

} // namespace test
} // namespace gaudi

#endif
