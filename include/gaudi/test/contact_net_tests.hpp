#ifndef GAUDI_TEST_CONTACT_NET_TESTS_HPP
#define GAUDI_TEST_CONTACT_NET_TESTS_HPP

/// The contact element net (kusama/contact_net.hpp; scripts/sympy/cr_poisson/
/// MATH.md 19.10 to 19.12): read-out, generated rows, fixed points, and
/// convergence on the shared suite (the Python `lie_compare.cases()` meshes).

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/stencil_walker.hpp"
#include "gaudi/hepworth/block/contact_net_gn.hpp"
#include "gaudi/hepworth/block/frame_rows.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/kusama/contact_net.hpp"
#include "gaudi/test/test.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace gaudi {
namespace test {
namespace contact_net_detail {

namespace cnet = kusama::contact_net;
using tri_list = std::vector<std::array<int, 3>>;

/// `cr_ops.torus`: vertex `i nv + j` at spine angle `2 pi i / nu`, tube angle
/// `2 pi j / nv`.
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

/// `cr_ops.strip_param` + `flat_from_param`: the grid in `(s, z)`, flat in
/// `y = 0`; `P` the parameters.
inline void strip(int ns, int nz, real w, real h, std::vector<vec3> &x, std::vector<std::array<real, 2>> &P,
                  tri_list &tris) {
  x.clear();
  P.clear();
  tris.clear();
  for (int i = 0; i < ns; ++i)
    for (int j = 0; j < nz; ++j) {
      const real s = -w / 2 + w * i / (ns - 1), z = -h / 2 + h * j / (nz - 1);
      P.push_back({s, z});
      x.push_back(vec3(s, 0, z));
    }
  for (int i = 0; i + 1 < ns; ++i)
    for (int j = 0; j + 1 < nz; ++j) {
      const int a = i * nz + j, b = (i + 1) * nz + j;
      tris.push_back({a, b, b + 1});
      tris.push_back({a, b + 1, a + 1});
    }
}

/// `closed_targets.cube`: `[-1, 1]^3`, `k x k` quads per side, outward.
inline void cube(int k, std::vector<vec3> &x, tri_list &tris) {
  x.clear();
  tris.clear();
  std::map<std::array<long, 3>, int> index;
  auto vid = [&](const vec3 &p) {
    const std::array<long, 3> key{std::lround(p[0] * 1e6), std::lround(p[1] * 1e6), std::lround(p[2] * 1e6)};
    auto it = index.find(key);
    if (it != index.end())
      return it->second;
    index[key] = static_cast<int>(x.size());
    x.push_back(p);
    return static_cast<int>(x.size()) - 1;
  };
  for (int axis = 0; axis < 3; ++axis)
    for (real side : {real(-1), real(1)}) {
      const int u = (axis + 1) % 3, v = (axis + 2) % 3;
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) {
          std::array<int, 4> q{};
          const int di[4] = {0, 1, 1, 0}, dj[4] = {0, 0, 1, 1};
          for (int c = 0; c < 4; ++c) {
            vec3 p = vec3::Zero();
            p[axis] = side;
            p[u] = real(-1) + real(2) * (i + di[c]) / k;
            p[v] = real(-1) + real(2) * (j + dj[c]) / k;
            q[static_cast<size_t>(c)] = vid(p);
          }
          tris.push_back({q[0], q[1], q[2]});
          tris.push_back({q[0], q[2], q[3]});
        }
    }
  for (auto &f : tris) {
    const vec3 n = (x[static_cast<size_t>(f[1])] - x[static_cast<size_t>(f[0])])
                       .cross(x[static_cast<size_t>(f[2])] - x[static_cast<size_t>(f[0])]);
    if (n.dot(x[static_cast<size_t>(f[0])] + x[static_cast<size_t>(f[1])] + x[static_cast<size_t>(f[2])]) < 0)
      std::swap(f[1], f[2]);
  }
}

/// `closed_targets.radial_ellipsoid`: each point pushed radially onto the
/// ellipsoid with semi-axes `abc`.
inline std::vector<vec3> radial_ellipsoid(const std::vector<vec3> &x, const vec3 &abc) {
  std::vector<vec3> y;
  for (const vec3 &p : x) {
    const real s = std::sqrt(std::pow(p[0] / abc[0], 2) + std::pow(p[1] / abc[1], 2) + std::pow(p[2] / abc[2], 2));
    y.push_back(p / s);
  }
  return y;
}

inline kusama::stencil_table table_of(int n, const tri_list &tris) {
  return kusama::build_stencil_table(asawa::shell::triangle_walker(n, tris));
}

/// Max distance after the best rigid alignment (Kabsch), relative to the
/// target's size.
inline real shape_error(const std::vector<vec3> &x, const std::vector<vec3> &y) {
  vec3 cx = vec3::Zero(), cy = vec3::Zero();
  for (size_t i = 0; i < x.size(); ++i) {
    cx += x[i];
    cy += y[i];
  }
  cx /= real(x.size());
  cy /= real(y.size());
  mat3 C = mat3::Zero();
  vec3 lo = y[0], hi = y[0];
  for (size_t i = 0; i < x.size(); ++i) {
    C += (y[i] - cy) * (x[i] - cx).transpose();
    lo = lo.cwiseMin(y[i]);
    hi = hi.cwiseMax(y[i]);
  }
  Eigen::JacobiSVD<mat3> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
  mat3 U = svd.matrixU();
  if ((U * svd.matrixV().transpose()).determinant() < 0)
    U.col(2) *= -1;
  const mat3 Rm = U * svd.matrixV().transpose();
  real e = 0;
  for (size_t i = 0; i < x.size(); ++i)
    e = std::max(e, (Rm * (x[i] - cx) + cy - y[i]).norm());
  return e / (hi - lo).maxCoeff();
}

struct test_case {
  std::string name;
  std::vector<vec3> x0, xt;
  tri_list tris;
};

/// The Python shared suite (`lie_compare.cases()`), minus the isometric
/// cylinder (its target needs the run_goal helper) and plus it as the plain
/// isometric wrap.
inline std::vector<test_case> suite() {
  std::vector<test_case> out;
  {
    test_case c{"torus tube 0.7->0.9"};
    torus(2.0, 0.7, 24, 12, c.x0, c.tris);
    tri_list t;
    torus(2.0, 0.9, 24, 12, c.xt, t);
    out.push_back(c);
  }
  {
    test_case c{"torus R 2->2.5"};
    torus(2.0, 0.7, 24, 12, c.x0, c.tris);
    tri_list t;
    torus(2.5, 0.7, 24, 12, c.xt, t);
    out.push_back(c);
  }
  {
    test_case c{"cube->sphere"};
    cube(6, c.x0, c.tris);
    const real R = std::sqrt(real(24) / (real(4) * real(M_PI)));
    c.xt = radial_ellipsoid(c.x0, vec3(R, R, R));
    out.push_back(c);
  }
  {
    test_case c{"cube->ellipsoid"};
    cube(6, c.x0, c.tris);
    c.xt = radial_ellipsoid(c.x0, vec3(1.8, 1.2, 0.9));
    out.push_back(c);
  }
  {
    test_case c{"sphere cap"};
    std::vector<std::array<real, 2>> P;
    strip(13, 9, 1.6, 1.2, c.x0, P, c.tris);
    for (const auto &p : P) {
      const real rr = std::hypot(p[0], p[1]);
      const real d0 = rr > 1e-12 ? p[0] / rr : 0, d1 = rr > 1e-12 ? p[1] / rr : 0, a = rr / 2;
      c.xt.push_back(vec3(2 * std::sin(a) * d0, 2 * (1 - std::cos(a)), 2 * std::sin(a) * d1));
    }
    out.push_back(c);
  }
  {
    test_case c{"saddle a 0.8"};
    std::vector<std::array<real, 2>> P;
    strip(13, 13, 1.6, 1.6, c.x0, P, c.tris);
    for (size_t i = 0; i < P.size(); ++i)
      c.xt.push_back(vec3(P[i][0], -real(0.4) * (P[i][0] * P[i][0] - P[i][1] * P[i][1]), P[i][1]));
    out.push_back(c);
  }
  {
    test_case c{"cylinder r 1.5"};
    std::vector<std::array<real, 2>> P;
    strip(17, 6, 2.0, 1.0, c.x0, P, c.tris);
    for (const auto &p : P) {
      const real t = p[0] / real(1.5) - real(M_PI) / 2;
      c.xt.push_back(vec3(1.5 * std::cos(t), 1.5 * std::sin(t) + 1.5, p[1]));
    }
    out.push_back(c);
  }
  return out;
}

struct run_stats {
  int iterates_to = -1; // first iterate with merit / start below 1e-9
  int iterations = 0, factorizations = 0, backsolves = 0;
  real shape = 0, seconds = 0;
};

inline run_stats run_case(const test_case &c, const cnet::config &cfg) {
  const kusama::stencil_table t = table_of(static_cast<int>(c.x0.size()), c.tris);
  const cnet::topology tp = cnet::build_topology(t);
  const cnet::targets T = cnet::read(tp, c.xt, cnet::vertex_normals(tp, c.xt));
  const auto t0 = std::chrono::steady_clock::now();
  const cnet::result r = cnet::solve(tp, c.x0, cnet::vertex_normals(tp, c.x0), T, cfg);
  run_stats s;
  s.seconds = std::chrono::duration<real>(std::chrono::steady_clock::now() - t0).count();
  for (size_t k = 0; k < r.merit.size(); ++k)
    if (r.merit[k] < real(1e-9) * r.merit[0]) {
      s.iterates_to = static_cast<int>(k);
      break;
    }
  s.iterations = r.iterations;
  s.factorizations = r.factorizations;
  s.backsolves = r.backsolves;
  s.shape = shape_error(r.x, c.xt);
  return s;
}

} // namespace contact_net_detail

GAUDI_TEST(contact_net_readout_is_exact) {
  using namespace contact_net_detail;
  std::vector<vec3> x;
  tri_list tris;
  torus(2.0, 0.7, 24, 12, x, tris);
  std::mt19937 rng(3);
  std::normal_distribution<real> nd(0, 0.03);
  for (vec3 &p : x)
    p += vec3(nd(rng), nd(rng), nd(rng));
  const cnet::topology tp = cnet::build_topology(table_of(static_cast<int>(x.size()), tris));
  const std::vector<vec3> N = cnet::vertex_normals(tp, x);
  const cnet::targets T = cnet::read(tp, x, N);
  real worst = 0;
  for (size_t e = 0; e < tp.edges.size(); ++e) {
    const int i = tp.edges[e][0], j = tp.edges[e][1];
    const vec3 dx = x[static_cast<size_t>(j)] - x[static_cast<size_t>(i)];
    vec3 nb = N[static_cast<size_t>(i)] + N[static_cast<size_t>(j)];
    nb.normalize();
    const vec3 rec = -T.k[e] * (dx - dx.dot(nb) * nb) - T.tau[e] * nb.cross(dx);
    worst = std::max(worst, (rec - (N[static_cast<size_t>(j)] - N[static_cast<size_t>(i)])).cwiseAbs().maxCoeff());
  }
  GAUDI_EXPECT(worst < 1e-13);
}

GAUDI_TEST(contact_net_twist_vanishes_on_curvature_lines) {
  using namespace contact_net_detail;
  const real R = 2.0, r = 0.7;
  const int nu = 48, nv = 24;
  std::vector<vec3> x;
  tri_list tris;
  torus(R, r, nu, nv, x, tris);
  const cnet::topology tp = cnet::build_topology(table_of(static_cast<int>(x.size()), tris));
  std::vector<vec3> N;
  for (const vec3 &p : x) {
    const real rho = std::hypot(p[0], p[1]);
    N.push_back((p - vec3(R * p[0] / rho, R * p[1] / rho, 0)).normalized());
  }
  const cnet::targets T = cnet::read(tp, x, N);
  real grid = 0, diag = 1e9;
  for (size_t e = 0; e < tp.edges.size(); ++e) {
    const int i = tp.edges[e][0], j = tp.edges[e][1];
    const bool same_u = i / nv == j / nv, same_v = i % nv == j % nv;
    if (same_u || same_v)
      grid = std::max(grid, std::abs(T.tau[e]));
    else
      diag = std::min(diag, std::abs(T.tau[e]));
  }
  GAUDI_EXPECT(grid < 1e-12);
  GAUDI_EXPECT(diag > 1e-3);
}

GAUDI_TEST(contact_net_generated_rows_match_differences) {
  std::mt19937 rng(11);
  std::normal_distribution<real> nd(0, 1);
  real worst = 0;
  for (int trial = 0; trial < 10; ++trial) {
    Eigen::Matrix<real, 12, 1> x;
    for (int i = 0; i < 12; ++i)
      x[i] = nd(rng);
    x.segment<3>(6).normalize();
    x.segment<3>(9) = x.segment<3>(6) + 0.2 * vec3(nd(rng), nd(rng), nd(rng));
    const Eigen::Matrix<real, 2, 1> c(nd(rng), nd(rng));
    Eigen::Matrix<real, 3, 1> y, yp, ym;
    Eigen::Matrix<real, 3, 12> J, Jd;
    kusama::contact_generated::contact_rodrigues(c, x, y, J);
    for (int i = 0; i < 12; ++i) {
      Eigen::Matrix<real, 12, 1> xp = x, xm = x;
      xp[i] += 1e-6;
      xm[i] -= 1e-6;
      kusama::contact_generated::contact_rodrigues(c, xp, yp, Jd);
      kusama::contact_generated::contact_rodrigues(c, xm, ym, Jd);
      worst = std::max(worst, ((yp - ym) / 2e-6 - J.col(i)).cwiseAbs().maxCoeff());
    }
  }
  GAUDI_EXPECT(worst < 1e-7);
}

GAUDI_TEST(contact_net_targets_are_fixed_points) {
  using namespace contact_net_detail;
  for (const auto &c : suite()) {
    const cnet::topology tp = cnet::build_topology(table_of(static_cast<int>(c.xt.size()), c.tris));
    const std::vector<vec3> N = cnet::vertex_normals(tp, c.xt);
    const cnet::targets T = cnet::read(tp, c.xt, N);
    for (cnet::kind k : {cnet::kind::gn, cnet::kind::pd}) {
      cnet::config cfg;
      cfg.solver = k;
      cfg.iterations = 5;
      const cnet::result r = cnet::solve(tp, c.xt, N, T, cfg);
      real d = 0;
      for (size_t v = 0; v < c.xt.size(); ++v)
        d = std::max(d, (r.x[v] - c.xt[v]).norm());
      GAUDI_EXPECT(d < 1e-9);
    }
  }
}

/// Gauss-Newton with the fading memory (defaults): Python reaches merit /
/// start 1e-9 in 7 to 15 iterates on every case (MATH.md 19.12).
GAUDI_TEST(contact_net_gn_converges_on_the_suite) {
  using namespace contact_net_detail;
  cnet::config cfg;
  cfg.iterations = 40;
  for (const auto &c : suite()) {
    const run_stats s = run_case(c, cfg);
    std::printf("  contact_net gn  %-22s 1e-9 at %3d  its %3d  factors %3d  solves %3d  shape %.1e  %.3fs\n",
                c.name.c_str(), s.iterates_to, s.iterations, s.factorizations, s.backsolves, s.shape, s.seconds);
    GAUDI_EXPECT(s.iterates_to > 0 && s.iterates_to <= 25);
    GAUDI_EXPECT(s.shape < 1e-9);
  }
}

/// The constant-factor projection solve with Anderson (default) and with
/// Chebyshev: Anderson reaches 1e-9 on the closed shapes in a few hundred
/// iterates; Chebyshev is slower.
GAUDI_TEST(contact_net_pd_accelerators) {
  using namespace contact_net_detail;
  const auto cases = suite();
  for (cnet::accel a : {cnet::accel::anderson, cnet::accel::chebyshev}) {
    cnet::config cfg;
    cfg.solver = cnet::kind::pd;
    cfg.acceleration = a;
    cfg.iterations = 3000;
    for (size_t k = 0; k < 4; ++k) { // the closed shapes
      const run_stats s = run_case(cases[k], cfg);
      std::printf("  contact_net pd %-9s %-22s 1e-9 at %4d  its %4d  solves %5d  shape %.1e  %.3fs\n",
                  a == cnet::accel::anderson ? "anderson" : "chebyshev", cases[k].name.c_str(), s.iterates_to,
                  s.iterations, s.backsolves, s.shape, s.seconds);
      GAUDI_EXPECT(s.shape < (a == cnet::accel::anderson ? 1e-8 : 1e-6));
    }
  }
}

/// The native rows' position stiffness against the cotan `bend` row at equal
/// weights (as `exponent_coupling_matches_bend_at_equal_weight`), the mean
/// squared column norm of their Jacobian on positions.
GAUDI_TEST(contact_net_gn_matches_bend_at_equal_weight) {
  namespace hb = hepworth::block;
  for (auto [R, r, nu, nv] : {std::tuple{1.0, 0.33, 48, 24}, std::tuple{2.0, 0.7, 24, 12}}) {
    asawa::shell::shell::ptr M = asawa::shell::load_torus(R, r, nu, nv);
    std::vector<vec3> &x = asawa::get_vec_data(*M, 0);
    const int n = static_cast<int>(x.size());
    const kusama::stencil_table t = kusama::build_stencil_table(asawa::shell::shell_walker(*M));
    std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Ones());
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    const kusama::contact_net::topology tp = kusama::contact_net::build_topology(t);
    std::vector<vec3> N = kusama::contact_net::vertex_normals(tp, x), nv0(n, vec3::Zero()), nf = nv0,
                      nm(n, vec3::Constant(real(1e-8)));
    auto nrm = hepworth::vec3_block::create(nm, N, nv0, nf);
    std::vector<vec3> hs(static_cast<size_t>(t.n_edges), vec3::Zero()), hv = hs, hf = hs,
                      hm(hs.size(), vec3::Constant(real(1e-8)));
    auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
    vecX q;
    pos->map_to_x(q);
    nrm->map_to_x(q);
    store->map_to_x(q);
    auto mean_pos_diag = [&](const std::vector<trip> &trips) {
      std::map<index_t, real> d;
      for (const trip &tr : trips)
        if (tr.col() < 3 * n)
          d[tr.col()] += tr.value() * tr.value();
      real sum = 0;
      for (const auto &[c, val] : d)
        sum += val;
      return d.empty() ? real(0) : sum / real(d.size());
    };
    std::vector<hepworth::projection_constraint::ptr> bend;
    hb::init_bending(*M, bend, x, real(1), {pos});
    std::vector<trip> tb;
    index_t id0 = 0;
    for (auto &c : bend)
      c->fill_A(id0, tb);
    auto net = hb::contact_net_gn::create(&t, x, N, hs, real(1),
                                          std::vector<hepworth::sim_block::ptr>{pos, nrm, store});
    std::vector<trip> tc;
    id0 = 0;
    net->fill_J(q, id0, tc);
    const real kb = mean_pos_diag(tb), kc = mean_pos_diag(tc);
    std::printf("  contact_net_gn w = 1 position stiffness on torus (%g, %g, %d x %d): bend %g net %g ratio %g\n", R,
                r, nu, nv, kb, kc, kc / kb);
    GAUDI_EXPECT(kc / kb > real(1) / 3 && kc / kb < real(3));
  }
}

/// The native rows through `projection_solver` as the demo uses them: every
/// frame the rows are rebuilt at the current state (normals reseeded from the
/// mesh, the store's increment the remaining change toward the target read-out,
/// held by a large mass), positions nearly massless.
///
/// Asserted: with the target's tangency offsets (`CNET_EXACT_C`, default on
/// here) the closed shapes reach machine precision (the rows and the solve are
/// right). Reported, not asserted (MATH.md 19.13), with `CNET_EXACT_C=0` (what
/// the demo has; `CNET_TANGENCY=vertex|edge`, `CNET_STORE=cr|net`,
/// `CNET_MEMORY=w`, `CNET_MASS=m`): the vertex-tensor target follows large
/// edits where the edge target holds the old shape; the rows alone have a
/// collapsed zero (every point and normal equal), so without the demo's memory
/// the cubes shrink to a point; positions need mass or memory as the damping
/// `projection_solver` lacks; CR store increments leave a floor of a few
/// percent (the net cannot express all of them).
GAUDI_TEST(contact_net_gn_converges_in_projection_solver) {
  using namespace contact_net_detail;
  namespace hb = hepworth::block;
  const char *ef = std::getenv("CNET_FRAMES");
  const int frames = ef ? std::atoi(ef) : 12;
  if (!std::getenv("CNET_EXACT_C"))
    setenv("CNET_EXACT_C", "1", 0);
  const bool exact = std::atoi(std::getenv("CNET_EXACT_C")) != 0;
  for (const auto &c : suite()) {
    std::vector<vec3> x = c.x0;
    const int n = static_cast<int>(x.size());
    const kusama::stencil_table t = table_of(n, c.tris);
    const cnet::topology tp = cnet::build_topology(t);
    const cnet::targets T1 = cnet::read(tp, c.xt, cnet::vertex_normals(tp, c.xt));
    const char *em = std::getenv("CNET_MASS");
    std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Constant(em ? real(std::atof(em)) : real(1e-6)));
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    const size_t nE = static_cast<size_t>(t.n_edges);
    hb::projection_solver solver;
    solver.gn_refresh_every = 5;
    solver.gn_backtracks = 8;
    int its = 0, refactors = 0;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<real> err_f;
    for (int k = 0; k < frames; ++k) {
      std::vector<vec3> N = cnet::vertex_normals(tp, x), nv0(n, vec3::Zero()), nf = nv0,
                        nm(n, vec3::Constant(real(1e-8)));
      auto nrm = hepworth::vec3_block::create(nm, N, nv0, nf);
      const cnet::targets T0 = cnet::read(tp, x, N);
      std::vector<vec3> h0(nE, vec3::Zero()), hs(nE), hv(nE, vec3::Zero()), hf = hv,
          hm(nE, vec3::Constant(real(1e12)));
      // CNET_RAMP=1: an edit spread over the frames as the demo's store edits
      // are (this frame takes 1 / frames-left of the remaining change).
      const char *er = std::getenv("CNET_RAMP");
      const real a = er && std::atoi(er) ? real(1) / real(frames - k) : real(1);
      // CNET_STORE=cr: the store moved by the CR edge readings' change toward
      // the target (all three components), as the demo's store edits are;
      // otherwise the net's own (k, tau) increments.
      const char *estore = std::getenv("CNET_STORE");
      const bool cr_store = estore && std::string(estore) == "cr";
      const std::vector<vec3> R0 = cr_store ? kusama::weingarten::edge_readings(t, x) : std::vector<vec3>{};
      const std::vector<vec3> R1 = cr_store ? kusama::weingarten::edge_readings(t, c.xt) : std::vector<vec3>{};
      for (size_t e = 0; e < nE; ++e)
        hs[e] = cr_store ? vec3(a * (R1[e] - R0[e])) : vec3(a * vec3(T1.k[e] - T0.k[e], T1.tau[e] - T0.tau[e], 0));
      auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
      std::vector<hepworth::sim_block::ptr> blocks{pos, nrm, store};
      // CNET_EXACT_C=1: hold the target mesh's tangency offsets instead of the
      // frame start's (separates the offset's approximation from the solve).
      std::vector<real> coff;
      const char *ex = std::getenv("CNET_EXACT_C");
      if (ex && std::atoi(ex)) {
        coff.resize(tp.spokes.size());
        std::map<std::pair<int, int>, int> eo;
        for (size_t e = 0; e < nE; ++e)
          eo[{tp.edges[e][0], tp.edges[e][1]}] = static_cast<int>(e);
        for (size_t s = 0; s < tp.spokes.size(); ++s) {
          const int i = tp.spokes[s][0], j = tp.spokes[s][1];
          const int e = eo[{std::min(i, j), std::max(i, j)}];
          coff[s] = T1.c[s] - real(0.5) * T1.k[static_cast<size_t>(e)] *
                                  (c.xt[static_cast<size_t>(j)] - c.xt[static_cast<size_t>(i)]).norm();
        }
      }
      auto row = hb::contact_net_gn::create(&t, x, N, h0, real(1), blocks, hb::k_contact_gn,
                                            coff.empty() ? nullptr : &coff);
      // CNET_TANGENCY=edge|vertex (default: edge with exact offsets, vertex otherwise).
      const char *etan = std::getenv("CNET_TANGENCY");
      const bool edge_tan = etan ? std::string(etan) == "edge" : !coff.empty();
      row->tangency_target = edge_tan ? hb::contact_net_gn::tangency::edge : hb::contact_net_gn::tangency::vertex;
      // CNET_MEMORY=w: the demo's frame-start Dirichlet memory beside the rows.
      const char *emem = std::getenv("CNET_MEMORY");
      const real w_mem = emem ? real(std::atof(emem)) : real(0);
      std::vector<hepworth::projection_constraint::ptr> cs{row};
      if (w_mem > 0)
        cs.push_back(hb::dirichlet_memory_gn::create(&t, x, w_mem, std::vector<hepworth::sim_block::ptr>{pos}));
      solver.set_constraints(cs);
      for (vec3 &vv : v)
        vv.setZero();
      solver.step(blocks, real(0.1), real(1), 30);
      its += solver.last_step.its;
      refactors += solver.last_step.refreshes;
      err_f.push_back(shape_error(x, c.xt));
    }
    const real secs = std::chrono::duration<real>(std::chrono::steady_clock::now() - t0).count();
    std::printf("  contact_net_gn in projection_solver %-22s %d frames: its %3d refactors %3d  shape %.1e -> %.1e "
                "(frame 1 %.1e, 3 %.1e)  %.3fs\n",
                c.name.c_str(), frames, its, refactors, shape_error(c.x0, c.xt), err_f.back(), err_f[0],
                err_f[std::min<size_t>(2, err_f.size() - 1)], secs);
    const bool closed = c.name.rfind("torus", 0) == 0 || c.name.rfind("cube", 0) == 0;
    if (exact && closed)
      GAUDI_EXPECT(err_f.back() < 1e-6);
  }
}

/// Frames mode (per-vertex quaternions in place of the normals block) lands
/// where the normals mode does on the closed suite cases: exact target
/// offsets, edge tangency, the net's own `(k, tau)` increments.
GAUDI_TEST(contact_net_gn_frames_match_normals) {
  using namespace contact_net_detail;
  namespace hb = hepworth::block;
  const int frames = 12;
  auto run = [&](const auto &c, bool quat) {
    std::vector<vec3> x = c.x0;
    const int n = static_cast<int>(x.size());
    const kusama::stencil_table t = table_of(n, c.tris);
    const cnet::topology tp = cnet::build_topology(t);
    const cnet::targets T1 = cnet::read(tp, c.xt, cnet::vertex_normals(tp, c.xt));
    std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Constant(real(1e-6)));
    auto pos = hepworth::vec3_block::create(m, x, v, f);
    const size_t nE = static_cast<size_t>(t.n_edges);
    std::map<std::pair<int, int>, int> eo;
    for (size_t e = 0; e < nE; ++e)
      eo[{tp.edges[e][0], tp.edges[e][1]}] = static_cast<int>(e);
    std::vector<real> coff(tp.spokes.size());
    for (size_t s = 0; s < tp.spokes.size(); ++s) {
      const int i = tp.spokes[s][0], j = tp.spokes[s][1];
      const int e = eo[{std::min(i, j), std::max(i, j)}];
      coff[s] = T1.c[s] - real(0.5) * T1.k[static_cast<size_t>(e)] *
                              (c.xt[static_cast<size_t>(j)] - c.xt[static_cast<size_t>(i)]).norm();
    }
    hb::projection_solver solver;
    solver.gn_refresh_every = 5;
    solver.gn_backtracks = 8;
    int its = 0;
    for (int k = 0; k < frames; ++k) {
      std::vector<vec3> N = cnet::vertex_normals(tp, x), nv0(n, vec3::Zero()), nf = nv0,
                        nm(n, vec3::Constant(real(1e-8)));
      std::vector<vec4> Y(n), ym(n, vec4::Constant(real(1e-8)));
      for (int i = 0; i < n; ++i) {
        const vec3 nn = N[static_cast<size_t>(i)];
        vec3 tt = std::abs(nn[0]) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
        tt = (tt - nn * nn.dot(tt)).normalized();
        mat3 R;
        R.col(0) = tt;
        R.col(1) = nn.cross(tt);
        R.col(2) = nn;
        Y[static_cast<size_t>(i)] = hb::quat_ops::quat_from_matrix(R);
      }
      hepworth::sim_block::ptr nb = quat ? hepworth::sim_block::ptr(hepworth::vec4_block::create(ym, Y))
                                         : hepworth::sim_block::ptr(hepworth::vec3_block::create(nm, N, nv0, nf));
      const cnet::targets T0 = cnet::read(tp, x, N);
      std::vector<vec3> h0(nE, vec3::Zero()), hs(nE), hv(nE, vec3::Zero()), hf = hv,
          hm(nE, vec3::Constant(real(1e12)));
      for (size_t e = 0; e < nE; ++e)
        hs[e] = vec3(T1.k[e] - T0.k[e], T1.tau[e] - T0.tau[e], 0);
      auto store = hepworth::vec3_block::create(hm, hs, hv, hf);
      std::vector<hepworth::sim_block::ptr> blocks{pos, nb, store};
      auto row = hb::contact_net_gn::create(&t, x, N, h0, real(1), blocks, hb::k_contact_gn, &coff);
      row->tangency_target = hb::contact_net_gn::tangency::edge;
      if (const char *eg = std::getenv("CNET_GAUGE"))
        row->w_frame_gauge = std::atof(eg);
      solver.set_constraints({row});
      for (vec3 &vv : v)
        vv.setZero();
      solver.step(blocks, real(0.1), real(1), 30);
      its += solver.last_step.its;
      if (std::getenv("CNET_DEBUG"))
        std::printf("    %s frame %d its %d shape %.2e E %.3e -> %.3e\n", quat ? "quat" : "nrm ", k,
                    solver.last_step.its, shape_error(x, c.xt), solver.last_step.E_first, solver.last_step.E_last);
    }
    return std::make_pair(shape_error(x, c.xt), its);
  };
  for (const auto &c : suite()) {
    const bool closed = c.name.rfind("torus", 0) == 0 || c.name.rfind("cube", 0) == 0;
    if (!closed)
      continue;
    const auto a = run(c, false), b = run(c, true);
    std::printf("  frames vs normals %-22s shape %.1e (its %d) against %.1e (its %d)\n", c.name.c_str(), b.first,
                b.second, a.first, a.second);
    GAUDI_EXPECT(b.first < 1e-6);
  }
}

/// Store-free frames mode: positions and quaternion frames only. The frame
/// start is an exact fixed point, and a growth drive through `frame_rows`
/// (vertex tensors as their own block) beside the metric memory tightens the
/// tube.
GAUDI_TEST(contact_net_gn_store_free_frames) {
  using namespace contact_net_detail;
  namespace hb = hepworth::block;
  std::vector<vec3> x;
  tri_list tris;
  torus(1.0, 0.35, 32, 16, x, tris);
  const int n = static_cast<int>(x.size());
  const kusama::stencil_table t = table_of(n, tris);
  const cnet::topology tp = cnet::build_topology(t);
  auto frames_of = [&](const std::vector<vec3> &N) {
    std::vector<vec4> Y(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
      const vec3 nn = N[static_cast<size_t>(i)];
      vec3 tt = std::abs(nn[0]) < real(0.9) ? vec3::UnitX() : vec3::UnitY();
      tt = (tt - nn * nn.dot(tt)).normalized();
      mat3 R;
      R.col(0) = tt;
      R.col(1) = nn.cross(tt);
      R.col(2) = nn;
      Y[static_cast<size_t>(i)] = hb::quat_ops::quat_from_matrix(R);
    }
    return Y;
  };
  auto tube = [&]() {
    real rmin = 1e30, rmax = 0;
    for (const vec3 &p : x) {
      rmin = std::min(rmin, std::hypot(p[0], p[1]));
      rmax = std::max(rmax, std::hypot(p[0], p[1]));
    }
    const real R = real(0.5) * (rmin + rmax);
    real s = 0;
    for (const vec3 &p : x)
      s += std::hypot(std::hypot(p[0], p[1]) - R, p[2]);
    return s / real(n);
  };
  std::vector<vec3> v(n, vec3::Zero()), f = v, m(n, vec3::Constant(real(1e-6)));
  auto pos = hepworth::vec3_block::create(m, x, v, f);
  hb::projection_solver solver;
  solver.gn_refresh_every = 5;
  solver.gn_backtracks = 8;
  const real tube0 = tube();
  for (int k = 0; k < 3; ++k) {
    const std::vector<vec3> N = cnet::vertex_normals(tp, x);
    std::vector<vec4> Y = frames_of(N), ym(static_cast<size_t>(n), vec4::Constant(real(1e-8)));
    auto fb = hepworth::vec4_block::create(ym, Y);
    std::vector<vec3> S(static_cast<size_t>(n), vec3::Zero()), sv = S, sf = S,
        sm(static_cast<size_t>(n), vec3::Constant(real(1e-8)));
    auto sb = hepworth::vec3_block::create(sm, S, sv, sf);
    std::vector<hepworth::sim_block::ptr> blocks{pos, fb, sb};
    auto net = hb::contact_net_gn::create(&t, x, N, std::vector<vec3>{}, real(1),
                                          std::vector<hepworth::sim_block::ptr>{pos, fb});
    hb::frame_rows_config fc;
    fc.w_growth = 0.5;
    fc.growth_x = 0.2;
    auto rows = hb::frame_rows::create(&t, x, Y, fc, blocks);
    S = rows->start_tensors();
    vecX q;
    for (auto &b : blocks)
      b->map_to_x(q);
    if (k == 0) {
      const real r0 = net->residual_sq(q);
      std::printf("  store-free net at the frame start: |r|^2 %.2e\n", r0);
      GAUDI_EXPECT(r0 < 1e-20);
    }
    // The demo's metric memory: tangency alone leaves slides and scale free.
    auto mem = hb::dirichlet_memory_gn::create(&t, x, real(0.5), std::vector<hepworth::sim_block::ptr>{pos});
    solver.set_constraints({net, rows, mem});
    for (vec3 &vv : v)
      vv.setZero();
    solver.step(blocks, real(0.1), real(1), 20);
    std::printf("  store-free growth frame %d: tube %.4f (start %.4f) its %d\n", k, tube(), tube0,
                solver.last_step.its);
  }
  GAUDI_EXPECT(tube() < tube0);
  for (const vec3 &p : x)
    GAUDI_EXPECT(p.allFinite());
}

} // namespace test
} // namespace gaudi

#endif
