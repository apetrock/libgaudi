#ifndef GAUDI_TP_CONSTRAINT_SOLVER_TESTS_HPP
#define GAUDI_TP_CONSTRAINT_SOLVER_TESTS_HPP

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/tp_config.hpp"
#include "gaudi/duchamp/modules/tp_constraint_solver.hpp"
#include "gaudi/duchamp/modules/tp_sobolev_jacobi.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

namespace {

inline std::vector<vec3> constraint_ngon(int n) {
  std::vector<vec3> pts;
  pts.reserve(static_cast<size_t>(n));
  const vec3 shift(2.0, -1.0, 0.4);
  for (int i = 0; i < n; ++i) {
    const real t = real(2) * real(3.14159265358979323846) * real(i) / real(n);
    pts.emplace_back(std::cos(t), std::sin(t), real(0));
    pts.back() += shift;
  }
  return pts;
}

inline duchamp::sobolev_tangent_point_force_config constraint_metric() {
  duchamp::sobolev_tangent_point_force_config cfg;
  cfg.alpha_B = 1.0;
  cfg.alpha_0 = 1.0;
  cfg.p1 = 2.0;
  cfg.bh_eps = 0.0;
  cfg.cone_eps = 0.0;
  cfg.solver.type = duchamp::sobolev_solver::cg;
  cfg.solver.iter = 40;
  cfg.num_v_cycles = 0;
  cfg.max_levels = 1;
  return cfg;
}

inline real length_row_max(asawa::rod::rod &rod, const std::vector<vec3> &g) {
  const auto E = duchamp::sobolev_detail::build_edges(rod);
  real m = 0.0;
  for (size_t e = 0; e < E.n; ++e) {
    const vec3 dg = g[static_cast<size_t>(E.ep1[e])] -
                    g[static_cast<size_t>(E.ep0[e])];
    m = std::max(m, std::abs(E.T[e].dot(dg)));
  }
  return m;
}

inline vec3 field_sum(const std::vector<vec3> &g) {
  vec3 s = vec3::Zero();
  for (const vec3 &a : g)
    s += a;
  return s;
}

inline real max_length_drift(asawa::rod::rod &rod,
                             const std::vector<real> &rest) {
  const auto E = duchamp::sobolev_detail::build_edges(rod);
  real m = 0.0;
  const size_t n = std::min(E.n, rest.size());
  for (size_t e = 0; e < n; ++e)
    m = std::max(m, std::abs(E.ell[e] - rest[e]));
  return m;
}

inline real c2_norm_at(const asawa::rod::rod &rod, int i) {
  const auto &x = rod.x();
  const auto ci = asawa::rod::corner_id(i);
  const auto cim1 = rod.prev(ci);
  const auto cip1 = rod.next(ci);
  const auto cim2 = rod.prev(cim1);
  const auto cip2 = rod.next(cip1);
  const vec3 d = -x[static_cast<size_t>(static_cast<int>(cim2))] +
                 real(2) * x[static_cast<size_t>(static_cast<int>(cim1))] -
                 real(2) * x[static_cast<size_t>(static_cast<int>(cip1))] +
                 x[static_cast<size_t>(static_cast<int>(cip2))];
  return d.norm();
}

inline real c1_norm_at(const asawa::rod::rod &rod, int i) {
  const auto &x = rod.x();
  const auto ci = asawa::rod::corner_id(i);
  const int ip = static_cast<int>(rod.prev(ci));
  const int in = static_cast<int>(rod.next(ci));
  const vec3 d = x[static_cast<size_t>(ip)] - real(2) * x[static_cast<size_t>(i)] +
                 x[static_cast<size_t>(in)];
  return d.norm();
}

} // namespace

/// Hard length rows kill a dilation; hard barycenter kills a translation.
GAUDI_TEST(sobolev_saddle_kills_dilation_and_translation) {
  auto R = asawa::rod::rod::create(constraint_ngon(8), true);
  const auto cfg = constraint_metric();

  const auto &x = R->x();
  vec3 mean = vec3::Zero();
  for (const vec3 &p : x)
    mean += p;
  mean /= real(x.size());

  std::vector<vec3> f(x.size(), vec3::Zero());
  for (size_t i = 0; i < x.size(); ++i)
    f[i] = x[i] - mean;

  const std::vector<vec3> g0 =
      duchamp::sobolev_mg_solve(*R, f, cfg, real(1));
  GAUDI_ASSERT(length_row_max(*R, g0) > real(1e-4));

  duchamp::sobolev_constraint_solver solver;
  solver.max_cg = 32;
  solver.add(duchamp::sobolev_edge_length_constraint::create());
  solver.add(duchamp::sobolev_barycenter_constraint::create());
  solver.init_rest(*R);
  const std::vector<vec3> g = solver.project_gradient(*R, f, cfg, real(1));

  GAUDI_ASSERT(length_row_max(*R, g) < real(1e-6));
  GAUDI_ASSERT(field_sum(g).norm() < real(1e-6));
}

/// The config list builds the same hard rows as add().
GAUDI_TEST(sobolev_constraint_list_kills_dilation) {
  auto R = asawa::rod::rod::create(constraint_ngon(8), true);
  const auto metric = constraint_metric();

  const auto &x = R->x();
  vec3 mean = vec3::Zero();
  for (const vec3 &p : x)
    mean += p;
  mean /= real(x.size());
  std::vector<vec3> f(x.size(), vec3::Zero());
  for (size_t i = 0; i < x.size(); ++i)
    f[i] = x[i] - mean;

  duchamp::tangent_point_solver_config cfg;
  cfg.sobolev = metric;
  cfg.constraint_cg = 64;
  cfg.hard_constraints = {
      {.type = duchamp::sobolev_constraint_type::edge_length},
      {.type = duchamp::sobolev_constraint_type::barycenter},
  };
  cfg.proj_constraints = {
      {.type = duchamp::sobolev_constraint_type::c1, .w = 1.0},
  };
  duchamp::sobolev_constraint_solver solver =
      duchamp::make_sobolev_constraint_solver(cfg);
  solver.init_rest(*R);
  const std::vector<vec3> g =
      solver.project_gradient(*R, f, metric, real(1));

  GAUDI_ASSERT(length_row_max(*R, g) < real(1e-6));
  GAUDI_ASSERT(field_sum(g).norm() < real(1e-6));
}

/// A spiked vertex is pulled back toward the rest polygon. Lengths stay put.
GAUDI_TEST(sobolev_saddle_c1_damps_spike) {
  auto R = asawa::rod::rod::create(constraint_ngon(8), true);
  const auto cfg = constraint_metric();
  const auto E0 = duchamp::sobolev_detail::build_edges(*R);
  const std::vector<real> rest(E0.ell.begin(), E0.ell.end());
  const vec3 mean0 = [&] {
    vec3 m = vec3::Zero();
    for (const vec3 &p : R->x())
      m += p;
    return m / real(R->x().size());
  }();

  duchamp::sobolev_constraint_solver solver;
  solver.max_cg = 64;
  solver.add(duchamp::sobolev_edge_length_constraint::create());
  solver.add(duchamp::sobolev_barycenter_constraint::create());
  solver.add(duchamp::sobolev_c1_constraint::create(real(4)));
  solver.init_rest(*R);

  R->x()[0] += real(0.15) * (R->x()[0] - mean0);
  const real c1_before = c1_norm_at(*R, 0);
  const real len_before = max_length_drift(*R, rest);
  GAUDI_ASSERT(len_before > real(1e-3));
  GAUDI_ASSERT(c1_before > real(1e-3));

  for (int k = 0; k < 6; ++k)
    solver.project_positions(*R, cfg, real(1));

  GAUDI_ASSERT(max_length_drift(*R, rest) < real(1e-4) * (real(1) + len_before));
  GAUDI_ASSERT(c1_norm_at(*R, 0) < c1_before);
  const vec3 mean1 = [&] {
    vec3 m = vec3::Zero();
    for (const vec3 &p : R->x())
      m += p;
    return m / real(R->x().size());
  }();
  GAUDI_ASSERT((mean1 - mean0).norm() < real(1e-4));
}

/// w == 0 leaves C2 out of the saddle. A positive weight damps the third difference.
GAUDI_TEST(sobolev_c2_gated_and_damps) {
  duchamp::tangent_point_solver_config off;
  off.hard_constraints = {
      {.type = duchamp::sobolev_constraint_type::c2, .w = 0.0}};
  duchamp::sobolev_constraint_solver empty =
      duchamp::make_sobolev_constraint_solver(off);
  GAUDI_ASSERT(empty.empty());

  auto R = asawa::rod::rod::create(constraint_ngon(8), true);
  const auto cfg = constraint_metric();
  const auto E0 = duchamp::sobolev_detail::build_edges(*R);
  const std::vector<real> rest(E0.ell.begin(), E0.ell.end());
  vec3 mean0 = vec3::Zero();
  for (const vec3 &p : R->x())
    mean0 += p;
  mean0 /= real(R->x().size());

  duchamp::tangent_point_solver_config spec;
  spec.constraint_cg = 64;
  spec.hard_constraints = {
      {.type = duchamp::sobolev_constraint_type::edge_length},
      {.type = duchamp::sobolev_constraint_type::c2, .w = 4.0},
  };
  duchamp::sobolev_constraint_solver solver =
      duchamp::make_sobolev_constraint_solver(spec);
  solver.init_rest(*R);

  R->x()[0] += real(0.15) * (R->x()[0] - mean0);
  const real c2_before = c2_norm_at(*R, 0);
  GAUDI_ASSERT(c2_before > real(1e-3));
  for (int k = 0; k < 6; ++k)
    solver.project_positions(*R, cfg, real(1));
  GAUDI_ASSERT(c2_norm_at(*R, 0) < c2_before);
  GAUDI_ASSERT(max_length_drift(*R, rest) < real(1e-4));
}

/// Rest lengths and the barycenter target survive a split and a collapse.
GAUDI_TEST(sobolev_rest_survives_split_and_collapse) {
  auto R = asawa::rod::rod::create(constraint_ngon(8), true);
  const auto cfg = constraint_metric();

  duchamp::sobolev_constraint_solver solver;
  solver.max_cg = 64;
  solver.add(duchamp::sobolev_edge_length_constraint::create());
  solver.add(duchamp::sobolev_barycenter_constraint::create());
  solver.init_rest(*R);

  auto rest_sum = [&]() {
    const auto E = duchamp::sobolev_detail::build_edges(*R);
    const auto &L = R->rest_l();
    real s = 0.0;
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      s += L[static_cast<size_t>(i0)];
    }
    return s;
  };
  auto mean_of = [&]() {
    vec3 m = vec3::Zero();
    for (const vec3 &p : R->x())
      m += p;
    return m / real(R->x().size());
  };
  auto length_drift = [&]() {
    const auto E = duchamp::sobolev_detail::build_edges(*R);
    const auto &L = R->rest_l();
    real m = 0.0;
    for (size_t e = 0; e < E.n; ++e) {
      const int i0 = E.ep0[e];
      m = std::max(m, std::abs(E.ell[e] - L[static_cast<size_t>(i0)]));
    }
    return m;
  };

  const real sum0 = rest_sum();
  const vec3 mean0 = mean_of();
  GAUDI_ASSERT(sum0 > real(0));

  auto Rd = asawa::rod::dynamic::create(R, real(0), real(1e9), real(0));
  Rd->split_edge(asawa::rod::corner_id(0));
  const real sum1 = rest_sum();
  GAUDI_ASSERT(sum1 > real(0.5) * sum0);
  GAUDI_ASSERT(std::abs(sum1 - sum0) < real(0.25) * sum0);

  const size_t n_split = R->corner_count();
  Rd->collapse_edge(asawa::rod::corner_id(3));
  R->pack();
  GAUDI_ASSERT(R->corner_count() + 1 == n_split);
  const real sum2 = rest_sum();
  GAUDI_ASSERT(std::abs(sum2 - sum1) < real(1e-8) * (real(1) + sum1));

  const vec3 slide(0.3, -0.2, 0.1);
  for (vec3 &p : R->x())
    p += slide;
  GAUDI_ASSERT((mean_of() - mean0).norm() > real(0.2));
  for (int k = 0; k < 8; ++k)
    solver.project_positions(*R, cfg, real(1));
  GAUDI_ASSERT((mean_of() - mean0).norm() < real(1e-3));

  const int i_stretch = 1;
  const auto ci = asawa::rod::corner_id(i_stretch);
  const int in = static_cast<int>(R->next(ci));
  R->x()[static_cast<size_t>(i_stretch)] +=
      real(0.2) * (R->x()[static_cast<size_t>(in)] - R->x()[static_cast<size_t>(i_stretch)]);
  const real drift_stretch = length_drift();
  GAUDI_ASSERT(drift_stretch > real(1e-3));
  for (int k = 0; k < 8; ++k)
    solver.project_positions(*R, cfg, real(1));
  GAUDI_ASSERT(length_drift() < drift_stretch);
  GAUDI_ASSERT((mean_of() - mean0).norm() < real(1e-3));
}

GAUDI_TEST(sobolev_hard_projection_reduces_phi) {
  auto R = asawa::rod::rod::create(constraint_ngon(12), true);
  auto len = duchamp::sobolev_edge_length_constraint::create();
  auto bar = duchamp::sobolev_barycenter_constraint::create();
  duchamp::sobolev_constraint_solver solver;
  solver.add(len);
  solver.add(bar);
  solver.init_rest(*R);

  vec3 mean = vec3::Zero();
  for (const vec3 &p : R->x())
    mean += p;
  mean /= real(R->x().size());
  for (vec3 &p : R->x())
    p = mean + real(1.08) * (p - mean) + vec3(real(0.05), 0, 0);

  auto phi_norm = [&]() {
    std::vector<real> phi;
    len->residual(*R, phi);
    bar->residual(*R, phi);
    real s = 0;
    for (real a : phi)
      s += a * a;
    return std::sqrt(s);
  };
  const real before = phi_norm();
  GAUDI_ASSERT(before > real(1e-3));
  const auto cfg = constraint_metric();
  solver.project_positions(*R, cfg, real(1));
  const real after = phi_norm();
  GAUDI_ASSERT(std::isfinite(after));
  GAUDI_ASSERT(after < before);
}

} // namespace test
} // namespace gaudi

#endif
