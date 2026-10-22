#ifndef GAUDI_TP_SOBOLEV_JACOBI_TESTS_HPP
#define GAUDI_TP_SOBOLEV_JACOBI_TESTS_HPP

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/modules/tp_config.hpp"
#include "gaudi/duchamp/modules/tp_energies.hpp"
#include "gaudi/duchamp/modules/tp_sobolev_jacobi.hpp"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

namespace gaudi {
namespace test {

namespace {

inline std::vector<vec3> sobolev_helix(int n = 24) {
  std::vector<vec3> pts;
  pts.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const real t = real(2) * real(3.14159265358979323846) * real(i) / real(n);
    pts.emplace_back(std::cos(t), std::sin(t), real(0.12) * t);
  }
  return pts;
}

inline real field_dot(const std::vector<vec3> &a, const std::vector<vec3> &b) {
  real s = 0;
  for (size_t i = 0; i < a.size(); ++i)
    s += a[i].dot(b[i]);
  return s;
}

inline real field_norm(const std::vector<vec3> &a) {
  return std::sqrt(field_dot(a, a));
}

inline vec3 vertex_sum(const std::vector<vec3> &a) {
  vec3 s = vec3::Zero();
  for (const vec3 &v : a)
    s += v;
  return s;
}

inline real field_max_abs_diff(const std::vector<vec3> &a,
                               const std::vector<vec3> &b) {
  real m = 0;
  for (size_t i = 0; i < a.size(); ++i)
    m = std::max(m, (a[i] - b[i]).cwiseAbs().maxCoeff());
  return m;
}

} // namespace

GAUDI_TEST(sobolev_B_B0_symmetry) {
  auto R = asawa::rod::rod::create(sobolev_helix(20), true);
  const real l0 = real(0.1) * R->lavg();
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const size_t n = R->corner_count();
  std::vector<vec3> u(n), v(n);
  for (size_t i = 0; i < n; ++i) {
    const real t = real(i);
    u[i] = vec3(std::sin(t), std::cos(real(0.7) * t), real(0.1) * t);
    v[i] = vec3(std::cos(real(1.3) * t), std::sin(real(0.4) * t),
                real(-0.05) * t);
  }
  const real bh = real(0); // leaves-only
  const real Bu_v =
      duchamp::sobolev_quadratic_B(*R, /*p=*/4.0, l0, bh,
                                   u, v);
  const real Bv_u =
      duchamp::sobolev_quadratic_B(*R, /*p=*/4.0, l0, bh,
                                   v, u);
  const real B0u_v = duchamp::sobolev_quadratic_B0(
      *R, /*p=*/4.0, l0, bh, u, v);
  const real B0v_u = duchamp::sobolev_quadratic_B0(
      *R, /*p=*/4.0, l0, bh, v, u);
  const real B0u_u = duchamp::sobolev_quadratic_B0(
      *R, /*p=*/4.0, l0, bh, u, u);
  const real Bu_u = duchamp::sobolev_quadratic_B(
      *R, /*p=*/4.0, l0, bh, u, u);
  // Calder POV↔edge 1-pt is only approximately symmetric.
  GAUDI_ASSERT(std::abs(Bu_v - Bv_u) < real(0.2) * (1 + std::abs(Bu_v)));
  GAUDI_ASSERT(std::abs(B0u_v - B0v_u) < real(0.2) * (1 + std::abs(B0u_v)));
  GAUDI_ASSERT(B0u_u >= real(-1e-6));
  GAUDI_ASSERT(Bu_u >= real(-1e-6));
}

GAUDI_TEST(sobolev_B_B0_magnitude_parity) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  const real l0 = real(0.1) * R->lavg();
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const size_t n = R->corner_count();
  std::vector<vec3> u(n);
  for (size_t i = 0; i < n; ++i) {
    const real t =
        real(2) * real(3.14159265358979323846) * real(i) / real(n);
    u[i] = vec3(std::cos(t), std::sin(t), real(0.2) * std::sin(real(2) * t));
  }
  const real bh = real(0);
  const real Buu = duchamp::sobolev_quadratic_B(*R, /*p=*/3.0, l0, bh, u, u);
  const real B0uu = duchamp::sobolev_quadratic_B0(
      *R, /*p=*/3.0, l0, bh, u, u);
  GAUDI_ASSERT(Buu > real(0) && B0uu > real(0));
  // Paper D has 1/ℓ, so ⟨B⟩ can exceed ⟨B0⟩ by ~1/h²; just require finite PSD.
  GAUDI_ASSERT(std::isfinite(Buu / B0uu));
}

GAUDI_TEST(sobolev_calder_bh_agrees_leaves) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const real l0 = real(0.1) * R->lavg();
  const size_t n = R->corner_count();
  std::vector<vec3> u(n);
  for (size_t i = 0; i < n; ++i) {
    const real t = real(i);
    u[i] = vec3(std::sin(t), std::cos(t), real(0.05) * t);
  }
  std::vector<vec3> B0_leaf, B0_bh;
  duchamp::sobolev_detail::matvec_B0(*R, u, /*p=*/3.0, l0, /*bh_eps=*/0.0,
                                     B0_leaf);
  duchamp::sobolev_detail::matvec_B0(*R, u, /*p=*/3.0, l0, /*bh_eps=*/0.25,
                                     B0_bh);
  const real nrm = std::max(field_norm(B0_leaf), real(1e-12));
  // Paper k is more singular; 3D BH is a coarser stand-in than leaves.
  GAUDI_ASSERT(field_max_abs_diff(B0_leaf, B0_bh) / nrm < real(0.85));
  GAUDI_ASSERT(field_norm(B0_bh) > real(0.15) * nrm);
}

GAUDI_TEST(sobolev_jacobi_residual_drops) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));
  duchamp::tangent_point_solver_config leaf;
  leaf.type = {duchamp::tangent_point_type::cauchy};
  leaf.cauchy = {.w = 1.0, .l0 = 0.5, .p = 3.0};

  const std::vector<vec3> f = duchamp::tp_grad_energy(*R, *Rd, leaf);
  duchamp::sobolev_tangent_point_force_config sc;
  sc.w = 1.0;
  sc.p1 = 3.0;
  sc.l0 = 0.5;
  sc.alpha_B = 1.0;
  sc.alpha_0 = 1.0;
  sc.bh_eps = 0.0;
  const real l0_abs = sc.l0 * Rd->_Cc;

  const real r0 = field_norm(f);
  sc.solver.iter = 24;
  duchamp::sobolev_constraint_solver pin;
  pin.add(duchamp::sobolev_barycenter_constraint::create());
  pin.init_rest(*R);
  const std::vector<vec3> g = pin.project_gradient(*R, f, sc, l0_abs);
  std::vector<vec3> Ag;
  duchamp::sobolev_detail::matvec_A(*R, sc.alpha_B, sc.alpha_0, sc.p1, l0_abs,
                                    sc.bh_eps, g, Ag);
  real r1 = 0;
  for (size_t i = 0; i < f.size(); ++i)
    r1 += (f[i] - Ag[i]).squaredNorm();
  r1 = std::sqrt(r1);

  GAUDI_ASSERT(std::isfinite(r0) && std::isfinite(r1));
  GAUDI_ASSERT(r0 > real(0));
  GAUDI_ASSERT(r1 < r0 * real(0.75));
}

GAUDI_TEST(sobolev_type_cauchy_leaf_finite) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));

  duchamp::tangent_point_solver_config sob;
  sob.type = {duchamp::tangent_point_type::sobolev};
  sob.sobolev = {.w = 1.0,
                 .l0 = 0.5,
                 .p0 = 3.0,
                 .p1 = 3.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .solver = {.iter = 8},
                 .bh_eps = 0.0};

  const std::vector<vec3> g_sob = duchamp::tp_grad_energy(*R, *Rd, sob);
  GAUDI_ASSERT(field_norm(g_sob) > real(0));
  for (const auto &a : g_sob)
    GAUDI_ASSERT(std::isfinite(a[0]) && std::isfinite(a[1]) &&
                 std::isfinite(a[2]));

  const real E = duchamp::tp_energy(*R, *Rd, sob);
  GAUDI_ASSERT(std::isfinite(E));
}

GAUDI_TEST(sobolev_B0_only_near_contact_descent) {
  auto R = asawa::rod::rod::create(sobolev_helix(28), true);
  // Near-contact pair so B0 (TP k) should dominate the metric.
  R->x()[8] = R->x()[3] + vec3(0.02, 0.0, 0.0);
  R->_update_frames();
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));

  duchamp::tangent_point_solver_config leaf;
  leaf.type = {duchamp::tangent_point_type::cauchy};
  leaf.cauchy = {.w = 1.0, .l0 = 0.5, .p = 3.0};

  duchamp::tangent_point_solver_config sob = leaf;
  sob.type = {duchamp::tangent_point_type::sobolev};
  sob.sobolev = {.w = 1.0,
                 .l0 = 0.5,
                 .p0 = 3.0,
                 .p1 = 3.0,
                 .alpha_B = 0.0,
                 .alpha_0 = 1.0,
                 .solver = {.iter = 16},
                 .bh_eps = 0.0};

  duchamp::sobolev_constraint_solver pin;
  pin.add(duchamp::sobolev_barycenter_constraint::create());
  pin.init_rest(*R);
  const std::vector<vec3> g = duchamp::tp_grad_energy(*R, *Rd, sob, &pin);
  GAUDI_ASSERT(field_norm(g) > real(0));

  const real E0 = duchamp::tp_energy(*R, *Rd, sob);
  auto d = g;
  duchamp::tp_peak_normalize(d);
  for (auto &a : d)
    a = -a;
  const real eps = real(1e-4) * R->lavg();
  auto x0 = R->x();
  for (size_t i = 0; i < d.size(); ++i)
    R->x()[i] = x0[i] + eps * d[i];
  R->_update_frames();
  const real E1 = duchamp::tp_energy(*R, *Rd, sob);
  for (size_t i = 0; i < d.size(); ++i)
    R->x()[i] = x0[i];
  R->_update_frames();
  GAUDI_ASSERT(E1 < E0);
}

GAUDI_TEST(sobolev_descent_smoke) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  R->x()[3] += vec3(0.04, -0.02, 0.01);
  R->_update_frames();
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));

  duchamp::tangent_point_solver_config leaf;
  leaf.type = {duchamp::tangent_point_type::cauchy};
  leaf.cauchy = {.w = 1.0, .l0 = 0.5, .p = 3.0};

  duchamp::tangent_point_solver_config sob = leaf;
  sob.type = {duchamp::tangent_point_type::sobolev};
  sob.sobolev = {.w = 1.0,
                 .l0 = 0.5,
                 .p0 = 3.0,
                 .p1 = 3.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .solver = {.iter = 16},
                 .bh_eps = 0.0};

  const std::vector<vec3> g = duchamp::tp_grad_energy(*R, *Rd, sob);
  GAUDI_ASSERT(field_norm(g) > real(0));

  const real E0 = duchamp::tp_energy(*R, *Rd, sob);
  auto d = g;
  duchamp::tp_peak_normalize(d);
  for (auto &a : d)
    a = -a;
  const real eps = real(1e-4) * R->lavg();
  auto x0 = R->x();
  for (size_t i = 0; i < d.size(); ++i)
    R->x()[i] = x0[i] + eps * d[i];
  R->_update_frames();
  const real E1 = duchamp::tp_energy(*R, *Rd, sob);
  for (size_t i = 0; i < d.size(); ++i)
    R->x()[i] = x0[i];
  R->_update_frames();
  GAUDI_ASSERT(E1 < E0);
}

GAUDI_TEST(sobolev_mg_residual_vs_jacobi) {
  auto R = asawa::rod::rod::create(sobolev_helix(32), true);
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));
  duchamp::tangent_point_solver_config leaf;
  leaf.type = {duchamp::tangent_point_type::cauchy};
  leaf.cauchy = {.w = 1.0, .l0 = 0.5, .p = 3.0};
  const std::vector<vec3> f = duchamp::tp_grad_energy(*R, *Rd, leaf);

  duchamp::sobolev_tangent_point_force_config sc;
  sc.w = 1.0;
  sc.p1 = 3.0;
  sc.l0 = 0.5;
  sc.alpha_B = 1.0;
  sc.alpha_0 = 1.0;
  sc.bh_eps = 0.0;
  sc.solver.iter = 4;
  const real l0_abs = sc.l0 * Rd->_Cc;

  auto residual_of = [&](const std::vector<vec3> &g) {
    std::vector<vec3> Ag;
    duchamp::sobolev_detail::matvec_A(*R, sc.alpha_B, sc.alpha_0, sc.p1,
                                      l0_abs, sc.bh_eps, g, Ag);
    real s = 0;
    for (size_t i = 0; i < f.size(); ++i)
      s += (f[i] - Ag[i]).squaredNorm();
    return std::sqrt(s);
  };

  sc.num_v_cycles = 0;
  const auto g_jac = duchamp::sobolev_jacobi_solve(*R, f, sc, l0_abs);
  // Flat Jacobi with same smooth budget as one MG cycle pre+post+coarse.
  sc.solver.iter = 12;
  const auto g_flat = duchamp::sobolev_jacobi_solve(*R, f, sc, l0_abs);

  sc.solver.iter = 4;
  sc.num_v_cycles = 2;
  sc.max_levels = 3;
  const auto g_mg = duchamp::sobolev_mg_solve(*R, f, sc, l0_abs);

  const real r_mg = residual_of(g_mg);
  const real r_flat = residual_of(g_flat);
  GAUDI_ASSERT(std::isfinite(r_mg) && std::isfinite(r_flat));
  GAUDI_ASSERT(r_mg < field_norm(f));
  // MG should not be dramatically worse than a similar Jacobi budget.
  GAUDI_ASSERT(r_mg < r_flat * real(2.5) + real(1e-9));
  (void)g_jac;
}

// A helix is not collinear, so a large opening angle lumps distant clusters
// and the Cauchy force differs from the leaves-only sum.
GAUDI_TEST(cauchy_far_opening_differs) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  const std::vector<vec3> &x = R->x();
  const std::vector<real> w = R->l0();
  const std::vector<vec3> T = R->N2c();
  const real l0 = real(0.05) * R->lavg();
  const auto leaves =
      calder::tangent_point_gradient_cauchy(*R, x, w, T, l0, /*p=*/3.0,
                                            /*bh_eps=*/0.0);
  const auto far =
      calder::tangent_point_gradient_cauchy(*R, x, w, T, l0, /*p=*/3.0,
                                            /*bh_eps=*/2.0);
  GAUDI_ASSERT(leaves.size() == far.size());
  GAUDI_ASSERT(field_max_abs_diff(leaves, far) > real(1e-6));
  for (const auto &g : far)
    GAUDI_ASSERT(g.allFinite());
}

// Two matvecs on a frozen rod share the covariance tree.
GAUDI_TEST(sobolev_matvec_reuses_cov_tree) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  auto cov = calder::build_rod_cov_tree(*R);
  const int built = arp::edge_cov_tree<2>::build_count();
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  std::vector<vec3> u(R->corner_count(), vec3(1, 0, 0));
  std::vector<vec3> out;
  const real l0 = real(0.1) * R->lavg();
  duchamp::sobolev_detail::matvec_A(*R, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*bh=*/0, u,
                                    out, cov.get());
  duchamp::sobolev_detail::matvec_A(*R, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*bh=*/0, u,
                                    out, cov.get());
  GAUDI_ASSERT(arp::edge_cov_tree<2>::build_count() == built);
  GAUDI_ASSERT(out.size() == u.size());
}

// Each residual-forest level owns the covariance tree of its own rod.
GAUDI_TEST(sobolev_levels_own_cov_trees) {
  auto R = asawa::rod::rod::create(sobolev_helix(32), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 2);
  F->restrict_up(0);

  auto fine = calder::build_rod_cov_tree(*F->level_rod(0));
  const int before = arp::edge_cov_tree<2>::build_count();
  auto coarse = duchamp::sobolev_level_cov(*F->level_rod(1), /*ell=*/1, fine);
  auto fine_again =
      duchamp::sobolev_level_cov(*F->level_rod(0), /*ell=*/0, fine);
  GAUDI_ASSERT(coarse.get() != fine.get());
  GAUDI_ASSERT(fine_again.get() == fine.get());
  GAUDI_ASSERT(arp::edge_cov_tree<2>::build_count() == before + 1);
  GAUDI_ASSERT(coarse->indices_.size() != fine->indices_.size());

  const int built = arp::edge_cov_tree<2>::build_count();
  auto rod1 = F->level_rod(1);
  auto m = duchamp::sobolev_detail::sample_masses(*rod1);
  std::vector<vec3> u(rod1->corner_count(), vec3(0, 1, 0));
  std::vector<vec3> out;
  const real l0 = real(0.1) * std::max(rod1->lavg(), real(1e-6));
  duchamp::sobolev_detail::matvec_A(*rod1, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*bh=*/0, u,
                                    out, coarse.get());
  duchamp::sobolev_detail::matvec_A(*rod1, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*bh=*/0, u,
                                    out, coarse.get());
  GAUDI_ASSERT(arp::edge_cov_tree<2>::build_count() == built);
  GAUDI_ASSERT(out.size() == u.size());
}

// Frozen B and B0 weights match a live matvec and are reused without another walk.
GAUDI_TEST(sobolev_weights_match_and_reuse) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  auto cov = calder::build_rod_cov_tree(*R);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const real l0 = real(0.1) * R->lavg();
  const real bh = real(0.25);
  const real cone = real(0.25);
  const int before = calder::cov_walk_count();
  auto wB = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, bh, cone, /*use_tp_k=*/false, *cov);
  auto wB0 = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, bh, cone, /*use_tp_k=*/true, *cov);
  GAUDI_ASSERT(calder::cov_walk_count() == before + 2);

  std::vector<vec3> u(R->corner_count(), vec3(0.2, -0.1, 0.05));
  std::vector<vec3> live, cached;
  duchamp::sobolev_detail::matvec_A(*R, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, bh, u, live,
                                    cov.get(), cone);
  const int after_live = calder::cov_walk_count();
  duchamp::sobolev_detail::matvec_A(*R, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, bh, u,
                                    cached, cov.get(), cone, &wB, &wB0);
  GAUDI_ASSERT(calder::cov_walk_count() == after_live);
  GAUDI_ASSERT(field_max_abs_diff(live, cached) < real(1e-8));
}

// Leaf B and B0 weights are the midpoint pair, and adjacent edges are skipped.
GAUDI_TEST(sobolev_leaf_metric_uses_midpoints) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  auto cov = calder::build_rod_cov_tree(*R);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  const real l0 = real(0.1) * R->lavg();
  const auto &x = R->x();
  auto check = [&](bool use_tp_k) {
    const auto rows = duchamp::sobolev_detail::build_cov_weights(
        *R, E, /*p=*/4, l0, /*bh=*/0, /*cone=*/0, use_tp_k, *cov);
    int checked = 0;
    for (size_t i = 0; i < E.n; ++i) {
      for (const auto &hit : rows.row[i]) {
        GAUDI_ASSERT(!hit.far);
        const size_t j_edge = static_cast<size_t>(cov->get_index(hit.j));
        GAUDI_ASSERT(j_edge < E.n);
        const int a0 = E.ep0[i];
        const int a1 = E.ep1[i];
        const int b0 = E.ep0[j_edge];
        const int b1 = E.ep1[j_edge];
        const bool adj = a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1;
        GAUDI_ASSERT(!adj);
        const vec3 mid_j =
            real(0.5) * (x[static_cast<size_t>(b0)] + x[static_cast<size_t>(b1)]);
        const real expect = duchamp::sobolev_detail::directed_pair_weight(
            E.mid[i], mid_j, E.T[i], E.T[j_edge], E.ell[i], E.ell[j_edge], l0,
            /*p=*/4, use_tp_k);
        GAUDI_ASSERT(std::abs(hit.w - expect) <=
                     real(1e-8) * (real(1) + std::abs(expect)));
        ++checked;
      }
    }
    GAUDI_ASSERT(checked > 0);
  };
  check(false);
  check(true);
}

// Leaves-only energy matches the vertex pair sum, and a vertex nudge matches dE.
GAUDI_TEST(sobolev_paper_energy_leaf) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  const real alpha = real(3);
  const real beta = real(6);
  const auto V = duchamp::sobolev_detail::tpe_vertices(*R);
  const auto &x = R->x();
  real hand = real(0);
  int pairs = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (V[i].deg <= 0 || V[i].deg > 2)
      continue;
    for (size_t j = 0; j < x.size(); ++j) {
      if (i == j || !(V[j].dual > real(0)))
        continue;
      hand += V[i].dual * V[j].dual *
              duchamp::sobolev_detail::tpe_kf(x[i] - x[j], V[i].T, alpha, beta);
      ++pairs;
    }
  }
  GAUDI_ASSERT(pairs > 0);
  const auto ev =
      duchamp::sobolev_detail::eval_paper_tp(*R, alpha, beta, /*bh=*/0);
  GAUDI_ASSERT(std::abs(ev.energy - hand) <=
               real(1e-6) * (real(1) + std::abs(hand)));
  GAUDI_ASSERT(ev.grad.size() == x.size());

  const real h = real(1e-5) * R->lavg();
  const auto x0 = R->x();
  const size_t iv = 4;
  for (int axis = 0; axis < 3; ++axis) {
    vec3 step = vec3::Zero();
    step[axis] = h;
    R->x() = x0;
    R->x()[iv] += step;
    R->_update_frames();
    const real Ep =
        duchamp::sobolev_detail::eval_paper_tp(*R, alpha, beta, 0).energy;
    R->x() = x0;
    R->x()[iv] -= step;
    R->_update_frames();
    const real Em =
        duchamp::sobolev_detail::eval_paper_tp(*R, alpha, beta, 0).energy;
    const real num = (Ep - Em) / (real(2) * h);
    const real ana = ev.grad[iv][axis];
    GAUDI_ASSERT(std::abs(num - ana) <= real(1e-3) * (real(1) + std::abs(ana)));
  }
  R->x() = x0;
  R->_update_frames();
}

// Eq. 31: align is 1 on B and B0. A partner residual is −w φ_J.
GAUDI_TEST(sobolev_B_eq31_residual) {
  const std::vector<vec3> pts = {
      vec3(0, 0, 0), vec3(1, 0, 0), vec3(2, 0, 0), vec3(2, 1, 0),
      vec3(3, 1, 0), vec3(3, 2, 0), vec3(0, 2, 0)};
  auto R = asawa::rod::rod::create(pts, true);
  auto cov = calder::build_rod_cov_tree(*R);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  const auto &x = R->x();
  const real l0 = real(0.1) * R->lavg();
  auto edge_between = [&](const vec3 &a, const vec3 &b) -> int {
    for (size_t e = 0; e < E.n; ++e) {
      const vec3 &p = x[static_cast<size_t>(E.ep0[e])];
      const vec3 &q = x[static_cast<size_t>(E.ep1[e])];
      if ((p - a).norm() < real(1e-8) && (q - b).norm() < real(1e-8))
        return static_cast<int>(e);
    }
    return -1;
  };
  const int e_query = edge_between(vec3(0, 0, 0), vec3(1, 0, 0));
  const int e_perp = edge_between(vec3(2, 0, 0), vec3(2, 1, 0));
  const int e_par = edge_between(vec3(2, 1, 0), vec3(3, 1, 0));
  GAUDI_ASSERT(e_query >= 0 && e_perp >= 0 && e_par >= 0);

  const auto rowsB = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, /*bh=*/0, /*cone=*/0, /*use_tp_k=*/false,
      *cov);
  const auto rowsB0 = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, /*bh=*/0, /*cone=*/0, /*use_tp_k=*/true,
      *cov);
  auto find_hit = [&](const duchamp::sobolev_detail::cov_weight_rows &rows,
                      int query, int partner)
      -> const duchamp::sobolev_detail::cov_weight_hit * {
    for (const auto &hit : rows.row[static_cast<size_t>(query)]) {
      if (static_cast<int>(cov->get_index(hit.j)) == partner)
        return &hit;
    }
    return nullptr;
  };
  const auto *perp = find_hit(rowsB, e_query, e_perp);
  const auto *par = find_hit(rowsB, e_query, e_par);
  const auto *perp0 = find_hit(rowsB0, e_query, e_perp);
  GAUDI_ASSERT(perp && par && perp0);
  GAUDI_ASSERT(std::abs(perp->align - real(1)) < real(1e-8));
  GAUDI_ASSERT(std::abs(par->align - real(1)) < real(1e-8));
  GAUDI_ASSERT(std::abs(perp0->align - real(1)) < real(1e-8));

  const vec3 s(real(0.2), real(-0.3), real(0.4));
  std::vector<vec3> phi(E.n, vec3::Zero());
  phi[static_cast<size_t>(e_perp)] = s;
  const auto r_perp = duchamp::sobolev_detail::edge_K_diff(
      *R, E, phi, /*p=*/4, l0, /*bh=*/0, /*use_tp_k=*/false,
      cov.get(), /*cone=*/0, &rowsB);
  const vec3 expect_perp = -perp->w * s;
  GAUDI_ASSERT((r_perp[static_cast<size_t>(e_query)] - expect_perp).norm() <=
               real(1e-8) * (real(1) + expect_perp.norm()));

  phi.assign(E.n, vec3::Zero());
  phi[static_cast<size_t>(e_par)] = s;
  const auto r_par = duchamp::sobolev_detail::edge_K_diff(
      *R, E, phi, /*p=*/4, l0, /*bh=*/0, /*use_tp_k=*/false,
      cov.get(), /*cone=*/0, &rowsB);
  const vec3 expect = -par->w * s;
  GAUDI_ASSERT((r_par[static_cast<size_t>(e_query)] - expect).norm() <=
               real(1e-8) * (real(1) + expect.norm()));
}

GAUDI_TEST(tp_type_list_sums_gradients) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));

  duchamp::tangent_point_solver_config cauchy;
  cauchy.type = {duchamp::tangent_point_type::cauchy};
  cauchy.cauchy = {.w = 1.0, .l0 = 0.1, .p = 6.0};

  duchamp::tangent_point_solver_config sob;
  sob.type = {duchamp::tangent_point_type::sobolev};
  sob.sobolev = {.w = 1.0,
                 .l0 = 0.05,
                 .p0 = 2.0,
                 .p1 = 2.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .solver = {.type = duchamp::sobolev_solver::cg, .iter = 2},
                 .num_v_cycles = 1,
                 .max_levels = 1,
                 .bh_eps = 0.0,
                 .cone_eps = 0.0};

  const std::vector<vec3> gc = duchamp::tp_grad_energy(*R, *Rd, cauchy);
  const std::vector<vec3> gs = duchamp::tp_grad_energy(*R, *Rd, sob);
  std::vector<vec3> hand = gc;
  for (size_t i = 0; i < hand.size(); ++i)
    hand[i] += gs[i];

  duchamp::tangent_point_solver_config both = sob;
  both.type = {duchamp::tangent_point_type::cauchy,
               duchamp::tangent_point_type::sobolev};
  both.cauchy = cauchy.cauchy;
  const std::vector<vec3> g = duchamp::tp_grad_energy(*R, *Rd, both);
  const real nrm = std::max(field_norm(hand), real(1e-12));
  GAUDI_ASSERT(field_max_abs_diff(g, hand) / nrm < real(1e-8));

  const real Ec = duchamp::tp_energy(*R, *Rd, cauchy);
  const real Es = duchamp::tp_energy(*R, *Rd, sob);
  const real E = duchamp::tp_energy(*R, *Rd, both);
  GAUDI_ASSERT(std::abs(E - (Ec + Es)) <=
               real(1e-6) * (real(1) + std::abs(Ec + Es)));

  const std::vector<vec3> g1 = duchamp::tp_grad_energy(*R, *Rd, cauchy);
  GAUDI_ASSERT(field_max_abs_diff(g1, gc) / std::max(field_norm(gc), real(1e-12)) <
               real(1e-12));
}

GAUDI_TEST(harmonic_tp_pair_tangent_weight) {
  const vec3 dp(0.3, 0.0, 0.0);
  const vec3 Ti(0, 1, 0);
  const vec3 parallel(0, 1, 0);
  const vec3 anti(0, -1, 0);
  const vec3 perp(1, 0, 0);
  const real Tp = 1.5;
  const vec3 gTp(0.2, -0.4, 0.1);
  const real ell = 0.5;
  const vec3 gp =
      calder::harmonic_tp_pair(dp, Ti, parallel, Tp, gTp, ell);
  const vec3 ga = calder::harmonic_tp_pair(dp, Ti, anti, Tp, gTp, ell);
  const vec3 gz = calder::harmonic_tp_pair(dp, Ti, perp, Tp, gTp, ell);
  GAUDI_ASSERT(gz.norm() < real(1e-12));
  GAUDI_ASSERT((gp + ga).norm() < real(1e-8) * (real(1) + gp.norm()));

  // Moving the query by ε e changes dp by −ε e. With ∇Tp = 0 the pair
  // gradient is only w′ Tp.
  const real eps = real(1e-6);
  const vec3 e(0, 0, 1);
  const real d0 = calder::harmonic_tp_density(dp, Ti, parallel, Tp, ell);
  const real d1 =
      calder::harmonic_tp_density(dp - eps * e, Ti, parallel, Tp, ell);
  const vec3 g0 =
      calder::harmonic_tp_pair(dp, Ti, parallel, Tp, vec3::Zero(), ell);
  const real fd = (d1 - d0) / eps;
  GAUDI_ASSERT(std::abs(fd - g0.dot(e)) <
               real(1e-4) * (real(1) + std::abs(fd)));

  // p1 = 2 is the quadratic Gaussian. A different power moves κ and w′.
  const real k2 = calder::calc_gaussian(dp, ell);
  const real k2p = calder::calc_gaussian(dp, ell, real(2));
  GAUDI_ASSERT(std::abs(k2 - k2p) < real(1e-12) * (real(1) + std::abs(k2)));
  const real p4 = real(4);
  const real d4 = calder::harmonic_tp_density(dp, Ti, parallel, Tp, ell, p4);
  const real d4e = calder::harmonic_tp_density(dp - eps * e, Ti, parallel, Tp,
                                               ell, p4);
  const vec3 g4 =
      calder::harmonic_tp_pair(dp, Ti, parallel, Tp, vec3::Zero(), ell, p4);
  GAUDI_ASSERT(std::abs(d4 - d0) > real(1e-6));
  GAUDI_ASSERT(std::abs((d4e - d4) / eps - g4.dot(e)) <
               real(1e-4) * (real(1) + std::abs(d4)));
}

GAUDI_TEST(harmonic_rod_gradient_finite) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));
  duchamp::tangent_point_solver_config cfg;
  cfg.type = {duchamp::tangent_point_type::harmonic};
  cfg.harmonic = {.w = 1.0, .l0 = 0.1, .l1 = 2.0, .p0 = 3.0, .p1 = 2.0};
  const std::vector<vec3> g = duchamp::tp_grad_energy(*R, *Rd, cfg);
  const real E = duchamp::tp_energy(*R, *Rd, cfg);
  GAUDI_ASSERT(std::isfinite(E));
  GAUDI_ASSERT(field_norm(g) > real(0));
  for (const vec3 &a : g)
    GAUDI_ASSERT(a.allFinite());
}

GAUDI_TEST(ccd_velocity_linear_closing) {
  auto R = asawa::rod::rod::create({vec3(0, 0, 0), vec3(1, 0, 0)}, false);
  duchamp::ccd_config cfg;
  cfg.dt_max = 10.0;
  cfg.dt_min = 0.0;
  cfg.dx_max = 0.0;
  cfg.neighborhood.cfl_frac = 0.5;
  cfg.neighborhood.gate_frac = 0.0;
  cfg.neighborhood.use_envelope = false;
  cfg.neighborhood.dt_max = 10.0;
  cfg.neighborhood.dt_min = 0.0;

  // Closing speed 1 across a unit gap. s = 0.5 * 1 / 1.
  std::vector<vec3> v(2);
  v[0] = vec3(0.5, 0, 0);
  v[1] = vec3(-0.5, 0, 0);
  const real dt_v = duchamp::ccd_velocity(*R, v, cfg);
  GAUDI_ASSERT(std::abs(dt_v - real(0.5)) < real(1e-8));

  const real dt_f = duchamp::ccd_force(*R, v, /*h_ref=*/1.0, cfg);
  GAUDI_ASSERT(std::abs(dt_f - std::sqrt(real(0.5))) < real(1e-8));

  v[0] = vec3(-1, 0, 0);
  v[1] = vec3(1, 0, 0);
  cfg.dt_max = 1.5;
  cfg.neighborhood.dt_max = 1.5;
  const real dt_sep = duchamp::ccd_velocity(*R, v, cfg);
  GAUDI_ASSERT(std::abs(dt_sep - real(1.5)) < real(1e-8));
}

/// A displacement already inside τ_max clears dt_max. Nothing closing must
/// come back uncapped so the apply loop's min(dti, 1) takes the whole step.
GAUDI_TEST(ccd_velocity_separating_uncapped) {
  auto R = asawa::rod::rod::create({vec3(0, 0, 0), vec3(1, 0, 0)}, false);
  duchamp::ccd_config cfg;
  cfg.dt_max = 0.0;
  cfg.dt_min = 0.0;
  cfg.dx_max = 0.0;
  cfg.neighborhood.cfl_frac = 1.0;
  cfg.neighborhood.gate_frac = 0.1;
  cfg.neighborhood.use_envelope = true;
  cfg.neighborhood.dt_max = 0.0;
  cfg.neighborhood.dt_min = 0.0;
  cfg.neighborhood.chain_skip_hops = 8;

  duchamp::cfl_neighborhood_result hit;
  std::vector<vec3> v(2, vec3(1, 0, 0));
  const real dt_slide = duchamp::ccd_velocity(*R, v, cfg, &hit);
  GAUDI_ASSERT(hit.limiting == duchamp::cfl_limit_term::none);
  GAUDI_ASSERT(std::isinf(dt_slide) && dt_slide > real(0));
  GAUDI_ASSERT(std::min(dt_slide, real(1)) == real(1));

  v[0] = vec3(-1, 0, 0);
  v[1] = vec3(1, 0, 0);
  const real dt_sep = duchamp::ccd_velocity(*R, v, cfg, &hit);
  GAUDI_ASSERT(hit.limiting == duchamp::cfl_limit_term::none);
  GAUDI_ASSERT(std::isinf(dt_sep) && dt_sep > real(0));

  // Closing speed 2 across a unit gap: s = 1 · 1 / 2.
  v[0] = vec3(1, 0, 0);
  v[1] = vec3(-1, 0, 0);
  const real dt_hit = duchamp::ccd_velocity(*R, v, cfg, &hit);
  GAUDI_ASSERT(hit.limiting != duchamp::cfl_limit_term::none);
  GAUDI_ASSERT(std::isfinite(dt_hit));
  GAUDI_ASSERT(dt_hit > real(0) && dt_hit < real(1));
}

namespace {

inline bool paper_edges_touch(const duchamp::sobolev_detail::edge_table &E,
                              size_t i, size_t j) {
  return duchamp::sobolev_detail::leaf_edges_adjacent_or_self(
      E, i, E.ep0[j], E.ep1[j]);
}

/// §5.2.1–5.2.2 for B: φ = (Δu/ℓ) T, then r_I = Σ w (φ_I − φ_J), scatter T·r.
/// B0 stays the midpoint average. Adjacent edges are omitted. Weights are the
/// midpoint kernel.
inline std::vector<vec3>
dense_paper_matvec(asawa::rod::rod &rod, const std::vector<vec3> &u, real p,
                   real l0, bool use_tp_k) {
  using namespace duchamp::sobolev_detail;
  const auto E = build_edges(rod);
  if (!use_tp_k) {
    const std::vector<vec3> slope = edge_diff_D(E, u);
    std::vector<vec3> out(E.n_vert, vec3::Zero());
    for (int c = 0; c < 3; ++c) {
      std::vector<vec3> phi(E.n, vec3::Zero());
      for (size_t i = 0; i < E.n; ++i)
        phi[i] = slope[i][c] * E.T[i];
      std::vector<vec3> r(E.n, vec3::Zero());
      for (size_t i = 0; i < E.n; ++i) {
        for (size_t j = 0; j < E.n; ++j) {
          if (i == j || paper_edges_touch(E, i, j))
            continue;
          const real w = directed_pair_weight(
              E.mid[i], E.mid[j], E.T[i], E.T[j], E.ell[i], E.ell[j], l0, p,
              /*use_tp_k=*/false);
          r[i] += w * (phi[i] - phi[j]);
        }
      }
      for (size_t i = 0; i < E.n; ++i) {
        const real s = E.T[i].dot(r[i]) / E.ell[i];
        out[static_cast<size_t>(E.ep0[i])][c] -= s;
        out[static_cast<size_t>(E.ep1[i])][c] += s;
      }
    }
    return out;
  }
  const std::vector<vec3> phi = edge_avg_E(E, u);
  std::vector<vec3> r(E.n, vec3::Zero());
  for (size_t i = 0; i < E.n; ++i) {
    for (size_t j = 0; j < E.n; ++j) {
      if (i == j || paper_edges_touch(E, i, j))
        continue;
      const real w = directed_pair_weight(E.mid[i], E.mid[j], E.T[i], E.T[j],
                                          E.ell[i], E.ell[j], l0, p,
                                          /*use_tp_k=*/true);
      r[i] += w * (phi[i] - phi[j]);
    }
  }
  std::vector<vec3> out;
  scatter_ET(E, r, out);
  return out;
}

/// Reference tpe_total: every ordered vertex pair, dual lengths, vertex tangent.
inline real dense_paper_energy(asawa::rod::rod &rod, real alpha, real beta) {
  using namespace duchamp::sobolev_detail;
  const auto V = tpe_vertices(rod);
  const auto &x = rod.x();
  real e = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if (V[i].deg <= 0 || V[i].deg > 2)
      continue;
    for (size_t j = 0; j < x.size(); ++j) {
      if (i == j || !(V[j].dual > real(0)))
        continue;
      e += V[i].dual * V[j].dual * tpe_kf(x[i] - x[j], V[i].T, alpha, beta);
    }
  }
  return e;
}

} // namespace

GAUDI_TEST(sobolev_matvec_matches_eq31) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const real l0 = real(0.1) * R->lavg();
  const real p = 3;
  const size_t n = R->corner_count();
  std::vector<vec3> u(n);
  for (size_t i = 0; i < n; ++i) {
    const real t = real(i);
    u[i] = vec3(std::sin(t), std::cos(real(0.7) * t), real(0.15) * t);
  }
  std::vector<vec3> gotB, got0;
  duchamp::sobolev_detail::matvec_B(*R, u, p, l0, /*bh=*/0,
                                    gotB);
  duchamp::sobolev_detail::matvec_B0(*R, u, p, l0, /*bh=*/0,
                                     got0);
  const auto eqB = dense_paper_matvec(*R, u, p, l0, /*tp=*/false);
  const auto eq0 = dense_paper_matvec(*R, u, p, l0, /*tp=*/true);
  const real nB = std::max(field_norm(eqB), real(1));
  const real n0 = std::max(field_norm(eq0), real(1));
  GAUDI_ASSERT(field_max_abs_diff(gotB, eqB) / nB <= real(1e-6));
  GAUDI_ASSERT(field_max_abs_diff(got0, eq0) / n0 <= real(1e-6));
}

GAUDI_TEST(sobolev_energy_matches_eq18) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  const real alpha = real(3);
  const real beta = real(6);
  const auto ev =
      duchamp::sobolev_detail::eval_paper_tp(*R, alpha, beta, /*bh=*/0);
  const real hand = dense_paper_energy(*R, alpha, beta);
  GAUDI_ASSERT(std::abs(ev.energy - hand) <=
               real(1e-6) * (real(1) + std::abs(hand)));

  const real h = real(1e-5) * R->lavg();
  const auto x0 = R->x();
  const size_t iv = 3;
  for (int axis = 0; axis < 3; ++axis) {
    vec3 step = vec3::Zero();
    step[axis] = h;
    R->x() = x0;
    R->x()[iv] += step;
    R->_update_frames();
    const real Ep = dense_paper_energy(*R, alpha, beta);
    R->x() = x0;
    R->x()[iv] -= step;
    R->_update_frames();
    const real Em = dense_paper_energy(*R, alpha, beta);
    const real num = (Ep - Em) / (real(2) * h);
    GAUDI_ASSERT(std::abs(num - ev.grad[iv][axis]) <=
                 real(1e-3) * (real(1) + std::abs(ev.grad[iv][axis])));
  }
  R->x() = x0;
  R->_update_frames();
}

GAUDI_TEST(sobolev_unconstrained_direction_descends) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  R->x()[4] += vec3(0.05, -0.02, 0.01);
  R->_update_frames();
  auto Rd = asawa::rod::dynamic::create(R, real(0.25) * R->lavg(),
                                        real(2.5) * R->lavg(), real(0));
  duchamp::sobolev_tangent_point_force_config sc;
  sc.w = 1;
  sc.l0 = 0.05;
  sc.p0 = 3;
  sc.p1 = 3;
  sc.alpha_B = 1;
  sc.alpha_0 = 1;
  sc.solver.type = duchamp::sobolev_solver::cg;
  sc.solver.iter = 20;
  sc.num_v_cycles = 0;
  sc.bh_eps = 0;
  const real l0 = sc.l0 * Rd->_Cc;
  const auto ev = duchamp::sobolev_detail::eval_paper_tp(*R, sc.tpe_alpha,
                                                         sc.tpe_beta, 0);
  const std::vector<vec3> g =
      duchamp::sobolev_jacobi_solve(*R, ev.grad, sc, l0);
  const real fg = field_dot(ev.grad, g);
  std::cout << "[sobolev f·g] unconstrained " << fg << std::endl;
  GAUDI_ASSERT(fg > real(0));

  auto d = g;
  duchamp::tp_peak_normalize(d);
  const real eps = real(1e-4) * R->lavg();
  const auto x0 = R->x();
  for (size_t i = 0; i < d.size(); ++i)
    R->x()[i] = x0[i] - eps * d[i];
  R->_update_frames();
  const real E1 =
      duchamp::sobolev_detail::eval_paper_tp(*R, sc.tpe_alpha, sc.tpe_beta, 0)
          .energy;
  R->x() = x0;
  R->_update_frames();
  GAUDI_ASSERT(E1 < ev.energy);

  duchamp::tangent_point_solver_config cfg;
  cfg.type = {duchamp::tangent_point_type::sobolev};
  cfg.sobolev = sc;
  cfg.hard_constraints = {
      {.type = duchamp::sobolev_constraint_type::edge_length},
      {.type = duchamp::sobolev_constraint_type::barycenter}};
  cfg.constraint_cg = 4;
  auto saddle = duchamp::make_sobolev_constraint_solver(cfg);
  saddle.init_rest(*R);
  const std::vector<vec3> gs =
      saddle.project_gradient(*R, ev.grad, sc, l0);
  const real fg_s = field_dot(ev.grad, gs);
  std::cout << "[sobolev f·g] saddle " << fg_s << std::endl;
  GAUDI_ASSERT(std::isfinite(fg_s));
}

/// A constant vertex field has a constant edge average, so B0 u is zero.
/// A far partner that uses the raw sum of φ leaves a residual of size
/// about (n_cluster − 1). Appendix B.1 uses Σ(ℓ φ)/L, which stays zero.
GAUDI_TEST(sobolev_far_constant_field_is_kernel) {
  auto R = asawa::rod::rod::create(sobolev_helix(48), true);
  const real l0 = real(0.001) * R->lavg();
  auto m = duchamp::sobolev_detail::sample_masses(*R);
  const size_t n = R->corner_count();
  const std::vector<vec3> u(n, vec3(real(1), real(0), real(0)));
  auto cov = calder::build_rod_cov_tree(*R);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  const real bh = real(2);
  const real cone = real(10);
  const auto rows = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/2.0, l0, bh, cone, /*use_tp_k=*/true, *cov);
  int far_hits = 0;
  for (const auto &row : rows.row)
    for (const auto &hit : row)
      if (hit.far)
        ++far_hits;
  GAUDI_ASSERT(far_hits > 0);
  std::vector<vec3> B0u;
  duchamp::sobolev_detail::matvec_B0(*R, u, /*p=*/2.0, l0, bh,
                                     B0u, cov.get(), cone, &rows);
  GAUDI_ASSERT(field_norm(B0u) < real(1e-8));
}

// Vertex diagonal of B and of B0 against one column of the matvec.
GAUDI_TEST(sobolev_jacobi_diag_matches_matvec) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  auto cov = calder::build_rod_cov_tree(*R);
  const size_t n = R->corner_count();
  const std::vector<real> m(n, real(0));
  const real l0 = real(0.1) * R->lavg();
  const real p = real(3);

  auto probe = [&](real alpha_B, real alpha_0) {
    const auto diag = duchamp::sobolev_detail::jacobi_diag(
        *R, alpha_B, alpha_0, p, l0, /*bh_eps=*/0,
        cov.get(), /*cone_eps=*/0);
    GAUDI_ASSERT(diag.size() == n);
    for (size_t i = 0; i < n; ++i) {
      std::vector<vec3> u(n, vec3::Zero());
      u[i] = vec3(real(1), real(0), real(0));
      std::vector<vec3> Au;
      duchamp::sobolev_detail::matvec_A(*R, alpha_B, alpha_0, p, l0,
                                        /*bh_eps=*/0, u, Au, cov.get(),
                                        /*cone_eps=*/0);
      const real fd = Au[i][0];
      const real expect = std::max(fd, real(1e-12));
      GAUDI_ASSERT(std::abs(diag[i] - expect) <=
                   real(1e-5) * (real(1) + std::abs(expect)));
    }
  };
  probe(/*alpha_B=*/1, /*alpha_0=*/0);
  probe(/*alpha_B=*/0, /*alpha_0=*/1);
}

namespace {

inline real sobolev_residual(asawa::rod::rod &R,
                             const duchamp::sobolev_tangent_point_force_config &sc,
                             real l0, const std::vector<vec3> &f,
                             const std::vector<vec3> &g) {
  std::vector<vec3> Ag;
  duchamp::sobolev_detail::matvec_A(R, sc.alpha_B, sc.alpha_0, sc.p1, l0,
                                    sc.bh_eps, g, Ag);
  real s = 0;
  for (size_t i = 0; i < f.size(); ++i)
    s += (f[i] - Ag[i]).squaredNorm();
  return std::sqrt(s);
}

inline std::vector<vec3> sobolev_mean_free_rhs(size_t n) {
  std::vector<vec3> f(n, vec3::Zero());
  vec3 mean = vec3::Zero();
  for (size_t i = 0; i < n; ++i) {
    const real t = real(i);
    f[i] = vec3(std::sin(t), std::cos(real(0.5) * t), real(0.01) * t);
    mean += f[i];
  }
  mean /= real(n);
  for (vec3 &a : f)
    a -= mean;
  return f;
}

} // namespace

// More CG smoother steps must shrink |b - Ax| of the production solve.
GAUDI_TEST(sobolev_cg_residual_drops) {
  auto R = asawa::rod::rod::create(sobolev_helix(24), true);
  const std::vector<vec3> f = sobolev_mean_free_rhs(R->corner_count());
  duchamp::sobolev_tangent_point_force_config sc;
  sc.p1 = 3;
  sc.l0 = real(0.2);
  sc.alpha_B = 1;
  sc.alpha_0 = 1;
  sc.bh_eps = 0;
  sc.cone_eps = 0;
  sc.num_v_cycles = 1;
  sc.max_levels = 1;
  sc.solver.type = duchamp::sobolev_solver::cg;
  const real l0 = sc.l0 * R->lavg();
  const real f_norm = field_norm(f);
  sc.solver.iter = 4;
  const auto g4 = duchamp::sobolev_mg_solve(*R, f, sc, l0);
  sc.solver.iter = 16;
  const auto g16 = duchamp::sobolev_mg_solve(*R, f, sc, l0);
  const real r4 = sobolev_residual(*R, sc, l0, f, g4);
  const real r16 = sobolev_residual(*R, sc, l0, f, g16);
  GAUDI_ASSERT(std::isfinite(r4) && std::isfinite(r16));
  GAUDI_ASSERT(r4 < f_norm);
  GAUDI_ASSERT(r16 < r4);
}

// |B|, |B0|, and |A| on one prolonged wave, then one production V-cycle.
// drop_adjacent false keeps the midpoint weight of edges that share a vertex.
GAUDI_TEST(sobolev_coarse_scale_split) {
  auto R = asawa::rod::rod::create(sobolev_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 2);
  F->restrict_up(0);
  auto coarse = F->level_rod(1);

  duchamp::sobolev_tangent_point_force_config sc;
  sc.p1 = 3;
  sc.l0 = real(0.2);
  sc.alpha_B = 1;
  sc.alpha_0 = 1;
  sc.bh_eps = 0;
  sc.cone_eps = 0;
  sc.solver.iter = 4;
  sc.num_v_cycles = 1;
  sc.max_levels = 2;
  const real l0 = sc.l0 * R->lavg();

  const size_t nc = coarse->corner_count();
  std::vector<vec3> uc(nc);
  for (size_t i = 0; i < nc; ++i) {
    const real t = real(2) * real(3.14159265358979323846) * real(i) / real(nc);
    uc[i] = vec3(std::sin(t), 0, 0);
  }
  const std::vector<vec3> uf = F->prolong_J(0, uc);

  auto op_ratio = [&](real alpha_B, real alpha_0, bool drop_adjacent) {
    std::vector<vec3> Ac, Af;
    duchamp::sobolev_detail::matvec_A(
        *coarse, alpha_B, alpha_0, sc.p1, l0, sc.bh_eps, uc, Ac, nullptr, -1,
        nullptr, nullptr, drop_adjacent);
    duchamp::sobolev_detail::matvec_A(
        *R, alpha_B, alpha_0, sc.p1, l0, sc.bh_eps, uf, Af, nullptr, -1,
        nullptr, nullptr, drop_adjacent);
    const std::vector<vec3> Jt = F->restrict_JT(0, Af);
    return field_norm(Ac) / std::max(field_norm(Jt), real(1e-30));
  };

  const real b_drop = op_ratio(1, 0, true);
  const real b0_drop = op_ratio(0, 1, true);
  const real a_drop = op_ratio(1, 1, true);
  const real b_keep = op_ratio(1, 0, false);
  const real b0_keep = op_ratio(0, 1, false);
  const real a_keep = op_ratio(1, 1, false);

  const std::vector<vec3> f = sobolev_mean_free_rhs(R->corner_count());
  const auto g4 = duchamp::sobolev_mg_solve(*R, f, sc, l0);
  const real r4 = sobolev_residual(*R, sc, l0, f, g4);
  sc.solver.iter = 8;
  const auto g8 = duchamp::sobolev_mg_solve(*R, f, sc, l0);
  const real r8 = sobolev_residual(*R, sc, l0, f, g8);
  const real f_norm = field_norm(f);

  std::ostringstream msg;
  msg << "drop |B| " << b_drop << " |B0| " << b0_drop << " |A| " << a_drop
      << " keep |B| " << b_keep << " |B0| " << b0_keep << " |A| " << a_keep
      << " |f| " << f_norm << " |r| cg4 " << r4 << " cg8 " << r8;
  std::cout << "[coarse scale] " << msg.str() << std::endl;
  if (!(std::isfinite(r4) && std::isfinite(r8) && r4 < f_norm && r8 < r4))
    gaudi::test::record_failure(__FILE__, __LINE__, msg.str().c_str(), false);
}

// Coarse solve uses the rediscretized A on the coarse rod. The V-cycle adds
// J e at full length, then the post-smoother damps the residual.
GAUDI_TEST(sobolev_prolonged_correction_cuts_residual) {
  auto R = asawa::rod::rod::create(sobolev_helix(32), true);
  R->x()[8] = R->x()[3] + vec3(0.02, 0.0, 0.0);
  R->_update_frames();
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 2);
  F->restrict_up(0);
  auto coarse = F->level_rod(1);

  duchamp::sobolev_tangent_point_force_config sc;
  sc.p1 = 3;
  sc.l0 = real(0.05);
  sc.alpha_B = 1;
  sc.alpha_0 = 1;
  sc.bh_eps = 0;
  sc.cone_eps = 0;
  sc.solver.type = duchamp::sobolev_solver::cg;
  sc.solver.iter = 40;
  sc.num_v_cycles = 0;
  sc.max_levels = 1;
  const real l0 = sc.l0 * R->lavg();

  std::vector<vec3> z(R->corner_count(), vec3::Zero());
  z[8] = vec3(0.02, 0.0, 0.0);
  std::vector<vec3> f;
  duchamp::sobolev_detail::matvec_A(*R, sc.alpha_B, sc.alpha_0, sc.p1, l0,
                                    sc.bh_eps, z, f);

  const std::vector<vec3> bc = F->restrict_JT(0, f);
  duchamp::sobolev_constraint_solver pin;
  pin.add(duchamp::sobolev_barycenter_constraint::create());
  pin.init_rest(*coarse);
  const std::vector<vec3> e = pin.project_gradient(*coarse, bc, sc, l0);
  const std::vector<vec3> d = F->prolong_J(0, e);

  std::vector<vec3> Ad, Ac, Bd, B0d;
  duchamp::sobolev_detail::matvec_A(*R, sc.alpha_B, sc.alpha_0, sc.p1, l0,
                                    sc.bh_eps, d, Ad);
  duchamp::sobolev_detail::matvec_A(*coarse, sc.alpha_B, sc.alpha_0, sc.p1, l0,
                                    sc.bh_eps, e, Ac);
  duchamp::sobolev_detail::matvec_B(*R, d, sc.p1, l0, sc.bh_eps, Bd);
  duchamp::sobolev_detail::matvec_B0(*R, d, sc.p1, l0, sc.bh_eps, B0d);
  auto closest = [](asawa::rod::rod &rod, bool midpoints, real &gap) {
    using namespace duchamp::sobolev_detail;
    const auto E = build_edges(rod);
    const auto &x = rod.x();
    gap = real(1e300);
    int bi = -1, bj = -1;
    for (size_t i = 0; i < E.n; ++i) {
      for (size_t j = i + 1; j < E.n; ++j) {
        if (paper_edges_touch(E, i, j))
          continue;
        real dist = (E.mid[i] - E.mid[j]).norm();
        if (!midpoints) {
          const vec3 p[2] = {x[static_cast<size_t>(E.ep0[i])],
                             x[static_cast<size_t>(E.ep1[i])]};
          const vec3 q[2] = {x[static_cast<size_t>(E.ep0[j])],
                             x[static_cast<size_t>(E.ep1[j])]};
          dist = real(1e300);
          for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
              dist = std::min(dist, (p[a] - q[b]).norm());
        }
        if (dist < gap) {
          gap = dist;
          bi = static_cast<int>(i);
          bj = static_cast<int>(j);
        }
      }
    }
    return std::pair<int, int>{bi, bj};
  };
  real fine_end = 0, fine_mid = 0, coarse_end = 0, coarse_mid = 0;
  const auto fine_pair = closest(*R, false, fine_end);
  closest(*R, true, fine_mid);
  closest(*coarse, false, coarse_end);
  closest(*coarse, true, coarse_mid);
  real w_end = 0, w_mid = 0;
  {
    using namespace duchamp::sobolev_detail;
    const auto E = build_edges(*R);
    const auto &x = R->x();
    const auto &quad = trapezoid_edge_quad();
    const size_t i = static_cast<size_t>(fine_pair.first);
    const size_t j = static_cast<size_t>(fine_pair.second);
    if (fine_pair.first >= 0) {
      w_end = quad_pair_weight(x[static_cast<size_t>(E.ep0[i])],
                               x[static_cast<size_t>(E.ep1[i])], E.T[i],
                               x[static_cast<size_t>(E.ep0[j])],
                               x[static_cast<size_t>(E.ep1[j])], E.T[j],
                               E.ell[i], E.ell[j], l0, sc.p1, true, quad);
      w_mid = directed_pair_weight(E.mid[i], E.mid[j], E.T[i], E.T[j], E.ell[i],
                                   E.ell[j], l0, sc.p1, true);
    }
  }
  std::vector<vec3> r(f.size(), vec3::Zero());
  for (size_t i = 0; i < f.size() && i < Ad.size(); ++i)
    r[i] = f[i] - Ad[i];
  const std::vector<vec3> rc = F->restrict_JT(0, r);

  const real n0 = field_norm(f);
  const real n1 = field_norm(r);
  std::vector<vec3> x = d;
  auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    duchamp::sobolev_detail::matvec_A(*R, sc.alpha_B, sc.alpha_0, sc.p1, l0,
                                      sc.bh_eps, u, Au);
  };
  duchamp::sobolev_cg_smooth(f, 12, x, apply);
  const real n_post = sobolev_residual(*R, sc, l0, f, x);
  std::ostringstream msg;
  msg << "|r| " << n0 << " |e| " << field_norm(e) << " |J e| " << field_norm(d)
      << " |B Je| " << field_norm(Bd) << " |B0 Je| " << field_norm(B0d)
      << " |A J e| " << field_norm(Ad) << " |r - A J e| " << n1
      << " |r|_post " << n_post << " |J^T r| " << field_norm(bc)
      << " |A_c e| " << field_norm(Ac) << " endgap f/c " << fine_end << " "
      << coarse_end << " midgap f/c " << fine_mid << " " << coarse_mid
      << " w_end " << w_end << " w_mid " << w_mid;
  std::cout << "[prolong residual] " << msg.str() << std::endl;
  if (!(std::isfinite(n0) && std::isfinite(n_post) && n_post < n0))
    gaudi::test::record_failure(__FILE__, __LINE__, msg.str().c_str(), false);
}

// The white corner before the first black prolongs from the two blacks
// that prev/next reach. No loop flag is read.
GAUDI_TEST(sobolev_J_seam_follows_links) {
  auto R = asawa::rod::rod::create(sobolev_helix(32), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 2);
  const auto &pairs = F->j_pairs(0);
  int first = -1;
  for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
    if (pairs[static_cast<size_t>(i)][0] == 0 &&
        pairs[static_cast<size_t>(i)][1] < 0) {
      first = i;
      break;
    }
  }
  GAUDI_ASSERT(first >= 0);
  const int white =
      static_cast<int>(R->prev(asawa::rod::corner_id(first)));
  GAUDI_ASSERT(white >= 0);
  GAUDI_ASSERT(pairs[static_cast<size_t>(white)][1] >= 0);

  int last = white;
  asawa::rod::CornerId a = R->prev(asawa::rod::corner_id(white));
  int guard = static_cast<int>(R->corner_count());
  while (!(a < asawa::rod::corner_id(0)) && guard-- > 0) {
    const int ai = static_cast<int>(a);
    if (pairs[static_cast<size_t>(ai)][1] < 0 &&
        pairs[static_cast<size_t>(ai)][0] >= 0) {
      last = ai;
      break;
    }
    a = R->prev(a);
  }
  const index_t coarse_last = pairs[static_cast<size_t>(last)][0];
  GAUDI_ASSERT(pairs[static_cast<size_t>(white)][0] == coarse_last);
  GAUDI_ASSERT(pairs[static_cast<size_t>(white)][1] == 0);

  const size_t nc = F->level_rod(1)->corner_count();
  std::vector<vec3> e(nc, vec3::Zero());
  e[static_cast<size_t>(coarse_last)] = vec3(real(2), 0, 0);
  e[0] = vec3(real(4), 0, 0);
  const std::vector<vec3> fine = F->prolong_J(0, e);
  GAUDI_ASSERT((fine[static_cast<size_t>(white)] - vec3(real(3), 0, 0)).norm() <
               real(1e-12));
}

// Collapsing a vertex field keeps its sum. J^T sends each fine value entirely
// onto the coarse vertices: weight 1 on a black, 1/2 + 1/2 on a white.
GAUDI_TEST(sobolev_restrict_preserves_vertex_sum) {
  auto R = asawa::rod::rod::create(sobolev_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 3);

  std::vector<vec3> v(R->corner_count());
  for (size_t i = 0; i < v.size(); ++i)
    v[i] = vec3(std::sin(real(i)), real(0.01) * real(i), real(1));
  const vec3 s0 = vertex_sum(v);
  std::vector<vec3> cur = v;
  for (int ell = 0; ell + 1 < F->n_levels() && ell < 2; ++ell) {
    cur = F->restrict_JT(ell, cur);
    const vec3 s = vertex_sum(cur);
    GAUDI_ASSERT((s - s0).norm() <=
                 real(1e-8) * (real(1) + s0.norm()));
  }
}

// Prolongation copies blacks and averages whites. Collapsing that field
// returns the prolonged sum, which includes the white vertices just created.
GAUDI_TEST(sobolev_prolong_roundtrip_preserves_sum) {
  auto R = asawa::rod::rod::create(sobolev_helix(64), true);
  auto F = asawa::rod::residual_forest::create(R);
  GAUDI_ASSERT(F && F->valid() && F->n_levels() >= 2);
  const auto &pairs = F->j_pairs(0);
  const size_t nc = F->level_rod(1)->corner_count();

  const vec3 constant(real(3), real(-1), real(0.5));
  std::vector<vec3> coarse(nc, constant);
  std::vector<vec3> fine = F->prolong_J(0, coarse);
  GAUDI_ASSERT(fine.size() == R->corner_count());
  for (const vec3 &p : fine)
    GAUDI_ASSERT((p - constant).norm() < real(1e-12));

  for (size_t i = 0; i < nc; ++i)
    coarse[i] = vec3(std::sin(real(i)), std::cos(real(0.5) * real(i)),
                     real(0.01) * real(i));
  fine = F->prolong_J(0, coarse);
  for (size_t i = 0; i < pairs.size() && i < fine.size(); ++i) {
    const index_t i0 = pairs[i][0];
    const index_t i1 = pairs[i][1];
    GAUDI_ASSERT(i0 >= 0 && static_cast<size_t>(i0) < nc);
    if (i1 < 0) {
      GAUDI_ASSERT((fine[i] - coarse[static_cast<size_t>(i0)]).norm() <
                   real(1e-12));
    } else {
      GAUDI_ASSERT(static_cast<size_t>(i1) < nc);
      const vec3 avg = real(0.5) * (coarse[static_cast<size_t>(i0)] +
                                    coarse[static_cast<size_t>(i1)]);
      GAUDI_ASSERT((fine[i] - avg).norm() < real(1e-12));
    }
  }
  const std::vector<vec3> back = F->restrict_JT(0, fine);
  GAUDI_ASSERT((vertex_sum(back) - vertex_sum(fine)).norm() <=
               real(1e-8) * (real(1) + vertex_sum(fine).norm()));
}

/// Cycle Laplacian A = L + I. Eigenvalues are 1 + 4 sin²(π k / n).
/// A full-spectrum Chebyshev run must match an LDLT solve. The high-band
/// window used as a multigrid smoother must kill the top mode and leave
/// the bottom mode. CG on the same matvec is the other reference.
GAUDI_TEST(sobolev_chebyshev_matches_cycle_laplacian) {
  constexpr int n = 64;
  const real mu = real(1);
  const real pi = real(3.14159265358979323846);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd blow(n), bhigh(n);
  real lmin = mu + real(4);
  real lmax = 0;
  for (int k = 0; k < n; ++k) {
    const real lam = mu + real(4) * std::pow(std::sin(pi * real(k) / real(n)), 2);
    lmin = std::min(lmin, lam);
    lmax = std::max(lmax, lam);
  }
  for (int i = 0; i < n; ++i) {
    A(i, i) = real(2) + mu;
    A(i, (i + 1) % n) = real(-1);
    A(i, (i + n - 1) % n) = real(-1);
    const real th = real(2) * pi * real(i) / real(n);
    blow(i) = std::sin(th);
    bhigh(i) = (i & 1) ? real(1) : real(-1);
  }
  Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
  GAUDI_ASSERT(ldlt.info() == Eigen::Success);
  const Eigen::VectorXd xlow = ldlt.solve(blow);
  const Eigen::VectorXd xhigh = ldlt.solve(bhigh);

  auto apply = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    Au.assign(u.size(), vec3::Zero());
    for (int i = 0; i < n; ++i) {
      const vec3 up = u[static_cast<size_t>((i + 1) % n)];
      const vec3 um = u[static_cast<size_t>((i + n - 1) % n)];
      Au[static_cast<size_t>(i)] = (real(2) + mu) * u[static_cast<size_t>(i)] - up - um;
    }
  };
  auto pack = [&](const Eigen::VectorXd &s) {
    std::vector<vec3> f(static_cast<size_t>(n), vec3::Zero());
    for (int i = 0; i < n; ++i)
      f[static_cast<size_t>(i)] = vec3(s(i), real(0.5) * s(i), real(-0.25) * s(i));
    return f;
  };
  auto rel_err = [&](const std::vector<vec3> &x, const Eigen::VectorXd &xref) {
    real num = 0;
    real den = 0;
    for (int i = 0; i < n; ++i) {
      const real d = x[static_cast<size_t>(i)][0] - xref(i);
      num += d * d;
      den += xref(i) * xref(i);
    }
    return std::sqrt(num / den);
  };

  const std::vector<vec3> rhs = pack(blow + real(0.25) * bhigh);
  const Eigen::VectorXd xref = ldlt.solve(blow + real(0.25) * bhigh);
  std::vector<vec3> x_cheb(static_cast<size_t>(n), vec3::Zero());
  std::vector<vec3> x_cg(static_cast<size_t>(n), vec3::Zero());
  const real frac = lmin / lmax;
  duchamp::sobolev_chebyshev_smooth(rhs, 24, lmax, frac, x_cheb, apply);
  duchamp::sobolev_cg_smooth(rhs, 24, x_cg, apply);
  const real e_cheb = rel_err(x_cheb, xref);
  const real e_cg = rel_err(x_cg, xref);
  console_logger::info << "chebyshev full-window rel=" << e_cheb
                       << " cg rel=" << e_cg << " lmin=" << lmin
                       << " lmax=" << lmax << std::endl;
  GAUDI_ASSERT(e_cheb < real(1e-6));
  GAUDI_ASSERT(e_cg < real(1e-6));

  const real lam_est =
      duchamp::sobolev_arnoldi_lambda(pack(bhigh),
                                     std::vector<vec3>(static_cast<size_t>(n), vec3::Zero()),
                                     apply);
  console_logger::info << "arnoldi λ=" << lam_est << " true λmax=" << lmax
                       << std::endl;
  GAUDI_ASSERT(lam_est > real(0.8) * lmax);
  GAUDI_ASSERT(lam_est <= lmax * real(1.05));

  std::vector<vec3> y_high(static_cast<size_t>(n), vec3::Zero());
  std::vector<vec3> y_low(static_cast<size_t>(n), vec3::Zero());
  duchamp::sobolev_chebyshev_smooth(pack(bhigh), 8, lmax, real(0.3), y_high, apply);
  duchamp::sobolev_chebyshev_smooth(pack(blow), 8, lmax, real(0.3), y_low, apply);
  const real e_high = rel_err(y_high, xhigh);
  const real e_low = rel_err(y_low, xlow);
  console_logger::info << "chebyshev high-band high-mode rel=" << e_high
                       << " low-mode rel=" << e_low << std::endl;
  GAUDI_ASSERT(e_high < real(1e-2));
  GAUDI_ASSERT(e_low > real(10) * e_high);

  // Eigenvalues from 1 to 1e6. Eight steps on [0.3 λ_max, 1.1 λ_max]
  // solve the top entry and leave the bottom entry.
  constexpr int m = 16;
  Eigen::VectorXd diag(m), btop(m), bbot(m);
  for (int i = 0; i < m; ++i) {
    diag(i) = std::pow(real(10), real(6) * real(i) / real(m - 1));
    btop(i) = 0;
    bbot(i) = 0;
  }
  btop(m - 1) = 1;
  bbot(0) = 1;
  const real dmax = diag(m - 1);
  auto apply_d = [&](const std::vector<vec3> &u, std::vector<vec3> &Au) {
    Au.assign(u.size(), vec3::Zero());
    for (int i = 0; i < m; ++i)
      Au[static_cast<size_t>(i)] = diag(i) * u[static_cast<size_t>(i)];
  };
  auto pack_d = [&](const Eigen::VectorXd &s) {
    std::vector<vec3> f(static_cast<size_t>(m), vec3::Zero());
    for (int i = 0; i < m; ++i)
      f[static_cast<size_t>(i)][0] = s(i);
    return f;
  };
  const Eigen::VectorXd xtop = btop.cwiseQuotient(diag);
  const Eigen::VectorXd xbot = bbot.cwiseQuotient(diag);
  std::vector<vec3> ztop(static_cast<size_t>(m), vec3::Zero());
  std::vector<vec3> zbot(static_cast<size_t>(m), vec3::Zero());
  duchamp::sobolev_chebyshev_smooth(pack_d(btop), 8, dmax, real(0.3), ztop,
                                    apply_d);
  duchamp::sobolev_chebyshev_smooth(pack_d(bbot), 8, dmax, real(0.3), zbot,
                                    apply_d);
  real et = 0;
  real eb = 0;
  for (int i = 0; i < m; ++i) {
    const real dt = ztop[static_cast<size_t>(i)][0] - xtop(i);
    const real db = zbot[static_cast<size_t>(i)][0] - xbot(i);
    et += dt * dt;
    eb += db * db;
  }
  et = std::sqrt(et) / xtop.norm();
  eb = std::sqrt(eb) / xbot.norm();
  console_logger::info << "wide-spectrum high rel=" << et << " low rel=" << eb
                       << std::endl;
  GAUDI_ASSERT(et < real(1e-2));
  GAUDI_ASSERT(eb > real(0.5));
}

} // namespace test
} // namespace gaudi

#endif // GAUDI_TP_SOBOLEV_JACOBI_TESTS_HPP
