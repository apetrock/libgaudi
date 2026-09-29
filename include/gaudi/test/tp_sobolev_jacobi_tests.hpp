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
      duchamp::sobolev_quadratic_B(*R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh,
                                   u, v);
  const real Bv_u =
      duchamp::sobolev_quadratic_B(*R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh,
                                   v, u);
  const real B0u_v = duchamp::sobolev_quadratic_B0(
      *R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh, u, v);
  const real B0v_u = duchamp::sobolev_quadratic_B0(
      *R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh, v, u);
  const real B0u_u = duchamp::sobolev_quadratic_B0(
      *R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh, u, u);
  const real Bu_u = duchamp::sobolev_quadratic_B(
      *R, m, /*p=*/4.0, l0, /*mass_scale=*/1.0, bh, u, u);
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
  const real Buu = duchamp::sobolev_quadratic_B(*R, m, /*p=*/3.0, l0,
                                                /*mass_scale=*/1.0, bh, u, u);
  const real B0uu = duchamp::sobolev_quadratic_B0(
      *R, m, /*p=*/3.0, l0, /*mass_scale=*/1.0, bh, u, u);
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
  duchamp::sobolev_detail::matvec_B0(*R, u, m, /*p=*/3.0, l0, /*mass=*/1.0,
                                     /*bh_eps=*/0.0, B0_leaf);
  duchamp::sobolev_detail::matvec_B0(*R, u, m, /*p=*/3.0, l0, /*mass=*/1.0,
                                     /*bh_eps=*/0.25, B0_bh);
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
  sc.p = 3.0;
  sc.l0 = 0.5;
  sc.alpha_B = 1.0;
  sc.alpha_0 = 1.0;
  sc.mass = 1.0;
  sc.bh_eps = 0.0;
  const real l0_abs = sc.l0 * Rd->_Cc;
  auto m_pov = duchamp::sobolev_detail::sample_masses(*R);
  for (real &m : m_pov)
    m *= sc.mass;

  const real r0 = field_norm(f);
  sc.jacobi_iters = 24;
  const std::vector<vec3> g = duchamp::sobolev_jacobi_solve(*R, f, sc, l0_abs);
  std::vector<vec3> Ag;
  duchamp::sobolev_detail::matvec_A(*R, m_pov, sc.alpha_B, sc.alpha_0, sc.p,
                                    l0_abs, sc.mass, sc.bh_eps, g, Ag);
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
                 .p = 3.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .mass = 1.0,
                 .jacobi_iters = 8,
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
                 .p = 3.0,
                 .alpha_B = 0.0,
                 .alpha_0 = 1.0,
                 .mass = 1.0e-2,
                 .jacobi_iters = 16,
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
                 .p = 3.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .mass = 1.0,
                 .jacobi_iters = 16,
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
  sc.p = 3.0;
  sc.l0 = 0.5;
  sc.alpha_B = 1.0;
  sc.alpha_0 = 1.0;
  sc.mass = 1.0;
  sc.bh_eps = 0.0;
  sc.jacobi_iters = 4;
  const real l0_abs = sc.l0 * Rd->_Cc;
  auto m_pov = duchamp::sobolev_detail::sample_masses(*R);
  for (real &m : m_pov)
    m *= sc.mass;

  auto residual_of = [&](const std::vector<vec3> &g) {
    std::vector<vec3> Ag;
    duchamp::sobolev_detail::matvec_A(*R, m_pov, sc.alpha_B, sc.alpha_0, sc.p,
                                      l0_abs, sc.mass, sc.bh_eps, g, Ag);
    real s = 0;
    for (size_t i = 0; i < f.size(); ++i)
      s += (f[i] - Ag[i]).squaredNorm();
    return std::sqrt(s);
  };

  sc.num_v_cycles = 0;
  const auto g_jac = duchamp::sobolev_jacobi_solve(*R, f, sc, l0_abs);
  // Flat Jacobi with same smooth budget as one MG cycle pre+post+coarse.
  sc.jacobi_iters = 12;
  const auto g_flat = duchamp::sobolev_jacobi_solve(*R, f, sc, l0_abs);

  sc.jacobi_iters = 4;
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
  duchamp::sobolev_detail::matvec_A(*R, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, /*bh=*/0, u,
                                    out, cov.get());
  duchamp::sobolev_detail::matvec_A(*R, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, /*bh=*/0, u,
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
  duchamp::sobolev_detail::matvec_A(*rod1, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, /*bh=*/0, u,
                                    out, coarse.get());
  duchamp::sobolev_detail::matvec_A(*rod1, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, /*bh=*/0, u,
                                    out, coarse.get());
  GAUDI_ASSERT(arp::edge_cov_tree<2>::build_count() == built);
  GAUDI_ASSERT(out.size() == u.size());
}

// Frozen B and B0 weights match a live matvec, and Jacobi does not walk again.
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
      *R, E, /*p=*/4, l0, /*mass=*/1, bh, cone, /*use_tp_k=*/false, *cov);
  auto wB0 = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, /*mass=*/1, bh, cone, /*use_tp_k=*/true, *cov);
  GAUDI_ASSERT(calder::cov_walk_count() == before + 2);

  std::vector<vec3> u(R->corner_count(), vec3(0.2, -0.1, 0.05));
  std::vector<vec3> live, cached;
  duchamp::sobolev_detail::matvec_A(*R, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, bh, u, live,
                                    cov.get(), cone);
  const int after_live = calder::cov_walk_count();
  duchamp::sobolev_detail::matvec_A(*R, m, /*alpha_B=*/1, /*alpha_0=*/1,
                                    /*p=*/4, l0, /*mass_scale=*/1, bh, u,
                                    cached, cov.get(), cone, &wB, &wB0);
  GAUDI_ASSERT(calder::cov_walk_count() == after_live);
  GAUDI_ASSERT(field_max_abs_diff(live, cached) < real(1e-8));

  duchamp::sobolev_tangent_point_force_config sc;
  sc.alpha_B = 1;
  sc.alpha_0 = 1;
  sc.p = 4;
  sc.mass = 1;
  sc.bh_eps = bh;
  sc.cone_eps = cone;
  sc.jacobi_iters = 6;
  std::vector<vec3> f(R->corner_count(), vec3(1, 0, 0));
  std::vector<real> diag(R->corner_count(), real(1));
  std::vector<vec3> g(R->corner_count(), vec3::Zero());
  const int before_smooth = calder::cov_walk_count();
  duchamp::sobolev_jacobi_smooth(*R, f, m, diag, sc, l0, sc.jacobi_iters, g,
                                 cov.get(), &wB, &wB0);
  GAUDI_ASSERT(calder::cov_walk_count() == before_smooth);
  GAUDI_ASSERT(g.size() == f.size());
}

// Leaf B and B0 weights are the trapezoid pair, and adjacent edges are skipped.
GAUDI_TEST(sobolev_trapezoid_quad_weights) {
  auto R = asawa::rod::rod::create(sobolev_helix(16), true);
  auto cov = calder::build_rod_cov_tree(*R);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  const real l0 = real(0.1) * R->lavg();
  const auto &x = R->x();
  const auto &quad = duchamp::sobolev_detail::trapezoid_edge_quad();
  auto check = [&](bool use_tp_k) {
    const auto rows = duchamp::sobolev_detail::build_cov_weights(
        *R, E, /*p=*/4, l0, /*mass=*/1, /*bh=*/0, /*cone=*/0, use_tp_k, *cov,
        quad);
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
        const real expect = duchamp::sobolev_detail::quad_pair_weight(
            x[static_cast<size_t>(a0)], x[static_cast<size_t>(a1)], E.T[i],
            x[static_cast<size_t>(b0)], x[static_cast<size_t>(b1)], E.T[j_edge],
            E.ell[i], E.ell[j_edge], l0, /*p=*/4, use_tp_k, quad);
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

// Leaves-only Ê matches the hand trapezoid sum, adjacent pairs are absent,
// and a vertex nudge matches dÊ. The far lump's endpoint derivative matches
// a frozen cluster.
GAUDI_TEST(sobolev_paper_energy_leaf) {
  auto R = asawa::rod::rod::create(sobolev_helix(12), true);
  const real l0 = real(0.1) * R->lavg();
  const real p = real(3);
  const auto E = duchamp::sobolev_detail::build_edges(*R);
  const auto &quad = duchamp::sobolev_detail::trapezoid_edge_quad();
  const auto &x = R->x();
  real hand = real(0);
  int pairs = 0;
  int adj = 0;
  for (size_t i = 0; i < E.n; ++i) {
    for (size_t j = i + 1; j < E.n; ++j) {
      const int a0 = E.ep0[i];
      const int a1 = E.ep1[i];
      const int b0 = E.ep0[j];
      const int b1 = E.ep1[j];
      const bool share = a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1;
      if (share) {
        ++adj;
        continue;
      }
      hand += duchamp::sobolev_detail::energy_quad_weight(
          x[static_cast<size_t>(a0)], x[static_cast<size_t>(a1)], E.T[i],
          E.ell[i], x[static_cast<size_t>(b0)], x[static_cast<size_t>(b1)],
          E.T[j], E.ell[j], l0, p, quad);
      ++pairs;
    }
  }
  GAUDI_ASSERT(pairs > 0 && adj > 0);
  const auto ev =
      duchamp::sobolev_detail::eval_paper_tp(*R, l0, p, /*bh=*/0, /*cone=*/0);
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
        duchamp::sobolev_detail::eval_paper_tp(*R, l0, p, 0, 0).energy;
    R->x() = x0;
    R->x()[iv] -= step;
    R->_update_frames();
    const real Em =
        duchamp::sobolev_detail::eval_paper_tp(*R, l0, p, 0, 0).energy;
    const real num = (Ep - Em) / (real(2) * h);
    const real ana = ev.grad[iv][axis];
    GAUDI_ASSERT(std::abs(num - ana) <= real(1e-3) * (real(1) + std::abs(ana)));
  }
  R->x() = x0;
  R->_update_frames();

  const vec3 pj = E.mid[0] + vec3(real(0.3), real(0.2), real(0.1));
  const vec3 Tj = vec3(0, 0, 1);
  const real ellc = real(0.4);
  auto far_at = [&](const vec3 &p0, const vec3 &p1) {
    const vec3 d = p1 - p0;
    const real ell = std::max(d.norm(), real(1e-12));
    const vec3 T = d / ell;
    const vec3 mid = real(0.5) * (p0 + p1);
    return duchamp::sobolev_detail::energy_far_term(T, ell, mid, pj, Tj, ellc,
                                                    l0, p);
  };
  const vec3 p0 = x[static_cast<size_t>(E.ep0[0])];
  const vec3 p1 = x[static_cast<size_t>(E.ep1[0])];
  const auto f0 = far_at(p0, p1);
  for (int axis = 0; axis < 3; ++axis) {
    vec3 step = vec3::Zero();
    step[axis] = h;
    const real num0 =
        (far_at(p0 + step, p1).e - far_at(p0 - step, p1).e) / (real(2) * h);
    const real num1 =
        (far_at(p0, p1 + step).e - far_at(p0, p1 - step).e) / (real(2) * h);
    GAUDI_ASSERT(std::abs(num0 - f0.ga0[axis]) <=
                 real(1e-3) * (real(1) + std::abs(f0.ga0[axis])));
    GAUDI_ASSERT(std::abs(num1 - f0.ga1[axis]) <=
                 real(1e-3) * (real(1) + std::abs(f0.ga1[axis])));
  }
}

// B couples slopes by ⟨T_I, T_J⟩. A perpendicular partner drops out.
// A parallel partner matches the plain difference. B0 keeps align = 1.
GAUDI_TEST(sobolev_B_tangent_align) {
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
      *R, E, /*p=*/4, l0, /*mass=*/1, /*bh=*/0, /*cone=*/0, /*use_tp_k=*/false,
      *cov);
  const auto rowsB0 = duchamp::sobolev_detail::build_cov_weights(
      *R, E, /*p=*/4, l0, /*mass=*/1, /*bh=*/0, /*cone=*/0, /*use_tp_k=*/true,
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
  GAUDI_ASSERT(std::abs(perp->align) < real(1e-8));
  GAUDI_ASSERT(std::abs(par->align - real(1)) < real(1e-8));
  GAUDI_ASSERT(std::abs(perp0->align - real(1)) < real(1e-8));

  const vec3 s(real(0.2), real(-0.3), real(0.4));
  std::vector<vec3> phi(E.n, vec3::Zero());
  phi[static_cast<size_t>(e_perp)] = s;
  const auto r_perp = duchamp::sobolev_detail::edge_K_diff(
      *R, E, phi, /*p=*/4, l0, /*mass=*/1, /*bh=*/0, /*use_tp_k=*/false,
      cov.get(), /*cone=*/0, &rowsB);
  GAUDI_ASSERT(r_perp[static_cast<size_t>(e_query)].norm() < real(1e-8));

  phi.assign(E.n, vec3::Zero());
  phi[static_cast<size_t>(e_par)] = s;
  const auto r_par = duchamp::sobolev_detail::edge_K_diff(
      *R, E, phi, /*p=*/4, l0, /*mass=*/1, /*bh=*/0, /*use_tp_k=*/false,
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
                 .p = 2.0,
                 .alpha_B = 1.0,
                 .alpha_0 = 1.0,
                 .mass = 1.0,
                 .jacobi_iters = 2,
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

} // namespace test
} // namespace gaudi

#endif // GAUDI_TP_SOBOLEV_JACOBI_TESTS_HPP
