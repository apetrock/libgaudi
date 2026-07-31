#ifndef __GAUDI_MLS_JET_BOOTSTRAP_TESTS_HPP__
#define __GAUDI_MLS_JET_BOOTSTRAP_TESTS_HPP__

#include "gaudi/albers/osculating_torus.hpp"
#include "gaudi/albers/quadric.hpp"
#include "gaudi/albers/shape_operator_constrained_fit.hpp"
#include "gaudi/asawa/shell/shape_operator.hpp"
#include "gaudi/calder/mls_jet_bootstrap.hpp"
#include "gaudi/calder/shape_operator_weights.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"
#include "gaudi/test/test.hpp"

#include "shape_operator_fit_generated.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(albers_mk_quad_W_rows_matches_vech_PHP) {
  const vec3 n = vec3(0.2, -0.5, 0.8).normalized();
  albers::vec10 Q = albers::vec10::Zero();
  Q[0] = 1.1;
  Q[1] = 0.7;
  Q[2] = 0.4;
  Q[3] = -0.3;
  Q[4] = 0.2;
  Q[5] = -0.1;
  Q[6] = 0.05;
  Q[7] = -0.02;
  Q[8] = 0.01;
  Q[9] = -1.0;

  const mat3 H =
      (mat3() << 2.0 * Q[0], Q[3], Q[4], Q[3], 2.0 * Q[1], Q[5], Q[4],
       Q[5], 2.0 * Q[2])
          .finished();
  const mat3 P = mat3::Identity() - n * n.transpose();
  const mat3 Php = P * H * P;
  const vec6 expected = albers::vech3(Php);

  const Eigen::Matrix<real, 6, 10> M =
      albers::medial_generated::mk_quad_W_rows(n);
  const vec6 got = M * Q;
  const real err = (got - expected).norm();
  GAUDI_EXPECT(err < 1e-9);
}

GAUDI_TEST(albers_osculating_torus_from_sphere_like_W) {
  // Sphere of radius r=2, outward n → κ=1/2, κ=1/2.
  const real r = 2.0;
  const vec3 foot = vec3::Zero();
  const vec3 n = vec3::UnitZ();
  const mat3 W = (1.0 / r) * (mat3::Identity() - n * n.transpose());
  const albers::osculating_torus T =
      albers::osculating_torus_from_shape_operator(foot, n, W);
  GAUDI_EXPECT(T.valid);
  GAUDI_EXPECT(std::abs(T.r - r) < 1e-6);
  // Degenerate major ≈ r for equal curvatures with s_prod=+1 → R≈0;
  // still produces a finite tube prior.
  GAUDI_EXPECT(std::isfinite(T.R));
  GAUDI_EXPECT(std::abs(albers::torus_sdf(foot, T)) < 1e-6);
  // Sphere limit: geometric center opposite outward n → foot - r*n.
  GAUDI_EXPECT((T.center + T.r * n).norm() < 1e-5);
  GAUDI_EXPECT(T.sign == 1);
}

GAUDI_TEST(asawa_face_shape_operator_finite_on_torus) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  TorusMesh mesh = make_offset_torus_shell(
      24, 12, major_radius, minor_radius,
      make_torus_frame(vec3::Zero(), vec3::UnitZ()));
  asawa::shell::shell &M = *mesh.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<mat3> Wf = asawa::shell::face_shape_operators(M, x);
  int ok = 0;
  for (const mat3 &W : Wf) {
    if (W.allFinite() && W.norm() > 1e-8) {
      ++ok;
    }
  }
  GAUDI_EXPECT(ok > static_cast<int>(Wf.size()) / 2);
}

GAUDI_TEST(calder_mls_jet_bootstrap_torus_smoke) {
  const real major_radius = 1.25;
  const real minor_radius = 0.35;
  TorusMesh mesh = make_offset_torus_shell(
      32, 16, major_radius, minor_radius,
      make_torus_frame(vec3::Zero(), vec3::UnitZ()));
  asawa::shell::shell &M = *mesh.shell;
  const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(M, x);
  const real avg_len = asawa::shell::avg_length(M, x);
  const real l0 = std::max(real(2.0) * avg_len, real(1e-12));

  const int vi = static_cast<int>(x.size() / 2);
  const std::vector<vec3> p_fit = {x[static_cast<size_t>(vi)]};
  const std::vector<vec3> n_fit = {N[static_cast<size_t>(vi)]};

  calder::mls_jet_bootstrap_params params;
  params.ablation =
      calder::mls_bootstrap_ablation::soft_darboux_torus_green_dp;
  params.torus_sign = albers::torus_sign_mode::pat_eq22;
  params.fit_p = 3.0;
  params.fit_w0 = 1.0;
  params.radius_scale = 1.0;

  const std::vector<albers::vec14> Q =
      calder::darboux_fit_bootstrapped(M, p_fit, n_fit, l0, params);
  GAUDI_ASSERT(Q.size() == 1);
  GAUDI_EXPECT(Q[0].allFinite());

  const vec3 g = albers::darboux_grad(Q[0], vec3::Zero());
  GAUDI_EXPECT(g.allFinite());
  GAUDI_EXPECT(g.norm() > 1e-12);
  const real align = std::abs(g.normalized().dot(n_fit[0].normalized()));
  GAUDI_EXPECT(align > 0.35);

  std::cerr << "\n[mls_jet_bootstrap softDarbouxTorusGreen] |g|=" << g.norm()
            << " align=" << align << "\n";
}

GAUDI_TEST(calder_aniso_mahalanobis_weight_positive) {
  albers::principal_curvature_frame pc;
  pc.n = vec3::UnitZ();
  pc.e_min = vec3::UnitX();
  pc.e_max = vec3::UnitY();
  pc.k_min = 0.5;
  pc.k_max = 2.0;
  pc.valid = true;
  const real w =
      calder::aniso_inv_dist_weight(vec3(0.1, 0.0, 0.0), 0.05, 3.0, pc);
  GAUDI_EXPECT(std::isfinite(w));
  GAUDI_EXPECT(w > 0.0);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_MLS_JET_BOOTSTRAP_TESTS_HPP__
