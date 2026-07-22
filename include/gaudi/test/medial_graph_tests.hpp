#ifndef __GAUDI_TEST_MEDIAL_GRAPH_TESTS_HPP__
#define __GAUDI_TEST_MEDIAL_GRAPH_TESTS_HPP__

#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/albers/darboux_medial_geometry.hpp"
#include "gaudi/albers/medial_shape_search.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/calder_graph_nodes.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/duchamp/medial_graph_nodes.hpp"
#include "gaudi/duchamp/medial_search_helpers.hpp"
#include "gaudi/test/test.hpp"
#include "medial_kernels.h"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

GAUDI_TEST(medial_kernels_smoke_include) {
  const vec3 g(0.0, 0.0, 1.0);
  const mat3 H = mat3::Identity();
  mat3 W;
  albers::medial_generated::shape_operator_from_GH(g, H, W);
  GAUDI_ASSERT(W.allFinite());
}

GAUDI_TEST(medial_shape_search_sphere_smoke) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real avg_len = asawa::shell::avg_length(*M, x);
  const real l0 = 0.2 * avg_len;
  const std::vector<albers::vec14> fits =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);
  GAUDI_ASSERT(!fits.empty());

  int accepted = 0;
  for (size_t i = 0; i < std::min<size_t>(fits.size(), 8); ++i) {
    const vec3 dir = albers::aligned_inward_ray_dir(fits[i], N[i]);
    if (dir.squaredNorm() < 1e-24) {
      continue;
    }
    const auto result = albers::search_medial_along_ray(
        fits[i], x[i], dir, 4.0 * avg_len);
    if (result.converged && result.travel > 1e-6) {
      ++accepted;
    }
  }
  GAUDI_ASSERT(accepted > 0);
}

GAUDI_TEST(medial_graph_cyclide_fit_node) {
  auto M = asawa::shell::load_sphere(1.0, 10, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real l0 = 0.15 * asawa::shell::avg_length(*M, x);
  const std::vector<albers::vec14> direct =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  liblombardi::GraphContext ctx;
  auto body = ctx.create_node<duchamp::body_constant_node>(
      duchamp::make_shell_body(M));
  auto positions = ctx.create_node<duchamp::position_snapshot_node>();
  auto normals = ctx.create_node<duchamp::vertex_normals_snapshot_node>();
  auto fit = ctx.create_node<duchamp::darboux_cyclide_fit_node>(l0, 3.0);
  ctx.link(body->output(), positions->body());
  ctx.link(body->output(), normals->body());
  ctx.link(body->output(), fit->body());
  ctx.link(positions->output(), fit->pov());
  ctx.link(normals->output(), fit->n_pov());
  ctx.run();

  const auto &got =
      fit->get_datum<duchamp::darboux_cyclide_fit_node::OutputPortDef>()->data();
  GAUDI_ASSERT(got.size() == direct.size());
  real max_err = 0.0;
  for (size_t i = 0; i < direct.size(); ++i) {
    max_err = std::max(max_err, (got[i] - direct[i]).norm());
  }
  GAUDI_ASSERT(max_err < 1e-10);
}

GAUDI_TEST(medial_graph_end_to_end_sphere) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const real avg_len = asawa::shell::avg_length(*M, x);
  const real l0 = 0.2 * avg_len;
  const real max_travel = 1000.0 * avg_len;

  liblombardi::GraphContext ctx;
  auto body = ctx.create_node<duchamp::body_constant_node>(
      duchamp::make_shell_body(M));
  auto positions = ctx.create_node<duchamp::position_snapshot_node>();
  auto normals = ctx.create_node<duchamp::vertex_normals_snapshot_node>();
  auto fit = ctx.create_node<duchamp::darboux_cyclide_fit_node>(l0, 3.0);
  auto search =
      ctx.create_node<duchamp::medial_shape_energy_search_node>(max_travel);
  ctx.link(body->output(), positions->body());
  ctx.link(body->output(), normals->body());
  ctx.link(body->output(), fit->body());
  ctx.link(positions->output(), fit->pov());
  ctx.link(normals->output(), fit->n_pov());
  ctx.link(positions->output(), search->pov());
  ctx.link(fit->output(), search->cyclide());
  ctx.link(normals->output(), search->n_pov());
  ctx.run();

  const auto &mask =
      search->get_datum<duchamp::medial_shape_energy_search_node::MaskPortDef>()
          ->data();
  GAUDI_ASSERT(duchamp::count_medial_accepted(mask) > 0);
}

GAUDI_TEST(medial_backend_compare_sphere) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real avg_len = asawa::shell::avg_length(*M, x);
  const real l0 = 0.2 * avg_len;
  const real max_travel = 1000.0 * avg_len;
  const real max_newton_step = 1000.0 * l0;
  const std::vector<albers::vec14> fits =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  std::vector<duchamp::medial_point> legacy_results;
  std::vector<duchamp::medial_point> shape_results;
  legacy_results.reserve(fits.size());
  shape_results.reserve(fits.size());
  for (size_t i = 0; i < fits.size(); ++i) {
    legacy_results.push_back(duchamp::search_medial_legacy_ridge(
        fits[i], x[i], N[i], max_travel, 100, 1e-8, max_newton_step, -0.5));
    shape_results.push_back(duchamp::search_medial_shape_energy(
        fits[i], x[i], N[i], max_travel, {}));
  }
  GAUDI_ASSERT(duchamp::count_medial_accepted(legacy_results) > 0);
  GAUDI_ASSERT(duchamp::count_medial_accepted(shape_results) > 0);
}

GAUDI_TEST(medial_hessian_frame_symmetry) {
  albers::vec14 Q = albers::vec14::Zero();
  Q[0] = 1.0;
  Q[1] = 1.0;
  Q[2] = 1.0;
  Q[10] = -1.0;
  const vec3 x(0.1, -0.05, 0.02);
  const mat3 W = albers::shape_operator_at(Q, x);
  const mat3 asym = W - W.transpose();
  GAUDI_ASSERT(asym.norm() < 1e-8);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_MEDIAL_GRAPH_TESTS_HPP__
