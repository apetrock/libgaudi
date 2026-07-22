#ifndef __GAUDI_TEST_FIELD_GRAPH_TESTS_HPP__
#define __GAUDI_TEST_FIELD_GRAPH_TESTS_HPP__

#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/test/test.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

namespace {

inline std::vector<vec3> make_rod_loop_points(int n = 32) {
  std::vector<vec3> points;
  points.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    const real t = 2.0 * M_PI * real(i) / real(n);
    points.emplace_back(cos(t), sin(t), 0.0);
  }
  return points;
}

} // namespace

GAUDI_TEST(field_graph_shell_position_snapshot) {
  auto M = asawa::shell::load_cube();
  const std::vector<vec3> expect = asawa::const_get_vec_data(*M, 0);

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_shell_body(M));
  auto snap = ctx.create_node<duchamp::position_snapshot_node>();
  ctx.link(body_src->output(), snap->body());
  ctx.run();

  const std::vector<vec3> &got =
      snap->get_datum<duchamp::position_snapshot_node::OutputPortDef>()->data();
  GAUDI_ASSERT(got.size() == expect.size());
  real max_err = 0.0;
  for (size_t i = 0; i < expect.size(); ++i) {
    max_err = std::max(max_err, (got[i] - expect[i]).norm());
  }
  GAUDI_ASSERT(max_err < 1e-12);
}

GAUDI_TEST(field_graph_rod_position_snapshot) {
  auto R = asawa::rod::rod::create(make_rod_loop_points());
  const std::vector<vec3> expect = R->xc();

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_rod_body(R));
  auto snap = ctx.create_node<duchamp::position_snapshot_node>();
  ctx.link(body_src->output(), snap->body());
  ctx.run();

  const std::vector<vec3> &got =
      snap->get_datum<duchamp::position_snapshot_node::OutputPortDef>()->data();
  GAUDI_ASSERT(got.size() == expect.size());
  real max_err = 0.0;
  for (size_t i = 0; i < expect.size(); ++i) {
    max_err = std::max(max_err, (got[i] - expect[i]).norm());
  }
  GAUDI_ASSERT(max_err < 1e-12);
}

GAUDI_TEST(field_graph_snapshot_to_pov_wiring) {
  auto M = asawa::shell::load_sphere(1.0, 16, 12);

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_shell_body(M));
  auto snap = ctx.create_node<duchamp::position_snapshot_node>();
  ctx.link(body_src->output(), snap->body());
  ctx.run();

  const std::vector<vec3> &pov =
      snap->get_datum<duchamp::position_snapshot_node::OutputPortDef>()->data();
  GAUDI_ASSERT(!pov.empty());
  GAUDI_ASSERT(pov.size() == static_cast<size_t>(M->vert_count()));
}

GAUDI_TEST(field_graph_shell_nbody_eval_matches_calder) {
  auto M = asawa::shell::load_sphere(1.0, 24, 16);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const real l0 = asawa::shell::avg_length(*M, x);

  std::vector<real> source(static_cast<size_t>(M->face_count()), 2.5);
  const std::vector<vec3> pov = {vec3(0.2, 0.1, 0.9), vec3(-0.3, 0.4, 0.7)};

  const std::vector<real> expect =
      calder::mls_avg(*M, source, pov, l0);

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_shell_body(M));
  auto pov_src = ctx.create_node<duchamp::pov_constant_node>(pov);
  auto eval = ctx.create_node<duchamp::nbody_eval_node<duchamp::shell_mls_avg_tag>>(
      source, l0);
  ctx.link(body_src->output(), eval->body());
  ctx.link(pov_src->output(), eval->pov());
  ctx.run();

  const std::vector<real> &got =
      eval->get_datum<duchamp::nbody_eval_node<duchamp::shell_mls_avg_tag>::OutputPortDef>()
          ->data();
  GAUDI_ASSERT(got.size() == expect.size());
  real max_err = 0.0;
  for (size_t i = 0; i < expect.size(); ++i) {
    max_err = std::max(max_err, std::abs(got[i] - expect[i]));
  }
  GAUDI_ASSERT(max_err < 1e-10);
}

GAUDI_TEST(field_graph_rod_nbody_eval_matches_calder) {
  auto R = asawa::rod::rod::create(make_rod_loop_points());
  const real l0 = 0.25 * R->lavg();
  const auto edge_range = R->get_vert_range();
  std::vector<vec3> source(edge_range.size(), vec3(1.75, -0.25, 0.5));
  const std::vector<vec3> pov = {vec3(0.5, 0.0, 0.0), vec3(-0.5, 0.0, 0.0)};

  const std::vector<vec3> expect = calder::mls_avg(*R, source, pov, l0);

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_rod_body(R));
  auto pov_src = ctx.create_node<duchamp::pov_constant_node>(pov);
  auto eval = ctx.create_node<duchamp::nbody_eval_node<duchamp::rod_mls_avg_tag>>(
      source, l0);
  ctx.link(body_src->output(), eval->body());
  ctx.link(pov_src->output(), eval->pov());
  ctx.run();

  const std::vector<vec3> &got =
      eval->get_datum<duchamp::nbody_eval_node<duchamp::rod_mls_avg_tag>::OutputPortDef>()
          ->data();
  GAUDI_ASSERT(got.size() == expect.size());
  real max_err = 0.0;
  for (size_t i = 0; i < expect.size(); ++i) {
    max_err = std::max(max_err, (got[i] - expect[i]).norm());
  }
  GAUDI_ASSERT(max_err < 1e-10);
}

GAUDI_TEST(field_graph_snapshot_pov_to_nbody_eval) {
  auto M = asawa::shell::load_sphere(1.0, 20, 14);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const real l0 = asawa::shell::avg_length(*M, x);
  std::vector<real> source(static_cast<size_t>(M->face_count()), 3.0);

  liblombardi::GraphContext ctx;
  auto body_src =
      ctx.create_node<duchamp::body_constant_node>(duchamp::make_shell_body(M));
  auto snap = ctx.create_node<duchamp::position_snapshot_node>();
  auto eval = ctx.create_node<duchamp::nbody_eval_node<duchamp::shell_mls_avg_tag>>(
      source, l0);

  ctx.link(body_src->output(), snap->body());
  ctx.link(body_src->output(), eval->body());
  ctx.link(snap->output(), eval->pov());
  ctx.run();

  const std::vector<vec3> &pov =
      snap->get_datum<duchamp::position_snapshot_node::OutputPortDef>()->data();
  const std::vector<real> expect = calder::mls_avg(*M, source, pov, l0);
  const std::vector<real> &got =
      eval->get_datum<duchamp::nbody_eval_node<duchamp::shell_mls_avg_tag>::OutputPortDef>()
          ->data();

  GAUDI_ASSERT(got.size() == expect.size());
  real max_err = 0.0;
  for (size_t i = 0; i < expect.size(); ++i) {
    max_err = std::max(max_err, std::abs(got[i] - expect[i]));
  }
  GAUDI_ASSERT(max_err < 1e-10);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_FIELD_GRAPH_TESTS_HPP__
