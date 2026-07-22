#ifndef __GAUDI_TEST_KUSAMA_CYCLIDE_SMOOTH_TESTS_HPP__
#define __GAUDI_TEST_KUSAMA_CYCLIDE_SMOOTH_TESTS_HPP__

#include <cmath>
#include <vector>

#include "darboux_generated.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/duchamp/body_datum.hpp"
#include "gaudi/duchamp/calder_graph_nodes.hpp"
#include "gaudi/duchamp/field_graph_nodes.hpp"
#include "gaudi/kusama/cyclide_jet_smooth.hpp"
#include "gaudi/test/test.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

namespace {

inline real jet_edge_disagreement(const albers::vec14 &Qi,
                                  const albers::vec14 &Qj,
                                  const vec3 &x_j_in_i, real alpha_G,
                                  real alpha_H) {
  const vec3 Gi =
      albers::medial_generated::darboux_grad_generated(Qi, x_j_in_i);
  const vec3 Gj =
      albers::medial_generated::darboux_grad_generated(Qj, vec3::Zero());
  const mat3 Hi =
      albers::medial_generated::darboux_hessian_generated(Qi, x_j_in_i);
  const mat3 Hj =
      albers::medial_generated::darboux_hessian_generated(Qj, vec3::Zero());
  return alpha_G * (Gi - Gj).squaredNorm() +
         alpha_H * (Hi - Hj).squaredNorm();
}

inline real total_jet_disagreement(asawa::shell::shell &M,
                                   const std::vector<vec3> &pov,
                                   const std::vector<albers::vec14> &Q,
                                   real alpha_G, real alpha_H) {
  const size_t n = pov.size();
  real total = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const asawa::shell::VertId vi =
        asawa::shell::vert_id(static_cast<asawa::shell::index_t>(i));
    M.const_for_each_vertex(vi, [&](asawa::shell::CornerId c,
                                    const asawa::shell::shell &Ms) {
      const int j = Ms.vert(Ms.next(c));
      if (j < 0 || static_cast<size_t>(j) >= n) {
        return;
      }
      const vec3 x_j_in_i = pov[static_cast<size_t>(j)] - pov[i];
      total += jet_edge_disagreement(Q[i], Q[static_cast<size_t>(j)], x_j_in_i,
                                     alpha_G, alpha_H);
    });
  }
  return total;
}

} // namespace

GAUDI_TEST(kusama_cyclide_jet_smooth_sphere_finite) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real l0 = 0.2 * asawa::shell::avg_length(*M, x);
  const std::vector<albers::vec14> Q_global =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  kusama::cyclide_jet_smooth_params params;
  const std::vector<albers::vec14> Q_smooth =
      kusama::cyclide_jet_smooth(*M, x, Q_global, params);
  GAUDI_ASSERT(Q_smooth.size() == Q_global.size());
  for (const albers::vec14 &q : Q_smooth) {
    GAUDI_ASSERT(q.allFinite());
  }
}

GAUDI_TEST(kusama_cyclide_jet_smooth_reduces_jet_disagreement) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real l0 = 0.2 * asawa::shell::avg_length(*M, x);
  const std::vector<albers::vec14> Q_global =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  kusama::cyclide_jet_smooth_params params;
  params.wi = 0.5;
  const std::vector<albers::vec14> Q_smooth =
      kusama::cyclide_jet_smooth(*M, x, Q_global, params);

  const real before = total_jet_disagreement(*M, x, Q_global, params.alpha_G,
                                             params.alpha_H);
  const real after = total_jet_disagreement(*M, x, Q_smooth, params.alpha_G,
                                            params.alpha_H);
  GAUDI_ASSERT(after <= before + real(1e-6));
}

GAUDI_TEST(kusama_cyclide_jet_smooth_anchor_preserved) {
  auto M = asawa::shell::load_sphere(1.0, 10, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real l0 = 0.2 * asawa::shell::avg_length(*M, x);
  const std::vector<albers::vec14> Q_global =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  kusama::cyclide_jet_smooth_params params;
  params.wi = 0.5;
  params.alpha_G = 0.0;
  params.alpha_H = 0.0;
  const std::vector<albers::vec14> Q_smooth =
      kusama::cyclide_jet_smooth(*M, x, Q_global, params);

  real max_delta = 0.0;
  for (size_t i = 0; i < Q_global.size(); ++i) {
    max_delta = std::max(max_delta, (Q_smooth[i] - Q_global[i]).norm());
  }
  GAUDI_ASSERT(max_delta < real(1e-10));
}

GAUDI_TEST(duchamp_cyclide_smooth_node_matches_direct) {
  auto M = asawa::shell::load_sphere(1.0, 10, 8);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
  const std::vector<vec3> N = asawa::shell::vertex_normals(*M, x);
  const real l0 = 0.15 * asawa::shell::avg_length(*M, x);
  const std::vector<albers::vec14> Q_fit =
      calder::darboux_cyclide_shell_fit(*M, x, N, l0, 3.0);

  kusama::cyclide_jet_smooth_params params;
  const std::vector<albers::vec14> direct =
      kusama::cyclide_jet_smooth(*M, x, Q_fit, params);

  liblombardi::GraphContext ctx;
  auto body = ctx.create_node<duchamp::body_constant_node>(
      duchamp::make_shell_body(M));
  auto positions = ctx.create_node<duchamp::position_snapshot_node>();
  auto normals = ctx.create_node<duchamp::vertex_normals_snapshot_node>();
  auto fit = ctx.create_node<duchamp::darboux_cyclide_fit_node>(l0, 3.0);
  auto smooth = ctx.create_node<duchamp::darboux_cyclide_smooth_node>(params);
  ctx.link(body->output(), positions->body());
  ctx.link(body->output(), normals->body());
  ctx.link(body->output(), fit->body());
  ctx.link(positions->output(), fit->pov());
  ctx.link(normals->output(), fit->n_pov());
  ctx.link(body->output(), smooth->body());
  ctx.link(positions->output(), smooth->pov());
  ctx.link(fit->output(), smooth->cyclide_in());
  ctx.run();

  const auto &got =
      smooth->get_datum<duchamp::darboux_cyclide_smooth_node::CyclideOutPortDef>()
          ->data();
  GAUDI_ASSERT(got.size() == direct.size());
  real max_err = 0.0;
  for (size_t i = 0; i < direct.size(); ++i) {
    max_err = std::max(max_err, (got[i] - direct[i]).norm());
  }
  GAUDI_ASSERT(max_err < real(1e-10));
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_KUSAMA_CYCLIDE_SMOOTH_TESTS_HPP__
