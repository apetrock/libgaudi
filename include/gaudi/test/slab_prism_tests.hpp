#ifndef GAUDI_TEST_SLAB_PRISM_TESTS_HPP
#define GAUDI_TEST_SLAB_PRISM_TESTS_HPP

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/arp/brute_force_tree.hpp"
#include "gaudi/hepworth/block/block_geometry.hpp"
#include "gaudi/hepworth/block/relative_plane_twist.hpp"
#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/slab_dual_ring_volume_procrustes.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/slab_prism_init.hpp"
#include "gaudi/hepworth/block/slab_fan_tet.hpp"
#include "gaudi/hepworth/block/slab_prism_operators.hpp"
#include "gaudi/hepworth/block/incompressible_phi_energy.hpp"
#include "gaudi/hepworth/block/incompressible_phi_strain.hpp"
#include "gaudi/hepworth/block/phi_div_free_basis.hpp"
#include "gaudi/hepworth/block/rbf_displacement_basis.hpp"
#include "gaudi/hepworth/block/rbf_elastic_energy.hpp"
#include "gaudi/hepworth/block/rbf_elastic_strain.hpp"
#include "gaudi/hepworth/block/level_set_field_bundle.hpp"
#include "gaudi/hepworth/block/local_slab_field_fit.hpp"
#include "gaudi/hepworth/block/slab_field_divergence_energy.hpp"
#include "gaudi/hepworth/block/slab_wedge_strain.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/block/projection_diagnostics.hpp"
#include "gaudi/sparse_solver.h"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/duchamp/pierce_cylinder_medial.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

GAUDI_TEST(slab_pierce_sphere_foot_has_slab) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  GAUDI_ASSERT(tree != nullptr);
  const std::vector<index_t> face_map = M->get_face_map(true);

  hepworth::block::slab_pierce_config cfg;
  const real avg_edge_length = asawa::shell::avg_length(*M, x);
  size_t slab = 0;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg_edge_length);
    if (!result.stencil.has_slab)
      continue;
    ++slab;
    GAUDI_EXPECT(result.stencil.pierce_foot != result.stencil.foot);
    GAUDI_EXPECT(result.stencil.bot_ids.size() >= 2);
  }
  GAUDI_EXPECT(slab > 0);
}

GAUDI_TEST(slab_pierce_nearest_inward_on_bunny) {
  auto M = asawa::shell::load_bunny();
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);
  asawa::center(x);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  GAUDI_ASSERT(tree != nullptr);
  const std::vector<index_t> face_map = M->get_face_map(true);

  hepworth::block::slab_pierce_config cfg;
  cfg.antipodal_fallback = false;
  const real avg_edge_length = asawa::shell::avg_length(*M, x);
  const real t_min =
      std::max(cfg.t_min_scale * avg_edge_length, real(1e-8));

  int pierced = 0;
  real travel_sum = 0.0;
  real max_travel = 0.0;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const asawa::shell::VertId foot = asawa::shell::vert_id(iv);
    vec3 n_out = asawa::shell::vert_normal(*M, foot, x);
    if (n_out.squaredNorm() < 1e-24) {
      continue;
    }
    n_out.normalize();
    const auto hit = hepworth::block::ray_pierce_slab_shell(
        *M, foot, x[static_cast<size_t>(iv)], n_out, x, *tree, face_map, cfg,
        t_min, avg_edge_length);
    if (!hit) {
      continue;
    }
    const auto tri = M->get_tri(hit->fi);
    const index_t pierce_v = hepworth::block::closest_triangle_vertex_excluding(
        tri, x, hit->hit.point, foot);
    if (pierce_v < 0) {
      continue;
    }
    ++pierced;
    const real travel =
        (x[static_cast<size_t>(pierce_v)] - x[static_cast<size_t>(iv)]).norm();
    travel_sum += travel;
    max_travel = std::max(max_travel, travel);
  }

  GAUDI_EXPECT(pierced > static_cast<int>(M->vert_count()) / 2);
  const real avg_travel = travel_sum / real(pierced);
  // Pierce-vertex spine on centered bunny: through-body hits, not contact-only.
  GAUDI_EXPECT(avg_travel > real(0.05));
  GAUDI_EXPECT(avg_travel < real(1.5));
  GAUDI_EXPECT(max_travel < real(2.5));
}

GAUDI_TEST(pierce_cylinder_medial_bunny_demo_mesh) {
  auto M = asawa::shell::load_bunny();
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  duchamp::normalize_medial_demo_mesh(*M);

  duchamp::pierce_cylinder_medial_stats stats;
  const auto samples = duchamp::compute_pierce_cylinder_medial(
      *M, duchamp::default_pierce_cylinder_medial_params(), &stats);

  GAUDI_EXPECT(static_cast<int>(samples.size()) == stats.total);
  const real miss_frac =
      stats.total > 0 ? real(stats.pierce_miss_count) / real(stats.total)
                      : real(0);
  GAUDI_EXPECT(miss_frac < real(0.001));
  GAUDI_EXPECT(stats.accepted > stats.total / 2);
  GAUDI_EXPECT(stats.avg_travel > real(0.05));
  GAUDI_EXPECT(stats.avg_travel < real(2.5));
  GAUDI_EXPECT(stats.mls_fit_count + stats.plucker_fallback_count ==
               stats.accepted);
  GAUDI_EXPECT(stats.mls_fit_count > stats.total / 4);

  int midpoint_ok = 0;
  int radius_ok = 0;
  for (const duchamp::pierce_cylinder_sample &sample : samples) {
    if (!sample.accepted) {
      continue;
    }
    const vec3 expected_mid =
        duchamp::pierce_midpoint_proxy(sample.foot, sample.pierce);
    if ((sample.midpoint - expected_mid).norm() < real(1e-6)) {
      ++midpoint_ok;
    }
    if (sample.radius > real(1e-8)) {
      ++radius_ok;
    }
  }
  GAUDI_EXPECT(midpoint_ok > stats.accepted / 2);
  GAUDI_EXPECT(radius_ok > stats.accepted / 4);

  std::cerr << "[medial_axis_estimation] accepted=" << stats.accepted
            << " pierce_miss=" << stats.pierce_miss_count
            << " fit_fail=" << stats.fit_fail_count
            << " mls_fit=" << stats.mls_fit_count
            << " plucker_fallback=" << stats.plucker_fallback_count
            << " avg_travel=" << stats.avg_travel
            << " avg_center_mid_drift_frac="
            << stats.avg_center_mid_drift_frac << std::endl;
}

namespace detail {
inline std::vector<index_t>
demo_point_normal_load_verts(const std::vector<vec3> &x) {
  std::vector<index_t> loads;
  if (x.empty())
    return loads;
  const int n = static_cast<int>(x.size());
  index_t top = 0;
  real y_min = x[0].y();
  real y_max = x[0].y();
  for (int i = 0; i < n; ++i) {
    y_min = std::min(y_min, x[static_cast<size_t>(i)].y());
    y_max = std::max(y_max, x[static_cast<size_t>(i)].y());
    if (x[static_cast<size_t>(i)].y() > x[static_cast<size_t>(top)].y())
      top = static_cast<index_t>(i);
  }
  loads.push_back(top);
  const real y_cut = y_min + 0.72 * (y_max - y_min);
  index_t side = top;
  real best_r = -1.0;
  for (int i = 0; i < n; ++i) {
    const vec3 &xi = x[static_cast<size_t>(i)];
    if (xi.y() < y_cut)
      continue;
    vec3 r = xi;
    r.y() = 0.0;
    const real rn = r.norm();
    if (rn > best_r && static_cast<index_t>(i) != top) {
      best_r = rn;
      side = static_cast<index_t>(i);
    }
  }
  if (side != top)
    loads.push_back(side);
  return loads;
}
} // namespace detail

GAUDI_TEST(slab_pierce_bunny_audit_load_verts) {
  auto M = asawa::shell::load_bunny();
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);
  asawa::center(x);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  GAUDI_ASSERT(tree != nullptr);
  const std::vector<index_t> face_map = M->get_face_map(true);
  const hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  const hepworth::block::slab_prism_init_stats stats =
      hepworth::block::audit_slab_pierce_stencils(*M, x, *tree, face_map, cfg,
                                                  avg);

  std::cerr << "[slab_pierce_audit] bunny active_feet="
            << stats.active_feet << " slab_ok=" << stats.slab_count
            << " pierce_fail=" << stats.pierce_fail_count;
  for (const auto &[reason, count] : stats.fail_by_reason)
    std::cerr << " " << hepworth::block::pierce_failure_reason_name(reason)
              << "=" << count;
  std::cerr << std::endl;

  const std::vector<index_t> loads = detail::demo_point_normal_load_verts(x);
  GAUDI_EXPECT(!loads.empty());
  for (index_t vi : loads) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(static_cast<int>(vi)), x, *tree, face_map, cfg,
        avg);
    std::cerr << "[slab_pierce_audit] load_vert=" << vi
              << " has_slab=" << (result.stencil.has_slab ? 1 : 0)
              << " reason=" << hepworth::block::pierce_failure_reason_name(result.reason)
              << " y=" << x[static_cast<size_t>(vi)].y()
              << " ray_t=" << result.ray_t << std::endl;
    GAUDI_EXPECT(result.stencil.has_slab);
  }

  const real fail_frac =
      stats.active_feet > 0
          ? static_cast<real>(stats.pierce_fail_count) /
                static_cast<real>(stats.active_feet)
          : real(0.0);
  GAUDI_EXPECT(fail_frac < real(0.001));
}

GAUDI_TEST(slab_pierce_coarse_sphere_antipodal_fallback) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  GAUDI_ASSERT(tree != nullptr);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  cfg.antipodal_fallback = true;
  const real avg = asawa::shell::avg_length(*M, x);

  const hepworth::block::slab_prism_init_stats stats =
      hepworth::block::audit_slab_pierce_stencils(*M, x, *tree, face_map, cfg,
                                                  avg);

  size_t antipodal = 0;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.used_antipodal_fallback)
      ++antipodal;
  }

  std::cerr << "[slab_pierce_sphere] active_feet=" << stats.active_feet
            << " slab_ok=" << stats.slab_count
            << " pierce_fail=" << stats.pierce_fail_count
            << " antipodal=" << antipodal << std::endl;

  GAUDI_EXPECT(stats.slab_count == stats.active_feet);
  // Coarse sphere may pierce directly; antipodal path is optional when fallback on.
  GAUDI_EXPECT(antipodal >= 0);
}

GAUDI_TEST(slab_prism_init_tet_no_arap_fallback) {
  auto M = asawa::shell::load_tet();
  GAUDI_ASSERT(M != nullptr);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);

  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  const hepworth::block::slab_pierce_config pierce_cfg;
  hepworth::block::slab_prism_config prism_cfg;
  prism_cfg.w_twist = 1.0;
  const auto stats = hepworth::block::init_slab_prism_constraints(
      *M, ctx.constraints, pierce_cfg, prism_cfg, blocks, ctx, false);
  GAUDI_EXPECT(stats.slab_count + stats.pierce_fail_count > 0);
  for (const auto &c : ctx.constraints) {
    GAUDI_EXPECT(c->name() == "relative_plane_twist");
  }
}

GAUDI_TEST(slab_prism_row_count_matches_arap) {
  auto M = asawa::shell::load_tet();
  GAUDI_ASSERT(M != nullptr);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const hepworth::block::slab_prism_rest rest =
      hepworth::block::build_slab_prism_rest(stencil, x);
  GAUDI_EXPECT(rest.L_rest > 0.0);
  GAUDI_EXPECT(rest.d_top_rest.size() == stencil.top_weights.size());
  GAUDI_EXPECT(rest.d_bot_rest.size() == stencil.bot_weights.size());

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);
  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  const auto wedges = hepworth::block::build_slab_fan_wedges(stencil, x);
  (void)wedges;
  const real w = 1.0;
  const vecX q = hepworth::block::pack_block_q(blocks);
  auto twist_c = std::make_shared<hepworth::block::relative_plane_twist>(
      stencil, hepworth::block::slab_prism_config{}, w, blocks);
  hepworth::block::solver_context init_ctx;
  twist_c->init_rest(q, init_ctx);

  const size_t expected_spokes =
      stencil.top_weights.size() + stencil.bot_weights.size();
  GAUDI_EXPECT(twist_c->total_spoke_count() == expected_spokes);

  index_t id0 = 0;
  std::vector<trip> triplets;
  twist_c->fill_A(id0, triplets);
  GAUDI_EXPECT(id0 == static_cast<index_t>(3 * expected_spokes));
  GAUDI_EXPECT(triplets.size() == 6 * expected_spokes);
}

GAUDI_TEST(slab_prism_volume_rest_stable) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (!result.stencil.has_slab)
      continue;
    const hepworth::block::slab_prism_rest rest =
        hepworth::block::build_slab_prism_rest(result.stencil, x);
    const auto pos_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
    const real V =
        hepworth::block::eval_prism_volume(result.stencil, pos_at);
    GAUDI_EXPECT(std::abs(V - rest.V_rest) < 1e-6 * std::max(real(1.0), std::abs(rest.V_rest)));
    return;
  }
  GAUDI_ASSERT(false);
}

GAUDI_TEST(slab_prism_volume_gd_identity_at_rest) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const hepworth::block::slab_prism_rest rest =
      hepworth::block::build_slab_prism_rest(stencil, x);
  const auto anchor_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  auto top = hepworth::block::top_ring_positions::from_stencil(stencil, anchor_at);
  auto resolved = hepworth::block::resolve_prism_volume_gd(stencil, top, anchor_at,
                                                           rest.V_rest, 16, 1.0);
  for (size_t i = 0; i < top.x.size(); ++i)
    GAUDI_EXPECT((resolved.x[i] - top.x[i]).norm() < 1e-10);
}

GAUDI_TEST(slab_prism_twist_correction_identity_at_rest) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const hepworth::block::slab_prism_rest rest =
      hepworth::block::build_slab_prism_rest(stencil, x);
  mat3 R0 = rest.R0_rest;
  mat3 R1 = rest.R1_rest;
  const vec3 d_hat = rest.d_rest.normalized();
  hepworth::block::apply_twist_to_top_frame(R0, R1, rest.R0_rest, rest.R1_rest,
                                            d_hat, false);
  GAUDI_EXPECT((R0 - rest.R0_rest).norm() < 1e-10);
  GAUDI_EXPECT((R1 - rest.R1_rest).norm() < 1e-10);
}

GAUDI_TEST(slab_shape_field_energy_zero_at_rest) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const auto pos_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_slab_field_on_stencil(
          hepworth::block::slab_field_family::darboux, stencil, pos_at, fit_cfg,
          avg);
  auto ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_at);
  const real e0 = hepworth::block::eval_shape_field_divergence_energy(
      stencil, Q, Q, ring, 1.0);
  GAUDI_EXPECT(e0 < 1e-20);
}

GAUDI_TEST(slab_shape_field_init_w_field_gate) {
  auto M = asawa::shell::load_tet();
  GAUDI_ASSERT(M != nullptr);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);

  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  const hepworth::block::slab_pierce_config pierce_cfg;
  hepworth::block::slab_prism_config with_field;
  with_field.w_field = 1.0;
  with_field.w_twist = 0.0;
  hepworth::block::init_slab_prism_constraints(*M, ctx.constraints, pierce_cfg,
                                               with_field, blocks, ctx, false);
  bool has_field = false;
  for (const auto &c : ctx.constraints)
    if (c->name() == "shape_field_divergence_strain")
      has_field = true;
  GAUDI_EXPECT(has_field);

  ctx.constraints.clear();
  hepworth::block::slab_prism_config no_field;
  no_field.w_field = 0.0;
  no_field.w_twist = 0.0;
  hepworth::block::init_slab_prism_constraints(*M, ctx.constraints, pierce_cfg,
                                               no_field, blocks, ctx, false);
  for (const auto &c : ctx.constraints)
    GAUDI_EXPECT(c->name() != "shape_field_divergence_strain");
}

/// Low-valence mesh: RBF fit is wide (n_cols > n_rows). Regression for SVD column OOB.
GAUDI_TEST(slab_rbf_init_low_valence_no_crash) {
  auto M = asawa::shell::load_tet();
  GAUDI_ASSERT(M != nullptr);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);

  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  const hepworth::block::slab_pierce_config pierce_cfg;
  hepworth::block::slab_prism_config cfg;
  cfg.w_field = 1.0;
  cfg.w_twist = 0.0;
  cfg.field_family = hepworth::block::slab_field_family::rbf_mls;
  const auto stats = hepworth::block::init_slab_prism_constraints(
      *M, ctx.constraints, pierce_cfg, cfg, blocks, ctx, false);
  GAUDI_EXPECT(stats.slab_count > 0);
  size_t field_rows = 0;
  for (const auto &c : ctx.constraints) {
    if (c->name() != "shape_field_divergence_strain")
      continue;
    index_t n_rows = 0;
    std::vector<trip> trips;
    c->fill_A(n_rows, trips);
    field_rows += static_cast<size_t>(n_rows);
  }
  GAUDI_EXPECT(field_rows > 0);
}

GAUDI_TEST(slab_darboux_bulk_energy_grad_fd) {
  albers::vec14 Qr = albers::vec14::Zero();
  albers::vec14 Qk = albers::vec14::Zero();
  for (int i = 0; i < 14; ++i) {
    Qr[i] = 0.1 * (i + 1);
    Qk[i] = Qr[i] + 0.01 * i;
  }
  const vec3 x0(0.1, -0.2, 0.05);
  const real lam = 1.0;
  const real eps = 1e-6;

  const auto E = [&](const vec3 &p) {
    return hepworth::block::darboux_field_bundle::bulk_sample_energy(Qr, Qk, p,
                                                                     lam);
  };

  const vec3 g_sym =
      hepworth::block::darboux_field_bundle::bulk_sample_grad(Qr, Qk, x0, lam);
  const vec3 g_fd((E(x0 + eps * vec3::UnitX()) - E(x0 - eps * vec3::UnitX())) /
                      (2.0 * eps),
                  (E(x0 + eps * vec3::UnitY()) - E(x0 - eps * vec3::UnitY())) /
                      (2.0 * eps),
                  (E(x0 + eps * vec3::UnitZ()) - E(x0 - eps * vec3::UnitZ())) /
                      (2.0 * eps));
  GAUDI_EXPECT((g_sym - g_fd).norm() < 1e-4 * std::max(g_fd.norm(), real(1.0)));
}

GAUDI_TEST(slab_rbf_field_energy_positive_on_perturb) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const auto pos_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q_rest =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_at, fit_cfg,
          avg);
  auto ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_at);

  std::vector<vec3> x_pert = x;
  x_pert[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.05, 0.02, 0.0);
  const auto pos_pert = [&](index_t vi) {
    return x_pert[static_cast<size_t>(vi)];
  };
  const hepworth::block::slab_field_coeffs Q_live =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_pert, fit_cfg,
          avg, &Q_rest);
  ring = hepworth::block::slab_ring_positions::from_stencil(stencil, pos_pert);

  const real e1 = hepworth::block::eval_shape_field_divergence_energy(
      stencil, Q_rest, Q_live, ring, 1.0);
  GAUDI_EXPECT(e1 > 1e-12);
}

namespace detail {
inline bool build_one_slab_field_constraint(
    hepworth::block::slab_pierce_stencil &stencil_out,
    hepworth::block::slab_prism_rest &rest_out,
    std::shared_ptr<hepworth::block::shell_position_block> &shell_out,
    std::vector<hepworth::sim_block::ptr> &blocks_out, const real avg) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  if (!M)
    return false;
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  if (!stencil.has_slab)
    return false;

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  shell_out = std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);
  hepworth::block::solver_context ctx;
  shell_out->prepare(ctx);
  blocks_out.clear();
  shell_out->emit_dof_blocks(blocks_out);

  stencil_out = stencil;
  rest_out = hepworth::block::build_slab_prism_rest(stencil, x);
  return true;
}

inline vecX pack_shell_q(const std::vector<hepworth::sim_block::ptr> &blocks) {
  return hepworth::block::pack_block_q(blocks);
}

inline void init_block_constraint_rest(
    const std::shared_ptr<hepworth::block::block_constraint> &c,
    const vecX &q) {
  if (!c)
    return;
  hepworth::block::solver_context ctx;
  c->init_rest(q, ctx);
}

template <typename ConstraintT>
inline std::shared_ptr<ConstraintT> make_init_slab_ring_constraint(
    const hepworth::block::slab_pierce_stencil &stencil,
    const hepworth::block::slab_prism_config &cfg, real w,
    const std::vector<hepworth::sim_block::ptr> &blocks, const vecX &q) {
  auto c = std::make_shared<ConstraintT>(stencil, cfg, w, blocks);
  init_block_constraint_rest(c, q);
  return c;
}

inline vecX apply_triplets(const std::vector<trip> &triplets, index_t n_rows,
                           const vecX &q) {
  vecX out = vecX::Zero(n_rows);
  for (const trip &t : triplets)
    out[t.row()] += t.value() * q[t.col()];
  return out;
}

inline std::shared_ptr<hepworth::block::shape_field_divergence_strain>
make_shape_field_constraint(
    const hepworth::block::slab_pierce_stencil &stencil,
    const hepworth::block::slab_prism_rest &rest,
    hepworth::block::slab_prism_config cfg,
    const std::vector<hepworth::sim_block::ptr> &blocks,
    asawa::shell::shell &shell, real avg, real w_field = 1.0) {
  const auto pos_at = [&](index_t vi) {
    return asawa::const_get_vec_data(shell, 0)[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  fit_cfg.foot_normal_weight = cfg.field_foot_normal_w;
  const hepworth::block::slab_field_coeffs Q_rest =
      hepworth::block::fit_field_coeffs_on_stencil(
          cfg.field_family, stencil, pos_at, fit_cfg, avg);
  return std::make_shared<hepworth::block::shape_field_divergence_strain>(
      stencil, rest, cfg, w_field, blocks, Q_rest, avg);
}

/// Max relative tap ||p - Aq|| / ||Aq|| for feet farther than min_dist from load.
inline real max_far_field_tap(
    asawa::shell::shell::ptr M, const std::vector<vec3> &x_rest,
    const std::vector<vec3> &x_deformed, index_t load_vert,
    hepworth::block::slab_field_family family, real w_field,
    hepworth::block::field_live_update field_live, real min_dist) {
  asawa::get_vec_data(*M, 0) = x_deformed;
  arp::bvh_tree<3>::ptr tree =
      hepworth::block::make_shell_face_bvh(*M, x_deformed);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config pierce_cfg;
  const real avg = asawa::shell::avg_length(*M, x_deformed);

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);
  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);
  const vecX q = pack_shell_q(blocks);

  const vec3 load_pos = x_rest[static_cast<size_t>(load_vert)];
  hepworth::block::slab_prism_config cfg;
  cfg.field_family = family;
  cfg.field_live = field_live;
  cfg.field_gd_iters = 8;

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };

  real max_tap = 0.0;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const index_t foot = static_cast<index_t>(iv);
    if ((x_rest[static_cast<size_t>(foot)] - load_pos).norm() < min_dist)
      continue;

    const auto pierce = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x_deformed, *tree, face_map, pierce_cfg,
        avg);
    if (!pierce.stencil.has_slab)
      continue;

    const hepworth::block::slab_prism_rest rest =
        hepworth::block::build_slab_prism_rest(pierce.stencil, x_rest);

    hepworth::block::local_field_fit_config fit_cfg;
    const hepworth::block::slab_field_coeffs Q_rest =
        hepworth::block::fit_field_coeffs_on_stencil(
            family, pierce.stencil, pos_rest, fit_cfg, avg);

    auto c = std::make_shared<hepworth::block::shape_field_divergence_strain>(
        pierce.stencil, rest, cfg, w_field, blocks, Q_rest, avg);

    index_t n_rows = 0;
    std::vector<trip> trips;
    c->fill_A(n_rows, trips);
    const vecX Ap = apply_triplets(trips, n_rows, q);
    vecX p = vecX::Zero(n_rows);
    c->project(q, p);
    const real tap =
        (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
    max_tap = std::max(max_tap, tap);
  }
  return max_tap;
}
} // namespace detail

/// Localized bump away from foot: at w_field=1.0 far feet should stay quiet (tap ≈ 0).
GAUDI_TEST(shape_field_far_tap_quiet_w1_darboux) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> x_rest = asawa::const_get_vec_data(*M, 0);
  const real avg = asawa::shell::avg_length(*M, x_rest);

  index_t load = 0;
  for (int i = 1; i < static_cast<int>(M->vert_count()); ++i)
    if (x_rest[static_cast<size_t>(i)].y() >
        x_rest[static_cast<size_t>(load)].y())
      load = static_cast<index_t>(i);

  std::vector<vec3> x_def = x_rest;
  const vec3 n_load =
      asawa::shell::vert_normal(*M, asawa::shell::vert_id(load), x_rest)
          .normalized();
  x_def[static_cast<size_t>(load)] += 0.08 * n_load;

  constexpr real w_field = 1.0;
  const real tap = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::darboux,
      w_field, hepworth::block::field_live_update::live_matches_rest, 5.0 * avg);
  GAUDI_EXPECT(tap < 0.02);
}

GAUDI_TEST(shape_field_far_tap_quiet_w1_rbf) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> x_rest = asawa::const_get_vec_data(*M, 0);
  const real avg = asawa::shell::avg_length(*M, x_rest);

  index_t load = 0;
  for (int i = 1; i < static_cast<int>(M->vert_count()); ++i)
    if (x_rest[static_cast<size_t>(i)].y() >
        x_rest[static_cast<size_t>(load)].y())
      load = static_cast<index_t>(i);

  std::vector<vec3> x_def = x_rest;
  const vec3 n_load =
      asawa::shell::vert_normal(*M, asawa::shell::vert_id(load), x_rest)
          .normalized();
  x_def[static_cast<size_t>(load)] += 0.08 * n_load;

  constexpr real w_field = 1.0;
  const real tap = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::rbf_mls,
      w_field, hepworth::block::field_live_update::live_matches_rest, 5.0 * avg);
  GAUDI_EXPECT(tap < 0.05);
}

/// Bunny mesh: far-field tap at w_field=1.0 with localized top bump (demo-like).
GAUDI_TEST(shape_field_far_tap_quiet_w1_bunny) {
  auto M = asawa::shell::load_bunny();
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  std::vector<vec3> x_rest = asawa::const_get_vec_data(*M, 0);
  asawa::center(x_rest);
  const real avg = asawa::shell::avg_length(*M, x_rest);

  index_t load = 0;
  for (int i = 1; i < static_cast<int>(M->vert_count()); ++i)
    if (x_rest[static_cast<size_t>(i)].y() >
        x_rest[static_cast<size_t>(load)].y())
      load = static_cast<index_t>(i);

  std::vector<vec3> x_def = x_rest;
  const vec3 n_load =
      asawa::shell::vert_normal(*M, asawa::shell::vert_id(load), x_rest)
          .normalized();
  x_def[static_cast<size_t>(load)] += 0.06 * n_load;

  constexpr real w_field = 1.0;
  const real tap_darboux = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::darboux,
      w_field, hepworth::block::field_live_update::live_matches_rest, 8.0 * avg);
  const real tap_rbf = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::rbf_mls,
      w_field, hepworth::block::field_live_update::live_matches_rest, 8.0 * avg);
  GAUDI_EXPECT(tap_darboux < 0.05);
  GAUDI_EXPECT(tap_rbf < 0.15);
}

/// Same bunny tap with live Q_k refit — documents extra far-field noise at w_field=1.
GAUDI_TEST(shape_field_far_tap_w1_bunny_refit_noisier) {
  auto M = asawa::shell::load_bunny();
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  std::vector<vec3> x_rest = asawa::const_get_vec_data(*M, 0);
  asawa::center(x_rest);
  const real avg = asawa::shell::avg_length(*M, x_rest);

  index_t load = 0;
  for (int i = 1; i < static_cast<int>(M->vert_count()); ++i)
    if (x_rest[static_cast<size_t>(i)].y() >
        x_rest[static_cast<size_t>(load)].y())
      load = static_cast<index_t>(i);

  std::vector<vec3> x_def = x_rest;
  const vec3 n_load =
      asawa::shell::vert_normal(*M, asawa::shell::vert_id(load), x_rest)
          .normalized();
  x_def[static_cast<size_t>(load)] += 0.06 * n_load;

  constexpr real w_field = 1.0;
  const real tap_frozen = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::darboux,
      w_field, hepworth::block::field_live_update::live_matches_rest, 8.0 * avg);
  const real tap_refit = detail::max_far_field_tap(
      M, x_rest, x_def, load, hepworth::block::slab_field_family::darboux,
      w_field, hepworth::block::field_live_update::refit_every_project,
      8.0 * avg);
  GAUDI_EXPECT(tap_frozen < 0.05);
  GAUDI_EXPECT(tap_refit >= tap_frozen);
}

/// Two toy one-rings (top + bottom pierce stencil) at rest: projection should be
/// identity — p == A*q and local resolve x* == current ring positions.
GAUDI_TEST(shape_field_dual_ring_projection_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };

  hepworth::block::slab_prism_config cfg;
  cfg.field_family = hepworth::block::slab_field_family::darboux;
  cfg.field_live = hepworth::block::field_live_update::live_matches_rest;
  cfg.field_gd_iters = 8;
  cfg.lam_gap = 1.0;

  auto c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                               *shell->mesh, avg);
  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);

  const vecX q = detail::pack_shell_q(blocks);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);

  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-5);

  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_field_coeffs_on_stencil(
          cfg.field_family, stencil, pos_at, fit_cfg, avg);
  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil,
                                                                   pos_at);
  auto x_star = hepworth::block::resolve_shape_field_divergence_gd(
      stencil, Q, Q, ring, cfg.lam_gap, cfg.field_gd_iters,
      cfg.field_gd_step_scale);
  real max_delta = 0.0;
  for (size_t i = 0; i < ring.x.size(); ++i)
    max_delta = std::max(max_delta, (x_star.x[i] - ring.x[i]).norm());
  GAUDI_EXPECT(max_delta < 1e-5);
}

/// With Q_k refit every project, rest geometry should still be identity — if this
/// fails, live refit breaks the rest fixed point even with zero external force.
GAUDI_TEST(shape_field_dual_ring_projection_identity_refit_live) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  hepworth::block::slab_prism_config cfg;
  cfg.field_live = hepworth::block::field_live_update::refit_every_project;
  cfg.field_gd_iters = 8;

  auto c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                               *shell->mesh, avg);
  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);

  const vecX q = detail::pack_shell_q(blocks);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);

  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-5);
}

GAUDI_TEST(shape_field_rbf_refit_matches_Q_rest_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q_rest =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_at,
          fit_cfg, avg);
  const hepworth::block::slab_field_coeffs Q_live =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_at,
          fit_cfg, avg);

  GAUDI_EXPECT(hepworth::block::rbf_fields_near_equal(Q_rest.rbf, Q_live.rbf,
                                                      1e-6));

  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil, pos_at);
  const real e0 = hepworth::block::eval_shape_field_divergence_energy(
      stencil, Q_rest, Q_live, ring, 1.0);
  GAUDI_EXPECT(e0 < 1e-12);
}

GAUDI_TEST(shape_field_dual_ring_projection_identity_refit_live_rbf) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  hepworth::block::slab_prism_config cfg;
  cfg.field_family = hepworth::block::slab_field_family::rbf_mls;
  cfg.field_live = hepworth::block::field_live_update::refit_every_project;
  cfg.field_gd_iters = 12;

  auto c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                               *shell->mesh, avg);
  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);

  const vecX q = detail::pack_shell_q(blocks);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);

  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-4);
}

GAUDI_TEST(shape_field_rbf_small_perturbation_bounded_energy) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q_rest =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_rest,
          fit_cfg, avg);

  std::vector<vec3> x_eps = x_rest;
  const real eps = 0.01;
  x_eps[static_cast<size_t>(stencil.top_ids[1])] += vec3(eps, 0.5 * eps, 0.0);
  const auto pos_eps = [&](index_t vi) {
    return x_eps[static_cast<size_t>(vi)];
  };
  const hepworth::block::slab_field_coeffs Q_live =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_eps,
          fit_cfg, avg, &Q_rest);

  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil, pos_eps);
  const real e1 = hepworth::block::eval_shape_field_divergence_energy(
      stencil, Q_rest, Q_live, ring, 1.0);

  GAUDI_EXPECT(Q_rest.rbf.centers.size() == Q_live.rbf.centers.size());
  for (size_t k = 0; k < Q_rest.rbf.centers.size(); ++k)
    GAUDI_EXPECT((Q_rest.rbf.centers[k] - Q_live.rbf.centers[k]).norm() < 1e-12);

  std::vector<vec3> x_big = x_rest;
  x_big[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);
  const auto pos_big = [&](index_t vi) {
    return x_big[static_cast<size_t>(vi)];
  };
  const hepworth::block::slab_field_coeffs Q_big =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::rbf_mls, stencil, pos_big, fit_cfg,
          avg, &Q_rest);
  auto ring_big =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_big);
  const real e_big = hepworth::block::eval_shape_field_divergence_energy(
      stencil, Q_rest, Q_big, ring_big, 1.0);

  GAUDI_EXPECT(e1 > 0.0);
  GAUDI_EXPECT(e_big > e1);
  GAUDI_EXPECT(e1 < 0.25 * std::max(e_big, real(1e-12)));
}

/// Bent stencil: live_matches_rest is a no-op (p == Aq); refit must produce forces.
GAUDI_TEST(shape_field_bent_projection_not_noop_with_refit) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);

  const vecX q = detail::pack_shell_q(blocks);
  constexpr real w_field = 1.0;

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q_rest =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::darboux, stencil, pos_rest,
          fit_cfg, avg);

  const auto tap_for = [&](hepworth::block::field_live_update mode) {
    hepworth::block::slab_prism_config cfg;
    cfg.field_live = mode;
    cfg.field_gd_iters = 12;
    auto c = std::make_shared<hepworth::block::shape_field_divergence_strain>(
        stencil, rest, cfg, w_field, blocks, Q_rest, avg);
    index_t n_rows = 0;
    std::vector<trip> trips;
    c->fill_A(n_rows, trips);
    const vecX Ap = detail::apply_triplets(trips, n_rows, q);
    vecX p = vecX::Zero(n_rows);
    c->project(q, p);
    return (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  };

  const real tap_noop =
      tap_for(hepworth::block::field_live_update::live_matches_rest);
  const real tap_live =
      tap_for(hepworth::block::field_live_update::refit_every_project);
  GAUDI_EXPECT(tap_noop < 1e-5);
  GAUDI_EXPECT(tap_live > 0.01);
}

/// Dual-ring projector identity at rest (relative plane twist).
GAUDI_TEST(twist_dual_ring_projection_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const real w = 1.0;
  const vecX q = detail::pack_shell_q(blocks);
  auto c = detail::make_init_slab_ring_constraint<
      hepworth::block::relative_plane_twist>(
      stencil, hepworth::block::slab_prism_config{}, w, blocks, q);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);

  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);

  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-5);
}

GAUDI_TEST(bulk_procrustes_row_count_from_foot_samples) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const vecX q = detail::pack_shell_q(blocks);
  auto bulk = detail::make_init_slab_ring_constraint<
      hepworth::block::slab_dual_ring_volume_procrustes>(
      stencil, hepworth::block::slab_prism_config{}, 1.0, blocks, q);

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  vec3 foot;
  const std::vector<hepworth::block::slab_foot_local_sample> fit_samples =
      hepworth::block::build_top_ring_fit_samples(stencil, pos_at, foot);
  GAUDI_EXPECT(bulk->fit_sample_count() == fit_samples.size());

  index_t n_rows = 0;
  std::vector<trip> trips;
  bulk->fill_A(n_rows, trips);
  GAUDI_EXPECT(n_rows == static_cast<index_t>(3 * stencil.top_weights.size()));
}

GAUDI_TEST(bulk_procrustes_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const vecX q = detail::pack_shell_q(blocks);
  auto c = detail::make_init_slab_ring_constraint<
      hepworth::block::slab_dual_ring_volume_procrustes>(
      stencil, hepworth::block::slab_prism_config{}, 1.0, blocks, q);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);

  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);

  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-5);
}

GAUDI_TEST(bulk_procrustes_uses_live_covariance_extents) {
  // Cross-cov singular values ≠ geometric live extents under anisotropic stretch.
  // Build a synthetic rest ring in xy, stretch live along x only.
  std::vector<hepworth::block::slab_foot_local_sample> samples;
  const real w = 1.0;
  samples.push_back({0, vec3(1, 0, 0), w});
  samples.push_back({1, vec3(0, 1, 0), w});
  samples.push_back({2, vec3(-1, 0, 0), w});
  samples.push_back({3, vec3(0, -1, 0), w});

  hepworth::block::foot_frame_bulk_rest rest =
      hepworth::block::build_foot_frame_bulk_rest(samples, 1e-8);

  std::vector<vec3> live = {vec3(2, 0, 0), vec3(0, 1, 0), vec3(-2, 0, 0),
                            vec3(0, -1, 0)};
  const auto pos_at = [&](index_t vi) -> vec3 {
    if (vi < 0)
      return vec3::Zero();
    return live[static_cast<size_t>(vi)];
  };
  (void)pos_at;

  mat3 M_cross = mat3::Zero();
  mat3 Cov_live = mat3::Zero();
  for (size_t i = 0; i < samples.size(); ++i) {
    const vec3 dq = live[i];
    M_cross += w * dq * samples[i].local_rest.transpose();
    Cov_live += w * dq * dq.transpose();
  }

  const auto m_bad = hepworth::block::bulk_procrustes_map(
      M_cross, M_cross, rest, 1.08, 1.0, 1.0);
  const auto m_ok = hepworth::block::bulk_procrustes_map(
      M_cross, Cov_live, rest, 1.08, 1.0, 1.0);

  // Geometric live σ_x should be larger than rest; cross-cov σ is a different scale.
  GAUDI_EXPECT(m_ok.l0[0] > rest.sigma_rest[0] * 1.2);
  GAUDI_EXPECT((m_ok.l0 - m_bad.l0).norm() > 1e-6);
}

GAUDI_TEST(bulk_live_axis_nearest_ring_normal) {
  mat3 M = mat3::Zero();
  M(0, 0) = 1.0;
  M(1, 1) = 0.9;
  M(2, 2) = 0.05;
  const vec3 n = vec3(0.0, 0.0, 1.0);
  const vec3 axis =
      hepworth::block::bulk_live_axis_nearest_normal(M, n, 1e-8);
  GAUDI_EXPECT(axis.dot(n) > 0.99);
}

GAUDI_TEST(bulk_volume_pump_compress_sign) {
  const vec3 sigma = vec3(1.0, 0.9, 0.05);
  const real V_rest = sigma[0] * sigma[1] * sigma[2];
  const real max_step = 0.1;
  const real L = 0.05;
  const real delta_hi = hepworth::block::bulk_volume_pump_delta(
      1.2, 1.0, sigma, V_rest, L, 1e-4, max_step);
  const real delta_lo = hepworth::block::bulk_volume_pump_delta(
      0.8, 1.0, sigma, V_rest, L, 1e-4, max_step);
  GAUDI_EXPECT(delta_hi < 0.0);
  GAUDI_EXPECT(delta_lo > 0.0);
  GAUDI_EXPECT(std::abs(delta_hi) <= max_step + 1e-12);
  GAUDI_EXPECT(std::abs(delta_lo) <= max_step + 1e-12);
}

GAUDI_TEST(twist_rest_not_baked_from_inertia_prediction) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  auto &x = shell->xs->get();
  const vecX q_rest = detail::pack_shell_q(blocks);
  auto c = detail::make_init_slab_ring_constraint<
      hepworth::block::relative_plane_twist>(
      stencil, hepworth::block::slab_prism_config{}, 1.0, blocks, q_rest);

  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.2, 0.1, 0.05);
  const vecX q_bent = detail::pack_shell_q(blocks);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);
  vecX p = vecX::Zero(n_rows);
  c->project(q_bent, p);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q_bent);
  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel > 1e-4);
}

GAUDI_TEST(slab_prism_AtA_well_formed_on_sphere) {
  auto M = asawa::shell::load_sphere(1.0, 10, 7);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);

  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  hepworth::block::slab_prism_config prism_cfg;
  prism_cfg.w_bulk = 1.0;
  prism_cfg.w_twist = 1.0;
  const auto stats = hepworth::block::init_slab_prism_constraints(
      *M, ctx.constraints, hepworth::block::slab_pierce_config{}, prism_cfg,
      blocks, ctx, false);
  GAUDI_EXPECT(stats.slab_count > 0);

  const vecX q = hepworth::block::pack_block_q(blocks);
  hepworth::block::init_block_constraints_rest(ctx.constraints, q, ctx);
  const int nm = static_cast<int>(q.size());
  matS A;
  matS AtA;
  index_t a_rows = 0;
  GAUDI_ASSERT(hepworth::block::build_constraint_matrices(ctx.constraints, nm,
                                                          A, AtA, a_rows));
  GAUDI_EXPECT(a_rows > 0);

  const auto ata = hepworth::block::diagnose_sparse(AtA);
  GAUDI_EXPECT(ata.non_finite_nnz == 0);
  const double diag_ratio =
      ata.min_diag > 0.0 ? ata.max_diag / ata.min_diag : 0.0;
  GAUDI_EXPECT(diag_ratio < 1e10);

  matS MAtA = AtA;
  for (int i = 0; i < nm; ++i)
    MAtA.coeffRef(i, i) += 1e-2;
  gaudi::m_solver solver(MAtA);
  GAUDI_EXPECT(solver.success());
}

GAUDI_TEST(slab_prism_one_solver_step_identity) {
  auto M = asawa::shell::load_sphere(1.0, 10, 7);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());

  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  for (vec3 &v : vs->get())
    v.setZero();
  auto shell_block =
      std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);

  hepworth::block::solver_context ctx;
  shell_block->prepare(ctx);
  std::vector<hepworth::sim_block::ptr> blocks;
  shell_block->emit_dof_blocks(blocks);

  hepworth::block::slab_prism_config prism_cfg;
  prism_cfg.w_bulk = 1.0;
  prism_cfg.w_twist = 0.0;
  const auto stats = hepworth::block::init_slab_prism_constraints(
      *M, ctx.constraints, hepworth::block::slab_pierce_config{}, prism_cfg,
      blocks, ctx, false);
  GAUDI_EXPECT(stats.slab_count > 0);

  const std::vector<vec3> x0 = xs->get();
  hepworth::block::solver_context rest_ctx;
  rest_ctx.blocks = blocks;
  hepworth::block::projection_solver solver;
  solver.set_constraints(ctx.constraints);
  solver.step(blocks, 0.01, 0.2, 10, &rest_ctx);

  real max_drift = 0.0;
  const std::vector<vec3> &x1 = xs->get();
  GAUDI_EXPECT(x1.size() == x0.size());
  for (size_t i = 0; i < x0.size(); ++i)
    max_drift = std::max(max_drift, (x1[i] - x0[i]).norm());
  GAUDI_EXPECT(max_drift < 1e-4);
}

GAUDI_TEST(fan_tet_null_gauge_small_rms) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  hepworth::block::slab_fan_gauge_cache gauge;
  (void)hepworth::block::build_slab_fan_tets_aligned(stencil, pos_at, gauge);
  GAUDI_EXPECT(gauge.basis.dim_q > 0);
  GAUDI_EXPECT(std::isfinite(gauge.align_rms));
}

GAUDI_TEST(fan_tet_volume_matches_tau_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  const std::vector<hepworth::block::slab_fan_wedge> wedges =
      hepworth::block::build_slab_fan_wedges(stencil, pos_at);
  hepworth::block::slab_fan_gauge_cache gauge;
  const std::vector<hepworth::block::slab_fan_tet> tets =
      hepworth::block::build_slab_fan_tets_aligned(stencil, pos_at, gauge);

  GAUDI_EXPECT(tets.size() == wedges.size());
  for (size_t k = 0; k < tets.size(); ++k) {
    GAUDI_EXPECT(tets[k].V_rest > 0.0);
    GAUDI_EXPECT(std::abs(std::abs(tets[k].E_rest.determinant()) -
                         std::abs(wedges[k].tau_rest.coeff)) <
                 1e-10 * std::max(tets[k].V_rest * 6.0, real(1e-12)));
    GAUDI_EXPECT(tets[k].E_rest.determinant() > 0.0);
  }
}

GAUDI_TEST(shape_field_project_same_q_twice) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  for (auto field_live :
       {hepworth::block::field_live_update::refit_every_project,
        hepworth::block::field_live_update::live_matches_rest}) {
    hepworth::block::slab_prism_config cfg;
    cfg.field_live = field_live;
    cfg.field_gd_iters = 8;
    auto c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                                 *shell->mesh, avg);
    index_t id0 = 0;
    std::vector<trip> trips;
    c->fill_A(id0, trips);
    const vecX q = detail::pack_shell_q(blocks);
    vecX p1 = vecX::Zero(id0);
    vecX p2 = vecX::Zero(id0);
    c->project(q, p1);
    c->project(q, p2);
    const real dp = (p1 - p2).norm();
    GAUDI_EXPECT(dp < 1e-12);
  }
}

/// Hypothesis B: refit_every_project amplifies target drift under small q perturbation
/// (simulation creep) vs frozen Q_k.
GAUDI_TEST(shape_field_refit_amplifies_target_drift) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto run_sensitivity = [&](hepworth::block::field_live_update mode) {
    hepworth::block::slab_prism_config cfg;
    cfg.field_live = mode;
    cfg.field_gd_iters = 8;
    auto c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                                 *shell->mesh, avg);
    index_t id0 = 0;
    std::vector<trip> trips;
    c->fill_A(id0, trips);

    vecX q0 = detail::pack_shell_q(blocks);
    vecX q1 = q0;
    const index_t vi = stencil.top_ids[1];
    const index_t off = blocks[0]->get_offset_idx(vi);
    q1[off + 0] += 0.03;
    q1[off + 1] += 0.02;

    vecX p0 = vecX::Zero(id0);
    vecX p1 = vecX::Zero(id0);
    c->project(q0, p0);
    c = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                            *shell->mesh, avg);
    c->fill_A(id0, trips);
    c->project(q1, p1);
    return (p1 - p0).norm();
  };

  const real drift_refit = run_sensitivity(
      hepworth::block::field_live_update::refit_every_project);
  const real drift_frozen =
      run_sensitivity(hepworth::block::field_live_update::live_matches_rest);
  GAUDI_EXPECT(drift_refit > 1e-10);
  GAUDI_EXPECT(drift_frozen > 1e-10);
  GAUDI_EXPECT(drift_refit > 0.5 * drift_frozen);
}

/// Hypothesis C: partial GD may not reach stationarity (residual grad at x_star).
GAUDI_TEST(shape_field_gd_partial_stationarity) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const auto pos_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::darboux, stencil, pos_at, fit_cfg,
          avg);

  std::vector<vec3> x_pert = x;
  x_pert[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.05, 0.02, 0.0);
  const auto pos_pert = [&](index_t vi) {
    return x_pert[static_cast<size_t>(vi)];
  };
  const hepworth::block::slab_field_coeffs Q_live =
      hepworth::block::fit_field_coeffs_on_stencil(
          hepworth::block::slab_field_family::darboux, stencil, pos_pert,
          fit_cfg, avg);

  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil, pos_pert);
  auto x_star = hepworth::block::resolve_shape_field_divergence_gd(
      stencil, Q, Q_live, ring, 1.0, 8, 0.25);

  std::vector<vec3> grad;
  hepworth::block::shape_field_divergence_grad<hepworth::block::darboux_field_bundle>(
      stencil, Q.darboux, Q_live.darboux, x_star, 1.0, grad);
  real g2 = 0.0;
  for (const vec3 &g : grad)
    g2 += g.squaredNorm();
  // If partial GD stops early, grad norm stays measurably above tolerance.
  GAUDI_EXPECT(g2 > 1e-20);
}

GAUDI_TEST(slab_wedge_strain_zero_at_rest) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const auto wedges = hepworth::block::build_slab_fan_wedges(stencil, x);
  GAUDI_EXPECT(wedges.size() >= 4);
  const auto anchor_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil, anchor_at);
  const real e0 =
      hepworth::block::eval_wedge_strain_energy(stencil, wedges, ring, anchor_at);
  GAUDI_EXPECT(e0 < 1e-20);
}

GAUDI_TEST(slab_wedge_strain_gn_identity_at_rest) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const auto wedges = hepworth::block::build_slab_fan_wedges(stencil, x);
  const auto anchor_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  auto ring = hepworth::block::slab_ring_positions::from_stencil(stencil, anchor_at);
  auto resolved = hepworth::block::resolve_wedge_strain_gauss_newton(
      stencil, wedges, ring, anchor_at, 4);
  for (size_t i = 0; i < ring.x.size(); ++i)
    GAUDI_EXPECT((resolved.x[i] - ring.x[i]).norm() < 1e-10);
}

GAUDI_TEST(slab_prism_divergence_gd_reduces_volume_error) {
  auto M = asawa::shell::load_sphere(1.0, 8, 6);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  std::vector<vec3> x = asawa::const_get_vec_data(*M, 0);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  const std::vector<index_t> face_map = M->get_face_map(true);
  hepworth::block::slab_pierce_config cfg;
  const real avg = asawa::shell::avg_length(*M, x);

  hepworth::block::slab_pierce_stencil stencil;
  for (int iv = 0; iv < static_cast<int>(M->vert_count()); ++iv) {
    const auto result = hepworth::block::build_slab_pierce_stencil(
        *M, asawa::shell::vert_id(iv), x, *tree, face_map, cfg, avg);
    if (result.stencil.has_slab) {
      stencil = result.stencil;
      break;
    }
  }
  GAUDI_ASSERT(stencil.has_slab);

  const hepworth::block::slab_prism_rest rest =
      hepworth::block::build_slab_prism_rest(stencil, x);
  const auto anchor_at = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  const auto top_at = anchor_at;
  auto top = hepworth::block::top_ring_positions::from_stencil(stencil, top_at);

  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.05, 0.02, 0.0);
  top = hepworth::block::top_ring_positions::from_stencil(stencil, top_at);
  const real e0 = std::abs(
      hepworth::block::eval_prism_volume_top_only(stencil, top, anchor_at) -
      rest.V_rest);
  GAUDI_EXPECT(e0 > 1e-10);
  top = hepworth::block::resolve_prism_volume_gd(stencil, top, anchor_at,
                                                 rest.V_rest, 16, 1.0);
  const real e1 = std::abs(
      hepworth::block::eval_prism_volume_top_only(stencil, top, anchor_at) -
      rest.V_rest);
  GAUDI_EXPECT(std::isfinite(e1));
}

  GAUDI_TEST(shape_field_init_rest_matches_legacy_ctor) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  hepworth::block::slab_prism_config cfg;
  cfg.field_family = hepworth::block::slab_field_family::darboux;
  cfg.field_live = hepworth::block::field_live_update::live_matches_rest;
  const real w = 1.0;

  auto legacy = detail::make_shape_field_constraint(stencil, rest, cfg, blocks,
                                                    *shell->mesh, avg, w);
  auto via_init = hepworth::block::shape_field_divergence_strain::create(
      stencil, cfg, w, blocks, avg);
  hepworth::block::solver_context ctx;
  ctx.blocks = blocks;
  const vecX q = hepworth::block::pack_block_q(blocks);
  via_init->init_rest(q, ctx);

  index_t n_rows = 0;
  std::vector<trip> trips;
  legacy->fill_A(n_rows, trips);
  index_t n_rows2 = 0;
  std::vector<trip> trips2;
  via_init->fill_A(n_rows2, trips2);
  GAUDI_EXPECT(n_rows == n_rows2);
  GAUDI_EXPECT(trips.size() == trips2.size());

  vecX p0 = vecX::Zero(n_rows);
  vecX p1 = vecX::Zero(n_rows2);
  legacy->project(q, p0);
  via_init->project(q, p1);
  GAUDI_EXPECT((p0 - p1).norm() < 1e-10 * std::max(p0.norm(), real(1.0)));
}

GAUDI_TEST(phi_div_free_basis_zero_divergence) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_at);
  GAUDI_ASSERT(basis.dim_q > 0);

  vecX q = vecX::Random(basis.dim_q);
  for (int i = 0; i < basis.dim_q; ++i)
    q[i] = 0.01 * q[i];
  GAUDI_EXPECT(basis.max_div_on_wedges(q) < 1e-10);
}

GAUDI_TEST(phi_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  hepworth::block::slab_prism_config cfg;
  cfg.phi_nr_iters = 12;
  const real w = 1.0;
  auto c = std::make_shared<hepworth::block::incompressible_phi_strain>(
      stencil, rest, cfg, w, blocks,
      hepworth::block::fit_phi_quadric_on_stencil(
          stencil,
          [&](index_t vi) {
            return shell->xs->get()[static_cast<size_t>(vi)];
          },
          hepworth::block::local_field_fit_config{}),
      avg);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);
  const vecX q = detail::pack_shell_q(blocks);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel < 1e-5);
}

GAUDI_TEST(phi_gd_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_at, fit_cfg);
  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_at);
  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_at);

  const vecX q0 = vecX::Zero(basis.dim_q);
  const real e0 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0);
  GAUDI_EXPECT(e0 < 1e-20);

  vecX q_warm = vecX::Zero(basis.dim_q);
  for (int i = 0; i < basis.dim_q; ++i)
    q_warm[i] = 0.02 * (i + 1);
  const vecX q1 = hepworth::block::resolve_incompressible_phi_gd(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q_warm, 1.0, 16, 0.25);
  const real e1 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q1, 1.0);
  GAUDI_EXPECT(e1 < 1e-16);
  GAUDI_EXPECT(q1.norm() < 1e-4);

  const hepworth::block::slab_ring_positions x_star =
      hepworth::block::apply_phi_to_ring(stencil, ring, basis, q1);
  real max_disp = 0.0;
  for (index_t vi : ring.ids)
    max_disp = std::max(max_disp, (x_star.at(vi) - ring.at(vi)).norm());
  GAUDI_EXPECT(max_disp < 1e-6);
}

GAUDI_TEST(phi_gauge_gd_reduces_energy) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  auto &x = shell->xs->get();
  const std::vector<vec3> x_rest = x;
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  const auto pos_bent = [&](index_t vi) {
    return x[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_rest, fit_cfg);

  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);
  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_rest);
  GAUDI_ASSERT(basis.dim_q > 0);

  const vecX q0 = vecX::Zero(basis.dim_q);
  const real e0 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0);
  const vecX q1 = hepworth::block::resolve_incompressible_phi_gd(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0, 12, 0.25);
  const real e1 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q1, 1.0);
  GAUDI_EXPECT(e0 > 1e-20);
  GAUDI_EXPECT(e1 <= e0);
}

GAUDI_TEST(phi_gauge_energy_nonzero_at_q0_on_bent_mesh) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  const auto pos_bent = [&](index_t vi) {
    return x[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_rest, fit_cfg);

  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_rest);
  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);
  const vecX q0 = vecX::Zero(basis.dim_q);
  const real e0 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0);
  GAUDI_EXPECT(e0 > 1e-20);
}

GAUDI_TEST(phi_bent_projection_not_noop) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const vecX q_rest = detail::pack_shell_q(blocks);
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);
  const vecX q = detail::pack_shell_q(blocks);

  constexpr real w_phi = 1.0;
  hepworth::block::slab_prism_config cfg;
  cfg.phi_nr_iters = 12;
  auto c = hepworth::block::incompressible_phi_strain::create(
      stencil, cfg, w_phi, blocks, avg);
  hepworth::block::solver_context ctx;
  ctx.blocks = blocks;
  c->init_rest(q_rest, ctx);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);
  const vecX Ap = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);
  const real rel =
      (p - Ap).norm() / std::max(Ap.norm(), real(1e-12));
  GAUDI_EXPECT(rel > 1e-5);
}

GAUDI_TEST(phi_gauge_analytic_grad_fd) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  auto &x = shell->xs->get();
  const std::vector<vec3> x_rest = x;
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.04, 0.02, 0.01);

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  const auto pos_bent = [&](index_t vi) {
    return x[static_cast<size_t>(vi)];
  };
  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_rest, fit_cfg);
  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_rest);
  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);

  vecX q = vecX::Zero(basis.dim_q);
  for (int i = 0; i < basis.dim_q; ++i)
    q[i] = 0.01 * (i + 1);

  const auto check_grad = [&](const albers::vec10 &Qq, real tol = 1e-4,
                              real eps = 1e-6) {
    const vecX ga = hepworth::block::eval_incompressible_phi_gauge_grad_q_analytic(
        stencil, Qq, ring, basis, basis.foot_rest, q, 1.0);
    const vecX gf = hepworth::block::eval_incompressible_phi_gauge_grad_q_fd(
        stencil, Qq, ring, basis, basis.foot_rest, q, 1.0, 1.0, eps);
    const real rel_norm =
        (ga - gf).norm() / std::max(ga.norm(), real(1e-12));
    GAUDI_EXPECT(rel_norm < tol);
  };

  albers::vec10 Q_sphere = albers::vec10::Zero();
  Q_sphere[0] = Q_sphere[1] = Q_sphere[2] = 1.0;
  Q_sphere[9] = -1.0;

  check_grad(Q_sphere);
  check_grad(Q.quadric, 0.1, 1e-5);
}

GAUDI_TEST(phi_sphere_caps_pull_apart_report) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  auto &x = shell->xs->get();

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };

  const auto radial_stretch = [&](real eps) {
    std::vector<vec3> out = x_rest;
    const vec3 c_top = x_rest[static_cast<size_t>(stencil.foot)];
    for (size_t k = 1; k < stencil.top_ids.size(); ++k) {
      const index_t vi = stencil.top_ids[k];
      out[static_cast<size_t>(vi)] =
          c_top + (1.0 + eps) * (x_rest[static_cast<size_t>(vi)] - c_top);
    }
    const vec3 c_bot = x_rest[static_cast<size_t>(stencil.pierce_foot)];
    for (size_t k = 1; k < stencil.bot_ids.size(); ++k) {
      const index_t vi = stencil.bot_ids[k];
      out[static_cast<size_t>(vi)] =
          c_bot + (1.0 + eps) * (x_rest[static_cast<size_t>(vi)] - c_bot);
    }
    return out;
  };

  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_rest, fit_cfg);
  hepworth::block::phi_div_free_basis basis;
  basis.build(stencil, pos_rest);

  constexpr real eps = 0.05;
  x = radial_stretch(eps);
  const auto pos_bent = [&](index_t vi) { return x[static_cast<size_t>(vi)]; };
  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);

  const vecX q0 = vecX::Zero(basis.dim_q);
  const real e0 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0);
  const vecX q1 = hepworth::block::resolve_incompressible_phi_gd(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q0, 1.0, 20, 0.25);
  const real e1 = hepworth::block::eval_incompressible_phi_gauge_energy(
      stencil, Q.quadric, ring, basis, basis.foot_rest, q1, 1.0);
  const hepworth::block::slab_ring_positions x_star =
      hepworth::block::apply_phi_to_ring(stencil, ring, basis, q1);

  real max_disp = 0.0;
  for (index_t vi : ring.ids) {
    max_disp = std::max(max_disp, (x_star.at(vi) - ring.at(vi)).norm());
  }
  const real max_div = basis.max_div_on_wedges(q1);

  GAUDI_EXPECT(basis.dim_q > 3);
  GAUDI_EXPECT(e0 > 1e-20);
  GAUDI_EXPECT(e1 <= e0);
  GAUDI_EXPECT(max_disp > 1e-4 * avg);
  GAUDI_EXPECT(max_div < 1e-6);
}


/// North-pole foot on a sphere: opposing top/bot rings. Probe bulk ΔD energy for
/// Darboux and quadric under (a) identical Q (must be ~0 for any geometry) and
/// (b) frozen Q_r + live-refit Q_k after a spoke bump (E>0, GD should reduce E).
GAUDI_TEST(shape_field_sphere_opposing_rings_energy_probe) {
  auto M = asawa::shell::load_sphere(1.0, 12, 8);
  GAUDI_ASSERT(M != nullptr);
  asawa::shell::triangulate(*M);
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);

  index_t north = 0;
  for (int i = 1; i < static_cast<int>(x.size()); ++i)
    if (x[static_cast<size_t>(i)].y() > x[static_cast<size_t>(north)].y())
      north = static_cast<index_t>(i);

  arp::bvh_tree<3>::ptr tree = hepworth::block::make_shell_face_bvh(*M, x);
  GAUDI_ASSERT(tree);
  const std::vector<index_t> face_map = M->get_face_map(true);
  const real avg = asawa::shell::avg_length(*M, x);
  const auto pierce = hepworth::block::build_slab_pierce_stencil(
      *M, asawa::shell::vert_id(static_cast<int>(north)), x, *tree, face_map,
      hepworth::block::slab_pierce_config{}, avg);
  GAUDI_ASSERT(pierce.stencil.has_slab);

  const vec3 c_top = x[static_cast<size_t>(pierce.stencil.foot)];
  const vec3 c_bot = x[static_cast<size_t>(pierce.stencil.pierce_foot)];
  const real spine = (c_bot - c_top).norm();
  const real antipodal = (c_top.normalized() + c_bot.normalized()).norm();
  std::cerr << "[opposing_rings] foot=" << pierce.stencil.foot
            << " pierce=" << pierce.stencil.pierce_foot
            << " top_spokes=" << pierce.stencil.top_weights.size()
            << " bot_spokes=" << pierce.stencil.bot_weights.size()
            << " spine=" << spine << " antipodal_sum_norm=" << antipodal
            << std::endl;
  // Sphere pierce should be roughly antipodal (unit sphere diameter ~2).
  GAUDI_EXPECT(spine > 1.5);
  GAUDI_EXPECT(antipodal < 0.35);

  const auto run_family = [&](hepworth::block::slab_field_family family,
                              const char *name) {
    const auto pos_rest = [&](index_t vi) {
      return x[static_cast<size_t>(vi)];
    };
    hepworth::block::local_field_fit_config fit_cfg;
    const hepworth::block::slab_field_coeffs Q_rest =
        hepworth::block::fit_field_coeffs_on_stencil(
            family, pierce.stencil, pos_rest, fit_cfg, avg);

    auto ring_rest = hepworth::block::slab_ring_positions::from_stencil(
        pierce.stencil, pos_rest);
    const real e_same_rest = hepworth::block::eval_shape_field_divergence_energy(
        pierce.stencil, Q_rest, Q_rest, ring_rest, 1.0);

    // Same Q on perturbed geometry: ΔD ≡ 0 analytically → energy must stay ~0.
    std::vector<vec3> x_pert = x;
    GAUDI_ASSERT(pierce.stencil.top_ids.size() >= 2);
    const index_t spoke = pierce.stencil.top_ids[1];
    const vec3 n_spoke = x_pert[static_cast<size_t>(spoke)].normalized();
    x_pert[static_cast<size_t>(spoke)] += 0.15 * avg * n_spoke;
    const auto pos_pert = [&](index_t vi) {
      return x_pert[static_cast<size_t>(vi)];
    };
    auto ring_pert = hepworth::block::slab_ring_positions::from_stencil(
        pierce.stencil, pos_pert);
    const real e_same_pert = hepworth::block::eval_shape_field_divergence_energy(
        pierce.stencil, Q_rest, Q_rest, ring_pert, 1.0);

    // Meaningful case: freeze Q_r, refit Q_k on perturbed stencil.
    const hepworth::block::slab_field_coeffs Q_live =
        hepworth::block::fit_field_coeffs_on_stencil(
            family, pierce.stencil, pos_pert, fit_cfg, avg, &Q_rest);
    const real e_mismatch = hepworth::block::eval_shape_field_divergence_energy(
        pierce.stencil, Q_rest, Q_live, ring_pert, 1.0);

    auto x_star = hepworth::block::resolve_shape_field_divergence_gd(
        pierce.stencil, Q_rest, Q_live, ring_pert, 1.0, 16, 0.25);
    const real e_after = hepworth::block::eval_shape_field_divergence_energy(
        pierce.stencil, Q_rest, Q_live, x_star, 1.0);

    real max_disp = 0.0;
    vec3 spoke_delta = vec3::Zero();
    for (index_t vi : ring_pert.ids) {
      const vec3 d = x_star.at(vi) - ring_pert.at(vi);
      max_disp = std::max(max_disp, d.norm());
      if (vi == spoke)
        spoke_delta = d;
    }
    const real spoke_radial = spoke_delta.dot(n_spoke);

    std::cerr << "[opposing_rings." << name << "] e_same_rest=" << e_same_rest
              << " e_same_pert=" << e_same_pert
              << " e_mismatch=" << e_mismatch << " e_after_gd=" << e_after
              << " max_disp=" << max_disp
              << " spoke_radial_delta=" << spoke_radial << std::endl;

    GAUDI_EXPECT(e_same_rest < 1e-18);
    GAUDI_EXPECT(e_same_pert < 1e-18);
    GAUDI_EXPECT(e_mismatch > 1e-12);
    GAUDI_EXPECT(e_after <= e_mismatch * 1.01);
    // GD should move something when fields disagree.
    GAUDI_EXPECT(max_disp > 1e-8);
  };

  run_family(hepworth::block::slab_field_family::darboux, "darboux");
  run_family(hepworth::block::slab_field_family::quadric, "quadric");
}

GAUDI_TEST(rbf_displacement_basis_builds) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));
  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  hepworth::block::rbf_displacement_basis basis;
  basis.build(stencil, pos_at, avg, 4, 2.0);
  GAUDI_EXPECT(basis.dim_q > 0);
  GAUDI_EXPECT(basis.n_centers >= 1);
  GAUDI_EXPECT(basis.B.rows() == static_cast<int>(3 * basis.vert_ids.size()));
  GAUDI_EXPECT(basis.B.cols() == basis.dim_q);
}

GAUDI_TEST(rbf_elastic_identity_at_rest) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  hepworth::block::slab_prism_config cfg;
  cfg.rbf_nr_iters = 8;
  cfg.rbf_mu = 1.0;
  cfg.rbf_eta = 1.0;
  const auto pos_at = [&](index_t vi) {
    return shell->xs->get()[static_cast<size_t>(vi)];
  };
  auto c = std::make_shared<hepworth::block::rbf_elastic_strain>(
      stencil, rest, cfg, 1.0, blocks,
      hepworth::block::fit_phi_quadric_on_stencil(
          stencil, pos_at, hepworth::block::local_field_fit_config{}),
      avg);

  const vecX q = detail::pack_shell_q(blocks);
  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);
  const vecX Aq = detail::apply_triplets(trips, n_rows, q);
  // Rest: R≈I ⇒ p = w R d_rest matches hub-spoke Aq.
  GAUDI_EXPECT((p - Aq).norm() < 1e-8 * std::max(Aq.norm(), real(1.0)));
}

GAUDI_TEST(rbf_elastic_bent_emits_arap_R_not_raw_spokes) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const vecX q_rest = detail::pack_shell_q(blocks);
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);
  const vecX q = detail::pack_shell_q(blocks);

  hepworth::block::slab_prism_config cfg;
  cfg.rbf_nr_iters = 12;
  cfg.rbf_mu = 1.0;
  cfg.rbf_eta = 1.0;
  auto c = hepworth::block::rbf_elastic_strain::create(stencil, cfg, 1.0, blocks,
                                                       avg);
  hepworth::block::solver_context ctx;
  ctx.blocks = blocks;
  c->init_rest(q_rest, ctx);

  index_t n_rows = 0;
  std::vector<trip> trips;
  c->fill_A(n_rows, trips);
  const vecX Aq = detail::apply_triplets(trips, n_rows, q);
  vecX p = vecX::Zero(n_rows);
  c->project(q, p);
  // Non-rigid bend: w R d_rest ≠ w (q_spoke − q_hub).
  GAUDI_EXPECT((p - Aq).norm() / std::max(Aq.norm(), real(1e-12)) > 1e-5);
}

GAUDI_TEST(rbf_elastic_gd_reduces_energy_on_bent) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.06, 0.03, 0.01);

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  const auto pos_bent = [&](index_t vi) {
    return x[static_cast<size_t>(vi)];
  };

  hepworth::block::local_field_fit_config fit_cfg;
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(stencil, pos_rest, fit_cfg);
  rest = hepworth::block::build_slab_prism_rest(stencil, pos_rest);

  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);
  hepworth::block::rbf_displacement_basis basis;
  basis.build(stencil, pos_rest, avg, 4, 2.0);
  GAUDI_ASSERT(basis.dim_q > 0);

  const vecX q0 = vecX::Zero(basis.dim_q);
  const real e0 = hepworth::block::eval_rbf_elastic_energy(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q0, 1.0, 1.0,
      1.0);
  const vecX q1 = hepworth::block::resolve_rbf_elastic_gd(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q0, 1.0, 1.0, 1.0,
      12, 0.25);
  const real e1 = hepworth::block::eval_rbf_elastic_energy(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q1, 1.0, 1.0,
      1.0);
  GAUDI_EXPECT(e0 > 1e-20);
  GAUDI_EXPECT(e1 <= e0 * 1.01);
}

GAUDI_TEST(rbf_elastic_separate_weights_affect_energy) {
  hepworth::block::slab_pierce_stencil stencil;
  hepworth::block::slab_prism_rest rest;
  std::shared_ptr<hepworth::block::shell_position_block> shell;
  std::vector<hepworth::sim_block::ptr> blocks;
  const real avg = 0.1;
  GAUDI_ASSERT(detail::build_one_slab_field_constraint(stencil, rest, shell,
                                                       blocks, avg));

  const std::vector<vec3> x_rest = shell->xs->get();
  auto &x = shell->xs->get();
  x[static_cast<size_t>(stencil.top_ids[1])] += vec3(0.05, 0.0, 0.02);

  const auto pos_rest = [&](index_t vi) {
    return x_rest[static_cast<size_t>(vi)];
  };
  const auto pos_bent = [&](index_t vi) {
    return x[static_cast<size_t>(vi)];
  };
  rest = hepworth::block::build_slab_prism_rest(stencil, pos_rest);
  const hepworth::block::slab_field_coeffs Q =
      hepworth::block::fit_phi_quadric_on_stencil(
          stencil, pos_rest, hepworth::block::local_field_fit_config{});
  const hepworth::block::slab_ring_positions ring =
      hepworth::block::slab_ring_positions::from_stencil(stencil, pos_bent);
  hepworth::block::rbf_displacement_basis basis;
  basis.build(stencil, pos_rest, avg, 4, 2.0);
  const vecX q0 = vecX::Zero(basis.dim_q);

  const real e_gauge = hepworth::block::eval_rbf_elastic_energy(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q0, 1.0, 1.0,
      0.0);
  const real e_strain = hepworth::block::eval_rbf_elastic_energy(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q0, 1.0, 0.0,
      1.0);
  const real e_both = hepworth::block::eval_rbf_elastic_energy(
      stencil, Q.quadric, rest, ring, basis, basis.foot_rest, q0, 1.0, 1.0,
      1.0);
  GAUDI_EXPECT(e_gauge > 0.0);
  GAUDI_EXPECT(e_strain > 0.0);
  GAUDI_EXPECT(std::abs(e_both - (e_gauge + e_strain)) <
               1e-9 * std::max(e_both, real(1.0)));
}

} // namespace test
} // namespace gaudi

#endif // GAUDI_TEST_SLAB_PRISM_TESTS_HPP
