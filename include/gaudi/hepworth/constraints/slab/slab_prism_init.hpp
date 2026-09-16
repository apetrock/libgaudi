#ifndef GAUDI_HEPWORTH_SLAB_PRISM_INIT_HPP
#define GAUDI_HEPWORTH_SLAB_PRISM_INIT_HPP

#include "incompressible_phi_strain.hpp"
#include "rbf_elastic_strain.hpp"
#include "slab_dual_ring_volume_procrustes.hpp"
#include "gaudi/hepworth/block/block_geometry.hpp"
#include "relative_plane_twist.hpp"
#include "shape_field_divergence_strain.hpp"
#include "slab_pierce_failure_dump.hpp"
#include "slab_pierce_stencil.hpp"
#include "slab_prism_config.hpp"

#include "gaudi/hepworth/block/projection_diagnostics.hpp"
#include "gaudi/hepworth/block/slab_prism_log.hpp"

#include <atomic>
#include "gaudi/hepworth/projection_constraint.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>

namespace gaudi {
namespace hepworth {
namespace block {

struct slab_prism_init_stats {
  size_t slab_count = 0;
  size_t pierce_fail_count = 0;
  size_t active_feet = 0;
  std::map<pierce_failure_reason, size_t> fail_by_reason;
  real pierce_spine_min = 0.0;
  real pierce_spine_max = 0.0;
  real pierce_spine_mean = 0.0;
};

inline slab_prism_init_stats init_slab_prism_constraints(
    const asawa::shell::shell &shell,
    std::vector<projection_constraint::ptr> &constraints,
    const arp::bvh_tree<3> &face_tree, const std::vector<index_t> &face_map,
    const slab_pierce_config &pierce_cfg, const slab_prism_config &prism_cfg,
    std::vector<sim_block::ptr> blocks, solver_context &ctx,
    bool dump_failures = true) {
  slab_prism_init_stats stats;
  real spine_len_sum = 0.0;
  stats.pierce_spine_min = std::numeric_limits<real>::infinity();

  using clock = std::chrono::steady_clock;
  const auto t0 = clock::now();

  if (slab_prism_log_enabled())
    slab_prism_log("init_slab_prism begin verts=" +
                   std::to_string(shell.vert_count()));

  if (blocks.empty())
    return stats;

  const vecX q_rest = pack_block_q(blocks);
  const auto pos_at = make_block_pos_at(blocks[0], q_rest);
  const std::vector<vec3> x_shell =
      gather_shell_positions(shell, blocks[0], q_rest);
  const real avg_edge_length = asawa::shell::avg_length(shell, x_shell);

  const std::string dump_path = pierce_failure_dump_path();
  if (dump_failures) {
    std::error_code ec;
    std::filesystem::create_directories(
        std::filesystem::path(dump_path).parent_path(), ec);
    std::ofstream clear(dump_path, std::ios::trunc);
  }

  if (projection_diag_enabled())
    projection_diag_log("init_slab_prism avg_edge_length=" +
                        std::to_string(avg_edge_length));

  const size_t n_before = constraints.size();

  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x_shell);
    if (A < 1e-8)
      continue;

    auto ring = shell.get_one_ring(vid);
    if (ring.empty())
      continue;
    ++stats.active_feet;

    const slab_pierce_result pierce = build_slab_pierce_stencil(
        shell, vid, x_shell, face_tree, face_map, pierce_cfg, avg_edge_length);

    if ((projection_diag_enabled() || slab_prism_log_enabled()) &&
        (iv + 1) % 5000 == 0) {
      const double ms =
          std::chrono::duration<double, std::milli>(clock::now() - t0).count();
      slab_prism_log("init progress " + std::to_string(iv + 1) + "/" +
                     std::to_string(shell.vert_count()) + " slab=" +
                     std::to_string(stats.slab_count) + " fail=" +
                     std::to_string(stats.pierce_fail_count) + " ms=" +
                     std::to_string(ms));
    }

    if (!pierce.stencil.has_slab) {
      ++stats.pierce_fail_count;
      ++stats.fail_by_reason[pierce.reason];
      if (dump_failures) {
        pierce_failure_record rec;
        rec.foot = static_cast<index_t>(vid);
        rec.reason = pierce.reason;
        rec.x = pos_at(static_cast<index_t>(vid));
        vec3 n = asawa::shell::vert_normal(shell, vid, x_shell);
        if (n.squaredNorm() > 1e-24)
          rec.n_in = n.normalized();
        rec.area = A;
        rec.valence = ring.size();
        rec.hit_face = pierce.hit_face;
        rec.ray_t = pierce.ray_t;
        rec.miss_detail = pierce.miss_detail;
        append_pierce_failure_record(rec, dump_path);
      }
      continue;
    }

    const real spine_len =
        (pos_at(pierce.stencil.pierce_foot) - pos_at(pierce.stencil.foot)).norm();
    stats.pierce_spine_min = std::min(stats.pierce_spine_min, spine_len);
    stats.pierce_spine_max = std::max(stats.pierce_spine_max, spine_len);
    spine_len_sum += spine_len;

    if (prism_cfg.w_field > 0.0) {
      constraints.push_back(shape_field_divergence_strain::create(
          pierce.stencil, prism_cfg, prism_cfg.w_field, blocks,
          avg_edge_length));
    }
    if (prism_cfg.w_phi > 0.0) {
      constraints.push_back(incompressible_phi_strain::create(
          pierce.stencil, prism_cfg, prism_cfg.w_phi, blocks,
          avg_edge_length));
    }
    if (prism_cfg.w_rbf > 0.0) {
      constraints.push_back(rbf_elastic_strain::create(
          pierce.stencil, prism_cfg, prism_cfg.w_rbf, blocks,
          avg_edge_length));
    }
    if (prism_cfg.w_bulk > 0.0) {
      constraints.push_back(slab_dual_ring_volume_procrustes::create(
          pierce.stencil, prism_cfg, prism_cfg.w_bulk, blocks));
    }
    if (prism_cfg.w_twist > 0.0) {
      constraints.push_back(relative_plane_twist::create(
          pierce.stencil, prism_cfg, prism_cfg.w_twist, blocks));
    }
    ++stats.slab_count;
  }

  init_block_constraints_rest(
      std::vector<projection_constraint::ptr>(constraints.begin() +
                                              static_cast<std::ptrdiff_t>(n_before),
                                              constraints.end()),
      q_rest, ctx);

  if (stats.slab_count > 0) {
    stats.pierce_spine_mean = spine_len_sum / static_cast<real>(stats.slab_count);
    if (!std::isfinite(stats.pierce_spine_min))
      stats.pierce_spine_min = 0.0;
  }

  {
    const double ms =
        std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    std::string reason_line;
    for (const auto &[reason, count] : stats.fail_by_reason) {
      reason_line += " " + std::string(pierce_failure_reason_name(reason)) +
                     "=" + std::to_string(count);
    }
    static std::atomic<bool> init_traced{false};
    if (!init_traced.exchange(true)) {
      slab_prism_trace("init_slab_prism done slab=" +
                       std::to_string(stats.slab_count) + " active_feet=" +
                       std::to_string(stats.active_feet) + " pierce_fail=" +
                       std::to_string(stats.pierce_fail_count) + reason_line +
                       " ms=" + std::to_string(ms));
    }
  }

  if (projection_diag_enabled() || slab_prism_log_enabled()) {
    const double ms =
        std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    std::string reason_line;
    for (const auto &[reason, count] : stats.fail_by_reason) {
      reason_line += " " + std::string(pierce_failure_reason_name(reason)) +
                     "=" + std::to_string(count);
    }
    slab_prism_log("init done slab=" + std::to_string(stats.slab_count) +
                   " active_feet=" + std::to_string(stats.active_feet) +
                   " pierce_fail=" + std::to_string(stats.pierce_fail_count) +
                   reason_line + " ms=" + std::to_string(ms));
    if (projection_diag_enabled()) {
      projection_diag_log("init_slab_prism done slab=" +
                          std::to_string(stats.slab_count) + " active_feet=" +
                          std::to_string(stats.active_feet) + " pierce_fail=" +
                          std::to_string(stats.pierce_fail_count) + reason_line +
                          " ms=" + std::to_string(ms));
    }
    if (dump_failures && stats.pierce_fail_count > 0)
      slab_prism_log("pierce failure dump: " + dump_path);
  }

  return stats;
}

inline slab_prism_init_stats
audit_slab_pierce_stencils(const asawa::shell::shell &shell,
                           const std::vector<vec3> &x,
                           const arp::bvh_tree<3> &face_tree,
                           const std::vector<index_t> &face_map,
                           const slab_pierce_config &cfg, real avg_edge_length) {
  slab_prism_init_stats stats;
  for (int iv = 0; iv < static_cast<int>(shell.vert_count()); ++iv) {
    const asawa::shell::VertId vid = asawa::shell::vert_id(iv);
    const real A = asawa::shell::vert_area(shell, vid, x);
    if (A < 1e-8)
      continue;
    if (shell.get_one_ring(vid).empty())
      continue;
    ++stats.active_feet;
    const slab_pierce_result pierce = build_slab_pierce_stencil(
        shell, vid, x, face_tree, face_map, cfg, avg_edge_length);
    if (pierce.stencil.has_slab)
      ++stats.slab_count;
    else {
      ++stats.pierce_fail_count;
      ++stats.fail_by_reason[pierce.reason];
    }
  }
  return stats;
}

inline slab_prism_init_stats init_slab_prism_constraints(
    asawa::shell::shell &shell,
    std::vector<projection_constraint::ptr> &constraints,
    const slab_pierce_config &pierce_cfg, const slab_prism_config &prism_cfg,
    std::vector<sim_block::ptr> blocks, solver_context &ctx,
    bool dump_failures = true) {
  if (blocks.empty())
    return {};
  const vecX q_rest = pack_block_q(blocks);
  const std::vector<vec3> x_shell =
      gather_shell_positions(shell, blocks[0], q_rest);
  arp::bvh_tree<3>::ptr face_tree = make_shell_face_bvh(shell, x_shell);
  if (!face_tree)
    return {};
  const std::vector<index_t> face_map = shell.get_face_map(true);
  return init_slab_prism_constraints(shell, constraints, *face_tree, face_map,
                                   pierce_cfg, prism_cfg, blocks, ctx,
                                   dump_failures);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PRISM_INIT_HPP
