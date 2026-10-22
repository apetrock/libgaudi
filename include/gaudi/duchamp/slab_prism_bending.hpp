#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/projection_diagnostics.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/constraints/slab/slab_pierce_stencil.hpp"
#include "gaudi/hepworth/constraints/slab/slab_prism_init.hpp"
#include "gaudi/hepworth/block/slab_prism_log.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/constraints/slab/slab_foot_projection_probe.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

namespace gaudi {
namespace duchamp {

/// Slab prism bending: pierce-stencil constraints + optional shell ARAP/bending.
/// Edit weights / `_probe_accel` in `reset_scene()`.
class slab_prism_bending {
public:
  using ptr = std::shared_ptr<slab_prism_bending>;

  static ptr create(const std::string &mesh_path = {}) {
    return std::make_shared<slab_prism_bending>(mesh_path);
  }

  explicit slab_prism_bending(const std::string &mesh_path = {})
      : _mesh_path(mesh_path) {
    reset_scene();
  }

  void reset_scene() {
    if (_mesh_path.empty())
      __M = asawa::shell::load_bunny();
    else
      __M = asawa::shell::load_obj(_mesh_path);

    asawa::shell::triangulate(*__M);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    auto xs = std::make_shared<shell_vert_positions>(__M, 0);
    auto vs = std::make_shared<shell_vert_velocities>(__M, 1);
    _shell = std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);
    _shell->with_force([this]() { return compute_external_forces(); });

    const real avg_edge = asawa::shell::avg_length(*__M, xs->get());
    const real l0 = 3.0 * avg_edge;
    _pin_w = 1.0e0;
    _solver_dt = 0.01;
    _solver_iters = 10;
    _probe_accel = 1000.0 / std::max(l0, real(1e-4));
    _force_draw_scale = 3.0 * l0 / std::max(_probe_accel, real(1.0));
    _pin_sphere_r = std::max(0.35 * l0, real(0.004));
    __surf = asawa::shell::dynamic::create(__M, l0, 1.25 * l0, 0.5 * l0);
    dynamic_warmup(l0);
    zero_shell_velocities();

    _x_slab_rest = _shell->xs->get();
    _i_pin_rest = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    asawa::get_vec_data(*__M, _i_pin_rest) = _x_slab_rest;
    capture_pin_y_bounds(_x_slab_rest);
    compute_rest_height_bounds(_x_slab_rest);
    pick_point_normal_loads(_x_slab_rest);

    _prism_cfg.w_field = 0.0;
    _prism_cfg.w_phi = 0.0;
    _prism_cfg.w_rbf = 20.0;
    _prism_cfg.rbf_mu = 1.0;
    _prism_cfg.rbf_eta = 1.0;
    _prism_cfg.rbf_nr_iters = 8;
    _prism_cfg.rbf_n_centers = 4;
    _prism_cfg.w_bulk = 0.0;
    _prism_cfg.w_twist = 0.0;
    _shell_arap_w = 0.0;
    _shell_bending_w = 0.0;
    _prism_cfg.field_live = hepworth::block::field_live_update::refit_every_n;
    _prism_cfg.field_refit_period = 4;
    _prism_cfg.field_gd_iters = 8;
    _prism_cfg.field_family = hepworth::block::slab_field_family::darboux;

    _use_external_forces = true;

    rebuild_solver_config();
    if (!_point_loads.empty() &&
        (hepworth::block::projection_diag_enabled() ||
         hepworth::block::slab_prism_log_enabled()))
      hepworth::block::set_slab_foot_probe(_point_loads[0].vert);
  }

  void step(int frame) {
    (void)frame;
    hepworth::block::run_solver_step(_config, _solver);
    __surf->step(true);
    // Pin targets live in a vert datum (pack keeps them aligned). Re-pick loads
    // after remesh; slab constraints rebuild every solver step via the bundle.
    compute_rest_height_bounds(_shell->xs->get());
    pick_point_normal_loads(_shell->xs->get());
  }

  void draw_debug() const {
    draw_pinned_verts();
    draw_external_forces();
  }

  index_t primary_load_vert() const {
    return _point_loads.empty() ? -1 : _point_loads[0].vert;
  }

  real probe_accel() const { return _probe_accel; }

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;

private:
  void dynamic_warmup(real l0) {
    const int steps = 4;
    for (int k = 0; k < steps; ++k)
      __surf->step(true);
    (void)l0;
  }

  void zero_shell_velocities() {
    for (vec3 &vi : _shell->vs->get())
      vi.setZero();
  }

  void rebuild_solver_config() {
    auto builder =
        hepworth::block::block_solver_builder<
            hepworth::block::shell_position_block>::create()
            .with_blocks(_shell)
            .with_bundle(hepworth::block::make_shell_slab_prism_bundle<0>(
                _shell, _pierce_cfg, _prism_cfg));
    if (_shell_arap_w > 0.0)
      builder.with_recompute(
          hepworth::block::make_shell_arap_recompute<0>(_shell, _shell_arap_w));
    if (_shell_bending_w > 0.0)
      builder.with_recompute(hepworth::block::make_shell_bending_recompute<0>(
          _shell, _shell_bending_w));
    _config = builder
                    .with_recompute([this](hepworth::block::solver_context &ctx) {
                      append_bottom_pins(ctx);
                    })
                    .dt(_solver_dt)
                    .damping(0.2)
                    .iterations(_solver_iters)
                    .build();
  }

  void compute_rest_height_bounds(const std::vector<vec3> &x) {
    _y_min = std::numeric_limits<real>::infinity();
    _y_max = -std::numeric_limits<real>::infinity();
    for (const vec3 &xi : x) {
      _y_min = std::min(_y_min, xi.y());
      _y_max = std::max(_y_max, xi.y());
    }
    if (!std::isfinite(_y_min) || !std::isfinite(_y_max) || _y_max <= _y_min) {
      _y_min = -1.0;
      _y_max = 1.0;
    }
  }

  void capture_pin_y_bounds(const std::vector<vec3> &pin_rest) {
    _pin_y_min = std::numeric_limits<real>::infinity();
    _pin_y_max = -std::numeric_limits<real>::infinity();
    for (const vec3 &xi : pin_rest) {
      _pin_y_min = std::min(_pin_y_min, xi.y());
      _pin_y_max = std::max(_pin_y_max, xi.y());
    }
    if (!std::isfinite(_pin_y_min) || !std::isfinite(_pin_y_max) ||
        _pin_y_max <= _pin_y_min) {
      _pin_y_min = -1.0;
      _pin_y_max = 1.0;
    }
  }

  const std::vector<vec3> &pin_rest_positions() const {
    return asawa::const_get_vec_data(*__M, _i_pin_rest);
  }

  void append_bottom_pins(hepworth::block::solver_context &ctx) {
    if (!(_pin_w > 0.0) || _i_pin_rest < 0)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> &pin_rest = pin_rest_positions();
    if (pin_rest.size() != x.size())
      return;
    const real y_cut =
        _pin_y_min + _bottom_pin_frac * (_pin_y_max - _pin_y_min);
    std::vector<hepworth::sim_block::ptr> blocks;
    _shell->emit_dof_blocks(blocks);
    for (int iv = 0; iv < __M->vert_count(); ++iv) {
      const size_t si = static_cast<size_t>(iv);
      if (x[si].y() > y_cut)
        continue;
      ctx.constraints.push_back(hepworth::block::pinned::create(
          std::vector<index_t>({static_cast<index_t>(iv)}), pin_rest[si],
          _pin_w, blocks));
    }
  }

  void draw_pinned_verts() const {
    if (!(_pin_w > 0.0) || _i_pin_rest < 0)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> &pin_rest = pin_rest_positions();
    if (pin_rest.size() != x.size())
      return;
    const real y_cut =
        _pin_y_min + _bottom_pin_frac * (_pin_y_max - _pin_y_min);
    const vec4 pin_col(0.15, 1.0, 0.35, 1.0);
    const vec4 pin_rest_col(0.15, 1.0, 0.35, 0.55);
    for (int iv = 0; iv < __M->vert_count(); ++iv) {
      const size_t si = static_cast<size_t>(iv);
      if (x[si].y() > y_cut)
        continue;
      const vec3 &p = x[si];
      geometry_logger::sphere(p, _pin_sphere_r, pin_col);
      geometry_logger::line(p, pin_rest[si], pin_rest_col, 0.003);
    }
  }

  void draw_external_forces() const {
    if (!_use_external_forces || _point_loads.empty())
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> N = asawa::shell::vertex_normals(*__M, x);
    const vec4 load_col(1.0, 0.45, 0.12, 1.0);
    const real arrow_r = 0.004;
    for (const point_normal_load &load : _point_loads) {
      const size_t vi = static_cast<size_t>(load.vert);
      if (vi >= x.size() || vi >= N.size())
        continue;
      const vec3 n_hat = N[vi].normalized();
      const vec3 f = load.accel * n_hat;
      geometry_logger::sphere(x[vi], _pin_sphere_r * 1.8, load_col);
      geometry_logger::line(x[vi], x[vi] + _force_draw_scale * f, load_col,
                            arrow_r);
    }
  }

  std::vector<vec3> compute_external_forces() const {
    const size_t n = static_cast<size_t>(__M->vert_count());
    std::vector<vec3> f(n, vec3::Zero());
    if (!_use_external_forces || _point_loads.empty())
      return f;

    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<vec3> N = asawa::shell::vertex_normals(*__M, x);
    for (const point_normal_load &load : _point_loads) {
      const size_t vi = static_cast<size_t>(load.vert);
      if (vi >= n || vi >= N.size())
        continue;
      const real nn = N[vi].norm();
      if (nn < 1e-12)
        continue;
      f[vi] += load.accel * (N[vi] / nn);
    }
    return f;
  }

  struct point_normal_load {
    index_t vert = -1;
    real accel = 0.0;
  };

  void pick_point_normal_loads(const std::vector<vec3> &x) {
    _point_loads.clear();
    if (x.empty())
      return;

    const int n = static_cast<int>(x.size());
    index_t top = 0;
    for (int i = 1; i < n; ++i) {
      if (x[static_cast<size_t>(i)].y() > x[static_cast<size_t>(top)].y())
        top = static_cast<index_t>(i);
    }
    _point_loads.push_back({top, _probe_accel});

    const real y_cut = _y_min + 0.72 * (_y_max - _y_min);
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
      _point_loads.push_back({side, 0.85 * _probe_accel});
  }

  std::string _mesh_path;
  hepworth::block::shell_position_block::ptr _shell;
  std::vector<vec3> _x_slab_rest;
  /// Pin targets live in a vert datum so subdivide/collapse/pack keep indices aligned.
  index_t _i_pin_rest = -1;
  real _pin_y_min = -1.0;
  real _pin_y_max = 1.0;
  hepworth::block::slab_pierce_config _pierce_cfg;
  hepworth::block::slab_prism_config _prism_cfg;

  real _shell_arap_w = 0.0;
  real _shell_bending_w = 0.0;

  real _y_min = -1.0;
  real _y_max = 1.0;

  real _pin_w = 1.0e0;
  real _solver_dt = 0.01;
  int _solver_iters = 6;

  bool _use_external_forces = false;
  real _probe_accel = 1.0;
  real _bottom_pin_frac = 0.08;
  real _force_draw_scale = 0.002;
  real _pin_sphere_r = 0.004;
  std::vector<point_normal_load> _point_loads;

  hepworth::block::block_solver_config<hepworth::block::shell_position_block>
      _config;
  hepworth::block::projection_solver _solver;
};

} // namespace duchamp
} // namespace gaudi
