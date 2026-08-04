#ifndef GAUDI_DUCHAMP_BRAID_PLANAR_DEMO_HPP
#define GAUDI_DUCHAMP_BRAID_PLANAR_DEMO_HPP

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/braid_planar_rod.hpp"
#include "gaudi/duchamp/braid_sphere_rod.hpp"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/modules/ccd_step.hpp"
#include "gaudi/duchamp/modules/rod_forces.hpp"
#include "gaudi/duchamp/modules/rod_savitzky_golay.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_velocity_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/windychien/braid.hpp"
#include "gaudi/windychien/catalog.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

struct braid_planar_demo_config {
  bool use_plain_weave = false;
  int weave_strands = 50;
  int weave_frames = 50;
  std::string knot_name = "K11a359";
  braid_planar_params planar{};
  braid_sphere_params sphere{.n_lat = 256};
  vec3 rod_color = vec3(0.95, 0.55, 0.15);
  vec3 planar_color = vec3(0.35, 0.75, 1.0);
  real rod_radius = 0.01;
  bool show_planar_graph = true;

  real edge_stretch_w = 1.0e-4;
  real smooth_w = 1.0e-3;
  real min_kink_w = 1.0e-3;
  real collision_w = 1.0;
  real dt = 0.05;
  real damping = 0.5;
  int iterations = 20;

  tangent_point_solver_config tp{
      .type = tangent_point_type::soft,
      //.type = tangent_point_type::regularized,
      //.type = tangent_point_type::harmonic,
      .regularized = {.w = 1.0, .l0 = 1.0, .p = 6.0},
      .soft = {.w = 1.0e-1,
               .p = 6.0,
               .R_min_frac = 0.1,
               .tau_frac = 1.0,
               //.mode = soft_tp_mode::gradient,
               .mode = soft_tp_mode::hessian_force,
               //.mode = soft_tp_mode::newton,
               .newton = {
                   .mode = soft_tp_newton_mode::full_step,
                   //.mode = soft_tp_newton_mode::single_iter,
                   .newton_iters = 3,
                   .verbose = true}},
      .harmonic = {.w = 1.0e-1,
                   .l0 = 1.0,  ///< sharp local TP (shell: 0.1·Cc)
                   .l1 = 1.0,  ///< harmonic smooth over neighbors (shell: 4·Cc)
                   .p0 = 6.0,
                   .p1 = 3.0}};

  bool use_velocity_block = true;

  /// Savitzky–Golay smooth of TP drive along each strand (before CCD / dx cap).
  rod_sg_filter_config sg_filter{
      .enable = false,
      .window = 8,
      .poly_order = 4};

  /// CCD substep sizing (geometry only — does not rescale TP drive).
  ccd_config ccd{.dx_max = 0.01, .dt_max = 0.05,  .max_substeps=16};
  /// dx_max ≤ 0 ⇒ ccd_dx_max_frac · lavg at init.
  real ccd_dx_max_frac = 0.01;
  ccd_neighborhood_config neighborhood{};
  bool ccd_enable = true;
  bool ccd_verbose = true;
  bool show_tp_gradients = false; ///< draw SG-smoothed drive before CCD
};

inline const braid_planar_demo_config k_braid_planar_demo_default{};

/// Graph: (soft|raw) TP → solver. Frame loop: ccd_step until dt is filled.
class braid_planar_demo {
public:
  using ptr = std::shared_ptr<braid_planar_demo>;

  static ptr create(const braid_planar_demo_config &cfg =
                        k_braid_planar_demo_default) {
    return std::make_shared<braid_planar_demo>(cfg);
  }

  explicit braid_planar_demo(
      const braid_planar_demo_config &cfg = k_braid_planar_demo_default)
      : _cfg(cfg) {
    const windychien::braid b =
        _cfg.use_plain_weave
            ? windychien::make_plain_weave(_cfg.weave_strands,
                                           _cfg.weave_frames)
            : windychien::load_braid(_cfg.knot_name);

    braid_planar_params planar = _cfg.planar;
    planar.close = true;
    planar.center = false;
    planar.eps_z = _cfg.rod_radius + real(0.1) * _cfg.rod_radius;
    if (planar.dx <= 0.0)
      planar.dx = 1.0;
    if (planar.dy <= 0.0)
      planar.dy = planar.dx;

    braid_planar_params planar_viz = planar;
    planar_viz.close = false;
    __R_planar = braid_to_planar_rod(b, planar_viz);

    __R = braid_to_sphere_rod(b, planar, _cfg.sphere);
    __R->_r = _cfg.rod_radius;
    __R_planar->_r = _cfg.rod_radius;

    if (_cfg.show_planar_graph) {
      real x_min = __R_planar->x()[0][0], x_max = x_min;
      real y_min = __R_planar->x()[0][1], y_max = y_min;
      for (const vec3 &q : __R_planar->x()) {
        x_min = std::min(x_min, q[0]);
        x_max = std::max(x_max, q[0]);
        y_min = std::min(y_min, q[1]);
        y_max = std::max(y_max, q[1]);
      }
      const vec3 cen(0.5 * (x_min + x_max), 0.5 * (y_min + y_max), 0.0);
      const vec3 shift(0.0, 0.0, -_cfg.sphere.radius - 0.75);
      for (vec3 &q : __R_planar->x()) {
        q[0] -= cen[0];
        q[1] -= cen[1];
        q += shift;
      }
      __R_planar->_init_params();
    }

    const real lavg = std::max(__R->lavg(), real(1e-6));
    __Rd = asawa::rod::dynamic::create(__R, 1.0 * lavg, 2.0 * lavg,
                                       0.25 * lavg);

    _rod_pos = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);
    _rod_vel = std::make_shared<hepworth::block::rod_velocity_block>(__R, __Rd);

    _ccd_dx_max_explicit = (_cfg.ccd.dx_max > 0.0);
    ccd_config ccd = _cfg.ccd;
    if (ccd.dt_max <= 0.0)
      ccd.dt_max = _cfg.dt;
    if (!_ccd_dx_max_explicit && _cfg.ccd_dx_max_frac > 0.0)
      ccd.dx_max = _cfg.ccd_dx_max_frac * lavg;
    ccd.neighborhood = _cfg.neighborhood;
    if (ccd.neighborhood.radius <= 0.0)
      ccd.neighborhood.radius = _cfg.rod_radius;
    if (ccd.neighborhood.dt_max <= 0.0)
      ccd.neighborhood.dt_max = _cfg.dt;
    _ccd = ccd;

    if (_cfg.use_velocity_block) {
      _solver_config_vel =
          hepworth::block::block_solver_builder<
              hepworth::block::rod_velocity_block>::create()
              .with_blocks(_rod_vel)
              .with_recompute(
                  hepworth::block::make_rod_edge_stretch_recompute<0>(
                      __R, _cfg.edge_stretch_w))
              .with_recompute(hepworth::block::make_rod_smooth_recompute<0>(
                  __R, _cfg.smooth_w))
              .with_recompute(hepworth::block::make_rod_min_kink_recompute<0>(
                  __R, _cfg.min_kink_w))
              .with_recompute(
                  hepworth::block::make_rod_collisions_recompute<0>(
                      __R, __Rd, _cfg.collision_w))
              .dt(_cfg.dt)
              .damping(_cfg.damping)
              .iterations(_cfg.iterations)
              .build();

      _solver_node_vel =
          _graph.create_node<hepworth::block::block_solver_node<
              hepworth::block::rod_velocity_block>>(_solver_config_vel);
    } else {
      _solver_config_pos =
          hepworth::block::block_solver_builder<
              hepworth::block::rod_position_block>::create()
              .with_blocks(_rod_pos)
              .with_recompute(
                  hepworth::block::make_rod_edge_stretch_recompute<0>(
                      __R, _cfg.edge_stretch_w))
              .with_recompute(hepworth::block::make_rod_smooth_recompute<0>(
                  __R, _cfg.smooth_w))
              .with_recompute(hepworth::block::make_rod_min_kink_recompute<0>(
                  __R, _cfg.min_kink_w))
              .with_recompute(
                  hepworth::block::make_rod_collisions_recompute<0>(
                      __R, __Rd, _cfg.collision_w))
              .dt(_cfg.dt)
              .damping(_cfg.damping)
              .iterations(_cfg.iterations)
              .build();

      _solver_node_pos =
          _graph.create_node<hepworth::block::block_solver_node<
              hepworth::block::rod_position_block>>(_solver_config_pos);
    }

    _tp = tangent_point_drive::create(_graph, __R, __Rd, _cfg.tp,
                                      _cfg.rod_radius);
    _tp.set_velocity_output(_cfg.use_velocity_block);
    if (_cfg.use_velocity_block)
      _tp.link_to(_graph, _solver_node_vel->input_at<0>());
    else
      _tp.link_to(_graph, _solver_node_pos->input_at<0>());
  }

  const braid_planar_demo_config &config() const { return _cfg; }
  const asawa::rod::rod &rod() const { return *__R; }
  const asawa::rod::rod &planar_rod() const { return *__R_planar; }
  asawa::rod::rod::ptr rod_ptr() const { return __R; }
  const ccd_step_stats &last_ccd_stats() const { return _ccd_stats; }
  const ccd_neighborhood_result &last_neighborhood() const {
    return _ccd_stats.last;
  }

  void step(int frame) {
    auto sync_dx_max = [this]() {
      if (_ccd_dx_max_explicit || !(_cfg.ccd_dx_max_frac > 0.0))
        return;
      const real lavg = std::max(__R->lavg(), real(1e-6));
      _ccd.dx_max = _cfg.ccd_dx_max_frac * lavg;
    };

    auto refresh_tp = [this](real h) {
      _tp.set_step_h(h);
      _tp.compute();
    };

    auto mutable_drive_field = [this]() -> std::vector<vec3> & {
      return _tp.data();
    };

    auto filter_drive = [this](std::vector<vec3> &drive) {
      rod_sg_filter_field(*__R, drive, _cfg.sg_filter);
    };

    auto draw_drive = [this](const std::vector<vec3> &drive) {
      if (_cfg.show_tp_gradients)
        draw_rod_drive_field(*__R, drive);
    };

    auto scale_drive_to_dx_max = [this](std::vector<vec3> &drive, real h) {
      if (_cfg.use_velocity_block) {
        _ccd_breakdown.drive_scale =
            ccd_scale_field_to_dx_max(drive, _ccd.dx_max, _ccd.eps);
      } else {
        _ccd_breakdown.drive_scale =
            ccd_scale_force_to_dx_max(drive, h, _ccd.dx_max, _ccd.eps);
      }
    };

    auto prepare_drive = [this, &filter_drive, &draw_drive,
                          &scale_drive_to_dx_max](std::vector<vec3> &drive,
                                                  real h) {
      filter_drive(drive);
      draw_drive(drive);
      scale_drive_to_dx_max(drive, h);
    };

    auto run_solver = [this](real dti) {
      if (_cfg.use_velocity_block) {
        _solver_node_vel->set_dt(dti);
        _solver_node_vel->compute();
      } else {
        _solver_node_pos->set_dt(dti);
        _solver_node_pos->compute();
      }
    };

    auto log_ccd_diag = [this](const char *tag, int frame, int sub,
                               real dti) {
      if (!_cfg.ccd_verbose)
        return;
      const ccd_field_stats &d = _ccd_breakdown.drive;
      const ccd_field_stats &x = _ccd_breakdown.rod_x;
      std::cout << "[braid CCD] " << tag << " frame=" << frame << " sub=" << sub
                << " dti=" << dti << " dt_geom=" << _ccd_breakdown.dt_geom
                << " dt_dx=" << _ccd_breakdown.dt_dx << " dx_max=" << _ccd.dx_max
                << " drive_scale=" << _ccd_breakdown.drive_scale
                << " limit=" << cfl_limit_term_name(_ccd_stats.last.limiting)
                << " max|u|/lavg=" << _ccd_stats.last.u_over_lavg << "\n"
                << "  drive n=" << d.n << " peak=" << d.peak
                << " rms=" << d.rms << " i_peak=" << d.i_peak
                << " nonfinite=" << d.n_nonfinite << "\n"
                << "  rod   n=" << x.n << " |x|_max=" << x.peak
                << " rms=" << x.rms << " i_max=" << x.i_peak
                << " nonfinite=" << x.n_nonfinite << std::endl;
    };

    if (!_cfg.ccd_enable) {
      refresh_tp(_cfg.dt);
      prepare_drive(mutable_drive_field(), _cfg.dt);
      run_solver(_cfg.dt);
      __Rd->step();
      __Rd->step();
      __Rd->step();
      __R->_update_frames();
    } else {
      std::cout << "[braid CCD] begin frame=" << frame << std::endl;
      int sub_i = 0;

      _ccd_stats = ccd_step(
          _cfg.dt,
          [this, &sync_dx_max, &refresh_tp, &mutable_drive_field, &prepare_drive,
           &run_solver, &sub_i, &log_ccd_diag, frame](real accum) -> real {
            sync_dx_max();
            const real h = (_ccd_stats.last_dti > 0.0) ? _ccd_stats.last_dti
                                                       : _cfg.dt;
            refresh_tp(h);
            std::vector<vec3> &drive = mutable_drive_field();
            prepare_drive(drive, h);
            real dti = 0.0;
            if (_cfg.use_velocity_block) {
              dti = ccd_velocity(*__R, drive, _ccd, &_ccd_stats.last,
                                 &_ccd_breakdown);
            } else {
              _ccd_breakdown.drive = ccd_field_stats_vec(drive);
              ccd_diagnose_rod(*__R, &_ccd_breakdown.rod_x);
              dti = ccd_force(*__R, drive, h, _ccd, &_ccd_stats.last);
              _ccd_breakdown.dt_geom = _ccd_stats.last.dt;
              _ccd_breakdown.dt = dti;
            }
            log_ccd_diag("pre", frame, sub_i, dti);

            if (!(dti > 0.0) || !std::isfinite(dti))
              dti = _cfg.dt - accum;
            dti = ccd_residual_dt(dti, _cfg.dt, accum);
            const bool at_floor =
                (_ccd.dt_min > 0.0 && dti <= _ccd.dt_min + _ccd.eps);
            if (at_floor)
              dti = _cfg.dt - accum;

            run_solver(dti);
            __Rd->step();
            __R->_update_frames();

            ccd_diagnose_rod(*__R, &_ccd_breakdown.rod_x);
            log_ccd_diag("post", frame, sub_i, dti);
            ++sub_i;
            return dti;
          },
          _ccd);

      std::cout << "[braid CCD] frame=" << frame
                << " substeps=" << _ccd_stats.substeps
                << " accum=" << _ccd_stats.dt_accum << "/" << _cfg.dt
                << " dti=" << _ccd_stats.last_dti << std::endl;
    }
  }

private:
  braid_planar_demo_config _cfg;
  asawa::rod::rod::ptr __R_planar;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;

  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_velocity_block::ptr _rod_vel;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block>
      _solver_config_pos;
  hepworth::block::block_solver_config<hepworth::block::rod_velocity_block>
      _solver_config_vel;
  liblombardi::GraphContext _graph;
  tangent_point_drive _tp{};
  hepworth::block::block_solver_node<hepworth::block::rod_position_block>::ptr
      _solver_node_pos;
  hepworth::block::block_solver_node<hepworth::block::rod_velocity_block>::ptr
      _solver_node_vel;
  ccd_config _ccd{};
  ccd_step_stats _ccd_stats{};
  ccd_dt_breakdown _ccd_breakdown{};
  bool _ccd_dx_max_explicit = false;
};

class braid_planar_demo_adapter : public demo_adapter<braid_planar_demo> {
public:
  static ptr create(const braid_planar_demo_config &cfg =
                        k_braid_planar_demo_default) {
    return std::make_shared<braid_planar_demo_adapter>(cfg);
  }

  explicit braid_planar_demo_adapter(
      const braid_planar_demo_config &cfg = k_braid_planar_demo_default)
      : demo_adapter<braid_planar_demo>(braid_planar_demo::create(cfg),
                                        "braid_sphere") {}

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), demo()->config().rod_color);
  }

  std::optional<rod_snapshot> rod_polyline() const override {
    if (!demo()->config().show_planar_graph)
      return std::nullopt;
    return make_rod_snapshot(demo()->planar_rod(),
                             demo()->config().planar_color);
  }
};

} // namespace duchamp
} // namespace gaudi

#endif
