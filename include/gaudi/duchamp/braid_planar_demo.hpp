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
#include "gaudi/duchamp/modules/tp_energies.hpp"
#include "gaudi/duchamp/modules/tp_line_search.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
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
  // windychien::knot_spec knot{.name = "K11a359"};
  windychien::knot_spec knot{
      .name = "plain_weave", .strands = 50, .frames = 50};
  braid_planar_params planar{};
  braid_sphere_params sphere{.n_lat = 256};
  vec3 rod_color = vec3(0.95, 0.55, 0.15);
  vec3 planar_color = vec3(0.35, 0.75, 1.0);
  real rod_radius = 0.01;
  bool show_planar_graph = true;

  real edge_stretch_w = 1.0e-5;
  /// >0: quaternion block + stretch_shear + four-frame straight. 0: edge_stretch only.
  real straight_w = 1.0;
  real smooth_w = 1.0e-4;
  real min_kink_w = 1.0e-3;
  real collision_w = 1.0;
  real dt = 0.2;
  real damping = 0.05;
  int iterations = 20;

  tangent_point_solver_config tp{
      .type = {tangent_point_type::cauchy, tangent_point_type::sobolev},
      //.type = {tangent_point_type::sobolev},

      .regularized = {.w = 1.0, .l0 = 1.0, .p = 3.0},
      .soft = {.w = 1.0, .p = 6.0, .R_min_frac = 0.1, .tau_frac = 1.0},
      .cauchy = {.w = 1e-6, .l0 = 0.1, .p = 6.0},
      .sobolev = {.w = 1.0,
                  .l0 = 0.001,
                  .p = 2.0,
                  .alpha_B = 1.0,
                  .alpha_0 = 1.0,
                  .mass = 0.0001,
                  .jacobi_iters = 16,
                  .num_v_cycles = 4,
                  .max_levels = 12,
                  .bh_eps = 0.25,
                  .cone_eps = 1.0},
      .line_search = {.alpha_min = 1.0e-6,
                      .alpha_max = 0.0,
                      .fallback_alpha = 0.0,
                      .trust_frac_lavg = 128.0,
                      .shrink = 0.5,
                      .max_backtracks = 32,
                      .verbose = true},
      .peak_scale = 1.0,
      .ccd_dx_max_frac = 0.1,
      .ccd = {.dx_max = 0.0, .dt_max = 1.0, .dt_min = 1.0e-8, .max_substeps = 32},
      .neighborhood = {.cfl_frac = 0.8,
                       .agg = cfl_agg_mode::hard_min,
                       .gate_frac = 0.1,
                       .use_envelope = true,
                       .chain_skip_hops = 8}};

  bool use_velocity_block = true;
  bool ccd_enable = true;
  bool ccd_verbose = true;
  bool show_tp_gradients = true;
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
    const windychien::braid b = windychien::make_knot(_cfg.knot);

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
    const bool use_frames = _cfg.straight_w > real(0);
    if (use_frames)
      _rod_quat =
          std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);

    _ccd_dx_max_explicit = (_cfg.tp.ccd.dx_max > 0.0);
    ccd_config ccd = _cfg.tp.ccd;
    if (ccd.dt_max <= 0.0)
      ccd.dt_max = _cfg.dt;
    if (!_cfg.use_velocity_block && !_ccd_dx_max_explicit &&
        _cfg.tp.ccd_dx_max_frac > 0.0)
      ccd.dx_max = _cfg.tp.ccd_dx_max_frac * lavg;
    ccd.neighborhood = _cfg.tp.neighborhood;
    if (ccd.neighborhood.radius <= 0.0)
      ccd.neighborhood.radius = _cfg.rod_radius;
    if (ccd.neighborhood.dt_max <= 0.0)
      ccd.neighborhood.dt_max = ccd.dt_max;
    _ccd = ccd;

    if (_cfg.use_velocity_block && use_frames) {
      _solver_config_vel_frame =
          hepworth::block::block_solver_builder<
              hepworth::block::rod_velocity_block,
              hepworth::block::rod_quaternion_block>::create()
              .with_blocks(_rod_vel, _rod_quat)
              .with_recompute(
                  hepworth::block::make_rod_stretch_shear_recompute<0, 1>(
                      __R, _cfg.edge_stretch_w))
              .with_recompute(
                  hepworth::block::make_rod_straight_recompute<1>(
                      __R, _cfg.straight_w))
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

      _solver_node_vel_frame =
          _graph.create_node<hepworth::block::block_solver_node<
              hepworth::block::rod_velocity_block,
              hepworth::block::rod_quaternion_block>>(_solver_config_vel_frame);
    } else if (_cfg.use_velocity_block) {
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
    } else if (use_frames) {
      _solver_config_pos_frame =
          hepworth::block::block_solver_builder<
              hepworth::block::rod_position_block,
              hepworth::block::rod_quaternion_block>::create()
              .with_blocks(_rod_pos, _rod_quat)
              .with_recompute(
                  hepworth::block::make_rod_stretch_shear_recompute<0, 1>(
                      __R, _cfg.edge_stretch_w))
              .with_recompute(
                  hepworth::block::make_rod_straight_recompute<1>(
                      __R, _cfg.straight_w))
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

      _solver_node_pos_frame =
          _graph.create_node<hepworth::block::block_solver_node<
              hepworth::block::rod_position_block,
              hepworth::block::rod_quaternion_block>>(_solver_config_pos_frame);
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
    if (_solver_node_vel_frame)
      _tp.link_to(_graph, _solver_node_vel_frame->input_at<0>());
    else if (_solver_node_vel)
      _tp.link_to(_graph, _solver_node_vel->input_at<0>());
    else if (_solver_node_pos_frame)
      _tp.link_to(_graph, _solver_node_pos_frame->input_at<0>());
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
      if (_cfg.use_velocity_block || _ccd_dx_max_explicit ||
          !(_cfg.tp.ccd_dx_max_frac > 0.0))
        return;
      const real lavg = std::max(__R->lavg(), real(1e-6));
      _ccd.dx_max = _cfg.tp.ccd_dx_max_frac * lavg;
    };

    auto refresh_tp = [this](real h) {
      _tp.set_step_h(h);
      _tp.compute();
    };

    auto mutable_drive_field = [this]() -> std::vector<vec3> & {
      return _tp.data();
    };

    auto draw_drive = [this](const std::vector<vec3> &drive) {
      if (_cfg.show_tp_gradients)
        draw_rod_drive_field(*__R, drive);
    };

    auto peak_scale_drive = [this](std::vector<vec3> &drive) {
      if (!(_cfg.tp.peak_scale > 0.0)) {
        _ccd_breakdown.drive_scale = 1.0;
        return;
      }
      const real peak = ccd_peak_norm(drive);
      if (!(peak > _ccd.eps) || !std::isfinite(peak)) {
        _ccd_breakdown.drive_scale = 1.0;
        return;
      }
      const real alpha = real(1.0) / (_cfg.tp.peak_scale * peak);
      for (vec3 &a : drive)
        a *= alpha;
      _ccd_breakdown.drive_scale = alpha;
    };

    auto prepare_drive = [&draw_drive, &peak_scale_drive](std::vector<vec3> &drive) {
      peak_scale_drive(drive);
      draw_drive(drive);
    };

    auto run_solver = [this](real dti) {
      if (_solver_node_vel_frame) {
        _solver_node_vel_frame->set_dt(dti);
        _solver_node_vel_frame->compute();
      } else if (_solver_node_vel) {
        _solver_node_vel->set_dt(dti);
        _solver_node_vel->compute();
      } else if (_solver_node_pos_frame) {
        _solver_node_pos_frame->set_dt(dti);
        _solver_node_pos_frame->compute();
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
      prepare_drive(mutable_drive_field());
      run_solver(_cfg.dt);
      __Rd->step();
      __Rd->step();
      __Rd->step();
      __R->_update_frames();
      return;
    }

    std::cout << "[braid CCD] begin frame=" << frame
              << " schedule=fill_alpha" << std::endl;

    real alpha_star = 0.0;
    real alpha_left = 0.0;
    real E0_frame = 0.0;
    real E_star = 0.0;
    std::vector<vec3> frozen_disp;
    bool have_frozen = false;

    {
      tp_curvature_propose prop =
          tp_propose_curvature_step(*__R, *__Rd, _cfg.tp);
      frozen_disp = std::move(prop.displacement);
      const tp_line_search_result &ls = prop.ls;
      alpha_star = ls.alpha;
      E0_frame = ls.E0;
      E_star = ls.E_accepted;
      alpha_left = ccd_peak_norm(frozen_disp);
      have_frozen = true;
      if (_cfg.ccd_verbose) {
        std::cout << "[braid α*] frame=" << frame << " alpha=" << alpha_star
                  << " |δ|_peak=" << alpha_left << " E0=" << ls.E0
                  << " E*=" << ls.E_accepted
                  << " dE_ray=" << (ls.E_accepted - ls.E0)
                  << " backtracks=" << ls.n_backtracks
                  << " fallback=" << (ls.used_fallback ? 1 : 0)
                  << " propose_ms=" << prop.propose_ms << std::endl;
      }
    }

    const real target = std::max(alpha_left, _ccd.eps);
    int sub_i = 0;

    _ccd_stats = ccd_step(
        target,
        [this, &sync_dx_max, &mutable_drive_field, &run_solver, &sub_i,
         &log_ccd_diag, frame, &alpha_left, &frozen_disp,
         &have_frozen](real /*accum*/) -> real {
          if (alpha_left <= _ccd.eps)
            return 0.0;

          sync_dx_max();
          const real h = std::max(_cfg.dt, real(1e-20));

          if (have_frozen) {
            std::vector<vec3> &drive = mutable_drive_field();
            drive = frozen_disp;
            if (_cfg.use_velocity_block) {
              // v = δ. ccd_velocity returns the safe fraction of δ.
              _ccd_breakdown.drive_scale = 1.0;
            } else {
              const real h2 = std::max(h * h, real(1e-20));
              for (auto &a : drive)
                a /= h2;
              real bite_cap = _ccd.dx_max;
              if (alpha_left > 0.0) {
                if (!(bite_cap > 0.0) || alpha_left < bite_cap)
                  bite_cap = alpha_left;
              }
              if (bite_cap > 0.0) {
                _ccd_breakdown.drive_scale = ccd_scale_force_to_dx_max(
                    drive, h, bite_cap, _ccd.eps);
              }
            }
            if (_cfg.show_tp_gradients)
              draw_rod_drive_field(*__R, drive);
          }

          std::vector<vec3> &drive = mutable_drive_field();
          real dti = 0.0;
          if (_cfg.use_velocity_block) {
            dti = ccd_velocity(*__R, drive, _ccd, &_ccd_stats.last,
                               &_ccd_breakdown);
            // v is the remaining δ, so h > 1 would step past it.
            if (have_frozen)
              dti = std::min(dti, real(1.0));
          } else {
            _ccd_breakdown.drive = ccd_field_stats_vec(drive);
            ccd_diagnose_rod(*__R, &_ccd_breakdown.rod_x);
            dti = ccd_force(*__R, drive, h, _ccd, &_ccd_stats.last);
            _ccd_breakdown.dt_geom = _ccd_stats.last.dt;
            _ccd_breakdown.dt = dti;
          }
          log_ccd_diag("pre", frame, sub_i, dti);

          if (!(dti > 0.0) || !std::isfinite(dti))
            return 0.0;

          real peak_disp = 0.0;
          if (_cfg.use_velocity_block)
            peak_disp = dti * ccd_peak_norm(drive);
          else
            peak_disp = (dti * dti) * ccd_peak_norm(drive);

          const bool null_bite =
              (peak_disp <= _ccd.eps ||
               (_ccd.dt_min > 0.0 && dti <= _ccd.dt_min + _ccd.eps &&
                _ccd_breakdown.dt_geom <= _ccd.eps));
          if (null_bite) {
            alpha_left = 0.0;
            if (_cfg.ccd_verbose) {
              std::cout << "[braid α*] abort null bite dti=" << dti
                        << " dt_geom=" << _ccd_breakdown.dt_geom
                        << " bite_peak=" << peak_disp << std::endl;
            }
            return 0.0;
          }

          run_solver(dti);
          __Rd->step();
          __R->_update_frames();

          if (have_frozen) {
            const real peak0 = ccd_peak_norm(frozen_disp);
            if (peak0 > _ccd.eps && peak_disp > _ccd.eps) {
              const real frac = std::min(real(1.0), peak_disp / peak0);
              for (auto &a : frozen_disp)
                a *= (real(1.0) - frac);
              alpha_left = ccd_peak_norm(frozen_disp);
            } else {
              alpha_left = 0.0;
            }
            if (_cfg.ccd_verbose) {
              std::cout << "[braid α*] alpha_left=" << alpha_left
                        << " bite_peak=" << peak_disp << std::endl;
            }
          }

          ccd_diagnose_rod(*__R, &_ccd_breakdown.rod_x);
          log_ccd_diag("post", frame, sub_i, dti);
          ++sub_i;

          if (!(peak_disp > _ccd.eps))
            return 0.0;
          return peak_disp;
        },
        _ccd);

    std::cout << "[braid CCD] frame=" << frame
              << " substeps=" << _ccd_stats.substeps
              << " accum=" << _ccd_stats.dt_accum << "/" << target
              << " dti=" << _ccd_stats.last_dti << " alpha*=" << alpha_star
              << " alpha_left=" << alpha_left << std::endl;

    if (_cfg.ccd_verbose && have_frozen) {
      const real E_tp_post = tp_energy(*__R, *__Rd, _cfg.tp);
      const real dE_frame = E_tp_post - E0_frame;
      const real dE_prev =
          _have_E_prev_post ? (E_tp_post - _E_prev_post) : real(0.0);
      std::cout << "[braid E] frame=" << frame << " E0=" << E0_frame
                << " E*=" << E_star << " Epost=" << E_tp_post
                << " dE_ray=" << (E_star - E0_frame)
                << " dE_frame=" << dE_frame << " dE_prev=" << dE_prev
                << (dE_frame > 0.0 ? " ASCENT" : "") << std::endl;
      _E_prev_post = E_tp_post;
      _have_E_prev_post = true;
    }
  }

private:
  braid_planar_demo_config _cfg;
  asawa::rod::rod::ptr __R_planar;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;

  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_velocity_block::ptr _rod_vel;
  hepworth::block::rod_quaternion_block::ptr _rod_quat;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block>
      _solver_config_pos;
  hepworth::block::block_solver_config<hepworth::block::rod_velocity_block>
      _solver_config_vel;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block,
                                       hepworth::block::rod_quaternion_block>
      _solver_config_pos_frame;
  hepworth::block::block_solver_config<hepworth::block::rod_velocity_block,
                                       hepworth::block::rod_quaternion_block>
      _solver_config_vel_frame;
  liblombardi::GraphContext _graph;
  tangent_point_drive _tp{};
  hepworth::block::block_solver_node<hepworth::block::rod_position_block>::ptr
      _solver_node_pos;
  hepworth::block::block_solver_node<hepworth::block::rod_velocity_block>::ptr
      _solver_node_vel;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                     hepworth::block::rod_quaternion_block>::ptr
      _solver_node_pos_frame;
  hepworth::block::block_solver_node<hepworth::block::rod_velocity_block,
                                     hepworth::block::rod_quaternion_block>::ptr
      _solver_node_vel_frame;
  ccd_config _ccd{};
  ccd_step_stats _ccd_stats{};
  ccd_dt_breakdown _ccd_breakdown{};
  bool _ccd_dx_max_explicit = false;
  real _E_prev_post = 0.0;
  bool _have_E_prev_post = false;
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
