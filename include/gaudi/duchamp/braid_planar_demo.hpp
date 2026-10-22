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
#include "gaudi/duchamp/demo_log.hpp"
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
  //windychien::knot_spec knot{.name = "K11a359"};
  //windychien::knot_spec knot{.name = "T15_14"};
  //windychien::knot_spec knot{.name = "10_137"};
  //windychien::knot_spec knot{.name = "K15a41275"};
  //windychien::knot_spec knot{.name = "K11n170"};
  //windychien::knot_spec knot{.name = "K11n101"};
  //windychien::knot_spec knot{.name = "K11n19"};

  //windychien::knot_spec knot{.name = "K11a359"};
  windychien::knot_spec knot{.name = "K11n42"};
  //windychien::knot_spec knot{.name = "random", .strands = 3, .frames = 97};

  //windychien::knot_spec knot{.name = "plain_weave", .strands = 10, .frames = 25};
  vec3 rod_color = vec3(0.95, 0.55, 0.15);
  vec3 planar_color = vec3(0.35, 0.75, 1.0);
  real rod_radius = 0.015;
  bool show_planar_graph = true;

  /// Hepworth stretch. 0 while the Sobolev saddle holds edge length.
  real edge_stretch_w = 0.0;
  /// >0: quaternion block + stretch_shear + four-frame straight. 0: edge_stretch only.
  real straight_w = 1.0;
  real smooth_w = 1.0e-4;
  real min_kink_w = 0.0e-3;
  real collision_w = 1.0;
  real dt = 0.5;
  real damping = 0.05;
  int iterations = 20;

  tangent_point_solver_config tp{
      .type = {tangent_point_type::sobolev},
      //.type = {tangent_point_type::cauchy, tangent_point_type::sobolev},
      //.type = {tangent_point_type::cauchy, tangent_point_type::harmonic},

      .regularized = {.w = 1.0, .l0 = 1.0, .p = 3.0},
      .soft = {.w = 1.0, .p = 6.0, .R_min_frac = 0.1, .tau_frac = 1.0},
      .cauchy = {.w = 1e-6, .l0 = 0.1, .p = 8.0},
      .sobolev = {.w = 1.0,
                  .l0 = 0.001,
                  .p0 = 6.0,
                  .p1 = 1.7,
                  .alpha_B = 1.0,
                  .alpha_0 = 1.0,
                  .num_v_cycles = 4,
                  #if 1
                  .solver = {.type = sobolev_solver::cg,
                             .iter = 12},
                  #else
                  .solver = {.type = sobolev_solver::chebyshev,
                             .iter = 8,
                             .chebyshev_lower = 0.3},
                  #endif
                  .dense_cg_threshold = 32,
                  .max_levels = 12,
                  .bh_eps = 0.25,
                  .cone_eps = 0.25,
                  .coulomb_w = 0.0,
                  .coulomb_p = 3,
                  .vortex_w = 0.0,
                  .vortex_p = 3,
                  .vortex_q = 1},
      .harmonic = {.w = 1.0, .l0 = 0.1, .l1 = 24.0, .p0 = 6.0, .p1 = 2.0},
      .line_search = {.trust_frac_lavg = 1.0,
                      .shrink = 0.5,
                      .max_backtracks = 32,
                      .max_refine = 8,
                      .verbose = true},
      .ccd = {.max_substeps = 32,
              .neighborhood = {.cfl_frac = 1.0,
                               .agg = cfl_agg_mode::hard_min,
                               .gate_frac = 0.1,
                               .use_envelope = true,
                               .chain_skip_hops = 8}},
      .hard_constraints = {{.type = sobolev_constraint_type::edge_length},
                           {.type = sobolev_constraint_type::barycenter}},
      .proj_constraints = { /*{.type = sobolev_constraint_type::c2, .w = 1.0}, */
                            {.type = sobolev_constraint_type::c4, .w = 5.0},
                            {.type = sobolev_constraint_type::c1, .w = 5.0},
                           /*{.type = sobolev_constraint_type::elastic, .w = 0.5}*/},
      .constraint_cg = 12};

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
    open_demo_log("dump/braid_demo.log");
    const windychien::braid b = windychien::make_knot(_cfg.knot);

    braid_planar_params planar{};
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

    const braid_sphere_params sphere{.n_lat = 256};
    __R = braid_to_sphere_rod(b, planar, sphere);
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
      const vec3 shift(0.0, 0.0, -sphere.radius - 0.75);
      for (vec3 &q : __R_planar->x()) {
        q[0] -= cen[0];
        q[1] -= cen[1];
        q += shift;
      }
      __R_planar->_init_params();
    }

    const real lavg = std::max(__R->lavg(), real(1e-6));
    real C = 2.0;
    __Rd = asawa::rod::dynamic::create(__R, C * lavg, 2.0 * C * lavg,
                                       0.25 * C * lavg);

    _tp_constraints = make_sobolev_constraint_solver(_cfg.tp);
    _tp_constraints.init_rest(*__R);

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

    auto run_solver = [this](real dti) { solve_rod(dti); };

    auto log_hepworth = [this](int frame, int sub,
                               const std::vector<vec3> &x0) {
      if (!_cfg.ccd_verbose || !__R)
        return;
      const ccd_field_stats xs = ccd_field_stats_vec(__R->x());
      real peak_dx = 0.0;
      const auto &xnow = __R->x();
      const size_t n = std::min(x0.size(), xnow.size());
      for (size_t i = 0; i < n; ++i) {
        const real d = (xnow[i] - x0[i]).norm();
        if (std::isfinite(d))
          peak_dx = std::max(peak_dx, d);
      }
      demo_log() << "[braid step] frame=" << frame << " sub=" << sub
                << " |x|_max=" << xs.peak << " rms=" << xs.rms
                << " peak|dx|=" << peak_dx << std::endl;
    };

    auto log_ccd_diag = [this](const char *tag, int frame, int sub,
                               real dti) {
      if (!_cfg.ccd_verbose)
        return;
      const ccd_field_stats &d = _ccd_breakdown.drive;
      const ccd_field_stats &x = _ccd_breakdown.rod_x;
      demo_log() << "[braid CCD] " << tag << " frame=" << frame << " sub=" << sub
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
      std::vector<vec3> x_bite;
      if (_cfg.ccd_verbose)
        x_bite = __R->x();
      flow_then_project(_cfg.dt);
      log_hepworth(frame, 0, x_bite);
      __Rd->step();
      __R->_update_frames();
      return;
    }

    demo_log() << "[braid CCD] begin frame=" << frame
              << " schedule=tau_then_search" << std::endl;

    real alpha_star = 0.0;
    real alpha_left = 0.0;
    real E0_frame = 0.0;
    real E_star = 0.0;
    real tau_max = 0.0;
    real propose_ms = 0.0;
    int n_backtracks = 0;
    bool used_fallback = false;
    std::vector<vec3> frozen_disp;
    bool have_frozen = false;

    {
      tp_curvature_propose prop =
          tp_propose_curvature_step(*__R, *__Rd, _cfg.tp, &_tp_constraints);
      frozen_disp = std::move(prop.displacement);
      const tp_line_search_result &ls = prop.ls;
      alpha_star = ls.alpha;
      E0_frame = ls.E0;
      E_star = ls.E_accepted;
      tau_max = prop.tau_max;
      propose_ms = prop.propose_ms;
      n_backtracks = ls.n_backtracks;
      used_fallback = ls.used_fallback;
      alpha_left = ccd_peak_norm(frozen_disp);
      have_frozen = true;
      if (_cfg.ccd_verbose) {
        demo_log() << "[braid α*] frame=" << frame << " tau_max=" << prop.tau_max
                  << " alpha0=" << ls.alpha_quad
                  << " alpha=" << alpha_star
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
         &log_ccd_diag, &log_hepworth, frame, &alpha_left, &frozen_disp,
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
            ccd_config probe = _ccd;
            if (have_frozen) {
              // δ is already inside τ_max. Re-query geometry only, so
              // dt_max / dt_min cannot chop or step through that displacement.
              probe.dt_max = 0.0;
              probe.dt_min = 0.0;
              probe.dx_max = 0.0;
              probe.neighborhood.dt_max = 0.0;
              probe.neighborhood.dt_min = 0.0;
            }
            dti = ccd_velocity(*__R, drive, probe, &_ccd_stats.last,
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
              demo_log() << "[braid α*] abort null bite dti=" << dti
                        << " dt_geom=" << _ccd_breakdown.dt_geom
                        << " bite_peak=" << peak_disp << std::endl;
            }
            return 0.0;
          }

          std::vector<vec3> x_bite;
          if (_cfg.ccd_verbose)
            x_bite = __R->x();
          const real dti_used = flow_then_project(dti);
          if (_cfg.use_velocity_block)
            peak_disp = dti_used * ccd_peak_norm(drive);
          else
            peak_disp = (dti_used * dti_used) * ccd_peak_norm(drive);
          log_hepworth(frame, sub_i, x_bite);

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
              demo_log() << "[braid α*] alpha_left=" << alpha_left
                        << " bite_peak=" << peak_disp << std::endl;
            }
          }

          ccd_diagnose_rod(*__R, &_ccd_breakdown.rod_x);
          log_ccd_diag("post", frame, sub_i, dti_used);
          ++sub_i;

          if (!(peak_disp > _ccd.eps))
            return 0.0;
          return peak_disp;
        },
        _ccd);

    demo_log() << "[braid CCD] frame=" << frame
              << " substeps=" << _ccd_stats.substeps
              << " accum=" << _ccd_stats.dt_accum << "/" << target
              << " dti=" << _ccd_stats.last_dti << " alpha*=" << alpha_star
              << " alpha_left=" << alpha_left << std::endl;

    if (_cfg.ccd_verbose && have_frozen) {
      const real E_tp_post = tp_energy(*__R, *__Rd, _cfg.tp);
      const real dE_frame = E_tp_post - E0_frame;
      const real dE_prev =
          _have_E_prev_post ? (E_tp_post - _E_prev_post) : real(0.0);
      demo_log() << "[braid E] frame=" << frame << " E0=" << E0_frame
                << " E*=" << E_star << " Epost=" << E_tp_post
                << " dE_ray=" << (E_star - E0_frame)
                << " dE_frame=" << dE_frame << " dE_prev=" << dE_prev
                << (dE_frame > 0.0 ? " ASCENT" : "") << std::endl;
      _E_prev_post = E_tp_post;
      _have_E_prev_post = true;
    }
    if (have_frozen) {
      const real dE = E_star - E0_frame;
      std::cout << "[braid] frame=" << frame << " n=" << __R->corner_count()
                << " E0=" << E0_frame << " E*=" << E_star << " dE=" << dE
                << " alpha=" << alpha_star << " tau_max=" << tau_max
                << " fallback=" << (used_fallback ? 1 : 0)
                << " backtracks=" << n_backtracks
                << " propose_ms=" << propose_ms
                << (used_fallback ? " STALL" : "") << std::endl;
    }
    __Rd->step();
    __R->_update_frames();
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
  sobolev_constraint_solver _tp_constraints;
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
  void solve_rod(real dti) {
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
  }

  real project_sobolev_constraints() {
    if (_tp_constraints.empty() || !__R || !__Rd)
      return real(0);
    const real l0 = _cfg.tp.sobolev.l0 * __Rd->_Cc;
    return _tp_constraints.project_positions(*__R, _cfg.tp.sobolev, l0,
                                             _cfg.ccd_verbose);
  }

  /// §5.3.2 line search: keep the flow step only when the Newton projection
  /// brings max |Φ| under 10^{-4}. Otherwise halve it and retake from the
  /// saved pose. Ten attempts matches LSBackprojectMultigrid.
  real flow_then_project(real dti) {
    if (!__R)
      return real(0);
    const real tol = sobolev_constraint_solver::backproj_tol;
    const std::vector<vec3> x_pre = __R->x();
    const std::vector<vec3> v_pre = __R->v();
    real dt_try = dti;
    real viol = real(0);
    int attempt = 0;
    for (; attempt < 10; ++attempt) {
      __R->x() = x_pre;
      __R->v() = v_pre;
      solve_rod(dt_try);
      __R->_update_frames();
      viol = project_sobolev_constraints();
      if (viol < tol)
        break;
      const real next = dt_try * real(0.5);
      if (!(next > real(0)) || !(next < dt_try))
        break;
      dt_try = next;
    }
    if (!(viol < tol)) {
      __R->x() = x_pre;
      __R->v() = v_pre;
      __R->_update_frames();
      viol = project_sobolev_constraints();
      dt_try = real(0);
    }
    if (_cfg.ccd_verbose) {
      demo_log() << "[sobolev backproj] attempts=" << (attempt + 1)
                << " dti=" << dt_try << " |Phi|_max=" << viol << std::endl;
    }
    return dt_try;
  }

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
