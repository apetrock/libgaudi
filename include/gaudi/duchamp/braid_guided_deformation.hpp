#ifndef GAUDI_DUCHAMP_BRAID_GUIDED_DEFORMATION_HPP
#define GAUDI_DUCHAMP_BRAID_GUIDED_DEFORMATION_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/braid_circle_bundle.hpp"
#include "gaudi/duchamp/modules/knotted_surface.hpp"
#include "gaudi/duchamp/modules/tp_constraint_solver.hpp"
#include "gaudi/duchamp/modules/tp_energies.hpp"

namespace gaudi {
namespace duchamp {

using namespace asawa;

/// Which field of `tangent_point_solver_config` is the rod velocity this frame.
enum class braid_rod_drive {
  cauchy,  ///< `tp.cauchy` gradient
  sobolev, ///< `tp.sobolev` gradient (constraint V-cycle when lists are set)
  ccd,     ///< accepted `tp.ccd` / line-search displacement
};

struct braid_guided_config {
  braid_circle_config braid{
    /*.knot = {.name = "plain_weave", .strands = 25, .frames = 25}, */
    //.knot = {.name = "K11a280"},
    .knot = {.name = "K11a333"},

    .sphere = {.lat_min = -1.4, .lat_max = 1.4, .n_lat = 256},

  };
  /// Dipole tunnel/weld cylinder radius = `dipole_radius_scale * rod->_r`.
  /// Applied after `knot`, so this wins over `knot.dipole_radius`.
  real dipole_radius_scale = 3.0;
  /// Frame step. `knotted_surface::step` uses this as the Hepworth dt (`x += h v`).
  real h = 0.5;
  /// Hepworth weights, edge-edge gap, and solver. Order matches `knotted_surface_config`.
  knotted_surface_config knot{
      // shell
      .w_willmore = 0.0,
      .w_edge_willmore = 0.0,
      .w_smooth = 0.0,
      .w_area = 0.0,
      .w_shell_strain = 1.0e-2,
      .w_shell_bending = 1.0e-2,
      // rod
      .w_rod_min_kink = 1.0e-1,
      .w_rod_strain = 1.0e-5,
      .w_rod_bend = 1.0e-3,
      .w_rod_twist = 1.0e-3,
      // weld, dipole, tunnel
      .w_rod_weld = 1.0,
      .w_shell_weld = 1.0,
      .w_tunnel_orientation = 0.01,
      .w_dipole_weld = 0.15,
      .w_dipole_weld_rod = 0.01,
      .w_tunnel_force = 0.75,
      /*
      .w_darboux_force = 0.0,
      .darboux_fit_scale = 0.1,
      .darboux_fit_p0 = 3.0,
      .darboux_fit_p1 = 8.0,
      .darboux_foot_normal_w = 1.0,
      .darboux_smooth = false,
      */
      // pin, edge-edge gap, pose blend, solver
      .w_rod_pin = 5.0,
      .rod_offset_min = 2.0, // enforced centerline gap >= rod_offset_min * rod._r
      .repel_rods = true,   // no rod–rod edge–edge; welds may release
      .position_mix = 1.0,
      .damping = 0.1,
      .iterations = 30,

  };
  braid_rod_drive drive = braid_rod_drive::ccd;
  /// Fraction of the accepted CCD displacement pinned this frame.
  /// 1 keeps the full step. 0 pins the rod where it already is.
  real ccd_fraction = 0.5;
  tangent_point_solver_config tp{
      .type = {tangent_point_type::sobolev},
      .cauchy = {.w = 1e-6, .l0 = 0.1, .p = 8.0},
      .sobolev = {.w = 1.0,
                  .l0 = 0.001,
                  .p0 = 6.0,
                  .p1 = 1.7,
                  .alpha_B = 1.0,
                  .alpha_0 = 1.0,
                  .solver = {.type = sobolev_solver::cg, .iter = 12},
                  .dense_cg_threshold = 32,
                  .num_v_cycles = 4,
                  .max_levels = 12,
                  .bh_eps = 0.25,
                  .cone_eps = 0.25,
                  .coulomb_w = 0.0,
                  .coulomb_p = 3,
                  .vortex_w = 0.0,
                  .vortex_p = 3,
                  .vortex_q = 1},
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
      .proj_constraints = {{.type = sobolev_constraint_type::c4, .w = 5.0},
                           {.type = sobolev_constraint_type::c1, .w = 5.0}},
      .constraint_cg = 12};
  /// Shell quaternion triad after remesh. Red tangent, green bitangent, blue normal.
  bool draw_shell_frames = true;
  /// Draw one triad per this many valid vertices.
  int shell_frame_stride = 1;
};

inline const braid_guided_config k_braid_guided_default{};

/// Sphere shell + braid knot. Rod load is one velocity per frame.
class braid_guided_deformation {
public:
  using ptr = std::shared_ptr<braid_guided_deformation>;

  static ptr create(const braid_guided_config &cfg = k_braid_guided_default) {
    return std::make_shared<braid_guided_deformation>(cfg);
  }

  static ptr create(const braid_circle_config &bundle_cfg) {
    braid_guided_config cfg;
    cfg.braid = bundle_cfg;
    return create(cfg);
  }

  explicit braid_guided_deformation(
      const braid_guided_config &cfg = k_braid_guided_default)
      : _cfg(cfg) {
    const braid_circle_bundle bundle = make_braid_circle_bundle(_cfg.braid);
    __M = bundle.shell;
    __R = bundle.rod;

    shell::triangulate(*__M);

    const std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    const real l0 = asawa::shell::avg_length(*__M, x);

    const real Ct = 2.0;
    const real C = 3.0 * Ct;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);

    const real lavg = std::max(__R->lavg(), real(1e-6));
    __Rd = rod::dynamic::create(__R, 0.5 * Ct * lavg, 2.0 * Ct * lavg, 0.25 * Ct * lavg);

    for (int i = 0; i < 5; i++) {
      __surf->step();
      __Rd->step();
    }

    _tp_constraints = make_sobolev_constraint_solver(_cfg.tp);
    _tp_constraints.init_rest(*__R);

    _knotted_surface = knotted_surface_module::create(__M, __surf, __R, __Rd);
    _knotted_surface->set_config(_cfg.knot);
    const real dipole_r = _cfg.dipole_radius_scale * __R->_r;
    _knotted_surface->set_dipole_radius(dipole_r);
  }

  const braid_guided_config &config() const { return _cfg; }

  void step(int /*frame*/) {
    _knotted_surface->init_step(_cfg.h);
    if (_cfg.drive == braid_rod_drive::ccd) {
      // Pin the rod to the CCD pose. The weld carries the shell; no velocity.
      const tp_curvature_propose prop = tp_propose_curvature_step(
          *__R, *__Rd, _cfg.tp, constraints());
      std::vector<vec3> dx = prop.displacement;
      const real f = _cfg.ccd_fraction;
      for (vec3 &d : dx)
        d *= f;
      _knotted_surface->pin_rod_to_delta(dx);
    } else {
      _knotted_surface->add_rod_velocity(rod_velocity());
    }
    _knotted_surface->step(_cfg.h);
    if (_cfg.drive != braid_rod_drive::cauchy && !_tp_constraints.empty()) {
      __R->_update_frames();
      const real l0 = _cfg.tp.sobolev.l0 * __Rd->_Cc;
      _tp_constraints.project_positions(*__R, _cfg.tp.sobolev, l0);
    }
    __surf->step(true);
    __Rd->step();
  }

  void draw_shell_frames() {
    if (_cfg.draw_shell_frames && _knotted_surface)
      _knotted_surface->draw_shell_frames(_cfg.shell_frame_stride);
  }

  knotted_surface_module::ptr _knotted_surface;
  braid_guided_config _cfg;
  sobolev_constraint_solver _tp_constraints;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;

private:
  sobolev_constraint_solver *constraints() {
    return _tp_constraints.empty() ? nullptr : &_tp_constraints;
  }

  static void negate(std::vector<vec3> &g) {
    for (vec3 &a : g)
      a *= real(-1);
  }

  void peak_scale(std::vector<vec3> &v) const {
    if (!(_cfg.tp.peak_scale > real(0)))
      return;
    real peak = 0.0;
    for (const vec3 &a : v) {
      const real n = a.norm();
      if (std::isfinite(n))
        peak = std::max(peak, n);
    }
    if (!(peak > real(0)) || !std::isfinite(peak))
      return;
    const real alpha = real(1) / (_cfg.tp.peak_scale * peak);
    for (vec3 &a : v)
      a *= alpha;
  }

  std::vector<vec3> rod_velocity() {
    if (_cfg.drive == braid_rod_drive::ccd) {
      tp_curvature_propose prop = tp_propose_curvature_step(
          *__R, *__Rd, _cfg.tp, constraints());
      return std::move(prop.displacement);
    }
    const tangent_point_type type = _cfg.drive == braid_rod_drive::cauchy
                                        ? tangent_point_type::cauchy
                                        : tangent_point_type::sobolev;
    std::vector<vec3> v =
        tp_grad_energy_one(*__R, *__Rd, _cfg.tp, type, constraints());
    negate(v);
    peak_scale(v);
    return v;
  }
};

} // namespace duchamp
} // namespace gaudi

#endif
