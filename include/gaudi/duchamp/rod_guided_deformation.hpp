#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

#include "gaudi/geometry_logger.hpp"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell_id.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/calder/least_squares_fit.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include "braid_circle_bundle.hpp"
#include "modules/knotted_surface.hpp"
#include "modules/rod_forces.hpp"
#include "utils/point_things.hpp"
#include "utils/string_things.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>

#ifndef __ROD_GUIDED__
#define __ROD_GUIDED__

namespace gaudi {
namespace duchamp {

using namespace asawa;

enum class rod_guided_scene {
  braid_circle, // closed braid on matching sphere shell
  bunny_walk,   // bunny + silly_walk rod (legacy)
  sphere_walk,  // sphere + silly_walk rod
};

struct rod_guided_config {
  rod_guided_scene scene = rod_guided_scene::braid_circle;
  braid_circle_config braid{
    .use_plain_weave = true,
    .weave_strands =25,
    .weave_frames = 25
  };
  silly_walk_config walk{
      .i0 = 0,
      .N_steps = 10000,
      .thet = M_PI / 2.0,
      .rotate = true,
      .twist_amp = 1.00,
      .twist_freq = 0.6,
      .align = true,
      .ca = vec3(0.05, 0.08, 2.3),
  };
  real rod_radius = 0.02;
  int sphere_n = 192; // shell tessellation for sphere_walk
  /// Dipole tunnel/weld cylinder radius = `dipole_radius_scale * rod->_r`.
  real dipole_radius_scale = 3.0;
  /// Rod tangent-point force. Set `w = 0` to disable.
  /// `l0` multiplies knotted-surface eps; `p` is the TP power.
  tangent_point_force_config tangent{.w = 0.0e-9, .l0 = 3.0, .p = 6.0};
};

///////////////////////////////////////
class rod_guided_deformation {
public:
  typedef std::shared_ptr<rod_guided_deformation> ptr;

  static ptr create(const rod_guided_config &cfg = {}) {
    return std::make_shared<rod_guided_deformation>(cfg);
  }

  // Convenience: braid+circle bundle only.
  static ptr create(const braid_circle_config &bundle_cfg) {
    rod_guided_config cfg;
    cfg.scene = rod_guided_scene::braid_circle;
    cfg.braid = bundle_cfg;
    return create(cfg);
  }

  explicit rod_guided_deformation(const rod_guided_config &cfg = {})
      : _tangent(cfg.tangent), _dipole_radius_scale(cfg.dipole_radius_scale) {
    if (cfg.scene == rod_guided_scene::braid_circle) {
      const braid_circle_bundle bundle = make_braid_circle_bundle(cfg.braid);
      __M = bundle.shell;
      __R = bundle.rod;
    } else if (cfg.scene == rod_guided_scene::bunny_walk) {
      __M = shell::load_bunny();
    } else {
      __M = shell::load_sphere(1.0, cfg.sphere_n, cfg.sphere_n / 2);
    }

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(asawa::shell::face_id(i)) > asawa::shell::corner_id(0)) {
        assert(__M->fsize(asawa::shell::face_id(i)) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    if (cfg.scene != rod_guided_scene::braid_circle) {
      // Walk scenes: normalize bounds, then walk on the normalized mesh.
      asawa::center(x, 2.0);
      const std::vector<vec3> x_w = cfg.walk.run(*__M);
      __R = rod::rod::create(x_w, false);
      __R->_r = cfg.rod_radius;
    }
    // Braid+circle: keep matched sphere radius (no recenter/rescale).

    const real l0 = asawa::shell::avg_length(*__M, x);
    const real C = 3.0;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);

    const real lavg = std::max(__R->lavg(), real(1e-6));
    __Rd = rod::dynamic::create(__R, 0.35 * lavg, 2.0 * lavg, 0.25 * lavg);

    for (int i = 0; i < 5; i++) {
      __surf->step();
      __Rd->step();
    }

    _knotted_surface = knotted_surface_module::create(__M, __surf, __R, __Rd);
    // Absolute dipole cylinder; knotted_surface treats <=0 as "unset" (legacy).
    const real dipole_r = _dipole_radius_scale * __R->_r;
    _knotted_surface->set_dipole_radius(dipole_r);
  };

  std::vector<vec3> calc_quadric_grad() {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<vec3> x_s_f = asawa::shell::face_centers(M, x_s);

    std::vector<vec3> N_s_f = asawa::shell::face_normals(*__M, x_s);
    const std::vector<vec3> &x_r = R.x();

    real eps = _knotted_surface->get_eps();
    std::vector<vec3> Nr = _knotted_surface->get_rod_normals(R, M, 1.0 * eps);

    std::vector<real> Q = calder::quadric_sdf(R, Nr, x_s_f, N_s_f, 0.25 * eps);

#if 1
    int i = 0;
    for (vec3 &N : N_s_f) {
      N *= -5.0 * Q[i];
      // gg::geometry_logger::line(x_s_f[i], x_s_f[i] + 1.0 * N,
      //                           vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    std::vector<vec3> Nss =
        calder::mls_avg<vec3>(*__M, N_s_f, x_s, 1.0 * eps, 2.0);
#if 0
    i = 0;
    for (vec3 &N : Nss) {
      geometry_logger::line(x_s[i], x_s[i] + 1.0 * N,
                                vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    return Nss;
  }

  std::vector<vec3> compute_tangent_point_gradient() {
    real eps = _knotted_surface->get_eps();
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();

    std::vector<vec3> g0 = calder::tangent_point_gradient(
        *__R, x, l, T, _tangent.l0 * eps, _tangent.p);
    for (vec3 &g : g0)
      g *= _tangent.w;
    return g0;
  }

  void step(int frame) {
    _knotted_surface->set_rod_offset(1.0 + 0.02 * real(frame));
    // walk(__surf->_Cc);
    if (frame < 1200) {
      _knotted_surface->init_step(_h);
      if (_tangent.w != 0.0)
        _knotted_surface->add_rod_force(compute_tangent_point_gradient());
      // _knotted_surface->add_shell_force(calc_quadric_grad());
      _knotted_surface->step(_h);
    }

    if (frame > 1200)
      exit(0);

    for (int k = 0; k < 1; k++) {
      __surf->step(true);
      __Rd->step();
    }
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  knotted_surface_module::ptr _knotted_surface;
  tangent_point_force_config _tangent;
  real _dipole_radius_scale = 2.0;

  real _h = 0.05;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif
