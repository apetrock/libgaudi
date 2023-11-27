#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/calder/quadric_fit.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include "modules/knotted_surface.hpp"
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
#include <zlib.h>

#ifndef __ROD_GUIDED__
#define __ROD_GUIDED__

namespace gaudi {
namespace duchamp {

using namespace asawa;

///////////////////////////////////////
class rod_guided_deformation {
public:
  typedef std::shared_ptr<rod_guided_deformation> ptr;

  static ptr create() { return std::make_shared<rod_guided_deformation>(); }

  rod_guided_deformation() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    /////////
    // dynamic surface
    /////////
    real l0 = asawa::shell::avg_length(*__M, x);
    real C = 0.6;
    // real C = 2.0;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);

    /////////////////////
    // Rod
    /////////////////////
    std::vector<vec3> x_w = walk(*__M, 0.0, 0, 4000);

    __R = rod::rod::create(x_w, false);
    //__R->_update_frames(normals);

    real lavg = 1.0 * l0;
    __R->_r = 0.020;
    __Rd = rod::dynamic::create(__R, 0.35 * lavg, 2.0 * lavg, 0.25 * lavg);

    for (int i = 0; i < 5; i++) {
      __surf->step();
      __Rd->step();
    }

    _knotted_surface = knotted_surface_module::create(__M, __surf, __R, __Rd);
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
      gg::geometry_logger::line(x_s[i], x_s[i] + 1.0 * N,
                                vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    return Nss;
  }

#if 1
  std::vector<vec3> compute_tangent_point_gradient() {

    real eps = _knotted_surface->get_eps();
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_gradient(*__R, x, l, T, 3.0 * eps, 6.0);
    for (int i = 0; i < g0.size(); i++) {
      // gg::geometry_logger::line(x[i], x[i] + 1.0e-7 * g0[i],
      //                           vec4(0.6, 0.0, 0.8, 1.0));
      g0[i] *= -2.0e-6;
    }
    return g0;
  }
#endif

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    std::cout << "  -curve" << std::endl;
    std::cout << "    -verts: " << __R->x().size() << std::endl;

    // walk(__surf->_Cc);
    if (frame < 1200) {
      _knotted_surface->init_step(_h);
      _knotted_surface->add_rod_force(compute_tangent_point_gradient());
      _knotted_surface->add_shell_force(calc_quadric_grad());
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

  real _h = 0.05;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif