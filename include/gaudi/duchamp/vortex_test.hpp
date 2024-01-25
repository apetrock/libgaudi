#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/duchamp/modules/vortex.hpp"
#include "modules/module_base.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __TEST_TEMPLATE__
#define __TEST_TEMPLATE__

namespace gaudi {
namespace duchamp {
using namespace asawa;

class vortex_test {
public:
  typedef std::shared_ptr<vortex_test> ptr;

  static ptr create() { return std::make_shared<vortex_test>(); }

  vortex_test() { load_shell(); };

  void load_shell() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    real l0 = asawa::shell::avg_length(*__M, x);

    /////////
    // dynamic surface
    /////////
    real C = 0.75;
    // real C = 1.5;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.5 * C * l0);

    _vortex_edge = std::dynamic_pointer_cast<module_base>( //
        vortex_edge::create(__M, 1e-5));

    _vortex_vert = std::dynamic_pointer_cast<module_base>( //
        vortex_vert::create(__M, 1e-5, 5e-2, 1e-4));

    for (int i = 0; i < 5; i++) {
      __surf->step();
    }

    _eps = asawa::shell::avg_length(*__M, x);
  }

  void step_dynamics(int frame) {

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
#if 0
    real mollifier0 = 8.0 * _eps;
    std::dynamic_pointer_cast<vortex_edge>(_vortex_edge)->step(_h);
    std::vector<vec3> w0 =
        std::dynamic_pointer_cast<vortex_edge>(_vortex_edge)->get_vorticity();
    std::cout << "calder" << std::endl;
    std::vector<vec3> f0 = calder::vortex_force(M, x, w0, mollifier0, 3.0);
#endif
#if 1
    real mollifier1 = 18.0 * _eps;
    std::dynamic_pointer_cast<vortex_vert>(_vortex_vert)->step(_h);
    std::vector<vec3> w1 =
        std::dynamic_pointer_cast<vortex_vert>(_vortex_vert)->get_vorticity();
    std::vector<vec3> f1 = calder::vortex_force(M, x, w1, mollifier1, 3.0);
#endif
    // std::vector<vec3> w = calc_edge_vorticity(M, 1e-4 * _h);

#if 1
    // std::vector<vec3> fc = asawa::shell::vert_to_face<vec3>(M, f);
    //  std::vector<vec3> fs = calder::mls_avg<vec3>(M, fc, x, 4.0 * _eps, 3.0);

    for (int i = 0; i < f1.size(); i++) {
      // x[i] += _h * (f0[i] + f1[i]);
      x[i] += _h * f1[i];
    }
#endif
    // asawa::center(x, 2.0);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;

    for (int k = 0; k < 2; k++) {
      __surf->step(true);
    }
    step_dynamics(frame);
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  real _h = 0.05;
  real _eps;
  index_t _ig0, _ig1;

  module_base::ptr _vortex_vert;
  module_base::ptr _vortex_edge;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif