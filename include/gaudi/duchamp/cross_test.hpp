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

#include "modules/cross.hpp"
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

class cross_test {
public:
  typedef std::shared_ptr<cross_test> ptr;

  static ptr create() { return std::make_shared<cross_test>(); }

  cross_test() { load_shell(); };

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
    // real C = 0.5;
    real C = 1.5;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.1 * C * l0);

    _cross = std::dynamic_pointer_cast<module_base>(cross::create(__M));

    for (int i = 0; i < 5; i++) {
      __surf->step();
    }

    _eps = asawa::shell::avg_length(*__M, x);
  }

  void step_dynamics(int frame) {

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);

#if 1
    real mollifier = 8.0 * _eps;
    std::dynamic_pointer_cast<cross>(_cross)->step(_h);
    std::vector<vec3> X =
        std::dynamic_pointer_cast<cross>(_cross)->get_smoothd_cross_grad();

    for (int i = 0; i < X.size(); i++) {
      x[i] += 1e-7 * _h * X[i];
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

  module_base::ptr _cross;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif