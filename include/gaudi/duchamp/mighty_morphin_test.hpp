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
#include "gaudi/duchamp/modules/mighty_morphin.hpp"
#include "modules/ccd.hpp"
#include "modules/mighty_morphin.hpp"

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

class mighty_morphin_test {
public:
  typedef std::shared_ptr<mighty_morphin_test> ptr;

  static ptr create() { return std::make_shared<mighty_morphin_test>(); }

  mighty_morphin_test() { load_shell(); };

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

    std::vector<index_t> face_vert_ids = __M->get_face_vert_ids();
    mighty_morphin::ptr mm4 = mighty_morphin::create(__M);
    real off = 0.55;
    real scale = 0.6;
    mm4->add_geometry(x, face_vert_ids, scale, vec3(off, off, 0.0));
    mm4->add_geometry(x, face_vert_ids, scale, vec3(-off, off, 0.0));
    mm4->add_geometry(x, face_vert_ids, scale, vec3(off, -off, 0.0));
    mm4->add_geometry(x, face_vert_ids, scale, vec3(-off, -off, 0.0));
    mighty_morphin::ptr mm1 = mighty_morphin::create(__M);
    mm1->add_geometry(x, face_vert_ids, 1.0, vec3(0, 0, 0.0));

    _mighty_morphin4 = std::dynamic_pointer_cast<module_base>(mm4);
    _mighty_morphin1 = std::dynamic_pointer_cast<module_base>(mm1);

    _ccd = std::dynamic_pointer_cast<module_base>( //
        continuous_collision_detection::create(__M, __surf, 0.1 * l0));

    for (int i = 0; i < 5; i++) {
      __surf->step();
    }

    _eps = asawa::shell::avg_length(*__M, x);
  }

  void step_dynamics(int frame) {

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    mighty_morphin::ptr mm =
        std::dynamic_pointer_cast<mighty_morphin>(_mighty_morphin4);
    if (frame > 400)
      mm = std::dynamic_pointer_cast<mighty_morphin>(_mighty_morphin1);
    // mm->debug_target();

    continuous_collision_detection::ptr ccd =
        std::dynamic_pointer_cast<continuous_collision_detection>(_ccd);

    real h_step = 0.0;
    ccd->init_step(_h);
    mm->step(_h);
    ccd->add_shell_force(mm->forces(), 1.0);
    ccd->step(_h);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    for (int k = 0; k < 1; k++) {
      __surf->step(true);
    }

    step_dynamics(frame);
    if (frame > 600)
      exit(0);

    //  step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  real _h = 0.01;
  real _eps;
  index_t _ig0, _ig1;

  module_base::ptr _mighty_morphin4;
  module_base::ptr _mighty_morphin1;

  module_base::ptr _ccd;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif