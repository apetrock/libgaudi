#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/hepworth/shell/constraints.hpp"
#include "gaudi/hepworth/shell/constraints_init.hpp"
#include "gaudi/hepworth/shell/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

class shell_constraints_test {
public:
  typedef std::shared_ptr<shell_constraints_test> ptr;

  static ptr create() { return std::make_shared<shell_constraints_test>(); }

  shell_constraints_test() {
    //__M = load_cube();
    __M = shell::load_bunny();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    /////////
    // dynamic surface
    /////////

    real l0 = 1.0 * asawa::shell::avg_length(*__M, x);
    __surf = shell::dynamic::create(__M, 1.0 * l0, 3.0 * l0, 1.0 * l0);

    /////////
    // constraints setup
    /////////
    std::vector<real> lr = asawa::shell::edge_lengths(*__M, x);

    datum_t<real>::ptr ldata = datum_t<real>::create(prim_type::EDGE, lr);
    _il0 = __M->insert_datum(ldata);
  };

  void smoothMesh(real C, int N) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    bontecou::laplacian3 M(__M, x, true);

    for (int k = 0; k < N; k++) {
      std::cout << "." << std::flush;
      M.init();
      x = M.smooth(x, C, C + 1e-6);
      int i = 0;
      for (auto xi : x) {
        if (!std::isfinite(xi.dot(xi))) {
          std::cout << xi.transpose() << std::endl;
          __M->vprintv(i);
          i++;
        }
      }
    }
    // x_datum->data() = x;
    std::cout << "done!" << std::endl;
  }

  void step_dynamics(int frame) {

    hepworth::shell::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> m = asawa::shell::vertex_areas(*__M, x);
    solver.set_mass(m);

    int i = 0;
    for (auto &l : l0)
      l *= 1.0 + 0.01;

    std::cout << "l[0]: " << l0[0] << std::endl;
    hepworth::shell::init_edge_strain(*__M, constraints, x, l0, 1.02, 1.0e-4);
    hepworth::shell::init_pinned(*__M, constraints, x, 1.0e-0);
    hepworth::shell::init_bending(*__M, constraints, x, 1.0e-3);
    hepworth::shell::init_laplacian(*__M, constraints, x, 3.0e-1, true);
    hepworth::shell::init_area(*__M, constraints, x, 5.0e-3);

    //  hepworth::shell::init_cross(*__M, constraints, 1.05, 0.1);
    hepworth::shell::init_edge_edge_collisions(*__M, *__surf, constraints, x,
                                               1.0 * __surf->_Cm, 1.2);
    hepworth::shell::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             1.0 * __surf->_Cm, 1.2);
    solver.set_constraints(constraints);

    real bnd = 2.0;

    std::vector<vec3> f(x.size(), vec3::Zero());
    for (int i = 0; i < x.size(); i++) {
      vec3 xi = x[i]; // points to origin
      real xn = xi.norm();
      if (xn > bnd) {
        f[i] = -0.1 * (xn - bnd) / _h / _h * xi;
      }
      // f[i] += 1e-2 * vec3::Random();
    }
    // f[0][1] = 10.0;
    solver.step(x, v, f, _h);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_dynamics(frame);
    __surf->step(false);
    smoothMesh(0.001, 1);
  }
  real _h = 0.1;
  index_t _il0;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif