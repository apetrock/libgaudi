#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

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

    __surf->step(true);
  };

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas_3(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);

    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, x);
    }

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());
    for (int i = 0; i < x.size(); i++) {
      vec3 N = asawa::shell::vert_normal(*__M, i, x);

      // f[i] += 1e-3 * N;
      f[i] += 1e-1 * vec3::Random();
    }
    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M, x, v, f);

    std::cout << "l[0]: " << l0[0] << std::endl;
    hepworth::block::init_edge_growth(*__M, constraints, x, l0, 1.01, 2.0e-5,
                                      {X});
    hepworth::block::init_edge_strain(*__M, constraints, x, l0, 3e-2, {X});
    hepworth::block::init_pinned(*__M, constraints, x, 5.0e-3, {X});
    hepworth::block::init_bending(*__M, constraints, x, 3.0e-1, {X});
    hepworth::block::init_laplacian(*__M, constraints, x, 1, 2.0e-1, {X});
    hepworth::block::init_area(*__M, constraints, x, 5.0e-2, {X});

    //  hepworth::shell::init_cross(*__M, constraints, 1.05, 0.1);
    real eps = 3.0 * __surf->_Cm;
    // hepworth::shell::init_edge_edge_collisions(*__M, *__surf, constraints,
    // x,
    //                                            eps, 0.5 * eps, 1.0e0);
    hepworth::block::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             1.0 * eps, 0.5 * eps, 1.0, {X, X});
    solver.set_constraints(constraints);

    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 30);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_dynamics(frame);
    __surf->step(false);
  }

  real _h = 0.1;
  index_t _il0;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif