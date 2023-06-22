#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"

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

  std::vector<vec3> parallel_transport(shell::shell &M,
                                       const std::vector<vec3> &x,
                                       const std::vector<vec3> &Nx, real l0,
                                       real thet) {

    std::vector<mat3> Us = calder::fast_frame(M, x, x, Nx, l0);
    std::vector<vec3> v(Us.size());
    for (int i = 0; i < Us.size(); i++) {

      // U.array().rowwise() *= s.transpose().array();
      mat3 U = Us[i];
      real cx = cos(thet), cy = sin(thet);
      vec3 pf = U * vec3(cx, cy, 0.0);
      v[i] = pf;
      gg::geometry_logger::line(x[i],             //
                                x[i] + 0.01 * pf, //
                                vec4(0.0, 0.7, 1.0, 1.0));
    }
    return v;
  }

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    std::vector<vec3> Nx = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> M = asawa::shell::vertex_areas(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);

    real al0 = 1.0 * asawa::shell::avg_length(*__M, x);

    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, x);
    }

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());
    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M, x, v, f);

    std::cout << "l[0]: " << l0[0] << std::endl;
    // hepworth::block::init_edge_growth(*__M, constraints, x, l0, 1.01, 2.0e-5,
    // {X});
    // hepworth::block::init_edge_strain(*__M, constraints, x, l0, 3e-2, {X});
    // hepworth::block::init_pinned(*__M, constraints, x, 7.0e-3, {X});
    // hepworth::block::init_bending(*__M, constraints, x, 1.0e-2, {X});
    hepworth::block::init_vec_align(*__M, constraints, x, 10.0e-0, {X});

    // hepworth::block::init_laplacian(*__M, constraints, x, 2, 4.0e-2, {X});
    // hepworth::block::init_area(*__M, constraints, x, 3.0e-2, {X});

    //  hepworth::shell::init_cross(*__M, constraints, 1.05, 0.1);
    real eps = 3.0 * __surf->_Cm;
    // hepworth::shell::init_edge_edge_collisions(*__M, *__surf, constraints,
    // x,
    //                                            eps, 0.5 * eps, 1.0e0);
    hepworth::block::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             1.0 * eps, 0.5 * eps, 1.0, {X, X});
    solver.set_constraints(constraints);

    real bnd = 2.0;

    // f[0][1] = 10.0;
    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 30);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;

    // step_dynamics(frame);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> Nx = asawa::shell::vertex_normals(*__M, x);

    real al0 = 1.0 * asawa::shell::avg_length(*__M, x);
    std::vector<vec3> v_p =
        parallel_transport(*__M, x, Nx, 4.0 * al0, 0.1 * M_PI);

    __surf->step(false);

    std::cout << "print_walk" << std::endl;
  }
  real _h = 0.1;
  index_t _il0;
  index_t _iwalk;
  std::vector<index_t> _walk;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif