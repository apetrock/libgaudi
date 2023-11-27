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

#include "modules/reaction_diffusion.hpp"
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

class growth_study {
public:
  typedef std::shared_ptr<growth_study> ptr;

  static ptr create() { return std::make_shared<growth_study>(); }

  growth_study() {
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

    real l0 = 1.5 * asawa::shell::avg_length(*__M, x);
    __surf = shell::dynamic::create(__M, 1.0 * l0, 2.5 * l0, 1.0 * l0);

    __surf->step(true);
    __surf->step(true);
    __surf->step(true);
    // real f = 0.075, k = 0.0615;
    real f = 0.031, k = 0.0585;
    // real f = 0.04, k = 0.065;
    //  real f = 0.049, k = 0.0629;
    //  real f = 0.025, k = 0.0535;

    real da = 2.00e-4, db = 0.5 * da;
    reaction_diffusion::ptr rx = reaction_diffusion::create(__M, f, k, da, db);
    _rx = std::dynamic_pointer_cast<module_base>(rx);

    init_origin();
  };

  void init_origin() {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);

    std::array<vec3, 2> ext = asawa::extents(x);
    _origin = vec3(                    //
        0.5 * (ext[0][0] + ext[1][0]), //
        ext[0][1],                     //
        0.5 * (ext[0][2] + ext[1][2]));
  }

  std::vector<vec4> get_mesh_colors() {
    std::vector<vec4> colors(__M->vert_count(), vec4(1.0, 0.0, 0.0, 1.0));
    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxb();

    vec4 col_a(0.75, 0.35, 0.0, 1.0);
    vec4 col_b(0.0, 2.0, 1.5, 1.0);

    for (int k = 0; k < __M->vert_count(); k++) {
      colors[k] = rxa[k] * col_a + 1.5 * rxb[k] * col_b;
      // if (k % 2 == 0)
      //   colors[k] = vec4(0.0, 0.0, 1.0, 1.0);
    }
    return colors;
  }

  std::vector<real> geodesic_weight(asawa::shell::shell &shell) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(shell, 0);

    index_t imin = 0.0;
    real min = std::numeric_limits<real>::max();
    for (int i = 0; i < x.size(); i++) {
      real d = (x[i] - _origin).norm();
      if (d < min) {
        min = d;
        imin = i;
      }
    }
    vec3 x_o = x[imin];
    _origin = x_o;

    std::vector<real> f(shell.vert_count(), 0.0);
    f[imin] = 1.0;

    bontecou::laplacian L(__M, x);
    std::vector<real> d = L.heatDist(f, 0.2);
#if 0
    for (int i = 0; i < shell.vert_count(); i++) {
      vec3 xi = x[i];
      vec3 N = vert_normal(*__M, i, x);
      logger::line(xi, xi + 0.1 * d[i] * N, vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
    return d;
  }

  std::vector<real> growth_weights(asawa::shell::shell &shell) {
    std::vector<real> g_vert_geodesic = geodesic_weight(shell);
    auto range = shell.get_edge_range();
    std::vector<real> g_edge(__M->edge_count(), 0.0);

    std::vector<real> &rxa =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxa();
    std::vector<real> &rxb =
        std::dynamic_pointer_cast<reaction_diffusion>(_rx)->get_rxb();

    for (auto c0 : range) {
      int i = shell.vert(c0);
      int j = shell.vert(shell.other(c0));
      real di = g_vert_geodesic[i];
      real dj = g_vert_geodesic[j];
      real d = 0.5 * (di + dj);
      real ra = rxa[i] + rxa[j];
      real rb = rxb[i] + rxb[j];
      real g = d * (2.0 * rb - 0.2 * ra);
      // real g = d * (6.0 * rb - 0.5 * ra);

      g_edge[c0 / 2] = 1.0 + 0.01 * g;
    }
    return g_edge;
  }

  void step_dynamics(int frame) {
    for (int i = 0; i < 10; i++)
      _rx->step(20.0 * _h);

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas_3(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());
    std::vector<real> g = growth_weights(*__M);

    hepworth::vec3_block::ptr X = hepworth::vec3_block::create(M, x, v, f);

    hepworth::block::init_edge_growth(*__M, constraints, x, li, g, 1.0e-4, {X});
    hepworth::block::init_edge_strain(*__M, constraints, x, li, 1e-1, {X});
    hepworth::block::init_bending(*__M, constraints, x, 2.0e-1, {X});
    hepworth::block::init_edge_willmore(*__M, constraints, 2.5e-1, {X});
    hepworth::block::init_area(*__M, constraints, x, 3.0e-1, {X}, false);

    // hepworth::block::init_pinned(*__M, constraints, x, 5.0e-3, {X});
    //  hepworth::block::init_laplacian(*__M, constraints, x, 1, 2.0e-1, {X});

    //  hepworth::shell::init_cross(*__M, constraints, 1.05, 0.1);
    real eps = 3.0 * __surf->_Cm;
    // hepworth::shell::init_edge_edge_collisions(*__M, *__surf, constraints,
    // x,
    //                                            eps, 0.5 * eps, 1.0e0);

    hepworth::block::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             4.0 * eps, 1.0 * eps, 1.0, {X, X});

    solver.set_constraints(constraints);

    std::vector<hepworth::sim_block::ptr> blocks = {X};
    solver.step(blocks, _h, 0.5, 30);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_dynamics(frame);
    __surf->step(false);
  }

  module_base::ptr _rx;

  vec3 _origin;
  real _h = 0.1;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif