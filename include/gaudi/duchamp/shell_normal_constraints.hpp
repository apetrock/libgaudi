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

class shell_normal_constraints {
public:
  typedef std::shared_ptr<shell_normal_constraints> ptr;

  static ptr create() { return std::make_shared<shell_normal_constraints>(); }

  shell_normal_constraints() {
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
    walk_and_subdivide();
    __surf->step(true);
  };

  void walk_and_subdivide() {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    index_t test = 1000;

    vec3 N = asawa::shell::edge_normal(*__M, test, x);
    vec3 T = asawa::shell::edge_tangent(*__M, test, x).normalized();
    vec3 B = N.cross(T).normalized();

    real thet = 1e-2 * 2.0 * M_PI * 0.2;

    vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
    std::vector<index_t> corners;
    std::vector<real> S;
    vec3 x0 = x[__M->vert(test)];
    asawa::shell::walk(*__M, x, test, dir, 0.5, 2000, 1e-2,
                       [&](const asawa::shell::shell &M,
                           const std::vector<vec3> &x, const index_t &corner,
                           const real &s, const real &accumulated_length,
                           vec3 &dir) {
                         asawa::shell::index_t ci = corner;
                         S.push_back(s);
                         corners.push_back(corner);
                         return true;
                       });

    std::vector<index_t> walk = shell::crack_edges(*__M, corners, S, 1e-2);
    _walk = walk;
    std::vector<index_t> walk_stitch = shell::stitch_walk(*__M, x, walk);
    std::vector<int> has_walk = std::vector<int>(__M->corner_count() / 2, 0);

    for (const index_t &c : walk_stitch) {
      has_walk[c / 2] = 1;
    }

    rod_datum::ptr wdata = rod_datum::create(prim_type::EDGE, has_walk);
    _iwalk = __M->insert_datum(wdata);

    __surf->set_flip_pred([this](shell::shell &M, const index_t &c0) {
      std::vector<int> &walk = asawa::get_rod_data(M, _iwalk);
      std::cout << "flip: " << c0 / 2 << " " << walk.size() << std::endl;
      return walk[c0 / 2] == 1;
    });

    __surf->set_collapse_pred([this](shell::shell &M, const index_t &c0) {
      std::vector<int> &walk = asawa::get_rod_data(M, _iwalk);
      index_t c1 = M.other(c0);
      index_t c0p = M.prev(c0);
      index_t c1p = M.prev(c1);
      bool val0 = walk[c0 / 2];
      if (val0)
        return false;
      bool val0p = walk[c0p / 2];
      bool val1p = walk[c1p / 2];
      return val0p == 1 || val1p == 1;
    });
  }

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

  void print_walk() {
    std::vector<int> &walk = asawa::get_rod_data(*__M, _iwalk);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    /*
    for (int i = 0; i < _walk.size() - 1; i++) {
      index_t c0 = _walk[i];
      index_t c1 = _walk[i + 1];
      index_t v0 = __M->vert(c0);
      index_t v1 = __M->vert(c1);

      shell::log_seg_v(*__M, x, v0, v1, 0.005, vec4(0.0, 1.0, 0.0, 1.0));
    }
    */

    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (walk[c0 / 2]) {
        log_edge(*__M, x, c0, 0.005, vec4(0.0, 0.0, 1.0, 1.0));
      }
    }
  }

  void step_dynamics(int frame) {

    hepworth::shell::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> m = asawa::shell::vertex_areas(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);
    solver.set_mass(m);
    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, x);
    }

    std::cout << "l[0]: " << l0[0] << std::endl;
    // hepworth::shell::init_edge_growth(*__M, constraints, x,
    //  l0, 1.01, 1.0e-3);
    hepworth::shell::init_edge_strain(*__M, constraints, x, l0, 1e-1);
    //  hepworth::shell::init_pinned(*__M, constraints, x, 1.0e-1);
    hepworth::shell::init_bending(*__M, constraints, x, 1.0e-3);
    hepworth::shell::init_laplacian(*__M, constraints, x, 1.0e-3, 0);
    hepworth::shell::init_area(*__M, constraints, x, 1.0e-4);

    //  hepworth::shell::init_cross(*__M, constraints, 1.05, 0.1);
    real eps = 0.5 * __surf->_Cm;
    // hepworth::shell::init_edge_edge_collisions(*__M, *__surf, constraints,
    // x,
    //                                            eps, 0.5 * eps, 1.0e0);
    hepworth::shell::init_pnt_tri_collisions(*__M, *__surf, constraints, x,
                                             1.0 * eps, 0.5 * eps, 1.0);
    solver.set_constraints(constraints);

    real bnd = 2.0;
    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());
    for (int i = 0; i < x.size(); i++) {
      vec3 N = asawa::shell::vert_normal(*__M, i, x);

      f[i] += 1e-2 * N;
      //  f[i] += 1e-4 * vec3::Random();
    }
    // f[0][1] = 10.0;
    solver.step(x, v, f, _h);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    step_dynamics(frame);
    __surf->step(true);

    std::cout << "print_walk" << std::endl;
    print_walk();
    // smoothMesh(0.001, 1);
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