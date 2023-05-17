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

class shell_walk_test {
public:
  typedef std::shared_ptr<shell_walk_test> ptr;

  static ptr create() { return std::make_shared<shell_walk_test>(); }

  shell_walk_test() {
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

  vec3 slerp_vec(const vec3 &veca, const vec3 &vecb, float t) {
    vec3 axis = veca.cross(vecb);
    double angle = atan2(axis.norm(), veca.dot(vecb));
    if (axis.norm() < 1e-8) {
      axis = vec3::UnitZ(); // if the vectors are parallel, set the
                            // axis to the z-axis
      angle = 0.0;
    }
    real interpAngle = (1.0 - t) * 0.0 + t * angle;
    vec3 interp = cos(interpAngle) * veca +
                  sin(interpAngle) * (axis.normalized().cross(veca));
    return interp;
  };

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    __surf->step(false);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    index_t test = 1000;
    while (__M->next(test) < 0)
      test++;
    vec3 N = asawa::shell::edge_normal(*__M, test, x);
    vec3 T = asawa::shell::edge_tangent(*__M, test, x).normalized();
    vec3 B = N.cross(T).normalized();
    real thet = 0.0000001 * 2.0 * M_PI * (frame + 10);
    vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
    std::vector<index_t> corners;
    std::vector<real> S;
    vec3 x0 = x[__M->vert(test)];
    vec3 N0 = asawa::shell::edge_normal(*__M, test, x);
    asawa::shell::walk(*__M, x, test, dir, 0.5, 5000, 1e-2,
                       [&](const asawa::shell::shell &M,
                           const std::vector<vec3> &x, const index_t &corner,
                           const real &s, const real &accumulated_length,
                           vec3 &dir) {
                         asawa::shell::index_t ci = corner;
                         real si = s;
                         vec3 xi = asawa::shell::edge_vert(M, ci, si, x);
                         real d = (xi - x0).norm();
                         vec3 Ni = asawa::shell::edge_normal(M, ci, x);
                         if (false) {
                           vec3 B = N0.cross(dir).normalized();
                           vec3 rdir = va::reject<real>(B, dir).normalized();
                           // dir = slerp_vec(dir, rdir, 0.75);
                           dir = rdir;
                           /*
                           gg::geometry_logger::line(xi + 0.005 * Ni,
                                                     xi + 0.005 * Ni + 0.01 * B,
                                                     vec4(0.5, 1.0, 0.5, 1.0));
                           gg::geometry_logger::line(xi + 0.005 * Ni,
                                                     xi + 0.005 * Ni + 0.01 *
                           dir, vec4(1.0, 0.5, 0.5, 1.0));
                           gg::geometry_logger::line(xi + 0.005 * Ni,
                                                     xi + 0.005 * Ni + 0.01 *
                           N0, vec4(0.5, 0.5, 1.0, 1.0));
                          */
                         }
                         // real t = 1.0 / (1.0 + d * d);
                         // std::cout << "t: " << t << " " << d << std::endl;
                         // dir = slerp_vec(dir, gdir, t);
                         S.push_back(si);
                         corners.push_back(ci);

                         return true;
                       });
    // shell::crack_edges(*__M, corners, S, 1e-2);
    // smoothMesh(0.01, 10);
  }
  real _h = 0.1;
  index_t _il0;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif