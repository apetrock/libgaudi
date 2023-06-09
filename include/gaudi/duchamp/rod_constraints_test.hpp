#include "Eigen/src/Geometry/Scaling.h"
#include "gaudi/common.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"

#include "gaudi/asawa/primitive_objects.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"

#include "gaudi/calder/tree_code.hpp"

#include <array>

#include <math.h>
#include <random>

#include <cmath>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

void center(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];
  real s = 2.0 / maxl;
  vec3 cen = (0.5 * (max + min) - min) * s;
  std::cout << " scale: " << s << std::endl;
  cen -= min;
  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= cen;
  }
}

class block_test {
public:
  typedef std::shared_ptr<block_test> ptr;

  static ptr create() { return std::make_shared<block_test>(); }

  block_test() {
    //__M = load_cube();

    __M = shell::load_bunny();
    load_loop_rod();

    triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    std::vector<vec3> &x = x_datum->data();
    center(x);
  };

  void load_loop_rod() {
    std::uniform_real_distribution<real> dist(0.5, 1.0);
    std::mt19937_64 re;
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    real r0 = p0.norm();
    real r1 = p1.norm();

    p0.normalize();
    p1.normalize();
    vec3 f2 = p0.cross(p1).normalized();
    vec3 f1 = p0.cross(f2).normalized();
    vec3 f0 = f1.cross(f2).normalized();

    std::cout << p0 << ::std::endl;
    std::cout << ::std::endl;
    std::cout << p1 << ::std::endl;

    std::vector<vec3> points;
    int N = 256;
    for (int i = 0; i < N; i++) {
      real thet = 2.0 * M_PI * real(i) / real(N);
      std::cout << cos(thet) << " " << sin(thet) << std::endl;
      vec3 pi = r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }
    __R = rod::rod::create(points);

    real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = __R->__l0;
    std::vector<real> lp(l0);
    real h = 0.1;
    real bnd = 2.0;
    real s_vol = 4.0 / 3.0 * M_PI * pow(bnd, 3.0);
    real r_vol = __R->get_total_volume();
    real k = 0.1;
    real t1 = 1.0;
    real t0 = 1.01;
    real att = t0 + (t1 - t0) / (1.0 + exp(-k * (r_vol)));
    
    std::vector<vec3> f(__R->__v.size(), vec3::Zero());
    for (int i = 0; i < __R->__x.size(); i++) {
      vec3 x = __R->__x[i];
      real xn = x.norm();
      if (xn > bnd) {
        f[i] = -0.01 * (xn - bnd) / h / h * x;
      }
      // f[i] += 1e-1 * vec3::Random();
    }
    
    hepworth::vec3_block::ptr x = hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr u = hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    std::cout << "vol:" << r_vol << "/" << s_vol << " attenuation: " << att
              << std::endl;

    for (auto &l : l0)
      l *= att;
    // hepworth::rod::init_smooth(*__R, constraints, 0.2);
    // hepworth::rod::init_cylinder(*__R, constraints, 0.1);
    hepworth::block::init_stretch_shear(*__R, constraints, l0, 0.25, {x, u});
    hepworth::block::init_bend_twist(*__R, constraints, 0.1, {u});
    //  hepworth::rod::init_smooth_bend(*__R, constraints, 0.01);
#if 1
     hepworth::block::init_angle(*__R, constraints, vec3(1.0, 0.0, 0.0),
                             0.25 * M_PI, 0.05, {u});
    //hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.1, 0.0),
    //                          0.33 * M_PI, 0.1, {u});
    hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.5, 1.0),
                              0.33 * M_PI, 0.1, {u});
#endif
    hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {x, x});
    solver.set_constraints(constraints);

    // f[0][0] = 1.0;
    std::vector<hepworth::sim_block::ptr> blocks = {x, u};
    solver.step(blocks, h);
  }

  void step(int frame) {

    _frame = frame;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    step_dynamics(frame);
    __Rd->step();
    //__R->debug();
  }

  int _frame;
  shell::shell::ptr __M;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif