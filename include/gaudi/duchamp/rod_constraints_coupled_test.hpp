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

#include "gaudi/hepworth/block/coupling_collisions_init.hpp"
#include "gaudi/hepworth/block/weld_constraints.hpp"

#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"

#include "gaudi/calder/tree_code.hpp"

#include <array>

#include <cassert>
#include <chrono>
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

class block_test {
public:
  typedef std::shared_ptr<block_test> ptr;

  static ptr create() { return std::make_shared<block_test>(); }

  block_test() {
    //__M = load_cube();

    load_loop_rod();
    load_loop_rod();
  };

  void load_loop_rod() {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    std::uniform_real_distribution<real> dist(-1.0, 1.0);
    std::mt19937_64 re(seed);
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    real r0 = p0.norm();
    real r1 = p1.norm();

    p0.normalize();
    p1.normalize();
    vec3 f2 = p0.cross(p1).normalized();
    vec3 f1 = p0.cross(f2).normalized();
    vec3 f0 = f1.cross(f2).normalized();

    std::cout << "p0 : " << p0 << ::std::endl;
    std::cout << ::std::endl;
    std::cout << "p1 : " << p1 << ::std::endl;

    std::vector<vec3> points;
    int N = 256;
    for (int i = 0; i < N; i++) {
      real thet = 2.0 * M_PI * real(i) / real(N);
      vec3 pi = r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }

    rod::rod::ptr R = rod::rod::create(points);
    real lavg = R->lavg();
    rod::dynamic::ptr Rd =
        rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
    __R.push_back(R);
    __Rd.push_back(Rd);
  }

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    real h = 0.1;
    std::vector<real> bnd = {2.0, 1.5};

    real tk = 0.1;
    real t1 = 1.0;
    real t0 = 1.01;

    std::vector<vec3> cens = {
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 0.0, 0.0),
    };

    std::vector<hepworth::sim_block::ptr> blocks;

    for (int k = 0; k < __R.size(); k++) {
      asawa::rod::rod::ptr R = __R[k];
      std::vector<real> &l0 = R->__l0;
      std::vector<real> lp(l0);

      real s_vol = 4.0 / 3.0 * M_PI * pow(bnd[k], 3.0);
      real r_vol = R->get_total_volume();
      real att = t0 + (t1 - t0) / (1.0 + exp(-tk * (r_vol)));
      std::vector<vec3> f(R->__v.size(), vec3::Zero());

      for (int i = 0; i < R->__x.size(); i++) {
        vec3 dx = R->__x[i] - cens[k];
        real xn = dx.norm();
        if (xn > bnd[k]) {
          f[i] = -0.02 * (xn - bnd[k]) / h / h * dx;
        }
        // f[i] += 1e-1 * vec3::Random();
      }

      hepworth::vec3_block::ptr x =
          hepworth::vec3_block::create(R->__M, R->__x, R->__v, f);
      hepworth::quat_block::ptr u =
          hepworth::quat_block::create(R->__J, R->__u, R->__o);

      blocks.push_back(x);
      blocks.push_back(u);

      std::cout << "vol:" << r_vol << "/" << s_vol << " attenuation: " << att
                << std::endl;

      for (auto &l : l0)
        l *= att;

      real str = 0.03;
      real twi = 0.01;
      if (k == 0) {
        str = 0.25;
        twi = 0.1;
      }
      // hepworth::rod::init_smooth(*__R, constraints, 0.2);
      hepworth::block::init_stretch_shear(*R, constraints, l0, str, {x, u});
      hepworth::block::init_bend_twist(*R, constraints, 0.1, {u});
#if 1
      if (k == 0) {
        hepworth::block::init_angle(*R, constraints, vec3(1.0, 0.0, 0.0),
                                    0.35 * M_PI, 0.05, {u});
        // hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.1, 0.0),
        //                           0.33 * M_PI, 0.1, {u});
        hepworth::block::init_angle(*R, constraints, vec3(0.0, 0.0, 1.0),
                                    0.13 * M_PI, 0.1, {u});
      } else {
        // hepworth::block::init_cylinder(*R, constraints, 0.01, {x});
      }
#endif
      hepworth::block::init_collisions(*R, *__Rd[k], constraints, 1.0, {x, x});
    }
    constraints.push_back(hepworth::block::pnt_pnt_weld::create(
        {0, 0}, 2.0, {blocks[0], blocks[2]}));
    hepworth::block::init_coupling_collisions(*__R[0], *__Rd[0], //
                                              *__R[1], *__Rd[1], constraints,
                                              1.0, {blocks[0], blocks[2]});

    solver.set_constraints(constraints);

    solver.step(blocks, h);
  }

  void step(int frame) {

    _frame = frame;

    step_dynamics(frame);
    for (auto Rd : __Rd)
      Rd->step();
    //__R->debug();
  }

  int _frame;
  std::vector<rod::rod::ptr> __R;
  std::vector<rod::dynamic::ptr> __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif