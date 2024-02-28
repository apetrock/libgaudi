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
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"

#include "gaudi/calder/tree_code.hpp"

#include "utils/sdf.hpp"
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

    load_loop_rod(1.0);
    load_loop_rod(1.0);
    load_loop_rod(1.0);
    _eps = __R[0]->lavg();
  };

  void load_loop_rod(real r = 1.0) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    std::uniform_real_distribution<real> dist(-1.0, 1.0);
    std::mt19937_64 re(seed);

    vec3 cen = r * vec3(dist(re), dist(re), dist(re));
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    real r0 = r;
    real r1 = r;

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
      vec3 pi = cen + r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }

    rod::rod::ptr R = rod::rod::create(points);
    real lavg = R->lavg();
    R->set_radius(3.5 * lavg);
    rod::dynamic::ptr Rd =
        rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
    __R.push_back(R);
    __Rd.push_back(Rd);
  }

  real sigmoid(real x, real k, real x0) {
    return 1.0 / (1.0 + exp(-(x - x0) / k));
  }

  real one_minus_sigmoid(real x, real k, real x0) {
    return 1.0 - sigmoid(x, k, x0);
  }

  real frame_ratio(int frame, int t0, int t1) {
    return real(frame - t0) / real(t1 - t0);
  }

  void load_sdf() {

    // shell::shell::ptr __M = shell::load_obj("assets/hand_trimmed.obj");
    shell::shell::ptr __M = shell::load_obj("assets/skull.obj");

    // shell::shell::ptr __M = shell::load_obj("assets/wolf_skull.obj");

    shell::triangulate(*__M);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);
    real l0 = asawa::shell::avg_length(*__M, x);

    std::vector<index_t> face_vert_ids = __M->get_face_vert_ids();
    geometry_sdf::ptr sdf = geometry_sdf::create();
    sdf->add_geometry(x, face_vert_ids, 2.5, vec3(0.0, 0.0, 0.0));
    __sdf = sdf;
  }

  std::vector<vec3> compute_boundary_gradients(const rod::rod &R) {
    if (!this->__sdf)
      load_sdf();

    std::vector<real> dists = __sdf->distance(R.__x);
    std::vector<vec3> gdists = __sdf->grad_distance(R.__x);
    std::vector<vec3> f(R.__x.size(), vec3::Zero());

    for (int i = 0; i < R.__x.size(); i++) {

      // lets see if we can hollow out the sdf...
      real d = dists[i];
      real k = -0.20;
      d = d < k ? -d + 2.0 * k : d;
      // vec3 fi = va::sgn(dists[i] * d) * d * gdists[i];
      vec3 fi = d * gdists[i];

#if 0
      vec3 x = R.__x[i];
      if (d > 0.0) {
        logger::line(x, x - fi, vec4(0.0, 1.0, 0.0, 1.0));
      } else {
        logger::line(x, x - fi, vec4(1.0, 0.0, 0.0, 1.0));
      }
#endif
      if (d > 0.0) {
        f[i] = fi;
      }
    }
    return f;
  }

#if 1
  std::vector<vec3> calc_tangent_point_gradient(asawa::rod::rod &R, real eps) {
    std::vector<vec3> &x = R.x();
    std::vector<real> l = R.l0();
    std::vector<vec3> T = R.N2c();
    std::vector<vec3> xc = R.xc();

    std::vector<vec3> g0 = calder::tangent_point_gradient(R, x, l, T, eps, 6.0);
    return std::move(g0);
  }
#endif

#if 1
  std::vector<vec3> calc_vortex_forces(int i0, std::vector<rod::rod::ptr> &R,
                                       real eps) {
    int Nr = R.size();
    int i1 = (i0 + 1) % Nr;
    int i2 = (i0 + 2) % Nr;
    rod::rod &R0 = *R[i0];
    rod::rod &R1 = *R[i1];
    rod::rod &R2 = *R[i2];
    std::vector<vec3> &x0 = R0.x();
    std::vector<vec3> &x1 = R1.x();
    std::vector<vec3> &x2 = R2.x();
    std::vector<real> l1 = R1.l0();
    std::vector<real> l2 = R2.l0();
    std::vector<vec3> g1 = calder::vortex_force(R1, x0, l1, eps, 3.0);
    std::vector<vec3> g2 = calder::vortex_force(R2, x0, l2, eps, 3.0);
    std::vector<vec3> g0 = g1 + g2; // does this work???
#if 0
    for (int i = 0; i < x0.size(); i++) {
      vec3 x = x0[i];
      logger::line(x, x + g0[i], vec4(0.0, 1.0, 0.0, 1.0));
    }
#endif
    return std::move(g0);
  }
#endif

  void step_dynamics(int frame) {
    std::cout << "frame: " << frame << std::endl;
    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    real h = 0.1;
    std::vector<real> bnd = {1.5, 1.25};

    std::vector<real> tk = {0.04, 0.125};
    real t1 = 1.0;
    real t0 = 1.01;

    std::vector<vec3> cens = {
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 0.0, 0.0),
    };

    std::vector<hepworth::sim_block::ptr> blocks;

    for (int k = 0; k < __R.size(); k++) {
      std::cout << "rod " << k << " size: " << __R[k]->__x.size() << std::endl;
      asawa::rod::rod::ptr R = __R[k];
      std::vector<real> &l0 = R->__l0;

      std::vector<vec3> f(R->__v.size(), vec3::Zero());

      int t_start = 1200;
      int t_end = 1500;
      if (frame > t_start && frame < t_end) {
        std::vector<vec3> ft =
            calc_tangent_point_gradient(*R, 2.0 * __Rd[k]->_Cc);
        std::vector<vec3> fv = calc_vortex_forces(k, __R, 6.0 * __Rd[k]->_Cc);
        real t = frame_ratio(frame, t_start, t_end);
        f = f - 2.0e-6 * ft;
        f = f - 5e-1 * fv;
      } else {
        std::vector<vec3> fb = compute_boundary_gradients(*R);
        std::vector<vec3> fv = calc_vortex_forces(k, __R, 6.0 * __Rd[k]->_Cc);
        f = 8.0 * fb - 5.0e-2 * one_minus_sigmoid(real(frame), 200.0, 600) * fv;
      }

      hepworth::vec3_block::ptr x =
          hepworth::vec3_block::create(R->__M, R->__x, R->__v, f);
      hepworth::quat_block::ptr u =
          hepworth::quat_block::create(R->__J, R->__u, R->__o);

      blocks.push_back(x);
      blocks.push_back(u);

      std::vector<real> dists = __sdf->distance(R->__x);
      std::vector<index_t> verts = R->get_vert_range();
      for (auto i : verts) {
        asawa::rod::consec_t c = R->consec(i);
        // if (dists[i] > 0.0)
        //   continue;
        // l0[i] *= 8.0 * pow(dists[i], 0.75) *
        //         one_minus_sigmoid(real(frame), 200.0, 1400);
        l0[i] *= 1.0 + 0.042 * one_minus_sigmoid(real(frame), 200.0, 1000);
        if (dists[i] > 0.0)
          l0[i] *= 1.0 - 0.02 * dists[i] * sigmoid(real(frame), 200.0, 1800);

        // l0[i] = std::max(l0[i], 1e-8);
      }

      // hepworth::rod::init_smooth(*__R, constraints, 0.2);
      hepworth::block::init_stretch_shear(*R, constraints, l0, 4.0e-2, {x, u});
      hepworth::block::init_straight(*R, constraints, 1.0e-2, {u}, true);

      hepworth::block::init_collisions(*R, *__Rd[k], constraints, 2.0, {x, x});
    }

    hepworth::block::init_coupling_collisions(*__R[0], *__Rd[0], //
                                              *__R[1], *__Rd[1], constraints,
                                              1.0, {blocks[0], blocks[2]});
    hepworth::block::init_coupling_collisions(*__R[0], *__Rd[0], //
                                              *__R[2], *__Rd[2], constraints,
                                              1.0, {blocks[0], blocks[4]});
    hepworth::block::init_coupling_collisions(*__R[1], *__Rd[1], //
                                              *__R[2], *__Rd[2], constraints,
                                              1.0, {blocks[2], blocks[4]});
    /*
    constraints.push_back(hepworth::block::pnt_pnt_weld::create(
        std::vector<index_t>({0, 0}), 1e-2,
        std::vector({blocks[0], blocks[2]})));
    constraints.push_back(hepworth::block::pnt_pnt_weld::create(
        std::vector<index_t>({0, 0}), 1e-2,
        std::vector({blocks[0], blocks[4]})));
    constraints.push_back(hepworth::block::pnt_pnt_weld::create(
        std::vector<index_t>({0, 0}), 1e-2,
        std::vector({blocks[2], blocks[4]})));
    */
    solver.set_constraints(constraints);

    solver.step(blocks, h, 0.5);
  }

  void step(int frame) {

    _frame = frame;

    if (frame < 2000) {
      step_dynamics(frame);
      for (auto Rd : __Rd)
        Rd->step();
      //__sdf->debug_target();
    }
    if (frame > 2400) {
      exit(0);
    }
    //__R->debug();
  }

  int _frame;
  real _eps = 1e-3;
  sdf_base::ptr __sdf;
  std::vector<rod::rod::ptr> __R;
  std::vector<rod::dynamic::ptr> __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif