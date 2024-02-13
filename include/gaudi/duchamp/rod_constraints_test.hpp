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
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/asawa/primitive_objects.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"

#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "utils/sdf.hpp"

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

class block_test {
public:
  typedef std::shared_ptr<block_test> ptr;

  static ptr create() { return std::make_shared<block_test>(); }

  block_test() {
    //load_fib_rod();
    load_loop_rod();
    int N = 13;
    real r1 = 1.5;
    real r11 = 0.5;
    real pi43 = 4.0 / 3.0 * M_PI;
    real v0 = real(13.0) * pi43 * pow(r11, 3.0);
    real r0 = std::pow(v0 / pi43, 1.0 / 3.0);
    std::cout << "radius 0"<<r0 << std::endl;
    __sdf0 = sdf_sphere::create(vec3(0.0, 0.0, 0.0), r0);
    __sdf1 = sdf_multi_sphere::create(get_fib(r1, 13), r11);
    // load_sdf();
  };

  std::vector<vec3> get_fib(real r0, int N = 13){
    
    real golden = 0.5 * (1.0 + sqrt(5));
    std::vector<vec3> cens(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      real theta = 2.0 * M_PI * i / golden;
      real phi = acos(1.0 - 2.0 * (i + 0.5) / real(N));
      vec3 ci =
          r0 * vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
      cens[i] = ci;
    }
    return cens;
  }

  void load_fib_rod() {
    __R = rod::rod::create(get_fib(1.5, 13));

    real lavg = 0.01 * __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

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
    //compute _l0 sum of distances
    _lt0 = 0.0;
    for (int i = 0; i < N; i++) {
      _lt0 += (points[(i + 1) % N] - points[i]).norm();
    }
    __R = rod::rod::create(points);

    real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

#if 1
  std::vector<vec3> compute_tangent_point_gradient() {
    real eps = __Rd->_Cc;
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_gradient(*__R, x, l, T, 1.0 * eps, 6.0);
    return g0;
  }
#endif
  
  sdf_base::ptr get_sdf(int frame) {
    if ((frame / 400) % 2 == 0) {
      return __sdf0;
    } else {
      return __sdf1;
    }
  }
  
    
  std::vector<real> compute_growth_weights(index_t frame){

    auto sdf = get_sdf(frame);
    
    std::vector<real> dists = sdf->distance(__R->__x);
    std::vector<real> w(__R->__x.size(), 0);
    std::vector<vec3> xc = __R->xc();
    std::vector<vec3> N = __R->N0c();
    std::vector<index_t> verts = __R->get_vert_range();
      for (auto i : verts) {
        asawa::rod::consec_t c = __R->consec(i);

        vec3 xi = xc[i];
        vec3 Ni = N[i];
        real di = dists[i];

        #if 0
          vec4 c0 = vec4(0.0, 1.0, 0.0, 1.0);
        vec4 c1 = vec4(1.0, 0.0, 0.0, 1.0);
        if(di > 0.0){
            gg::geometry_logger::line(xi, xi + 0.1*di * Ni, c0);
        } else{
            gg::geometry_logger::line(xi, xi + 0.1*di * Ni, c1);
        }
        #endif
        di = di < 0.0 ? -1.0 : 0.5*di;
        w[i] = (1.0 - 0.28 * di);
        w[i] = va::clamp(w[i], 0.0, 4.00);
        // l0[i] = std::max(l0[i], 1e-8);
      }
  
    return std::move(w);
  }

  std::vector<vec3> compute_boundary_gradients(index_t frame) {

    auto sdf = get_sdf(frame);
    
    std::vector<real> dists = sdf->distance(__R->__x);
    std::vector<vec3> gdists = sdf->grad_distance(__R->__x);
    std::vector<vec3> f(__R->__x.size(), vec3::Zero());
    std::vector<vec3> xc = __R->xc();
    for (int i = 0; i < __R->__x.size(); i++) {
      if (dists[i] > 0.0) {
        f[i] = -dists[i] * gdists[i];
      }
        #if 0
        vec4 c0 = vec4(0.0, 1.0, 0.0, 1.0);
        vec4 c1 = vec4(1.0, 0.0, 0.0, 1.0);
        if(dists[i] > 0.0){
            gg::geometry_logger::line(xc[i], xc[i] + f[i], c0);
        } else{
            gg::geometry_logger::line(xc[i], xc[i] + f[i], c1);   
        }
        #endif
    }

    
    return std::move(f);
  }

  void step_dynamics(int frame) {
    std::cout << "frame: " << frame << ", size: " << __R->__x.size()
              << std::endl;
    hepworth::block::projection_solver solver;

    std::vector<hepworth::projection_constraint::ptr> constraints;

    std::vector<real> &l0 = __R->__l0;
    std::vector<real> lp(l0);
    real h = 0.05;
    std::vector<vec3> f(__R->__v.size(), vec3::Zero());
    // std::vector<vec3> fr = compute_coulomb_gradient();
    // std::vector<vec3> fr = compute_null_coulomb_gradient();
    std::vector<vec3> fr = compute_tangent_point_gradient();
    std::vector<vec3> fb = compute_boundary_gradients(frame);
    real fwave = 0.5 + 0.5*cos(M_PI * real(frame)/250.0);
    //f = 16.0 * fb - 1e-6*fwave* fr;
    f = 8.0 * fb;

    hepworth::vec3_block::ptr x =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr u =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);



    // hepworth::rod::init_smooth(*__R, constraints, 0.2);


    //  hepworth::rod::init_smooth_bend(*__R, constraints, 0.01);
#if 1
    if(real(frame % 500) / 500.0 > 0.5 && 0){
      for(int i = 0; i < l0.size(); i++){
        l0[i] *= 1.035;
      }
      hepworth::block::init_stretch_shear(*__R, constraints, l0, 1e-1, {x, u});
      hepworth::block::init_bend_twist(*__R, constraints, 1e-1, {u}, true);
      hepworth::block::init_angle(*__R, constraints, vec3(1.0, 0.0, 0.0),
                                  0.28 * M_PI, 0.1, {u});
      //  hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.1, 0.0),
      hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.0, 1.0),
                                  0.25 * M_PI, 0.1, {u});
    }
    else{
      //real th = (1000.0 - real(frame))/1000.0;
      std::vector<real> w = compute_growth_weights(frame);
      real lt = 0.0;
      for(int i = 0; i < w.size(); i++){
        lt += l0[i];
      }
      real wl = pow(_lt0 / lt, 0.25);
      
      for(int i = 0; i < w.size(); i++){
        l0[i] *= w[i];
      }

      std::cout << " lt0 / lt: " << _lt0 / lt << " w: " << wl << std::endl;
      hepworth::block::init_helicity(*__R, constraints, 1e-0*wl, {x});
      
      hepworth::block::init_stretch_shear(*__R, constraints, l0, 2e-2, {x, u});
      hepworth::block::init_bend_twist(*__R, constraints, 2e-2, {u}, true);
    }
#endif
    hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {x, x});
    solver.set_constraints(constraints);

    // f[0][0] = 1.0;
    std::vector<hepworth::sim_block::ptr> blocks = {x, u};
    solver.step(blocks, h);
  }

  void step(int frame) {

    _frame = frame;


    step_dynamics(frame);
    __Rd->step();
    if(frame > 3000)
    exit(0);
    //__R->debug();
  }

  int _frame;
  real _lt0 = 1.0;
  sdf_base::ptr __sdf0;
  sdf_base::ptr __sdf1;
  
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif