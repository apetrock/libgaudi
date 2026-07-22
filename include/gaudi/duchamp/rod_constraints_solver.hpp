#ifndef __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__
#define __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/duchamp/rod_constraints_dynamics.hpp"
#include "gaudi/duchamp/rod_force_nodes.hpp"
#include "gaudi/duchamp/utils/sdf.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

using namespace asawa;

class rod_constraints_solver {
public:
  using ptr = std::shared_ptr<rod_constraints_solver>;

  static ptr create() { return std::make_shared<rod_constraints_solver>(); }

  rod_constraints_solver() {
    load_loop_rod();
    const int N = 13;
    const real r1 = 1.5;
    const real r11 = 0.5;
    const real pi43 = 4.0 / 3.0 * M_PI;
    const real v0 = real(N) * pi43 * pow(r11, 3.0);
    const real r0 = std::pow(v0 / pi43, 1.0 / 3.0);
    __sdf0 = sdf_sphere::create(vec3(0.0, 0.0, 0.0), r0);
    __sdf1 = sdf_multi_sphere::create(get_fib(r1, N), r11);

    _rod_pos = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);
    _rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);
    _config = hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                                    hepworth::block::rod_quaternion_block>::create()
                  .with_blocks(_rod_pos, _rod_quat)
                  .with_bundle(hepworth::block::make_rod_physics_bundle<0, 1>(
                      __R, __Rd, 1e-1, 2e-1, 1.0))
                  .dt(0.05)
                  .damping(0.01)
                  .build();

    _boundary = _graph.create_node<boundary_gradient_node>(__R, __sdf0, __sdf1);
    _tangent = _graph.create_node<tangent_point_gradient_node>(__R, __Rd, 0.0e-7);
    _vortex = _graph.create_node<vortex_force_node>(__R, __Rd, 1e-1, 4.0, 0.0);
    _add = _graph.create_node<vec3_junction_node<3>>();
    _solver_node = _graph.create_node<
        hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>>(_config);

    _graph.link(_boundary->output(), _add->template input<0>());
    _graph.link(_tangent->output(), _add->template input<1>());
    _graph.link(_vortex->output(), _add->template input<2>());
    _graph.link(_add->output(), _solver_node->input_at<0>());
  }

  const rod::rod &rod() const { return *__R; }
  rod::rod::ptr rod_ptr() const { return __R; }

  void step(int frame) {
    _frame = frame;
    grow_rod_rest_lengths(*__R, 1.0 - 0.0001);
    _boundary->set_frame(frame);
    _graph.run();
    __Rd->step();
  }

private:
  std::vector<vec3> get_fib(real r0, int N = 13) {
    const real golden = 0.5 * (1.0 + sqrt(5));
    std::vector<vec3> cens(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      const real theta = 2.0 * M_PI * i / golden;
      const real phi = acos(1.0 - 2.0 * (i + 0.5) / real(N));
      cens[i] = r0 * vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
    }
    return cens;
  }

  void load_loop_rod() {
    std::uniform_real_distribution<real> dist(0.5, 1.0);
    std::mt19937_64 re;
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    const real r0 = p0.norm();
    const real r1 = p1.norm();

    p0.normalize();
    p1.normalize();
    const vec3 f2 = p0.cross(p1).normalized();
    const vec3 f1 = p0.cross(f2).normalized();
    const vec3 f0 = f1.cross(f2).normalized();

    const int N = 256;
    auto make_ellipse_loop = [&](const vec3 &center, real er0, real er1, const vec3 &axis0,
                                 const vec3 &axis1) {
      std::vector<vec3> pts;
      pts.reserve(N);
      for (int i = 0; i < N; i++) {
        const real thet = 2.0 * M_PI * real(i) / real(N);
        pts.push_back(center + er0 * cos(thet) * axis0 + er1 * sin(thet) * axis1);
      }
      return pts;
    };

    const vec3 cen = vec3::Zero();
    const vec3 norm_axis = f2;
    const real C = 0.2;
    const vec3 c0 = cen - C * norm_axis;
    const vec3 c1 = cen + C * norm_axis;

    __R = rod::rod::create();
    __R->append_loop(make_ellipse_loop(c0, r0, r1, f0, f1));
    __R->append_loop(make_ellipse_loop(c1, r0, r1, f0, f1));
    const real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  int _frame = 0;
  sdf_base::ptr __sdf0;
  sdf_base::ptr __sdf1;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;

  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_quaternion_block::ptr _rod_quat;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block,
                                       hepworth::block::rod_quaternion_block> _config;
  liblombardi::GraphContext _graph;
  boundary_gradient_node::ptr _boundary;
  tangent_point_gradient_node::ptr _tangent;
  vortex_force_node::ptr _vortex;
  vec3_junction_node<3>::ptr _add;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                     hepworth::block::rod_quaternion_block>::ptr _solver_node;
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__
