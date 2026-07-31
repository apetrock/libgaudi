#ifndef __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__
#define __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/modules/rod_forces.hpp"
#include "gaudi/duchamp/utils/sdf.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/vec_addendum.h"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

using namespace asawa;

// Knobs aligned with legacy rod_constraints_test (step_dynamics).
// Flip k_rod_constraints_active to switch presets.
struct rod_constraints_config {
  rod_physics_config physics{
      .stretch_w = 1e-1,
      .bend_w = 1e-1,
      .twist_w = 1e-1,
      .collision_w = 1.0,
      .smooth_w = 0.0,
      .min_kink_w = 0.0,
      .squad_smooth_w = 0.0,
      .dt = 0.05,
      .damping = 0.1,
      .iterations = 1,
  };
  real grow = 1.04;

  boundary_force_config boundary{.w = 1.0};
  tangent_point_force_config tangent{.w = 1.0e-6, .l0 = 1.0, .p = 6.0};
  vortex_force_config vortex{};
  pca_curve_force_config curve{};

  world_space_moment_config world_space{};
  eigen_moment_config eigen_moment{};

  real frame_scale = 0.12;
};

inline const rod_constraints_config k_rod_constraints_default{
    .physics =
        {
            .stretch_w = 1e-1,
            .bend_w = 1e-1,
            .twist_w = 1e-1,
            .collision_w = 1.0,
            .smooth_w = 0.0,
            .min_kink_w = 0.0,
            .squad_smooth_w = 0.0,
            .dt = 0.05,
            .damping = 0.1,
            .iterations = 1,
        },
    .grow = 1.04,
    .boundary = {.w = 1.0},
    .tangent = {.w = 1.0e-6, .l0 = 1.0, .p = 6.0},
    .vortex = {},
    .curve = {},
    .world_space = {},
    .eigen_moment = {},
    .frame_scale = 0.12,
};

// Legacy baseline + PCA curve force / eigen moment for buckling experiments.
inline const rod_constraints_config k_rod_constraints_buckle{
    .physics =
        {
            .stretch_w = 1e-1,
            .bend_w = 1e-1,
            .twist_w = 2e-1,
            .collision_w = 1.0,
            .smooth_w = 0.0,
            .min_kink_w = 0.15,
            .squad_smooth_w = 8.0e-3,
            .dt = 0.05,
            .damping = 0.1,
            .iterations = 20,
        },
    .grow = 1.02,
    .boundary = {.w = 0.0},
    .tangent = {.w = 0.0e-8, .l0 = 1.0, .p = 6.0},
    .vortex = {},
    .curve = {.w = 1.0e-1, .half_width = 5},
    .world_space = {},
    .eigen_moment = {.w = 5.0e-5, .p = 1.0, .half_width = 5},
    .frame_scale = 0.12,
};

inline const rod_constraints_config k_rod_constraints_active =
    k_rod_constraints_buckle;

// Graph-based rod constraints demo: force nodes → PD solver (+ Vermeer adapter).
class rod_constraints_solver {
public:
  using ptr = std::shared_ptr<rod_constraints_solver>;

  static ptr create(const rod_constraints_config &cfg = k_rod_constraints_active) {
    return std::make_shared<rod_constraints_solver>(cfg);
  }

  explicit rod_constraints_solver(
      const rod_constraints_config &cfg = k_rod_constraints_active)
      : _cfg(cfg) {
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
    _rod_quat =
        std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);
    _rod_quat->with_torque_world([this]() {
      auto tau = compute_pca_curve_torque(*__R, _cfg.eigen_moment);
      return apply_world_space_moment(std::move(tau), _cfg.world_space);
    });

    const auto &phys = _cfg.physics;
    _config =
        hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                              hepworth::block::rod_quaternion_block>::
            create()
                .with_blocks(_rod_pos, _rod_quat)
                .with_bundle(hepworth::block::make_rod_physics_bundle<0, 1>(
                    __R, __Rd, phys.stretch_w, phys.bend_w, phys.twist_w,
                    phys.collision_w, /*free_hinge_i=*/-1,
                    /*free_twist_w=*/0.0, phys.edge_stretch_w))
                .with_recompute(
                    hepworth::block::make_rod_smooth_recompute<0>(
                        __R, phys.smooth_w))
                .with_recompute(
                    hepworth::block::make_rod_min_kink_recompute<0>(
                        __R, phys.min_kink_w))
                .with_recompute(
                    hepworth::block::make_rod_squad_smooth_recompute<1>(
                        __R, phys.squad_smooth_w))
                .dt(phys.dt)
                .damping(phys.damping)
                .iterations(phys.iterations)
                .build();

    _boundary = _graph.create_node<boundary_gradient_node>(
        __R, __sdf0, __sdf1, _cfg.boundary);
    _tangent = _graph.create_node<tangent_point_gradient_node>(
        __R, __Rd, _cfg.tangent);
    _vortex =
        _graph.create_node<vortex_force_node>(__R, __Rd, _cfg.vortex);
    _curve = _graph.create_node<pca_curve_force_node>(__R, _cfg.curve);
    _add = _graph.create_node<vec3_junction_node<4>>();
    _solver_node = _graph.create_node<
        hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>>(
        _config);

    _graph.link(_boundary->output(), _add->template input<0>());
    _graph.link(_tangent->output(), _add->template input<1>());
    _graph.link(_vortex->output(), _add->template input<2>());
    _graph.link(_curve->output(), _add->template input<3>());
    _graph.link(_add->output(), _solver_node->input_at<0>());
  }

  const rod_constraints_config &config() const { return _cfg; }
  const rod::rod &rod() const { return *__R; }
  rod::rod::ptr rod_ptr() const { return __R; }

  void step(int frame) {
    _frame = frame;
    grow_rod_rest_lengths(*__R, _cfg.grow);
    _boundary->set_frame(frame);
    _graph.run();
    __Rd->step();
#if 1 // draw_frames — flip to 0 to disable
    draw_frames();
#endif
  }

private:
  void draw_frames() const {
    const std::vector<vec3> &x = __R->x();
    const std::vector<quat> &u = __R->u();
    const size_t n = std::min(x.size(), u.size());
    for (size_t i = 0; i < n; ++i)
      geometry_logger::frame(u[i].toRotationMatrix(), x[i], _cfg.frame_scale);
  }

  std::vector<vec3> get_fib(real r0, int N = 13) {
    const real golden = 0.5 * (1.0 + sqrt(5));
    std::vector<vec3> cens(N, vec3::Zero());
    for (int i = 0; i < N; i++) {
      const real theta = 2.0 * M_PI * i / golden;
      const real phi = acos(1.0 - 2.0 * (i + 0.5) / real(N));
      cens[i] =
          r0 * vec3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
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
    auto make_ellipse_loop = [&](const vec3 &center, real er0, real er1,
                                 const vec3 &axis0, const vec3 &axis1) {
      std::vector<vec3> pts;
      pts.reserve(N);
      for (int i = 0; i < N; i++) {
        const real thet = 2.0 * M_PI * real(i) / real(N);
        pts.push_back(center + er0 * cos(thet) * axis0 +
                      er1 * sin(thet) * axis1);
      }
      return pts;
    };

    const vec3 cen = vec3::Zero();

    __R = rod::rod::create();
    __R->append_loop(make_ellipse_loop(cen, r0, r1, f0, f1));
    const real lavg = __R->lavg();
    __Rd = rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  rod_constraints_config _cfg;
  int _frame = 0;
  sdf_base::ptr __sdf0;
  sdf_base::ptr __sdf1;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;

  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_quaternion_block::ptr _rod_quat;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>
      _config;
  liblombardi::GraphContext _graph;
  boundary_gradient_node::ptr _boundary;
  tangent_point_gradient_node::ptr _tangent;
  vortex_force_node::ptr _vortex;
  pca_curve_force_node::ptr _curve;
  vec3_junction_node<4>::ptr _add;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                       hepworth::block::rod_quaternion_block>::ptr
      _solver_node;
};

class rod_constraints_solver_adapter
    : public demo_adapter<rod_constraints_solver> {
public:
  static ptr create() {
    return std::make_shared<rod_constraints_solver_adapter>();
  }

  rod_constraints_solver_adapter()
      : demo_adapter<rod_constraints_solver>(
            rod_constraints_solver::create(k_rod_constraints_active),
            "rod_constraints_solver") {}

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), vec3(1.0, 0.0, 0.7));
  }
};

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_ROD_CONSTRAINTS_SOLVER__
