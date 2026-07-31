#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/field_nodes.hpp"
#include "gaudi/duchamp/modules/rod_forces.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

// Knobs for torque_loop. Add new frozen presets below as combos prove useful.
struct torque_loop_config {
  rod_physics_config physics{
      .stretch_w = 5e-2,
      .bend_w = 1e-1,
      .twist_w = 1e-1,
      .collision_w = 1.0,
      .smooth_w = 0.0,
      .min_kink_w = 0.0,
      .squad_smooth_w = 0.0,
      .dt = 0.05,
      .damping = 0.01,
      .iterations = 10,
  };

  // Dual-scale tangent-point: near-field clump + far-field repel.
  tangent_point_force_config tangent_hf{.w = -1.0e-7, .l0 = 1.0, .p = 6.0};
  tangent_point_force_config tangent_lf{.w = 3.0e-4, .l0 = 6.0, .p = 3.0};
  pca_curve_force_config curve{};

  local_moment_config local{.Mbn = 2.0e-8, .Mn = 0.0, .Mt = 1.3e-8};
  // <0 → apply local to every element; >=0 → only that corner.
  index_t local_i = -1;
  world_space_moment_config world_space{};
  eigen_moment_config eigen_moment{};

  // Throwaway: soften twist on one hinge (corner id); <0 disables.
  // Bend stays at physics.bend_w (helix / writhe dump seam).
  index_t free_hinge_i = 0;
  real free_twist_w = 0.0;

  real frame_scale = 0.08;
};

// Frozen “this one is neat” — designated init so the combo is explicit/editable.
inline const torque_loop_config k_torque_loop_test{
    .physics =
        {
            .stretch_w = 5e-2,
            .bend_w = 8e-2,
            .twist_w = 4e-2,
            .collision_w = 1.0,
            .smooth_w = 0.0,
            .min_kink_w = 0.0,
            .squad_smooth_w = 0.0,
            .dt = 0.05,
            .damping = 0.01,
            .iterations = 10,
        },
    .tangent_hf = {.w = 0.0, .l0 = 1.0, .p = 6.0},
    .tangent_lf = {.w = 0.0, .l0 = 6.0, .p = 3.0},
    .curve = {},
    .local = {.Mbn = 2.0e-7, .Mn = 1.3e-6, .Mt = 1.3e-7},
    .local_i = -1,
    .world_space = {},
    .eigen_moment = {},
    .free_hinge_i = 0,
    .free_twist_w = 0.0,
    .frame_scale = 0.08,
};

// Frozen “this one is neat” — designated init so the combo is explicit/editable.
inline const torque_loop_config k_torque_loop_rope{
    .physics =
        {
            .stretch_w = 5e-2,
            .bend_w = 1e-3,
            .twist_w = 1e-3,
            .collision_w = 1.0,
            .smooth_w = 0.0,
            .min_kink_w = 0.3,
            .squad_smooth_w = 6.0e-3,
            .dt = 0.05,
            .damping = 0.01,
            .iterations = 10,
        },
    .tangent_hf = {.w = -1.0e-7, .l0 = 1.0, .p = 6.0},
    .tangent_lf = {.w = 3.0e-4, .l0 = 6.0, .p = 3.0},
    .curve = {},
    .local = {.Mbn = 5.0e-7, .Mn = 0.0, .Mt = 1.3e-7},
    .local_i = -1,
    .world_space = {},
    .eigen_moment = {},
    .free_hinge_i = 0,
    .free_twist_w = 0.0,
    .frame_scale = 0.08,
};

// Frozen “this one is neat” — designated init so the combo is explicit/editable.
inline const torque_loop_config k_torque_loop_buckle{
    .physics =
        {
            .stretch_w = 5e-2,
            .bend_w = 5e-2,
            .twist_w = 1e-1,
            .collision_w = 1.0,
            .smooth_w = 0.0,
            .min_kink_w = 1.0e-1,
            .squad_smooth_w = 1.0e-4,
            .dt = 0.05,
            .damping = 0.01,
            .iterations = 20,
        },
    .tangent_hf = {.w = 1.0e-8, .l0 = 1.0, .p = 6.0},
    .tangent_lf = {},
    .curve = {.w = 1.0e-1, .half_width = 4},
    .local = {.Mbn = 1.0e-8, .Mn = 0.0, .Mt = 1.0e-8},
    .local_i = -1,
    .world_space = {},
    .eigen_moment = {.w = 2.0e-7, .p = 1.0, .half_width = 5},
    .free_hinge_i = 1,
    .free_twist_w = 1.0e-9,
    .frame_scale = 0.08,
};

// Flip this to pick the active preset (POD copy — no .copy() needed).
// Or mutate a local: torque_loop_config cfg = k_torque_loop_rope; cfg.local.Mt *= 2;
inline const torque_loop_config k_torque_loop_active = k_torque_loop_buckle;

// Circular rod driven by per-element material-frame torques.
//   local (Mbn, Mn, Mt) about (N0, N1, N2) — mostly axial Mt, small off-axis.
//   Optional world_space / eigen_moment via with_torque_world.
//
// Path: with_torque_local → body τ → quat_block::integrate_inertia → projection
class torque_loop {
public:
  using ptr = std::shared_ptr<torque_loop>;
  static ptr create(const torque_loop_config &cfg = k_torque_loop_active) {
    return std::make_shared<torque_loop>(cfg);
  }

  explicit torque_loop(const torque_loop_config &cfg = k_torque_loop_active)
      : _cfg(cfg) {
    load_circle_rod();

    _rod_pos = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);
    _rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);
    _rod_quat->with_torque_local([this]() { return make_torques_local(); });
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
                    phys.collision_w, _cfg.free_hinge_i, _cfg.free_twist_w,
                    phys.edge_stretch_w))
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

    _tangent_hf = _graph.create_node<tangent_point_gradient_node>(
        __R, __Rd, _cfg.tangent_hf);
    _tangent_lf = _graph.create_node<tangent_point_gradient_node>(
        __R, __Rd, _cfg.tangent_lf);
    _curve = _graph.create_node<pca_curve_force_node>(__R, _cfg.curve);
    _tp_add = _graph.create_node<vec3_junction_node<3>>();
    _solver_node = _graph.create_node<
        hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>>(
        _config);
    _graph.link(_tangent_hf->output(), _tp_add->template input<0>());
    _graph.link(_tangent_lf->output(), _tp_add->template input<1>());
    _graph.link(_curve->output(), _tp_add->template input<2>());
    _graph.link(_tp_add->output(), _solver_node->input_at<0>());
  }

  const torque_loop_config &config() const { return _cfg; }
  const asawa::rod::rod &rod() const { return *__R; }

  void step(int frame) {
    (void)frame;
    _graph.run();
    __Rd->step();
  }

private:
  std::vector<vec3> make_torques_local() const {
    const size_t n = __R->u().size();
    const vec3 M(_cfg.local.Mbn, _cfg.local.Mn, _cfg.local.Mt);
    if (_cfg.local_i < 0)
      return std::vector<vec3>(n, M);
    std::vector<vec3> tau(n, vec3::Zero());
    const size_t i = static_cast<size_t>(_cfg.local_i);
    if (i < n)
      tau[i] = M;
    return tau;
  }

  void draw_debug() const {
    const std::vector<vec3> &x = __R->x();
    const std::vector<quat> &u = __R->u();
    const size_t n = std::min(x.size(), u.size());
    for (size_t i = 0; i < n; ++i)
      geometry_logger::frame(u[i].toRotationMatrix(), x[i], _cfg.frame_scale);
  }

  void load_circle_rod() {
    const int N = 256;
    std::vector<vec3> pts;
    pts.reserve(N);
    for (int i = 0; i < N; ++i) {
      const real thet = 2.0 * M_PI * real(i) / real(N);
      pts.emplace_back(std::cos(thet), std::sin(thet), 0.0);
    }
    __R = asawa::rod::rod::create(pts, true);
    const real lavg = __R->lavg();
    __Rd = asawa::rod::dynamic::create(__R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  }

  torque_loop_config _cfg;
  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;
  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_quaternion_block::ptr _rod_quat;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>
      _config;
  liblombardi::GraphContext _graph;
  tangent_point_gradient_node::ptr _tangent_hf;
  tangent_point_gradient_node::ptr _tangent_lf;
  pca_curve_force_node::ptr _curve;
  vec3_junction_node<3>::ptr _tp_add;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                       hepworth::block::rod_quaternion_block>::ptr
      _solver_node;
};

class torque_loop_adapter : public demo_adapter<torque_loop> {
public:
  static ptr create() { return std::make_shared<torque_loop_adapter>(); }
  torque_loop_adapter()
      : demo_adapter<torque_loop>(torque_loop::create(k_torque_loop_active),
                                  "torque_loop") {}
  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), vec3(0.2, 0.85, 0.95));
  }
};

} // namespace duchamp
} // namespace gaudi
