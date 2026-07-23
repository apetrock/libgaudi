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
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace duchamp {

// Circular rod driven by per-element material-frame torques.
//   (Mbn, Mn, Mt) about (N0, N1, N2) — mostly axial Mt, small off-axis perturbation.
//
// Path: with_torque_local → body τ → quat_block::integrate_inertia → projection
class torque_loop {
public:
  using ptr = std::shared_ptr<torque_loop>;
  static ptr create() { return std::make_shared<torque_loop>(); }

  real stretch_w = 1e-2;
  real bend_w = 3e-2;
  real twist_w = 1e-3;
  real smooth_w = 0.1;
  real sym_bezier_w = 1e-2;
  real collision_w = 1.0;
  real dt = 0.05;
  real damping = 0.01;

  // Material-frame moments about (N0=binormal, N1=normal, N2=tangent).
  real Mt = 0.0e-8;   // axial
  real Mn = 1.0e-7;   // small normal perturbation (~5% of Mt)
  real Mbn = 3.0e-7;   // binormal

  real frame_scale = 0.08;

  torque_loop() {
    load_circle_rod();

    _rod_pos = std::make_shared<hepworth::block::rod_position_block>(__R, __Rd);
    _rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(__R, __Rd);
    _rod_quat->with_torque_local([this]() { return make_torques_local(); });

    _config =
        hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                              hepworth::block::rod_quaternion_block>::
            create()
                .with_blocks(_rod_pos, _rod_quat)
                .with_bundle(hepworth::block::make_rod_physics_bundle<0, 1>(
                    __R, __Rd, stretch_w, bend_w, twist_w, collision_w))
                .with_recompute(
                    hepworth::block::make_rod_smooth_recompute<0>(__R, smooth_w))
                .with_recompute(
                    hepworth::block::make_rod_symmetric_bezier_recompute<0>(
                        __R, sym_bezier_w))
                .dt(dt)
                .damping(damping)
                .build();

    _solver_node = _graph.create_node<
        hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>>(
        _config);
  }

  const asawa::rod::rod &rod() const { return *__R; }

  void step(int frame) {
    (void)frame;
    _graph.run();
    __Rd->step();
    //draw_debug();
  }

private:
  // Body-frame components (M0, M1, M2) = (Mbn, Mn, Mt).
  std::vector<vec3> make_torques_local() const {
    const size_t n = __R->u().size();
    return std::vector<vec3>(n, vec3(Mbn, Mn, Mt));
  }

  void draw_debug() const {
    const std::vector<vec3> &x = __R->x();
    const std::vector<quat> &u = __R->u();
    const size_t n = std::min(x.size(), u.size());
    for (size_t i = 0; i < n; ++i)
      geometry_logger::frame(u[i].toRotationMatrix(), x[i], frame_scale);
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

  asawa::rod::rod::ptr __R;
  asawa::rod::dynamic::ptr __Rd;
  hepworth::block::rod_position_block::ptr _rod_pos;
  hepworth::block::rod_quaternion_block::ptr _rod_quat;
  hepworth::block::block_solver_config<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>
      _config;
  liblombardi::GraphContext _graph;
  hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                       hepworth::block::rod_quaternion_block>::ptr
      _solver_node;
};

class torque_loop_adapter : public demo_adapter<torque_loop> {
public:
  static ptr create() { return std::make_shared<torque_loop_adapter>(); }
  torque_loop_adapter()
      : demo_adapter<torque_loop>(torque_loop::create(), "torque_loop") {}
  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), vec3(0.2, 0.85, 0.95));
  }
};

} // namespace duchamp
} // namespace gaudi
