#ifndef __GAUDI_TEST_ROD_FORCE_GRAPH_TESTS_HPP__
#define __GAUDI_TEST_ROD_FORCE_GRAPH_TESTS_HPP__

#include <memory>
#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/duchamp/modules/rod_forces.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/test/test.hpp"
#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

namespace {

inline std::vector<vec3> make_force_loop_points() {
  std::vector<vec3> points;
  const int N = 32;
  for (int i = 0; i < N; ++i) {
    const real t = 2.0 * M_PI * real(i) / real(N);
    points.emplace_back(cos(t), sin(t), 0.0);
  }
  return points;
}

struct const_vec3_node : public liblombardi::Node {
  using ptr = std::shared_ptr<const_vec3_node>;
  enum class PortId { Output = 0 };
  using OutputPortDef = liblombardi::PortDef<duchamp::field_datum<vec3>, PortId::Output>;

  explicit const_vec3_node(std::vector<vec3> value) : _value(std::move(value)) {}

  void compute() override { get_datum<OutputPortDef>()->data() = _value; }
  unsigned int port_count() const override { return 1; }
  liblombardi::PortRef<const_vec3_node, OutputPortDef> output() { return {*this}; }

  std::vector<vec3> _value;
};

} // namespace

GAUDI_TEST(vec3_junction_node_two_inputs) {
  liblombardi::GraphContext ctx;
  auto a = ctx.create_node<const_vec3_node>(std::vector<vec3>{vec3(1, 2, 3), vec3(4, 5, 6)});
  auto b = ctx.create_node<const_vec3_node>(std::vector<vec3>{vec3(10, 20, 30), vec3(40, 50, 60)});
  auto add = ctx.create_node<duchamp::vec3_junction_node<2>>();

  ctx.link(a->output(), add->input<0>());
  ctx.link(b->output(), add->input<1>());
  ctx.run();

  auto out = add->get_datum<duchamp::vec3_junction_node<2>::OutputPortDef>();
  GAUDI_ASSERT(out->size() == 2);
  GAUDI_ASSERT((out->data()[0] - vec3(11, 22, 33)).norm() < 1e-12);
  GAUDI_ASSERT((out->data()[1] - vec3(44, 55, 66)).norm() < 1e-12);
}

GAUDI_TEST(rod_force_graph_solver_step) {
  auto R = asawa::rod::rod::create(make_force_loop_points());
  const real lavg = R->lavg();
  auto Rd = asawa::rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  auto rod_pos = std::make_shared<hepworth::block::rod_position_block>(R, Rd);
  auto rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(R, Rd);

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>::create()
          .with_blocks(rod_pos, rod_quat)
          .with_bundle(hepworth::block::make_rod_physics_bundle<0, 1>(R, Rd, 1e-1, 1e-1, 1.0))
          .dt(0.05)
          .damping(1.0)
          .iterations(2)
          .build();

  liblombardi::GraphContext ctx;
  auto f0 = ctx.create_node<const_vec3_node>(std::vector<vec3>(R->v().size(), vec3(0.01, 0.0, 0.0)));
  auto f1 = ctx.create_node<const_vec3_node>(std::vector<vec3>(R->v().size(), vec3(0.0, 0.02, 0.0)));
  auto add = ctx.create_node<duchamp::vec3_junction_node<2>>();
  auto solver = ctx.create_node<
      hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>>(config);

  ctx.link(f0->output(), add->input<0>());
  ctx.link(f1->output(), add->input<1>());
  ctx.link(add->output(), solver->input_at<0>());

  const vec3 x0 = R->x()[0];
  ctx.run();
  GAUDI_ASSERT((R->x()[0] - x0).norm() > 1e-12 || R->v()[0].norm() > 1e-12);
}

GAUDI_TEST(rod_force_graph_external_force_drives_motion) {
  auto R = asawa::rod::rod::create(make_force_loop_points());
  const real lavg = R->lavg();
  auto Rd = asawa::rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  auto rod_pos = std::make_shared<hepworth::block::rod_position_block>(R, Rd);
  auto rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(R, Rd);

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>::create()
          .with_blocks(rod_pos, rod_quat)
          .dt(0.05)
          .damping(1.0)
          .iterations(1)
          .build();

  liblombardi::GraphContext ctx;
  const vec3 force(1.0, 0.0, 0.0);
  auto forces =
      ctx.create_node<const_vec3_node>(std::vector<vec3>(R->v().size(), force));
  auto solver = ctx.create_node<
      hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>>(config);

  ctx.link(forces->output(), solver->input_at<0>());

  const vec3 x0 = R->x()[0];
  ctx.run();
  const vec3 delta = R->x()[0] - x0;
  GAUDI_ASSERT(delta.dot(force) > 1e-6);
}

GAUDI_TEST(rod_torque_graph_external_torque_drives_frame) {
  auto R = asawa::rod::rod::create(make_force_loop_points());
  const real lavg = R->lavg();
  auto Rd = asawa::rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  auto rod_pos = std::make_shared<hepworth::block::rod_position_block>(R, Rd);
  auto rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(R, Rd);

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>::create()
          .with_blocks(rod_pos, rod_quat)
          .dt(0.05)
          .damping(0.0)
          .iterations(1)
          .build();

  liblombardi::GraphContext ctx;
  // Graph torque port is world-space (mirrors forces).
  std::vector<vec3> tau(R->u().size(), vec3::Zero());
  tau[0] = R->N2c()[0] * 50.0;
  auto torques = ctx.create_node<const_vec3_node>(tau);
  auto solver = ctx.create_node<
      hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                         hepworth::block::rod_quaternion_block>>(config);

  ctx.link(torques->output(), solver->input_at<1>());

  const quat u0 = R->u()[0];
  ctx.run();
  GAUDI_ASSERT(u0.angularDistance(R->u()[0]) > 1e-6);
}

GAUDI_TEST(rod_torque_local_axial_drives_frame) {
  auto R = asawa::rod::rod::create(make_force_loop_points());
  const real lavg = R->lavg();
  auto Rd = asawa::rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  auto rod_pos = std::make_shared<hepworth::block::rod_position_block>(R, Rd);
  auto rod_quat = std::make_shared<hepworth::block::rod_quaternion_block>(R, Rd);

  // Local API: (M0, M1, M2) about (N0, N1, N2); M2 = axial.
  rod_quat->with_torque_local([&]() {
    std::vector<vec3> tau(R->u().size(), vec3::Zero());
    tau[0] = vec3(0.0, 0.0, 50.0);
    return tau;
  });

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>::create()
          .with_blocks(rod_pos, rod_quat)
          .dt(0.05)
          .damping(0.0)
          .iterations(1)
          .build();

  liblombardi::GraphContext ctx;
  ctx.create_node<hepworth::block::block_solver_node<hepworth::block::rod_position_block,
                                                    hepworth::block::rod_quaternion_block>>(
      config);

  const quat u0 = R->u()[0];
  ctx.run();
  GAUDI_ASSERT(u0.angularDistance(R->u()[0]) > 1e-6);
}

} // namespace test
} // namespace gaudi

#endif
