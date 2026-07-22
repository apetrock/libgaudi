#ifndef __GAUDI_TEST_SOLVER_NODE_TESTS_HPP__
#define __GAUDI_TEST_SOLVER_NODE_TESTS_HPP__

#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/block_solver_node.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/test/test.hpp"

#include "liblombardi/graph_context.hpp"

namespace gaudi {
namespace test {

namespace {

hepworth::block::shell_position_block::ptr make_shell_block_setup() {
  auto M = asawa::shell::load_cube();
  asawa::init_vert_datum<vec3>(*M, vec3::Zero());
  auto xs = std::make_shared<duchamp::shell_vert_positions>(M, 0);
  auto vs = std::make_shared<duchamp::shell_vert_velocities>(M, 1);
  return std::make_shared<hepworth::block::shell_position_block>(M, xs, vs);
}

} // namespace

GAUDI_TEST(solver_builder_smoke) {
  auto shell = make_shell_block_setup();

  auto config = hepworth::block::block_solver_builder<hepworth::block::shell_position_block>::create()
                    .with_blocks(shell)
                    .with_recompute(
                        hepworth::block::make_shell_bending_recompute<0>(shell, 1.0))
                    .dt(0.01)
                    .iterations(2)
                    .build();

  GAUDI_ASSERT(std::tuple_size<decltype(config.dof_blocks)>::value == 1);
  GAUDI_ASSERT(!config.recompute.empty());
}

GAUDI_TEST(solver_presolve_and_recompute_called) {
  auto shell = make_shell_block_setup();

  bool presolve_called = false;
  bool recompute_called = false;
  size_t constraint_count = 0;

  auto config = hepworth::block::block_solver_builder<hepworth::block::shell_position_block>::create()
                    .with_blocks(shell)
                    .with_presolve([&](hepworth::block::solver_context &) {
                      presolve_called = true;
                    })
                    .with_recompute([&](hepworth::block::solver_context &ctx) {
                      recompute_called = true;
                      hepworth::block::make_shell_bending_recompute<0>(shell, 1.0)(ctx);
                      constraint_count = ctx.constraints.size();
                    })
                    .iterations(1)
                    .build();

  hepworth::block::projection_solver solver;
  hepworth::block::run_solver_step(config, solver);

  GAUDI_ASSERT(presolve_called);
  GAUDI_ASSERT(recompute_called);
  GAUDI_ASSERT(constraint_count > 0);
  GAUDI_ASSERT(std::tuple_size<decltype(config.dof_blocks)>::value == 1);
}

GAUDI_TEST(solver_node_one_step) {
  auto shell = make_shell_block_setup();

  std::vector<vec3> &velocities = shell->vs->get();
  velocities[0] = vec3(0.1, 0.0, 0.0);

  auto config = hepworth::block::block_solver_builder<hepworth::block::shell_position_block>::create()
                    .with_blocks(shell)
                    .with_recompute(
                        hepworth::block::make_shell_bending_recompute<0>(shell, 10.0))
                    .dt(0.01)
                    .damping(0.5)
                    .iterations(5)
                    .build();

  liblombardi::GraphContext ctx;
  ctx.create_node<hepworth::block::block_solver_node<hepworth::block::shell_position_block>>(
      config);
  ctx.run();

  GAUDI_ASSERT(velocities[0].norm() > 1e-12 ||
               shell->xs->get()[0].norm() > 1e-12);
}

GAUDI_TEST(solver_shell_physics_bundle_size) {
  auto shell = make_shell_block_setup();
  auto bundle = hepworth::block::make_shell_physics_bundle<0>(shell, 1.0, 0.1);
  GAUDI_ASSERT(bundle.size() == 2);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_SOLVER_NODE_TESTS_HPP__
