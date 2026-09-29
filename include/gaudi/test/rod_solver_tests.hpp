#ifndef __GAUDI_TEST_ROD_SOLVER_TESTS_HPP__
#define __GAUDI_TEST_ROD_SOLVER_TESTS_HPP__

#include <cmath>
#include <memory>
#include <tuple>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/rod_quaternion_block.hpp"
#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/constraints/bundles.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace test {

namespace {

std::vector<vec3> make_loop_points() {
  std::vector<vec3> points;
  const int N = 32;
  for (int i = 0; i < N; ++i) {
    const real t = 2.0 * M_PI * real(i) / real(N);
    points.emplace_back(cos(t), sin(t), 0.0);
  }
  return points;
}

std::tuple<hepworth::block::rod_position_block::ptr, hepworth::block::rod_quaternion_block::ptr>
make_loop_rod_dof_blocks() {
  auto R = asawa::rod::rod::create(make_loop_points());
  const real lavg = R->lavg();
  auto Rd = asawa::rod::dynamic::create(R, 0.25 * lavg, 2.5 * lavg, 0.25 * lavg);
  return hepworth::block::make_rod_dof_blocks(R, Rd);
}

} // namespace

GAUDI_TEST(rod_block_emit) {
  auto [rod_pos, rod_quat] = make_loop_rod_dof_blocks();
  hepworth::block::solver_context ctx;
  rod_pos->prepare(ctx);
  rod_quat->prepare(ctx);
  rod_pos->emit_dof_blocks(ctx.blocks);
  rod_quat->emit_dof_blocks(ctx.blocks);
  GAUDI_ASSERT(ctx.blocks.size() == 2);
}

GAUDI_TEST(rod_solver_one_step) {
  auto [rod_pos, rod_quat] = make_loop_rod_dof_blocks();
  auto R = rod_pos->rod;

  rod_pos->forces.assign(R->v().size(), vec3::Zero());
  R->v()[0] = vec3(0.05, 0.0, 0.0);

  auto config =
      hepworth::block::block_solver_builder<hepworth::block::rod_position_block,
                                           hepworth::block::rod_quaternion_block>::create()
          .with_blocks(rod_pos, rod_quat)
          .with_bundle(
              hepworth::block::make_rod_physics_bundle<0, 1>(R, rod_pos->dynamic, 1e-1, 1e-1, 1.0))
          .dt(0.05)
          .damping(1.0)
          .iterations(2)
          .build();

  const vec3 x0 = R->x()[0];
  const vec3 v0 = R->v()[0];
  hepworth::block::projection_solver solver;
  hepworth::block::run_solver_step(config, solver);
  GAUDI_ASSERT((R->x()[0] - x0).norm() > 1e-12 || (R->v()[0] - v0).norm() > 1e-12);
}

GAUDI_TEST(straight_four_frame_sign_flip) {
  auto Rz = [](real a) {
    return quat(Eigen::AngleAxis<real>(a, vec3::UnitZ()));
  };
  quat Q[4] = {Rz(0.0), Rz(0.2), Rz(0.4), Rz(0.6)};
  Q[1].coeffs() *= real(-1);

  quat t1, t2;
  hepworth::block::straight::target_interiors(Q, t1, t2);

  quat e1 = va::slerp(Rz(0.0), Rz(0.6), static_cast<float>(1.0 / 3.0));
  quat e2 = va::slerp(Rz(0.0), Rz(0.6), static_cast<float>(2.0 / 3.0));
  if (e1.coeffs().dot(Q[1].coeffs()) < 0.0)
    e1.coeffs() *= real(-1);
  if (e2.coeffs().dot(Q[2].coeffs()) < 0.0)
    e2.coeffs() *= real(-1);

  GAUDI_ASSERT(t1.normalized().coeffs().dot(e1.normalized().coeffs()) >
               real(0.999));
  GAUDI_ASSERT(t2.normalized().coeffs().dot(e2.normalized().coeffs()) >
               real(0.999));
  GAUDI_ASSERT(t1.coeffs().dot(Q[1].coeffs()) > real(0));
  GAUDI_ASSERT(t2.coeffs().dot(Q[2].coeffs()) > real(0));
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_TEST_ROD_SOLVER_TESTS_HPP__
