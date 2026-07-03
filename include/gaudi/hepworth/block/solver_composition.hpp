#ifndef __GAUDI_HEPWORTH_SOLVER_COMPOSITION__
#define __GAUDI_HEPWORTH_SOLVER_COMPOSITION__

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

struct solver_context {
  std::vector<sim_block::ptr> blocks;
  std::vector<projection_constraint::ptr> constraints;

  real dt = 0.01;
  real damping = 0.5;
  int iterations = 10;
};

using presolve_fn = std::function<void(solver_context &)>;
using constraint_recompute_fn = std::function<void(solver_context &)>;
using constraint_bundle = std::vector<constraint_recompute_fn>;

template <typename... DofBlocks>
struct block_solver_config {
  std::tuple<std::shared_ptr<DofBlocks>...> dof_blocks;
  std::vector<presolve_fn> presolve;
  std::vector<constraint_recompute_fn> recompute;

  real dt = 0.01;
  real damping = 0.5;
  int iterations = 10;

  template <size_t I>
  auto block_at() const {
    return std::get<I>(dof_blocks);
  }
};

namespace detail {

template <typename... DofBlocks, typename Fn, size_t... Is>
inline void for_each_dof_block(block_solver_config<DofBlocks...> &config, Fn &&fn,
                               std::index_sequence<Is...>) {
  (fn(std::get<Is>(config.dof_blocks), Is), ...);
}

} // namespace detail

template <typename... DofBlocks>
inline void flush_external_inputs(block_solver_config<DofBlocks...> &config) {
  detail::for_each_dof_block(
      config,
      [&](auto &block, size_t) {
        if (block) {
          block->flush_external_inputs();
        }
      },
      std::index_sequence_for<DofBlocks...>{});
}

template <typename... DofBlocks>
inline void prepare_solver_step(block_solver_config<DofBlocks...> &config,
                                solver_context &ctx) {
  ctx.dt = config.dt;
  ctx.damping = config.damping;
  ctx.iterations = config.iterations;
  ctx.blocks.clear();

  detail::for_each_dof_block(
      config,
      [&](auto &block, size_t) {
        if (block) {
          block->prepare(ctx);
        }
      },
      std::index_sequence_for<DofBlocks...>{});

  for (const auto &fn : config.presolve) {
    if (fn) {
      fn(ctx);
    }
  }
}

template <typename... DofBlocks>
inline void solve_solver_step(block_solver_config<DofBlocks...> &config,
                              projection_solver &solver, solver_context &ctx) {
  detail::for_each_dof_block(
      config,
      [&](auto &block, size_t) {
        if (block) {
          block->emit_dof_blocks(ctx.blocks);
        }
      },
      std::index_sequence_for<DofBlocks...>{});

  ctx.constraints.clear();
  for (const auto &fn : config.recompute) {
    if (fn) {
      fn(ctx);
    }
  }

  solver.set_constraints(ctx.constraints);
  solver.step(ctx.blocks, ctx.dt, ctx.damping, ctx.iterations);
}

template <typename... DofBlocks>
inline void run_solver_step(block_solver_config<DofBlocks...> &config,
                            projection_solver &solver) {
  solver_context ctx;
  prepare_solver_step(config, ctx);
  solve_solver_step(config, solver, ctx);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_SOLVER_COMPOSITION__
