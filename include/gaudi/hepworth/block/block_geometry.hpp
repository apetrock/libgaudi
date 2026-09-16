#ifndef GAUDI_HEPWORTH_BLOCK_GEOMETRY_HPP
#define GAUDI_HEPWORTH_BLOCK_GEOMETRY_HPP

#include "block_constraint.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

#include <functional>
#include <vector>

namespace gaudi {
namespace hepworth {
namespace block {

inline vecX pack_block_q(const std::vector<sim_block::ptr> &blocks) {
  vecX q;
  if (!blocks.empty() && blocks[0])
    blocks[0]->map_to_x(q);
  return q;
}

inline std::function<vec3(index_t)>
make_block_pos_at(const sim_block::ptr &block, const vecX &q) {
  return [block, q](index_t vi) { return block->get_vec3(vi, q); };
}

inline std::vector<vec3>
gather_shell_positions(const asawa::shell::shell &shell,
                       const sim_block::ptr &block, const vecX &q) {
  std::vector<vec3> x(static_cast<size_t>(shell.vert_count()), vec3::Zero());
  if (!block)
    return x;
  for (int iv = 0; iv < shell.vert_count(); ++iv)
    x[static_cast<size_t>(iv)] = block->get_vec3(iv, q);
  return x;
}

inline void init_block_constraints_rest(
    const std::vector<projection_constraint::ptr> &constraints,
    const vecX &q, solver_context &ctx) {
  for (const auto &c : constraints) {
    if (auto bc = std::dynamic_pointer_cast<block_constraint>(c))
      bc->init_rest(q, ctx);
  }
}

inline void init_block_constraints_rest(
    const std::vector<projection_constraint::ptr> &constraints,
    const std::vector<sim_block::ptr> &blocks, solver_context &ctx) {
  if (blocks.empty())
    return;
  init_block_constraints_rest(constraints, pack_block_q(blocks), ctx);
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_BLOCK_GEOMETRY_HPP
