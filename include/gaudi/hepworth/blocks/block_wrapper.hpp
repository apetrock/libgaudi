#ifndef __GAUDI_HEPWORTH_BLOCK_WRAPPER__
#define __GAUDI_HEPWORTH_BLOCK_WRAPPER__

#include <memory>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/hepworth/block/sim_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

struct solver_context;

struct dof_block_wrapper {
  virtual ~dof_block_wrapper() = default;
  virtual void prepare(solver_context &ctx) = 0;
  virtual void emit_dof_blocks(std::vector<sim_block::ptr> &out) = 0;
  virtual void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) = 0;
  virtual void flush_external_inputs() = 0;
};

using dof_block_ptr = std::shared_ptr<dof_block_wrapper>;

using block_wrapper = dof_block_wrapper;
using block_ptr = dof_block_ptr;

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_BLOCK_WRAPPER__
