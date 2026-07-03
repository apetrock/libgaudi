#ifndef __GAUDI_HEPWORTH_SHELL_POSITION_BLOCK__
#define __GAUDI_HEPWORTH_SHELL_POSITION_BLOCK__

#include <memory>
#include <vector>

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/hepworth/blocks/dof_inputs.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

class shell_position_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<shell_position_block>;
  using input_port = vec3_forces_port;
  using positions_ptr = std::shared_ptr<duchamp::shell_vert_positions>;
  using velocities_ptr = std::shared_ptr<duchamp::shell_vert_velocities>;

  shell_position_block(asawa::shell::shell::ptr mesh, positions_ptr positions,
                       velocities_ptr velocities)
      : mesh(std::move(mesh)), xs(std::move(positions)), vs(std::move(velocities)) {}

  asawa::shell::shell::ptr mesh;
  positions_ptr xs;
  velocities_ptr vs;
  std::vector<vec3> mass;
  std::vector<vec3> forces;

  void prepare(solver_context &ctx) override {
    (void)ctx;
    mass = asawa::shell::vertex_areas_3(*mesh, xs->get());
    forces.assign(xs->size(), vec3::Zero());
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    if (forces.size() < input.size()) {
      forces.resize(input.size(), vec3::Zero());
    }
    const size_t n = std::min(forces.size(), input.size());
    for (size_t i = 0; i < n; ++i) {
      forces[i] += scale * input[i];
    }
  }

  void flush_external_inputs() override { forces.assign(xs->size(), vec3::Zero()); }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _vec3 = vec3_block::create(mass, xs->get(), vs->get(), forces);
    out.push_back(_vec3);
  }

  vec3_block::ptr dof_vec3() const { return _vec3; }

private:
  vec3_block::ptr _vec3;
};

using shell_block = shell_position_block;

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_SHELL_POSITION_BLOCK__
