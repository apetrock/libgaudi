#ifndef __GAUDI_HEPWORTH_ROD_VELOCITY_BLOCK__
#define __GAUDI_HEPWORTH_ROD_VELOCITY_BLOCK__

#include <functional>
#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/block_wrapper.hpp"
#include "gaudi/hepworth/blocks/dof_inputs.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

/// Position DOFs driven by ephemeral target velocity (GD / optimizer steps).
class rod_velocity_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<rod_velocity_block>;
  using vec3_velocity_fn = std::function<std::vector<vec3>()>;
  using input_port = vec3_velocity_port;

  rod_velocity_block(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic)
      : rod(std::move(rod)), dynamic(std::move(dynamic)) {}

  asawa::rod::rod::ptr rod;
  asawa::rod::dynamic::ptr dynamic;
  std::vector<vec3> v_drive;

  rod_velocity_block &with_velocity(vec3_velocity_fn fn) {
    if (fn)
      _velocity_fns.push_back(std::move(fn));
    return *this;
  }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    v_drive.assign(rod->x().size(), vec3::Zero());
    for (const auto &fn : _velocity_fns)
      apply_external_input(fn());
  }

  void apply_external_input(const std::vector<vec3> &input,
                            real scale = 1.0) override {
    if (v_drive.size() < input.size())
      v_drive.resize(input.size(), vec3::Zero());
    const size_t n = std::min(v_drive.size(), input.size());
    for (size_t i = 0; i < n; ++i)
      v_drive[i] += scale * input[i];
  }

  void flush_external_inputs() override {
    v_drive.assign(rod->x().size(), vec3::Zero());
  }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _vec3 = vec3_velocity_block::create(rod->M(), rod->x(), v_drive);
    out.push_back(_vec3);
  }

  vec3_velocity_block::ptr dof_vec3() const { return _vec3; }

private:
  vec3_velocity_block::ptr _vec3;
  std::vector<vec3_velocity_fn> _velocity_fns;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_ROD_VELOCITY_BLOCK__
