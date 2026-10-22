#ifndef __GAUDI_HEPWORTH_ROD_POSITION_BLOCK__
#define __GAUDI_HEPWORTH_ROD_POSITION_BLOCK__

#include <algorithm>
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

class rod_position_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<rod_position_block>;
  using vec3_force_fn = std::function<std::vector<vec3>()>;
  using input_port = vec3_forces_port;

  rod_position_block(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic)
      : rod(std::move(rod)), dynamic(std::move(dynamic)) {}

  asawa::rod::rod::ptr rod;
  asawa::rod::dynamic::ptr dynamic;
  std::vector<vec3> forces;

  using vec3_velocity_fn = std::function<std::vector<vec3>()>;

  rod_position_block &with_force(vec3_force_fn fn) {
    if (fn) {
      _force_fns.push_back(std::move(fn));
    }
    return *this;
  }

  /// When set, emit `vec3_velocity_block` (`x += h v`) instead of the force block.
  rod_position_block &with_velocity(vec3_velocity_fn fn) {
    if (fn)
      _velocity_fns.push_back(std::move(fn));
    return *this;
  }

  void set_integrate_velocity(bool on) { _integrate_velocity = on; }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    forces.assign(rod->v().size(), vec3::Zero());
    // Same scale as graph connectors (apply_external_input default = 1).
    for (const auto &fn : _force_fns) {
      apply_external_input(fn());
    }
    if (!_integrate_velocity)
      return;
    v_drive.assign(rod->x().size(), vec3::Zero());
    for (const auto &fn : _velocity_fns) {
      const std::vector<vec3> in = fn();
      if (v_drive.size() < in.size())
        v_drive.resize(in.size(), vec3::Zero());
      const size_t n = std::min(v_drive.size(), in.size());
      for (size_t i = 0; i < n; ++i)
        v_drive[i] += in[i];
    }
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

  void flush_external_inputs() override { forces.assign(rod->v().size(), vec3::Zero()); }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    if (_integrate_velocity) {
      _vel = vec3_velocity_block::create(rod->M(), rod->x(), v_drive);
      out.push_back(_vel);
      return;
    }
    _vec3 = vec3_block::create(rod->M(), rod->x(), rod->v(), forces);
    out.push_back(_vec3);
  }

  vec3_block::ptr dof_vec3() const { return _vec3; }

  /// Target velocity. Live only while `set_integrate_velocity(true)`.
  std::vector<vec3> v_drive;

private:
  vec3_block::ptr _vec3;
  vec3_velocity_block::ptr _vel;
  std::vector<vec3_force_fn> _force_fns;
  std::vector<vec3_velocity_fn> _velocity_fns;
  bool _integrate_velocity = false;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_ROD_POSITION_BLOCK__
