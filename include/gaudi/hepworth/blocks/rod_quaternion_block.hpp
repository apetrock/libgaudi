#ifndef __GAUDI_HEPWORTH_ROD_QUATERNION_BLOCK__
#define __GAUDI_HEPWORTH_ROD_QUATERNION_BLOCK__

#include <functional>
#include <tuple>
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

class rod_quaternion_block : public dof_block_wrapper {
public:
  using ptr = std::shared_ptr<rod_quaternion_block>;
  using vec3_torque_fn = std::function<std::vector<vec3>()>;
  using input_port = vec3_torques_port;

  rod_quaternion_block(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic)
      : rod(std::move(rod)), dynamic(std::move(dynamic)) {}

  asawa::rod::rod::ptr rod;
  asawa::rod::dynamic::ptr dynamic;

  // Canonical storage: world-space torques (same convention as forces).
  std::vector<vec3> torques;

  // World-frame API (and graph input_at<1>).
  rod_quaternion_block &with_torque_world(vec3_torque_fn fn) {
    if (fn)
      _world_fns.push_back(std::move(fn));
    return *this;
  }

  // Local / material frame: (M0, M1, M2) about (N0, N1, N2). Builds world and
  // forwards to the world API: τ_world = R * τ_local.
  rod_quaternion_block &with_torque_local(vec3_torque_fn fn) {
    if (!fn)
      return *this;
    _world_fns.push_back([this, fn = std::move(fn)]() {
      const std::vector<vec3> local = fn();
      const std::vector<quat> &u = rod->u();
      std::vector<vec3> world(local.size(), vec3::Zero());
      const size_t n = std::min(local.size(), u.size());
      for (size_t i = 0; i < n; ++i)
        world[i] = u[i] * local[i];
      return world;
    });
    return *this;
  }

  // Alias → local (per-element frame is the usual authoring API).
  rod_quaternion_block &with_torque(vec3_torque_fn fn) {
    return with_torque_local(std::move(fn));
  }

  void prepare(solver_context &ctx) override {
    (void)ctx;
    torques.assign(rod->u().size(), vec3::Zero());
    for (const auto &fn : _world_fns)
      apply_external_input(fn());
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    if (torques.size() < input.size())
      torques.resize(input.size(), vec3::Zero());
    const size_t n = std::min(torques.size(), input.size());
    for (size_t i = 0; i < n; ++i)
      torques[i] += scale * input[i];
  }

  void flush_external_inputs() override { torques.assign(rod->u().size(), vec3::Zero()); }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _quat = quat_block::create(rod->J(), rod->u(), rod->o(), torques);
    out.push_back(_quat);
  }

  quat_block::ptr dof_quat() const { return _quat; }

private:
  quat_block::ptr _quat;
  std::vector<vec3_torque_fn> _world_fns;
};

inline std::tuple<rod_position_block::ptr, rod_quaternion_block::ptr>
make_rod_dof_blocks(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic) {
  return {std::make_shared<rod_position_block>(rod, dynamic),
          std::make_shared<rod_quaternion_block>(rod, dynamic)};
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_ROD_QUATERNION_BLOCK__
