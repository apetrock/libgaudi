#ifndef __GAUDI_HEPWORTH_ROD_QUATERNION_BLOCK__
#define __GAUDI_HEPWORTH_ROD_QUATERNION_BLOCK__

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
  using input_port = vec3_torques_port;

  rod_quaternion_block(asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic)
      : rod(std::move(rod)), dynamic(std::move(dynamic)) {}

  asawa::rod::rod::ptr rod;
  asawa::rod::dynamic::ptr dynamic;
  std::vector<vec3> torques;

  void prepare(solver_context &ctx) override {
    (void)ctx;
    torques.assign(rod->u().size(), vec3::Zero());
  }

  void apply_external_input(const std::vector<vec3> &input, real scale = 1.0) override {
    if (torques.size() < input.size()) {
      torques.resize(input.size(), vec3::Zero());
    }
    const size_t n = std::min(torques.size(), input.size());
    for (size_t i = 0; i < n; ++i) {
      torques[i] += scale * input[i];
    }
  }

  void flush_external_inputs() override { torques.assign(rod->u().size(), vec3::Zero()); }

  void emit_dof_blocks(std::vector<sim_block::ptr> &out) override {
    _quat = quat_block::create(rod->J(), rod->u(), rod->o(), torques);
    out.push_back(_quat);
  }

  quat_block::ptr dof_quat() const { return _quat; }

private:
  quat_block::ptr _quat;
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
