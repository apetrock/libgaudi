#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/monge_fairing_demo.hpp"

namespace gaudi {
namespace duchamp {

class monge_fairing_demo_adapter : public demo_adapter<monge_fairing_demo> {
public:
  static ptr create(const monge_fairing_demo_config &cfg =
                        k_monge_fairing_demo_default) {
    return std::make_shared<monge_fairing_demo_adapter>(cfg);
  }

  explicit monge_fairing_demo_adapter(
      const monge_fairing_demo_config &cfg = k_monge_fairing_demo_default)
      : demo_adapter<monge_fairing_demo>(monge_fairing_demo::create(cfg),
                                         "monge_fairing") {}

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.62, 0.66, 0.72));
  }
};

} // namespace duchamp
} // namespace gaudi
