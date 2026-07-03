#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/rod_guided_deformation.hpp"

namespace gaudi {
namespace duchamp {

class rod_guided_demo_adapter : public demo_adapter<rod_guided_deformation> {
public:
  static ptr create() {
    return std::make_shared<rod_guided_demo_adapter>();
  }

  rod_guided_demo_adapter()
      : demo_adapter<rod_guided_deformation>(rod_guided_deformation::create(),
                                             "rod_guided_deformation") {}

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.62, 0.66, 0.72));
  }

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(*demo()->__R, vec3(1.0, 0.45, 0.15));
  }
};

} // namespace duchamp
} // namespace gaudi
