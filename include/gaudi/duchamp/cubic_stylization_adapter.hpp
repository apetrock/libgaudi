#pragma once

#include "gaudi/duchamp/cubic_stylization.hpp"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"

namespace gaudi {
namespace duchamp {

class cubic_stylization_adapter : public demo_adapter<cubic_stylization> {
public:
  static ptr create() {
    return std::make_shared<cubic_stylization_adapter>();
  }

  cubic_stylization_adapter()
      : demo_adapter<cubic_stylization>(cubic_stylization::create(),
                                        "cubic_stylization") {}

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.72, 0.62, 0.48));
  }
};

} // namespace duchamp
} // namespace gaudi
