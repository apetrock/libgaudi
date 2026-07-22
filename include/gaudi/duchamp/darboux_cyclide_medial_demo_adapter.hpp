#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/medial_axis_graph_demo.hpp"

namespace gaudi {
namespace duchamp {

class darboux_cyclide_medial_demo_adapter
    : public demo_adapter<medial_axis_graph_demo<>> {
public:
  static ptr create() {
    return std::make_shared<darboux_cyclide_medial_demo_adapter>();
  }

  darboux_cyclide_medial_demo_adapter()
      : demo_adapter<medial_axis_graph_demo<>>(
            medial_axis_graph_demo<>::create(), "darboux_cyclide_medial") {}

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->shell(), vec3(0.72, 0.74, 0.78));
  }
};

} // namespace duchamp
} // namespace gaudi
