#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/medial_axis_estimation_demo.hpp"

namespace gaudi {
namespace duchamp {

class medial_axis_estimation_adapter
    : public demo_adapter<medial_axis_estimation_demo> {
public:
  static ptr create() {
    return std::make_shared<medial_axis_estimation_adapter>();
  }

  medial_axis_estimation_adapter()
      : demo_adapter<medial_axis_estimation_demo>(
            medial_axis_estimation_demo::create(), "medial_axis_estimation") {}

  void debug_draw() override { demo()->draw_debug(); }

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->shell(), vec3(0.72, 0.74, 0.78));
  }
};

} // namespace duchamp
} // namespace gaudi
