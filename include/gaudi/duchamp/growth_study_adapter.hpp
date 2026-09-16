#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/growth_study.hpp"

namespace gaudi {
namespace duchamp {

class growth_study_adapter : public demo_adapter<growth_study> {
public:
  static ptr create() {
    return std::make_shared<growth_study_adapter>();
  }

  growth_study_adapter()
      : demo_adapter<growth_study>(
            growth_study::create(growth_rx_model::ginzburg_landau),
            "growth_study") {}

  void debug_draw() override { demo()->draw_debug(); }

  std::optional<mesh_snapshot> shell_mesh() const override {
    auto colors = demo()->get_mesh_colors();
    return make_shell_mesh_snapshot(*demo()->__M, colors);
  }
};

} // namespace duchamp
} // namespace gaudi
