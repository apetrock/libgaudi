#pragma once

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/engman_demo.hpp"

namespace gaudi {
namespace duchamp {

class engman_demo_adapter : public demo_adapter<engman_demo> {
public:
  static ptr create() { return std::make_shared<engman_demo_adapter>(); }

  engman_demo_adapter() : demo_adapter<engman_demo>(engman_demo::create(), "engman") {}

  void debug_draw() override { demo()->draw_debug(); }

  std::optional<mesh_snapshot> shell_mesh() const override { return demo()->snapshot(); }
};

} // namespace duchamp
} // namespace gaudi
