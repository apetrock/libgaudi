#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/slab_prism_bending.hpp"

namespace gaudi {
namespace duchamp {

class slab_prism_bending_adapter : public demo_adapter<slab_prism_bending> {
public:
  static ptr create() { return std::make_shared<slab_prism_bending_adapter>(); }

  slab_prism_bending_adapter()
      : demo_adapter<slab_prism_bending>(slab_prism_bending::create(),
                                         "slab_prism_bending") {}

  void debug_draw() override { demo()->draw_debug(); }

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.82, 0.72, 0.55));
  }
};

} // namespace duchamp
} // namespace gaudi
