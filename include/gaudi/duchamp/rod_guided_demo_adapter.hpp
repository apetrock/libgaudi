#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/rod_guided_deformation.hpp"

namespace gaudi {
namespace duchamp {

class rod_guided_demo_adapter : public demo_adapter<rod_guided_deformation> {
public:
  static ptr create(const rod_guided_config &cfg = {}) {
    return std::make_shared<rod_guided_demo_adapter>(cfg);
  }

  // Convenience overloads.
  static ptr create(const braid_circle_config &braid_cfg) {
    rod_guided_config cfg;
    cfg.scene = rod_guided_scene::braid_circle;
    cfg.braid = braid_cfg;
    return create(cfg);
  }

  static ptr create_bunny_walk() {
    rod_guided_config cfg;
    cfg.scene = rod_guided_scene::bunny_walk;
    return create(cfg);
  }

  explicit rod_guided_demo_adapter(const rod_guided_config &cfg = {})
      : demo_adapter<rod_guided_deformation>(rod_guided_deformation::create(cfg),
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
