#pragma once

#include "gaudi/duchamp/braid_guided_deformation.hpp"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"

namespace gaudi {
namespace duchamp {

class braid_guided_demo_adapter : public demo_adapter<braid_guided_deformation> {
public:
  static ptr create(const braid_guided_config &cfg = k_braid_guided_default) {
    return std::make_shared<braid_guided_demo_adapter>(cfg);
  }

  static ptr create(const braid_circle_config &braid_cfg) {
    braid_guided_config cfg;
    cfg.braid = braid_cfg;
    return create(cfg);
  }

  explicit braid_guided_demo_adapter(
      const braid_guided_config &cfg = k_braid_guided_default)
      : demo_adapter<braid_guided_deformation>(
            braid_guided_deformation::create(cfg), "braid_guided_deformation") {}

  void debug_draw() override { demo()->draw_shell_frames(); }

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.62, 0.66, 0.72));
  }

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(*demo()->__R, vec3(1.0, 0.45, 0.15));
  }
};

} // namespace duchamp
} // namespace gaudi
