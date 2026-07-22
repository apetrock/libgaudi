#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/dipole_tunneling_demo.hpp"

namespace gaudi {
namespace duchamp {

class dipole_tunneling_demo_adapter : public demo_adapter<dipole_tunneling_demo> {
public:
  static ptr create() {
    return std::make_shared<dipole_tunneling_demo_adapter>();
  }

  dipole_tunneling_demo_adapter()
      : demo_adapter<dipole_tunneling_demo>(dipole_tunneling_demo::create(),
                                            "dipole_tunneling") {}


  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.58, 0.64, 0.78));
  }

  std::optional<rod_snapshot> rod_polyline() const override {
    return make_rod_snapshot(*demo()->__R, vec3(1.0, 0.45, 0.15));
  }
};

} // namespace duchamp
} // namespace gaudi
