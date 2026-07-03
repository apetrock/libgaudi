#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/rod_constraints_solver.hpp"

namespace gaudi {
namespace duchamp {

class rod_constraints_solver_adapter
    : public demo_adapter<rod_constraints_solver> {
public:
  static ptr create() {
    return std::make_shared<rod_constraints_solver_adapter>();
  }

  rod_constraints_solver_adapter()
      : demo_adapter<rod_constraints_solver>(rod_constraints_solver::create(),
                                             "rod_constraints_solver") {}

  std::optional<mesh_snapshot> rod_mesh() const override {
    return make_rod_mesh_snapshot(demo()->rod(), vec3(1.0, 0.0, 0.7));
  }
};

} // namespace duchamp
} // namespace gaudi
