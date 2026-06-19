#pragma once

#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/duchamp/fast_summation_test.hpp"

namespace gaudi {
namespace duchamp {

class fast_summation_demo_adapter : public demo_adapter<fast_summation_test> {
public:
  static ptr create() {
    return std::make_shared<fast_summation_demo_adapter>();
  }

  fast_summation_demo_adapter()
      : demo_adapter<fast_summation_test>(fast_summation_test::create(),
                                          "fast_summation_test") {}

  std::optional<mesh_snapshot> shell_mesh() const override {
    return make_shell_mesh_snapshot(*demo()->__M, vec3(0.72, 0.74, 0.78));
  }
};

} // namespace duchamp
} // namespace gaudi
