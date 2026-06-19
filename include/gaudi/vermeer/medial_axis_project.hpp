#pragma once

#include <memory>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/duchamp_graph_project.hpp"
#include "lewitt/gbuffer_graph_helpers.hpp"

namespace gaudi {
namespace vermeer {

using medial_gbuffer_visualizer =
    duchamp_graph_project<lewitt::gbuffer_visualizer_pipeline>;

using medial_ssao_deferred = duchamp_graph_project<lewitt::ssao_deferred_pipeline>;

class medial_axis_project : public medial_ssao_deferred {
public:
  explicit medial_axis_project(duchamp::demo_trait::ptr demo)
      : medial_ssao_deferred(std::move(demo)) {}

  bool init(lewitt::gpu_context &ctx) override {
    if (!medial_ssao_deferred::init(ctx)) {
      return false;
    }
    pipeline().set_mode(lewitt::render_pipeline_mode::deferred_lit);
    return true;
  }

  void set_render_mode(lewitt::render_pipeline_mode mode) { pipeline().set_mode(mode); }
};

} // namespace vermeer
} // namespace gaudi
