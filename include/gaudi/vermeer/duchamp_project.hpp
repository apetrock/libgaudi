#pragma once

#include <memory>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/vermeer/debug_line_utils.hpp"
#include "gaudi/vermeer/duchamp_mesh_scene.hpp"
#include "gaudi/vermeer/vermeer_config.hpp"
#include "lewitt/debug_line_buffer.hpp"
#include "lewitt/gpu_session.hpp"
#include "lewitt/performance.hpp"
#include "lewitt/pipeline_graph.hpp"

namespace gaudi {
namespace vermeer {

class duchamp_project : public lewitt::project_renderer {
public:
  explicit duchamp_project(duchamp::demo_trait::ptr demo, vermeer_config config = {})
      : _demo(std::move(demo)), _config(config),
        _debug_lines(lewitt::debug_line_buffer::create()) {}

  bool init(lewitt::gpu_context &ctx) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::init");
    if (!ctx.scene || !_demo) {
      return false;
    }

    _scene = duchamp_mesh_scene::create(_demo);
    _demo->reset();
    _scene.update(ctx);
    sync_debug_lines(*_demo);
    _debug_lines->upload(ctx.device);

    lewitt::render_graph_config graph_config{};
    graph_config.type = _config.renderer_type;
    graph_config.record = _config.record;
    _graph = lewitt::make_pipeline_graph(ctx, _scene.mesh_refs(),
                                         debug_line_weak_refs({_debug_lines}), graph_config);
    return static_cast<bool>(_graph);
  }

  void update(lewitt::gpu_context &ctx, uint frame) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::update");
    geometry_logger::clear();
    _demo->step(static_cast<int>(frame));
    _scene.update(ctx);
    sync_debug_lines(*_demo);
    _debug_lines->upload(ctx.device);
  }

  void render(lewitt::gpu_context &ctx,
              const std::function<void(wgpu::RenderPassEncoder &)> &overlay =
                  nullptr) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::render");
    (void)overlay;
    ctx.scene->update();
    ctx.scene->update_uniforms(ctx.queue);
    _graph->set_meshes(_scene.mesh_refs());
    _graph->set_debug_lines(debug_line_weak_refs({_debug_lines}));
    _graph->compute(ctx);
  }

  void resize(lewitt::gpu_context &ctx) override { _graph->resize(ctx); }

private:
  void sync_debug_lines(const duchamp::demo_trait &demo) {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::sync_debug_lines");
    sync_debug_line_buffer(demo, *_debug_lines);
  }

  duchamp::demo_trait::ptr _demo;
  vermeer_config _config;
  duchamp_mesh_scene _scene;
  std::unique_ptr<lewitt::pipeline_graph> _graph;
  lewitt::debug_line_buffer::ptr _debug_lines;
};

} // namespace vermeer
} // namespace gaudi
