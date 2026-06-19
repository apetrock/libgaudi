#pragma once

#include <functional>
#include <memory>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/vermeer/debug_line_utils.hpp"
#include "gaudi/vermeer/duchamp_mesh_scene.hpp"
#include "lewitt/debug_line_buffer.hpp"
#include "lewitt/gpu_session.hpp"

namespace gaudi {
namespace vermeer {

template <typename PipelineT>
class duchamp_graph_project : public lewitt::project_renderer {
public:
  explicit duchamp_graph_project(duchamp::demo_trait::ptr demo)
      : _demo(std::move(demo)), _debug_lines(lewitt::debug_line_buffer::create()) {}

  bool init(lewitt::gpu_context &ctx) override {
    if (!ctx.scene || !_demo) {
      return false;
    }

    _scene = duchamp_mesh_scene::create(_demo);
    _demo->reset();
    _scene.update(ctx);
    sync_debug_lines(*_demo);
    _debug_lines->upload(ctx.device);

    _pipeline = PipelineT::create(ctx, _scene.mesh_refs(), debug_line_weak_refs({_debug_lines}));
    _pipeline.rebind_pool();
    return static_cast<bool>(_pipeline.gbuffer) && static_cast<bool>(_pipeline.sink);
  }

  void update(lewitt::gpu_context &ctx, uint frame) override {
    geometry_logger::clear();
    _demo->step(static_cast<int>(frame));
    _scene.update(ctx);
    sync_debug_lines(*_demo);
    _debug_lines->upload(ctx.device);
  }

  void render(lewitt::gpu_context &ctx,
              const std::function<void(wgpu::RenderPassEncoder &)> &overlay =
                  nullptr) override {
    ctx.scene->update();
    ctx.scene->update_uniforms(ctx.queue);
    _pipeline.set_meshes(_scene.mesh_refs());
    _pipeline.set_debug_lines(debug_line_weak_refs({_debug_lines}));
    _pipeline.render(ctx, overlay);
  }

  void resize(lewitt::gpu_context &ctx) override { _pipeline.resize(ctx); }

  PipelineT &pipeline() { return _pipeline; }
  const PipelineT &pipeline() const { return _pipeline; }

protected:
  void sync_debug_lines(const duchamp::demo_trait &demo) {
    sync_debug_line_buffer_from_rod(demo, *_debug_lines);
    append_debug_line_buffer_from_logger(*_debug_lines);
  }

  duchamp::demo_trait::ptr _demo;
  duchamp_mesh_scene _scene;
  PipelineT _pipeline{};
  lewitt::debug_line_buffer::ptr _debug_lines;
};

} // namespace vermeer
} // namespace gaudi
