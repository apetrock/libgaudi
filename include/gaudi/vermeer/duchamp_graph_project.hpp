#pragma once

#include <chrono>
#include <cstdint>
#include <functional>
#include <memory>
#include <thread>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/debug_line_utils.hpp"
#include "gaudi/vermeer/duchamp_mesh_scene.hpp"
#include "gaudi/vermeer/duchamp_playback.hpp"
#include "gaudi/vermeer/duchamp_sim_runtime.hpp"
#include "gaudi/vermeer/scene_frame.hpp"
#include "lewitt/debug_line_buffer.hpp"
#include "lewitt/gpu_session.hpp"
#include "lewitt/performance.hpp"

namespace gaudi {
namespace vermeer {

template <typename PipelineT>
class duchamp_graph_project : public lewitt::project_renderer {
public:
  explicit duchamp_graph_project(duchamp::demo_trait::ptr demo,
                                 std::shared_ptr<duchamp_playback> playback = nullptr)
      : _playback(std::move(playback)), _runtime(std::move(demo), _playback),
        _debug_lines(lewitt::debug_line_buffer::create()) {}

  ~duchamp_graph_project() override { _runtime.stop(); }

  bool init(lewitt::gpu_context &ctx) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_graph_project::init");
    if (!ctx.scene || !_runtime.demo()) {
      return false;
    }

    _scene = duchamp_mesh_scene::create(duchamp_mesh_scene::k_slot_count);
    _runtime.reset_demo();
    _runtime.start();

    for (int i = 0; i < 200 && !_last_frame; ++i) {
      if (auto frame = _runtime.channel().take_latest(_seen_generation)) {
        apply_frame(ctx, *frame);
        _last_frame = std::move(frame);
      } else {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
    }

    _pipeline = PipelineT::create(ctx, _scene.mesh_refs(), debug_line_weak_refs({_debug_lines}));
    _pipeline.rebind_pool();
    if (_last_frame)
      _pipeline.request_record_frame();
    return static_cast<bool>(_pipeline.gbuffer) && static_cast<bool>(_pipeline.sink);
  }

  void update(lewitt::gpu_context &ctx, uint frame) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_graph_project::update");
    (void)frame;
    if (_playback)
      poll_playback_input(*_playback, ctx.window);

    const bool live = ctx.scene && ctx.scene->camera_animating();
    if (auto next = _runtime.channel().take_latest(_seen_generation)) {
      apply_frame(ctx, *next);
      _last_frame = std::move(next);
      if (!live)
        _pipeline.request_record_frame();
    }
  }

  void render(lewitt::gpu_context &ctx,
              const std::function<void(wgpu::RenderPassEncoder &)> &overlay =
                  nullptr) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_graph_project::render");
    ctx.scene->update();
    ctx.scene->update_uniforms(ctx.queue);
    if (ctx.scene && ctx.scene->camera_animating())
      _pipeline.request_record_frame();
    _pipeline.set_meshes(_scene.mesh_refs());
    _pipeline.set_debug_lines(debug_line_weak_refs({_debug_lines}));
    _pipeline.render(ctx, overlay);
  }

  void resize(lewitt::gpu_context &ctx) override { _pipeline.resize(ctx); }

  PipelineT &pipeline() { return _pipeline; }
  const PipelineT &pipeline() const { return _pipeline; }

protected:
  void apply_frame(lewitt::gpu_context &ctx, const SceneFrame &frame) {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_graph_project::apply_frame");
    _scene.apply_snapshots(frame.shell, frame.rod, ctx.device);
    apply_scene_debug_lines(frame.debug_lines, *_debug_lines);
    _debug_lines->upload(ctx.device);
  }

  std::shared_ptr<duchamp_playback> _playback;
  duchamp_sim_runtime _runtime;
  duchamp_mesh_scene _scene;
  PipelineT _pipeline{};
  lewitt::debug_line_buffer::ptr _debug_lines;
  std::shared_ptr<const SceneFrame> _last_frame;
  uint64_t _seen_generation = 0;
};

} // namespace vermeer
} // namespace gaudi
