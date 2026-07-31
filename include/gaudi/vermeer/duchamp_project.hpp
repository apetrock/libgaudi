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
#include "gaudi/vermeer/bbox_framing.hpp"
#include "gaudi/vermeer/obj_dump.hpp"
#include "gaudi/vermeer/scene_frame.hpp"
#include "gaudi/vermeer/vermeer_config.hpp"
#include "lewitt/debug_line_buffer.hpp"
#include "lewitt/debug_sphere_buffer.hpp"
#include "lewitt/debug_torus_buffer.hpp"
#include "lewitt/gpu_session.hpp"
#include "lewitt/performance.hpp"
#include "lewitt/pipeline_graph.hpp"

namespace gaudi {
namespace vermeer {

class duchamp_project : public lewitt::project_renderer {
public:
  explicit duchamp_project(duchamp::demo_trait::ptr demo, vermeer_config config = {},
                           std::shared_ptr<duchamp_playback> playback = nullptr)
      : _config(config), _playback(std::move(playback)),
        _runtime(std::move(demo), _playback),
        _debug_lines(lewitt::debug_line_buffer::create()),
        _debug_spheres(lewitt::debug_sphere_buffer::create()),
        _debug_tori(lewitt::debug_torus_buffer::create()),
        _framing(bbox_framing_config{
            .auto_fit = config.fit_bbox_auto,
            .period_s = config.fit_period_s,
            .anisotropy_min = config.fit_anisotropy_min,
            .keep_alignment = config.fit_keep_alignment,
            .min_score_improvement = config.fit_min_score_improvement,
        }) {}

  ~duchamp_project() override { _runtime.stop(); }

  bool init(lewitt::gpu_context &ctx) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::init");
    if (!ctx.scene || !_runtime.demo()) {
      return false;
    }

    _scene = duchamp_mesh_scene::create(duchamp_mesh_scene::k_slot_count);
    _runtime.reset_demo();
    _runtime.start();

    // Wait briefly for the pre-step SceneFrame so the first present isn't empty.
    for (int i = 0; i < 200 && !_last_frame; ++i) {
      if (auto frame = _runtime.channel().take_latest(_seen_generation)) {
        apply_frame(ctx, *frame);
        _last_frame = std::move(frame);
      } else {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
    }

    lewitt::render_graph_config graph_config{};
    graph_config.type = _config.renderer_type;
    graph_config.record = _config.record;
    _graph = lewitt::make_pipeline_graph(ctx, _scene.mesh_refs(),
                                         debug_line_weak_refs({_debug_lines}), graph_config);
    if (_graph && _last_frame)
      _graph->request_record_frame();
    return static_cast<bool>(_graph);
  }

  void update(lewitt::gpu_context &ctx, uint frame) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::update");
    (void)frame;
    if (_playback)
      poll_playback_input(*_playback, ctx.window);

    if (_playback && _playback->toggle_camera_manual.exchange(false))
      _framing.toggle_manual();

    const bool live = ctx.scene && ctx.scene->camera_animating();
    if (auto next = _runtime.channel().take_latest(_seen_generation)) {
      apply_frame(ctx, *next);
      _last_frame = std::move(next);
      if (ctx.scene && ctx.scene->get_camera() && _last_frame) {
        const float sim_t =
            static_cast<float>(_last_frame->sim_frame) / 30.0f;
        _framing.on_capture(*ctx.scene->get_camera(), *_last_frame, sim_t);
        ctx.scene->sync_camera_uniforms(ctx.queue);
      }
      // Sim-synced capture when the camera is idle; live mode handles orbit/zoom.
      if (_graph && !live)
        _graph->request_record_frame();
    }

    if (_playback && _last_frame && _runtime.demo())
      maybe_dump_scene_obj(*_playback, *_last_frame, _runtime.demo()->name());
  }

  void render(lewitt::gpu_context &ctx,
              const std::function<void(wgpu::RenderPassEncoder &)> &overlay =
                  nullptr) override {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::render");
    (void)overlay;
    ctx.scene->update();
    ctx.scene->update_uniforms(ctx.queue);
    // Live capture while orbiting / coasting / scrolling so camera moves land in the dump.
    if (_graph && ctx.scene && ctx.scene->camera_animating())
      _graph->request_record_frame();
    _graph->set_meshes(_scene.mesh_refs());
    _graph->set_debug_lines(debug_line_weak_refs({_debug_lines}));
    _graph->set_debug_spheres(debug_sphere_weak_refs({_debug_spheres}));
    _graph->set_debug_tori(debug_torus_weak_refs({_debug_tori}));
    _graph->compute(ctx);
  }

  void resize(lewitt::gpu_context &ctx) override { _graph->resize(ctx); }

private:
  void apply_frame(lewitt::gpu_context &ctx, const SceneFrame &frame) {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_project::apply_frame");
    _scene.apply_snapshots(frame.shell, frame.rod, ctx.device);
    apply_scene_debug_lines(frame.debug_lines, *_debug_lines);
    apply_scene_debug_spheres(frame.debug_spheres, *_debug_spheres);
    apply_scene_debug_tori(frame.debug_tori, *_debug_tori);
    _debug_lines->upload(ctx.device);
    _debug_spheres->upload(ctx.device);
    _debug_tori->upload(ctx.device);
  }

  vermeer_config _config;
  std::shared_ptr<duchamp_playback> _playback;
  duchamp_sim_runtime _runtime;
  duchamp_mesh_scene _scene;
  std::unique_ptr<lewitt::pipeline_graph> _graph;
  lewitt::debug_line_buffer::ptr _debug_lines;
  lewitt::debug_sphere_buffer::ptr _debug_spheres;
  lewitt::debug_torus_buffer::ptr _debug_tori;
  std::shared_ptr<const SceneFrame> _last_frame;
  uint64_t _seen_generation = 0;
  bbox_framing_controller _framing;
};

} // namespace vermeer
} // namespace gaudi
