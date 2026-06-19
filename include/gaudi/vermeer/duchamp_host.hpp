#pragma once

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/vermeer/duchamp_render_helpers.hpp"
#include "gaudi/vermeer/snapshot_renderables.hpp"
#include "lewitt/application.h"
#include "lewitt/resource_handles.hpp"

namespace gaudi {
namespace vermeer {

// Thin Vermeer host: shared stepping + snapshot rendering for Duchamp demos.
class duchamp_host {
public:
  explicit duchamp_host(duchamp::demo_trait::ptr demo) : _demo(std::move(demo)) {}

  int run(int width = 1280, int height = 720) {
    lewitt::app_runner app;
    _frame = 0;

    app.set_init_callback([this, &app]() { return on_init(app); });
    app.set_frame_callback([this](uint frame) { on_frame(frame); });
    app.set_before_render_callback([this, &app](uint) { on_before_render(app); });

    if (!app.onInit())
      return 1;

    while (app.isRunning()) {
      app.onFrame(static_cast<uint>(_frame++));
    }

    app.onFinish();
    return 0;
  }

private:
  bool on_init(lewitt::app_runner &app) {
    if (!_demo)
      return false;
    install_snapshot_renderables(app);
    _demo->reset();
    return true;
  }

  void on_frame(uint frame) {
    geometry_logger::clear();
    _demo->step(static_cast<int>(frame));
  }

  void on_before_render(lewitt::app_runner &app) {
    sync_snapshots(app.device());
  }

  void install_snapshot_renderables(lewitt::app_runner &app) {
    _mesh_renderable = mesh_snapshot_renderable::create(app.device());
    _mesh_renderable->set_texture_format(app.color_format(), app.depth_format());
    app.render_scene()->renderables.push_back(_mesh_renderable);

    _mesh_gbuffer_renderable = mesh_gbuffer_snapshot_renderable::create(app.device());
    _mesh_gbuffer_renderable->set_texture_format(
        lewitt::resources::position_attachment::format(),
        lewitt::resources::depth_attachment::format());
    app.add_gbuffer_renderable(_mesh_gbuffer_renderable);

    _rod_renderable = rod_snapshot_renderable::create();
    _rod_renderable->set_texture_format(app.color_format(), app.depth_format());
    app.render_scene()->renderables.push_back(_rod_renderable);
  }

  void sync_snapshots(wgpu::Device device) {
    sync_shell_renderables(*_demo, device, _mesh_renderable, _mesh_gbuffer_renderable);
    sync_rod_renderable(*_demo, _rod_renderable);
  }

  duchamp::demo_trait::ptr _demo;
  mesh_snapshot_renderable::ptr _mesh_renderable;
  mesh_gbuffer_snapshot_renderable::ptr _mesh_gbuffer_renderable;
  rod_snapshot_renderable::ptr _rod_renderable;
  int _frame = 0;
};

} // namespace vermeer
} // namespace gaudi
