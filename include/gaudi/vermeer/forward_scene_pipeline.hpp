#pragma once

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/duchamp_render_helpers.hpp"
#include "gaudi/vermeer/snapshot_renderables.hpp"
#include "lewitt/application.h"
#include "lewitt/resource_handles.hpp"

namespace gaudi {
namespace vermeer {

// Installs snapshot renderables on a render_scene and syncs them from a Duchamp demo.
class forward_scene_pipeline {
public:
  void install(lewitt::app_runner &app, const duchamp::demo_trait::ptr &demo) {
    _demo = demo;
    _use_gbuffer = demo && demo->shell_mesh().has_value();

    _rod_renderable = mesh_snapshot_renderable::create(app.device());
    _rod_renderable->set_texture_format(app.color_format(), app.depth_format());
    app.render_scene()->renderables.push_back(_rod_renderable);

    if (!_use_gbuffer) {
      return;
    }

    _mesh_gbuffer_renderable = mesh_gbuffer_snapshot_renderable::create(app.device());
    _mesh_gbuffer_renderable->set_texture_format(
        lewitt::resources::position_attachment::format(),
        lewitt::resources::depth_attachment::format());
    app.add_gbuffer_renderable(_mesh_gbuffer_renderable);
  }

  bool use_gbuffer() const { return _use_gbuffer; }

  void sync(wgpu::Device device) {
    if (!_demo) {
      return;
    }
    if (_use_gbuffer) {
      sync_shell_renderables(*_demo, device, nullptr, _mesh_gbuffer_renderable);
    }
    sync_rod_mesh_renderable(*_demo, device, _rod_renderable);
  }

private:
  duchamp::demo_trait::ptr _demo;
  bool _use_gbuffer = false;
  mesh_gbuffer_snapshot_renderable::ptr _mesh_gbuffer_renderable;
  mesh_snapshot_renderable::ptr _rod_renderable;
};

} // namespace vermeer
} // namespace gaudi
