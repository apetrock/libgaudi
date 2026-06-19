#pragma once

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/snapshot_renderables.hpp"

namespace gaudi {
namespace vermeer {

inline void sync_shell_renderables(const duchamp::demo_trait &demo, wgpu::Device device,
                                   const mesh_snapshot_renderable::ptr &forward_mesh,
                                   const mesh_gbuffer_snapshot_renderable::ptr &gbuffer_mesh) {
  if (auto mesh = demo.shell_mesh()) {
    if (forward_mesh) {
      forward_mesh->update(*mesh, device);
    }
    if (gbuffer_mesh) {
      gbuffer_mesh->update(*mesh, device);
    }
  } else {
    if (forward_mesh) {
      forward_mesh->clear_snapshot();
    }
    if (gbuffer_mesh) {
      gbuffer_mesh->clear_snapshot();
    }
  }
}

inline void sync_rod_renderable(const duchamp::demo_trait &demo,
                                const rod_snapshot_renderable::ptr &rod) {
  if (!rod) {
    return;
  }
  if (auto polyline = demo.rod_polyline()) {
    rod->update(*polyline);
  } else {
    rod->clear_snapshot();
  }
}

} // namespace vermeer
} // namespace gaudi
