#pragma once

#include <memory>
#include <vector>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/mesh_snapshot_utils.hpp"
#include "lewitt/gpu_session.hpp"
#include "lewitt/mesh_buffer.hpp"

namespace gaudi {
namespace vermeer {

// Scene-owned mesh buffers updated from a Duchamp demo's shell snapshot.
class duchamp_mesh_scene {
public:
  static duchamp_mesh_scene create(duchamp::demo_trait::ptr demo,
                                   std::size_t mesh_count = 1) {
    duchamp_mesh_scene scene;
    scene._demo = std::move(demo);
    scene._meshes.reserve(mesh_count);
    for (std::size_t i = 0; i < mesh_count; ++i) {
      scene._meshes.push_back(lewitt::mesh_buffer::create());
    }
    return scene;
  }

  void update(lewitt::gpu_context &ctx,
              const glm::vec3 &default_color = glm::vec3(0.72f, 0.74f, 0.78f)) {
    if (!_demo || _meshes.empty()) {
      return;
    }

    auto mesh = _meshes.front();
    if (!mesh) {
      return;
    }

    if (auto snapshot = _demo->shell_mesh()) {
      update_mesh_buffer(*mesh, *snapshot, ctx.device, default_color);
    } else {
      clear_mesh_buffer(*mesh);
    }
  }

  const std::vector<lewitt::mesh_buffer::ptr> &meshes() const { return _meshes; }

  std::vector<std::weak_ptr<lewitt::mesh_buffer>> mesh_refs() const {
    return mesh_weak_refs(_meshes);
  }

private:
  duchamp::demo_trait::ptr _demo;
  std::vector<lewitt::mesh_buffer::ptr> _meshes;
};

} // namespace vermeer
} // namespace gaudi
