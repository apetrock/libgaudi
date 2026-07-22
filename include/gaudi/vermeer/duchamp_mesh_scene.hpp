#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/mesh_snapshot_utils.hpp"
#include "lewitt/gpu_session.hpp"
#include "lewitt/mesh_buffer.hpp"
#include "lewitt/performance.hpp"

namespace gaudi {
namespace vermeer {
namespace detail {

inline void hash_combine(std::size_t &seed, std::size_t value) {
  seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

inline void hash_vec3(std::size_t &seed, const vec3 &value) {
  hash_combine(seed, std::hash<double>{}(value.x()));
  hash_combine(seed, std::hash<double>{}(value.y()));
  hash_combine(seed, std::hash<double>{}(value.z()));
}

inline std::size_t mesh_snapshot_hash(const duchamp::mesh_snapshot &snapshot) {
  LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::mesh_snapshot_hash");
  std::size_t seed = snapshot.positions.size();
  hash_combine(seed, snapshot.indices.size());
  hash_combine(seed, snapshot.colors.size());
  for (const auto &position : snapshot.positions) {
    hash_vec3(seed, position);
  }
  for (const auto index : snapshot.indices) {
    hash_combine(seed, std::hash<uint32_t>{}(index));
  }
  for (const auto &color : snapshot.colors) {
    hash_vec3(seed, color);
  }
  return seed;
}

} // namespace detail

enum class mesh_slot : std::size_t { shell = 0, rod = 1 };

// Scene-owned mesh buffers updated from a Duchamp demo's shell and rod snapshots.
class duchamp_mesh_scene {
public:
  static constexpr std::size_t k_slot_count = 2;

  static duchamp_mesh_scene create(duchamp::demo_trait::ptr demo,
                                   std::size_t mesh_count = k_slot_count) {
    duchamp_mesh_scene scene = create(mesh_count);
    scene._demo = std::move(demo);
    return scene;
  }

  static duchamp_mesh_scene create(std::size_t mesh_count = k_slot_count) {
    duchamp_mesh_scene scene;
    scene._meshes.reserve(mesh_count);
    for (std::size_t i = 0; i < mesh_count; ++i) {
      scene._meshes.push_back(lewitt::mesh_buffer::create());
    }
    scene._last_snapshot_hashes.resize(scene._meshes.size());
    return scene;
  }

  void update(lewitt::gpu_context &ctx,
              const glm::vec3 &default_color = glm::vec3(0.72f, 0.74f, 0.78f)) {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_mesh_scene::update");
    if (!_demo || _meshes.empty()) {
      return;
    }

    if (_meshes.size() > static_cast<std::size_t>(mesh_slot::shell)) {
      sync_slot(mesh_slot::shell, [&]() { return _demo->shell_mesh(); }, ctx.device,
                default_color);
    }
    if (_meshes.size() > static_cast<std::size_t>(mesh_slot::rod)) {
      sync_slot(mesh_slot::rod, [&]() { return _demo->rod_mesh(); }, ctx.device,
                default_color);
    }
  }

  // Upload from a published SceneFrame (render thread; no demo access).
  // Always replace slot contents — each SceneFrame is intentionally new, so
  // skip the expensive full-mesh hash used by update().
  void apply_snapshots(const std::optional<duchamp::mesh_snapshot> &shell,
                       const std::optional<duchamp::mesh_snapshot> &rod,
                       wgpu::Device device,
                       const glm::vec3 &default_color =
                           glm::vec3(0.72f, 0.74f, 0.78f)) {
    LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_mesh_scene::apply_snapshots");
    force_slot(mesh_slot::shell, shell, device, default_color);
    force_slot(mesh_slot::rod, rod, device, default_color);
  }

  const std::vector<lewitt::mesh_buffer::ptr> &meshes() const { return _meshes; }

  std::vector<std::weak_ptr<lewitt::mesh_buffer>> mesh_refs() const {
    return mesh_weak_refs(_meshes);
  }

private:
  template <typename SnapshotFn>
  void sync_slot(mesh_slot slot, SnapshotFn &&snapshot_fn, wgpu::Device device,
                 const glm::vec3 &default_color) {
    const auto index = static_cast<std::size_t>(slot);
    if (index >= _meshes.size()) {
      return;
    }

    auto mesh = _meshes[index];
    if (!mesh) {
      return;
    }

    if (auto snapshot = snapshot_fn()) {
      const auto snapshot_hash = detail::mesh_snapshot_hash(*snapshot);
      if (_last_snapshot_hashes[index] == snapshot_hash) {
        return;
      }
      if (update_mesh_buffer(*mesh, *snapshot, device, default_color)) {
        _last_snapshot_hashes[index] = snapshot_hash;
      } else {
        _last_snapshot_hashes[index].reset();
      }
    } else {
      clear_mesh_buffer(*mesh);
      _last_snapshot_hashes[index].reset();
    }
  }

  void force_slot(mesh_slot slot,
                  const std::optional<duchamp::mesh_snapshot> &snapshot,
                  wgpu::Device device, const glm::vec3 &default_color) {
    const auto index = static_cast<std::size_t>(slot);
    if (index >= _meshes.size())
      return;
    auto mesh = _meshes[index];
    if (!mesh)
      return;
    if (snapshot)
      update_mesh_buffer(*mesh, *snapshot, device, default_color);
    else
      clear_mesh_buffer(*mesh);
    _last_snapshot_hashes[index].reset();
  }

  duchamp::demo_trait::ptr _demo;
  std::vector<lewitt::mesh_buffer::ptr> _meshes;
  std::vector<std::optional<std::size_t>> _last_snapshot_hashes;
};

} // namespace vermeer
} // namespace gaudi
