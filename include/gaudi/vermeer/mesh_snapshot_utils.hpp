#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include <glm/glm.hpp>
#include <webgpu/webgpu.hpp>

#include "gaudi/duchamp/demo_trait.hpp"
#include "lewitt/mesh_buffer.hpp"
#include "lewitt/performance.hpp"

namespace gaudi {
namespace vermeer {
namespace detail {

inline glm::vec3 to_glm_vec3(const vec3 &v) {
  return glm::vec3(static_cast<float>(v.x()), static_cast<float>(v.y()),
                   static_cast<float>(v.z()));
}

inline std::vector<glm::vec3>
compute_vertex_normals(const std::vector<glm::vec3> &positions,
                       const std::vector<uint32_t> &indices) {
  LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::compute_vertex_normals");
  std::vector<glm::vec3> normals(positions.size(), glm::vec3(0.0f));
  for (size_t i = 0; i + 2 < indices.size(); i += 3) {
    const uint32_t i0 = indices[i + 0];
    const uint32_t i1 = indices[i + 1];
    const uint32_t i2 = indices[i + 2];
    if (i0 >= positions.size() || i1 >= positions.size() ||
        i2 >= positions.size())
      continue;

    const glm::vec3 n =
        glm::cross(positions[i1] - positions[i0], positions[i2] - positions[i0]);
    normals[i0] += n;
    normals[i1] += n;
    normals[i2] += n;
  }

  for (glm::vec3 &n : normals) {
    const float len = glm::length(n);
    n = len > 1.0e-8f ? n / len : glm::vec3(0.0f, 0.0f, 1.0f);
  }
  return normals;
}

} // namespace detail

inline std::vector<glm::vec3>
positions_from_snapshot(const duchamp::mesh_snapshot &snapshot) {
  LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::positions_from_snapshot");
  std::vector<glm::vec3> positions;
  positions.reserve(snapshot.positions.size());
  std::transform(snapshot.positions.begin(), snapshot.positions.end(),
                 std::back_inserter(positions), detail::to_glm_vec3);
  return positions;
}

inline std::vector<lewitt::mesh_buffer::vertex>
mesh_vertices_from_snapshot(const duchamp::mesh_snapshot &snapshot,
                            const glm::vec3 &default_color = glm::vec3(0.72f, 0.74f,
                                                                         0.78f)) {
  LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::mesh_vertices_from_snapshot");
  if (snapshot.positions.empty() || snapshot.indices.empty()) {
    return {};
  }

  const std::vector<glm::vec3> positions = positions_from_snapshot(snapshot);
  const std::vector<glm::vec3> normals =
      detail::compute_vertex_normals(positions, snapshot.indices);

  std::vector<lewitt::mesh_buffer::vertex> vertices;
  vertices.reserve(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) {
    const glm::vec3 color =
        i < snapshot.colors.size() ? detail::to_glm_vec3(snapshot.colors[i])
                                   : default_color;
    vertices.push_back({positions[i], normals[i], color});
  }
  return vertices;
}

inline bool update_mesh_buffer(lewitt::mesh_buffer &mesh,
                               const duchamp::mesh_snapshot &snapshot,
                               wgpu::Device device,
                               const glm::vec3 &default_color = glm::vec3(0.72f, 0.74f,
                                                                          0.78f)) {
  LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::update_mesh_buffer");
  const auto vertices = mesh_vertices_from_snapshot(snapshot, default_color);
  if (vertices.empty()) {
    mesh.clear();
    return false;
  }
  mesh.set_data(vertices, snapshot.indices, device);
  return true;
}

inline void clear_mesh_buffer(lewitt::mesh_buffer &mesh) { mesh.clear(); }

inline std::vector<std::weak_ptr<lewitt::mesh_buffer>>
mesh_weak_refs(const std::vector<lewitt::mesh_buffer::ptr> &meshes) {
  std::vector<std::weak_ptr<lewitt::mesh_buffer>> refs;
  refs.reserve(meshes.size());
  for (const auto &mesh : meshes) {
    refs.push_back(mesh);
  }
  return refs;
}

} // namespace vermeer
} // namespace gaudi
