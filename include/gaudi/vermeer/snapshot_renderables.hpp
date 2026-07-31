#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/mesh_snapshot_utils.hpp"
#include "lewitt/mesh_buffer.hpp"
#include "lewitt/buffer_flags.h"
#include "lewitt/buffers.hpp"
#include "lewitt/doables.hpp"
#include "lewitt/geometry_logger.h"
#include "lewitt/shaders.hpp"

namespace gaudi {
namespace vermeer {

class mesh_snapshot_renderable : public lewitt::doables::forward_renderable {
public:
  using ptr = std::shared_ptr<mesh_snapshot_renderable>;
  using vertex = lewitt::mesh_buffer::vertex;

  static ptr create(wgpu::Device device) {
    return std::make_shared<mesh_snapshot_renderable>(device);
  }

  explicit mesh_snapshot_renderable(wgpu::Device device) {
    vertex_buffer = lewitt::buffers::buffer::create();
    vertex_buffer->set_usage(lewitt::flags::vertex::read);
    vertex_buffer->set_vertex_layout<glm::vec3, glm::vec3, glm::vec3>(
        wgpu::VertexStepMode::Vertex);

    index_buffer = lewitt::buffers::buffer::create();
    index_buffer->set_usage(lewitt::flags::index::read);

    set_shader(lewitt::shaders::render_shader::create_from_path(
        RESOURCE_DIR "/gaudi_snapshot.wgsl", device));
  }

  void update(const duchamp::mesh_snapshot &snapshot, wgpu::Device device) {
    if (snapshot.positions.empty() || snapshot.indices.empty()) {
      _active = false;
      return;
    }

    const auto vertices = mesh_vertices_from_snapshot(snapshot);
    vertex_buffer->write<vertex>(vertices, device);
    index_buffer->write<uint32_t>(snapshot.indices, device);
    _active = true;
  }

  void clear_snapshot() { _active = false; }

  void draw(wgpu::RenderPassEncoder renderpass, wgpu::Device device) override {
    if (!_active || !vertex_buffer || !index_buffer || index_buffer->count() == 0)
      return;
    if (!_layout_prepared) {
      prep_shader_vertex_format();
      _layout_prepared = true;
    }
    init_pipeline(device);
    renderpass.setPipeline(_shader->render_pipe_line());
    renderpass.setVertexBuffer(0, vertex_buffer->get_buffer(), 0,
                               vertex_buffer->get_buffer().getSize());
    renderpass.setBindGroup(0, _bindings->get_group(), 0, nullptr);
    renderpass.setIndexBuffer(index_buffer->get_buffer(),
                              wgpu::IndexFormat::Uint32, 0,
                              index_buffer->size());
    renderpass.drawIndexed(index_buffer->count(), 1, 0, 0, 0);
  }

private:
  bool _active = false;
  bool _layout_prepared = false;
};

class mesh_gbuffer_snapshot_renderable
    : public lewitt::doables::g_buffer_renderable {
public:
  using ptr = std::shared_ptr<mesh_gbuffer_snapshot_renderable>;
  using vertex = lewitt::mesh_buffer::vertex;

  static ptr create(wgpu::Device device) {
    return std::make_shared<mesh_gbuffer_snapshot_renderable>(device);
  }

  explicit mesh_gbuffer_snapshot_renderable(wgpu::Device device) {
    vertex_buffer = lewitt::buffers::buffer::create();
    vertex_buffer->set_usage(lewitt::flags::vertex::read);
    vertex_buffer->set_vertex_layout<glm::vec3, glm::vec3, glm::vec3>(
        wgpu::VertexStepMode::Vertex);

    index_buffer = lewitt::buffers::buffer::create();
    index_buffer->set_usage(lewitt::flags::index::read);

    set_shader(lewitt::shaders::render_shader::create_from_path(
        RESOURCE_DIR "/gaudi_gbuffer.wgsl", device));
  }

  void update(const duchamp::mesh_snapshot &snapshot, wgpu::Device device) {
    if (snapshot.positions.empty() || snapshot.indices.empty()) {
      _active = false;
      return;
    }

    const auto vertices = mesh_vertices_from_snapshot(snapshot);
    vertex_buffer->write<vertex>(vertices, device);
    index_buffer->write<uint32_t>(snapshot.indices, device);
    _active = true;
  }

  void clear_snapshot() { _active = false; }

  void draw(wgpu::RenderPassEncoder renderpass, wgpu::Device device) override {
    if (!_active || !vertex_buffer || !index_buffer || index_buffer->count() == 0)
      return;
    g_buffer_renderable::draw(renderpass, device);
  }

private:
  bool _active = false;
};

class rod_snapshot_renderable : public lewitt::doables::lineable {
public:
  using ptr = std::shared_ptr<rod_snapshot_renderable>;

  static ptr create() { return std::make_shared<rod_snapshot_renderable>(); }

  void update(const duchamp::rod_snapshot &snapshot) {
    clear();
    if (snapshot.positions.size() < 2) {
      _active = false;
      return;
    }

    const glm::vec3 default_color(1.0f, 0.45f, 0.15f);
    // Positions are edge pairs: (p0,p1), (p0,p1), ...
    for (size_t i = 0; i + 1 < snapshot.positions.size(); i += 2) {
      const glm::vec3 color =
          i < snapshot.colors.size() ? detail::to_glm_vec3(snapshot.colors[i])
                                     : default_color;
      add_line({detail::to_glm_vec3(snapshot.positions[i]),
                detail::to_glm_vec3(snapshot.positions[i + 1])},
               color, 0.018f);
    }
    _active = true;
  }

  void clear_snapshot() {
    clear();
    _active = false;
  }

  void draw(wgpu::RenderPassEncoder renderpass, wgpu::Device device) override {
    if (!_active || _p0.empty())
      return;
    lewitt::doables::lineable::draw(renderpass, device);
  }

private:
  bool _active = false;
};

} // namespace vermeer
} // namespace gaudi
