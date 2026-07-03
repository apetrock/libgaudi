#pragma once

#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "gaudi/duchamp/demo_trait.hpp"
#include "lewitt/debug_line_buffer.hpp"
#include "lewitt/geometry_logger.h"

namespace gaudi {
namespace vermeer {

inline std::vector<lewitt::debug_line_buffer::line>
lines_from_rod_snapshot(const duchamp::rod_snapshot &snapshot) {
  std::vector<lewitt::debug_line_buffer::line> lines;
  if (snapshot.positions.size() < 2) {
    return lines;
  }

  const glm::vec3 default_color(1.0f, 0.45f, 0.15f);
  for (size_t i = 1; i < snapshot.positions.size(); ++i) {
    lewitt::debug_line_buffer::line entry{};
    entry.p0 = glm::vec3(static_cast<float>(snapshot.positions[i - 1].x()),
                         static_cast<float>(snapshot.positions[i - 1].y()),
                         static_cast<float>(snapshot.positions[i - 1].z()));
    entry.p1 = glm::vec3(static_cast<float>(snapshot.positions[i].x()),
                         static_cast<float>(snapshot.positions[i].y()),
                         static_cast<float>(snapshot.positions[i].z()));
    if (i - 1 < snapshot.colors.size()) {
      entry.color = glm::vec3(static_cast<float>(snapshot.colors[i - 1].x()),
                              static_cast<float>(snapshot.colors[i - 1].y()),
                              static_cast<float>(snapshot.colors[i - 1].z()));
    } else {
      entry.color = default_color;
    }
    entry.radius = 0.018f;
    lines.push_back(entry);
  }
  return lines;
}

inline void sync_debug_line_buffer_from_rod(const duchamp::demo_trait &demo,
                                            lewitt::debug_line_buffer &buffer) {
  buffer.clear();
  if (auto polyline = demo.rod_polyline()) {
    buffer.set_lines(lines_from_rod_snapshot(*polyline));
  }
}

inline void append_debug_lines_from_logger(std::vector<lewitt::debug_line_buffer::line> &lines) {
  const auto &logger = lewitt::logger::geometry::get_instance();
  if (!logger.debugLines) {
    return;
  }

  for (const auto &line : logger.debugLines->exported_lines()) {
    lines.push_back({glm::vec3(line.p0.x, line.p0.y, line.p0.z),
                     glm::vec3(line.p1.x, line.p1.y, line.p1.z),
                     glm::vec3(line.color.x, line.color.y, line.color.z), line.radius});
  }
}

inline void append_debug_line_buffer_from_logger(lewitt::debug_line_buffer &buffer) {
  std::vector<lewitt::debug_line_buffer::line> lines;
  append_debug_lines_from_logger(lines);
  for (const auto &line : lines) {
    buffer.add_line(line.p0, line.p1, line.color, line.radius);
  }
}

inline void sync_debug_line_buffer(const duchamp::demo_trait &demo,
                                   lewitt::debug_line_buffer &buffer) {
  std::vector<lewitt::debug_line_buffer::line> lines;
  if (auto polyline = demo.rod_polyline()) {
    lines = lines_from_rod_snapshot(*polyline);
  }
  append_debug_lines_from_logger(lines);
  buffer.set_lines(lines);
}

inline std::vector<std::weak_ptr<lewitt::debug_line_buffer>>
debug_line_weak_refs(const std::vector<lewitt::debug_line_buffer::ptr> &buffers) {
  std::vector<std::weak_ptr<lewitt::debug_line_buffer>> refs;
  refs.reserve(buffers.size());
  for (const auto &buffer : buffers) {
    refs.push_back(buffer);
  }
  return refs;
}

} // namespace vermeer
} // namespace gaudi
