#include "gaudi/geometry_logger.hpp"

#include <array>
#include <vector>

#include "lewitt/geometry_logger.h"

namespace gaudi {
namespace geometry_logger {
namespace {

constexpr float kDebugLineRadius = 0.01f;
constexpr float kDebugPointRadius = 0.025f;

glm::vec3 to_glm_vec3(const vec3 &v) {
  return glm::vec3(static_cast<float>(v.x()), static_cast<float>(v.y()),
                   static_cast<float>(v.z()));
}

glm::vec3 to_glm_color(const vec4 &v) {
  return glm::vec3(static_cast<float>(v.x()), static_cast<float>(v.y()),
                   static_cast<float>(v.z()));
}

} // namespace

void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
  const std::array<glm::vec3, 2> segment = {to_glm_vec3(p0), to_glm_vec3(p1)};
  lewitt::logger::geometry::line(segment, to_glm_color(color), kDebugLineRadius);
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
  const vec3 mn = cen - h;
  const vec3 mx = cen + h;

  line(vec3(mn[0], mn[1], mn[2]), vec3(mx[0], mn[1], mn[2]), col);
  line(vec3(mx[0], mn[1], mn[2]), vec3(mx[0], mx[1], mn[2]), col);
  line(vec3(mx[0], mx[1], mn[2]), vec3(mn[0], mx[1], mn[2]), col);
  line(vec3(mn[0], mx[1], mn[2]), vec3(mn[0], mn[1], mn[2]), col);

  line(vec3(mn[0], mn[1], mx[2]), vec3(mx[0], mn[1], mx[2]), col);
  line(vec3(mx[0], mn[1], mx[2]), vec3(mx[0], mx[1], mx[2]), col);
  line(vec3(mx[0], mx[1], mx[2]), vec3(mn[0], mx[1], mx[2]), col);
  line(vec3(mn[0], mx[1], mx[2]), vec3(mn[0], mn[1], mx[2]), col);

  line(vec3(mn[0], mn[1], mn[2]), vec3(mn[0], mn[1], mx[2]), col);
  line(vec3(mx[0], mn[1], mn[2]), vec3(mx[0], mn[1], mx[2]), col);
  line(vec3(mx[0], mx[1], mn[2]), vec3(mx[0], mx[1], mx[2]), col);
  line(vec3(mn[0], mx[1], mn[2]), vec3(mn[0], mx[1], mx[2]), col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
  box(0.5 * (mn + mx), 0.5 * (mx - mn), col);
}

void frame(const mat3 &M, const vec3 &c, double C) {
  line(c, c + C * M.col(0), vec4(1.0, 0.0, 0.0, 1.0));
  line(c, c + C * M.col(1), vec4(0.0, 1.0, 0.0, 1.0));
  line(c, c + C * M.col(2), vec4(0.0, 0.0, 1.0, 1.0));
}

void point(const vec3 &p0, const vec4 &color) {
  lewitt::logger::geometry::point(to_glm_vec3(p0), to_glm_color(color),
                                  kDebugPointRadius);
}

void clear() {
  lewitt::logger::geometry::clear();
}

const std::vector<vec3> &get_lines() {
  static const std::vector<vec3> empty;
  return empty;
}

const std::vector<vec4> &get_line_colors() {
  static const std::vector<vec4> empty;
  return empty;
}

const std::vector<vec3> &get_points() {
  static const std::vector<vec3> empty;
  return empty;
}

const std::vector<vec4> &get_point_colors() {
  static const std::vector<vec4> empty;
  return empty;
}

} // namespace geometry_logger
} // namespace gaudi
