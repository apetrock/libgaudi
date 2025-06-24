#include "gaudi/logger.h"

#ifdef EMSCRIPTEN

#include "gaudi/wasm_geometry_logger.h"

namespace gaudi {
namespace logger {

void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
  wasm_geometry_logger::line(p0, p1, color);
}

void lines(const std::vector<vec3> &p0, const std::vector<vec3> &p1,
           const std::vector<vec4> &colors) {
  wasm_geometry_logger::lines(p0, p1, colors);
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
  wasm_geometry_logger::box(cen, h, col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
  wasm_geometry_logger::ext(mn, mx, col);
}

void field(const std::vector<vec3> &p, const std::vector<vec3> &dirs, double D) {
  wasm_geometry_logger::field(p, dirs, D);
}

void frame(const mat3 &M, const vec3 &c, double C) {
  wasm_geometry_logger::frame(M, c, C);
}

void point(const vec3 &p0, const vec4 &color) {
  wasm_geometry_logger::point(p0, color);
}

void clear() { wasm_geometry_logger::clear(); }

const std::vector<vec3> &get_lines() {
  return wasm_geometry_logger::get_instance().get_lines();
}

const std::vector<vec4> &get_line_colors() {
  return wasm_geometry_logger::get_instance().get_line_colors();
}

const std::vector<vec3> &get_points() {
  return wasm_geometry_logger::get_instance().get_points();
}

const std::vector<vec4> &get_point_colors() {
  return wasm_geometry_logger::get_instance().get_point_colors();
}

} // namespace logger
} // namespace gaudi

#endif // EMSCRIPTEN
