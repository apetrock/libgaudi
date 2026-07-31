#include "../include/wasm_geometry_logger.h"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/common.h"
#include <vector>

namespace gaudi {
namespace geometry_logger {

void line(const vec3 &p0, const vec3 &p1, const vec4 &color, real radius) {
    (void)radius;
    wasm_geometry_logger::line(p0, p1, color);
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
    wasm_geometry_logger::box(cen, h, col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
    wasm_geometry_logger::ext(mn, mx, col);
}

void frame(const mat3 &M, const vec3 &c, double C) {
    wasm_geometry_logger::frame(M, c, C);
}

void point(const vec3 &p0, const vec4 &color) {
    wasm_geometry_logger::point(p0, color);
}

void sphere(const vec3 & /*center*/, real /*radius*/, const vec4 & /*color*/) {}

void torus(const vec3 & /*center*/, const vec3 & /*axis*/, real /*major_radius*/,
           real /*minor_radius*/, const vec4 & /*color*/) {}

void clear() {
    wasm_geometry_logger::clear();
}

// Data access functions - these need to be implemented to return actual data from wasm_geometry_logger
const std::vector<vec3> &get_lines() {
    // TODO: Implement to return data from wasm_geometry_logger singleton
    static std::vector<vec3> empty;
    return empty;
}

const std::vector<vec4> &get_line_colors() {
    // TODO: Implement to return data from wasm_geometry_logger singleton
    static std::vector<vec4> empty;
    return empty;
}

const std::vector<vec3> &get_points() {
    // TODO: Implement to return data from wasm_geometry_logger singleton
    static std::vector<vec3> empty;
    return empty;
}

const std::vector<vec4> &get_point_colors() {
    // TODO: Implement to return data from wasm_geometry_logger singleton
    static std::vector<vec4> empty;
    return empty;
}

} // namespace geometry_logger
} // namespace gaudi 