#include "gaudi/logger.hpp"
#include <vector>
#include <cmath>

namespace {
    // Internal storage for the WASM logger implementation
    std::vector<gaudi::vec3> line_positions;
    std::vector<gaudi::vec4> line_colors;
    std::vector<gaudi::vec3> point_positions;
    std::vector<gaudi::vec4> point_colors;
}

namespace gaudi {
namespace logger {

void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
    line_positions.push_back(p0);
    line_positions.push_back(p1);
    line_colors.push_back(color);
    line_colors.push_back(color);
}

void lines(const std::vector<vec3> &p0, const std::vector<vec3> &p1,
           const std::vector<vec4> &colors) {
    size_t count = std::min({p0.size(), p1.size(), colors.size()});
    for (size_t i = 0; i < count; ++i) {
        line(p0[i], p1[i], colors[i]);
    }
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
    // Create a simple box using 12 lines (wireframe)
    vec3 corners[8] = {
        cen + vec3(-h.x(), -h.y(), -h.z()),
        cen + vec3( h.x(), -h.y(), -h.z()),
        cen + vec3( h.x(),  h.y(), -h.z()),
        cen + vec3(-h.x(),  h.y(), -h.z()),
        cen + vec3(-h.x(), -h.y(),  h.z()),
        cen + vec3( h.x(), -h.y(),  h.z()),
        cen + vec3( h.x(),  h.y(),  h.z()),
        cen + vec3(-h.x(),  h.y(),  h.z())
    };
    
    // Bottom face
    line(corners[0], corners[1], col);
    line(corners[1], corners[2], col);
    line(corners[2], corners[3], col);
    line(corners[3], corners[0], col);
    
    // Top face
    line(corners[4], corners[5], col);
    line(corners[5], corners[6], col);
    line(corners[6], corners[7], col);
    line(corners[7], corners[4], col);
    
    // Vertical edges
    line(corners[0], corners[4], col);
    line(corners[1], corners[5], col);
    line(corners[2], corners[6], col);
    line(corners[3], corners[7], col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
    vec3 center = (mn + mx) * 0.5;
    vec3 half_size = (mx - mn) * 0.5;
    box(center, half_size, col);
}

void field(const std::vector<vec3> &p, const std::vector<vec3> &dirs,
           double D) {
    size_t count = std::min(p.size(), dirs.size());
    vec4 field_color(0.8, 0.3, 0.3, 1.0);
    for (size_t i = 0; i < count; ++i) {
        vec3 end = p[i] + dirs[i] * D;
        line(p[i], end, field_color);
    }
}

void frame(const mat3 &M, const vec3 &c, double C) {
    vec4 x_color(1.0, 0.0, 0.0, 1.0);
    vec4 y_color(0.0, 1.0, 0.0, 1.0);
    vec4 z_color(0.0, 0.0, 1.0, 1.0);
    
    vec3 x_axis = c + M.col(0) * C;
    vec3 y_axis = c + M.col(1) * C;
    vec3 z_axis = c + M.col(2) * C;
    
    line(c, x_axis, x_color);
    line(c, y_axis, y_color);
    line(c, z_axis, z_color);
}

void point(const vec3 &p0, const vec4 &color) {
    point_positions.push_back(p0);
    point_colors.push_back(color);
}

void clear() {
    line_positions.clear();
    line_colors.clear();
    point_positions.clear();
    point_colors.clear();
}

const std::vector<vec3> &get_lines() {
    return line_positions;
}

const std::vector<vec4> &get_line_colors() {
    return line_colors;
}

const std::vector<vec3> &get_points() {
    return point_positions;
}

const std::vector<vec4> &get_point_colors() {
    return point_colors;
}

} // namespace logger
} // namespace gaudi
