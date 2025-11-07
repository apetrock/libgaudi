#include "gaudi/geometry_logger.hpp"
#include "gaudi/common.h"
#include <iostream>
#include <vector>

namespace gaudi {
namespace geometry_logger {

// Static data storage
static std::vector<vec3> s_lines;
static std::vector<vec4> s_line_colors;
static std::vector<vec3> s_points;
static std::vector<vec4> s_point_colors;

// Geometry visualization functions
void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
    s_lines.push_back(p0);
    s_lines.push_back(p1);
    s_line_colors.push_back(color);
    s_line_colors.push_back(color);
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
    // Create box edges
    vec3 min = cen - h;
    vec3 max = cen + h;
    
    // Bottom face
    line(vec3(min[0], min[1], min[2]), vec3(max[0], min[1], min[2]), col);
    line(vec3(max[0], min[1], min[2]), vec3(max[0], max[1], min[2]), col);
    line(vec3(max[0], max[1], min[2]), vec3(min[0], max[1], min[2]), col);
    line(vec3(min[0], max[1], min[2]), vec3(min[0], min[1], min[2]), col);
    
    // Top face
    line(vec3(min[0], min[1], max[2]), vec3(max[0], min[1], max[2]), col);
    line(vec3(max[0], min[1], max[2]), vec3(max[0], max[1], max[2]), col);
    line(vec3(max[0], max[1], max[2]), vec3(min[0], max[1], max[2]), col);
    line(vec3(min[0], max[1], max[2]), vec3(min[0], min[1], max[2]), col);
    
    // Vertical edges
    line(vec3(min[0], min[1], min[2]), vec3(min[0], min[1], max[2]), col);
    line(vec3(max[0], min[1], min[2]), vec3(max[0], min[1], max[2]), col);
    line(vec3(max[0], max[1], min[2]), vec3(max[0], max[1], max[2]), col);
    line(vec3(min[0], max[1], min[2]), vec3(min[0], max[1], max[2]), col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
    box(0.5 * (mn + mx), 0.5 * (mx - mn), col);
}

void frame(const mat3 &M, const vec3 &c, double C) {
    vec3 x_axis = C * M.col(0);
    vec3 y_axis = C * M.col(1);
    vec3 z_axis = C * M.col(2);
    
    line(c, c + x_axis, vec4(1.0, 0.0, 0.0, 1.0)); // X axis - red
    line(c, c + y_axis, vec4(0.0, 1.0, 0.0, 1.0)); // Y axis - green
    line(c, c + z_axis, vec4(0.0, 0.0, 1.0, 1.0)); // Z axis - blue
}

void point(const vec3 &p0, const vec4 &color) {
    s_points.push_back(p0);
    s_point_colors.push_back(color);
}

void clear() {
    s_lines.clear();
    s_line_colors.clear();
    s_points.clear();
    s_point_colors.clear();
}

// Data access functions
const std::vector<vec3> &get_lines() {
    return s_lines;
}

const std::vector<vec4> &get_line_colors() {
    return s_line_colors;
}

const std::vector<vec3> &get_points() {
    return s_points;
}

const std::vector<vec4> &get_point_colors() {
    return s_point_colors;
}

} // namespace geometry_logger
} // namespace gaudi 