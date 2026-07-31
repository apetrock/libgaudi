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
void line(const vec3 &p0, const vec3 &p1, const vec4 &color, real radius) {
    (void)radius;
    s_lines.push_back(p0);
    s_lines.push_back(p1);
    s_line_colors.push_back(color);
    s_line_colors.push_back(color);
}

void ext(const vec3 &min, const vec3 &max, const vec4 &col) {
    // Bottom face
    line(vec3(min[0], min[1], min[2]), vec3(max[0], min[1], min[2]), col);
    line(vec3(max[0], min[1], min[2]), vec3(max[0], max[1], min[2]), col);
    line(vec3(max[0], max[1], min[2]), vec3(min[0], max[1], min[2]), col);
    line(vec3(min[0], max[1], min[2]), vec3(min[0], min[1], min[2]), col);

    // Top face
    /*
    line(vec3(min[0], min[1], max[2]), vec3(max[0], min[1], max[2]), col);
    line(vec3(max[0], min[1], max[2]), vec3(max[0], max[1], max[2]), col);
    line(vec3(max[0], max[1], max[2]), vec3(min[0], max[1], max[2]), col);
    line(vec3(min[0], max[1], max[2]), vec3(min[0], min[1], max[2]), col);
    */
    // Vertical edges
    line(vec3(min[0], min[1], min[2]), vec3(min[0], min[1], max[2]), col);
    line(vec3(max[0], min[1], min[2]), vec3(max[0], min[1], max[2]), col);
    line(vec3(max[0], max[1], min[2]), vec3(max[0], max[1], max[2]), col);
    line(vec3(min[0], max[1], min[2]), vec3(min[0], max[1], max[2]), col);
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

void sphere(const vec3 & /*center*/, real /*radius*/, const vec4 & /*color*/) {}

void torus(const vec3 & /*center*/, const vec3 & /*axis*/, real /*major_radius*/,
           real /*minor_radius*/, const vec4 & /*color*/) {}

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
