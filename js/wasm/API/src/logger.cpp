#include "../include/wasm_logger.h"

namespace gaudi {
namespace logger {

void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
    wasm_logger::line(p0, p1, color);
}

void box(const vec3 &cen, const vec3 &h, const vec4 &col) {
    wasm_logger::box(cen, h, col);
}

void ext(const vec3 &mn, const vec3 &mx, const vec4 &col) {
    wasm_logger::ext(mn, mx, col);
}

void frame(const mat3 &M, const vec3 &c, double C) {
    wasm_logger::frame(M, c, C);
}

void point(const vec3 &p0, const vec4 &color) {
    wasm_logger::point(p0, color);
}

void clear() {
    wasm_logger::clear();
}

const std::vector<vec3> &get_lines() {
    // This would need to be implemented to return the actual line data
    // For now, return an empty vector
    static std::vector<vec3> empty;
    return empty;
}

const std::vector<vec4> &get_line_colors() {
    // This would need to be implemented to return the actual line color data
    // For now, return an empty vector
    static std::vector<vec4> empty;
    return empty;
}

const std::vector<vec3> &get_points() {
    // This would need to be implemented to return the actual point data
    // For now, return an empty vector
    static std::vector<vec3> empty;
    return empty;
}

const std::vector<vec4> &get_point_colors() {
    // This would need to be implemented to return the actual point color data
    // For now, return an empty vector
    static std::vector<vec4> empty;
    return empty;
}

} // namespace logger
} // namespace gaudi
