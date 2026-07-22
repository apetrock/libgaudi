#ifndef __LIBGAUDI_GEOMETRY_LOGGER__
#define __LIBGAUDI_GEOMETRY_LOGGER__

#include "common.h"
#include <string>
#include <vector>

namespace gaudi {
namespace geometry_logger {

inline vec4 sdf4(double d) {
  const vec4 inside(0.0, 1.0, 0.0, 1.0);
  const vec4 outside(1.0, 0.0, 0.0, 1.0);
  if (d < 0.0)
    return std::abs(d) * inside;
  return std::abs(d) * outside;
}

constexpr real k_default_line_radius = 0.01;

// Geometry visualization functions
void line(const vec3 &p0, const vec3 &p1, const vec4 &color,
          real radius = k_default_line_radius);
void ext(const vec3 &mn, const vec3 &mx, const vec4 &col);
void frame(const mat3 &M, const vec3 &c, double C);
void point(const vec3 &p0, const vec4 &color);
void clear();

// Data access functions
const std::vector<vec3> &get_lines();
const std::vector<vec4> &get_line_colors();
const std::vector<vec3> &get_points();
const std::vector<vec4> &get_point_colors();

} // namespace geometry_logger
} // namespace gaudi

#endif
