#ifndef __LIBGAUDI_LOGGER_H__
#define __LIBGAUDI_LOGGER_H__

#include "common.h"
#include <vector>

namespace gaudi {
namespace logger {

// Interface for the geometry logger.
// The implementation is provided in a separate .cpp file.

void line(const vec3 &p0, const vec3 &p1, const vec4 &color);
void lines(const std::vector<vec3> &p0, const std::vector<vec3> &p1, const std::vector<vec4> &colors);
void box(const vec3 &cen, const vec3 &h, const vec4 &col);
void ext(const vec3 &mn, const vec3 &mx, const vec4 &col);
void field(const std::vector<vec3> &p, const std::vector<vec3> &dirs, double D = 0.1);
void frame(const mat3 &M, const vec3 &c, double C);
void point(const vec3 &p0, const vec4 &color);
void clear();

} // namespace logger
} // namespace gaudi

#endif // __LIBGAUDI_LOGGER_H__
