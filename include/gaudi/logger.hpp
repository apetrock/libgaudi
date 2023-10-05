

#ifndef __LIBGAUDI_LOGGER__
#define __LIBGAUDI_LOGGER__
#include "GaudiGraphics/geometry_logger.h"
#include "common.h"

// loggers, visual debugger... perhaps build better screen logger
// wrap gaudi graphics logger so that it can be abstracted away
namespace gaudi {
namespace logger {
static void line(const vec3 &p0, const vec3 &p1, const vec4 &color) {
  gg::geometry_logger::line(p0, p1, color);
};

} // namespace logger
} // namespace gaudi
#endif