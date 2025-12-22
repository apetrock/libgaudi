#ifndef __GAUDI_GEOMETRY_LOGGER_DISPATCH__
#define __GAUDI_GEOMETRY_LOGGER_DISPATCH__

// Automatically dispatch to the correct geometry logger based on build target
#ifdef __EMSCRIPTEN__
  // WASM build - use the lightweight WASM logger
  #include "gaudi/geometry_logger.hpp"
#else
  // Desktop build - use the OpenGL-based logger
  #include "GaudiGraphics/geometry_logger.h"
#endif

#endif // __GAUDI_GEOMETRY_LOGGER_DISPATCH__

