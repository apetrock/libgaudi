# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library.

# Optionally glob, but only for CMake 3.12 or later:
# file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${ModernCMakeExample_SOURCE_DIR}/include/modern/*.hpp")

set(HEADER_LIST 
    "${Gaudi_SOURCE_DIR}/include/SIMPLE_PARSER.h"
    "${Gaudi_SOURCE_DIR}/include/TIMER.h"
    "${Gaudi_SOURCE_DIR}/include/GaudiMath/typedefs.hpp"
    "${Gaudi_SOURCE_DIR}/include/GaudiGraphics/buffers.hpp"
    "${Gaudi_SOURCE_DIR}/include/GaudiGraphics/geometry_logger.h"
)

set(SOURCES
    "${Gaudi_SOURCE_DIR}/src/SIMPLE_PARSER.cpp"
    "${Gaudi_SOURCE_DIR}/src/TIMER.cpp"
    "${Gaudi_SOURCE_DIR}/src/GaudiGraphics/geometry_logger.cpp"
)

# Make an automatic library - will be static or dynamic based on user setting
add_library(libGaudi ${SOURCES}
                     ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(libGaudi PUBLIC "${Gaudi_SOURCE_DIR}/include")

# This depends on (header only) boost
#target_link_libraries(libGaudi PRIVATE Boost::boost)

# All users of this library will need at least C++11
target_compile_features(libGaudi PUBLIC cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
