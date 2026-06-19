# cmake/gaudi_headless.cmake
# STATIC library providing headless logger backends for native builds.
# Links gaudi::core. No GL/nanogui dependency.
#
# Logger backend pattern: same headers (geometry_logger.hpp, logger.hpp),
# different .cpp implementations per environment:
#   - Native headless: src/geometry_logger.cpp, src/terminal_logger_impl.cpp
#   - Native Vermeer: src/vermeer/geometry_logger_impl.cpp
#   - WASM: js/wasm/API/src/geometry_logger_impl.cpp, terminal_logger_impl.cpp
#   - GL: legacy projects may still call gg::geometry_logger directly

if(TARGET gaudi_headless)
    return()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/gaudi_core.cmake")

get_filename_component(LIBGAUDI_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

add_library(gaudi_headless STATIC
    "${LIBGAUDI_ROOT}/src/geometry_logger.cpp"
    "${LIBGAUDI_ROOT}/src/terminal_logger_impl.cpp"
)
add_library(gaudi::headless ALIAS gaudi_headless)

target_link_libraries(gaudi_headless PUBLIC gaudi::core)
gaudi_msvc_disable_parallel_cl_writers(gaudi_headless)

message(STATUS "gaudi::headless configured (headless logger backends)")
