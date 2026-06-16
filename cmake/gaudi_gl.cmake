# cmake/gaudi_gl.cmake
# Builds nanogui and provides gaudi::gl target + add_gaudi_gl_project() helper.
# Only include this when GAUDI_WITH_GL is ON.

if(TARGET gaudi_gl)
    return()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/gaudi_core.cmake")

get_filename_component(LIBGAUDI_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

# --- nanogui setup ---
set(SAVED_BINARY_DIR ${CMAKE_BINARY_DIR})
set(CMAKE_BINARY_DIR "${LIBGAUDI_ROOT}/ext/nanogui/build")

set(NANOGUI_BUILD_SHARED  OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)
set(NANOGUI_EIGEN_INCLUDE_DIR "${LIBGAUDI_ROOT}/ext/eigen" CACHE FILEPATH " " FORCE)

add_subdirectory("${LIBGAUDI_ROOT}/ext/nanogui" "${LIBGAUDI_ROOT}/ext/nanogui/build")
set(CMAKE_BINARY_DIR ${SAVED_BINARY_DIR})

set_property(TARGET glfw glfw_objects nanogui PROPERTY FOLDER "dependencies")

# --- gaudi::gl INTERFACE target ---
add_library(gaudi_gl INTERFACE)
add_library(gaudi::gl ALIAS gaudi_gl)

target_link_libraries(gaudi_gl INTERFACE
    gaudi::core
    nanogui
    ${NANOGUI_EXTRA_LIBS}
)

target_include_directories(gaudi_gl INTERFACE
    "${NANOGUI_EXTRA_INCS}"
    "${LIBGAUDI_ROOT}/ext/nanogui/include"
)

target_compile_definitions(gaudi_gl INTERFACE
    ${NANOGUI_EXTRA_DEFS}
)

# --- libGaudi legacy static library (SIMPLE_PARSER, TIMER, GL geometry_logger) ---
add_library(gaudi_legacy STATIC
    "${LIBGAUDI_ROOT}/src/SIMPLE_PARSER.cpp"
    "${LIBGAUDI_ROOT}/src/TIMER.cpp"
    "${LIBGAUDI_ROOT}/src/geometry_logger.cpp"
    "${LIBGAUDI_ROOT}/src/GaudiGraphics/geometry_logger.cpp"
    "${LIBGAUDI_ROOT}/src/terminal_logger_impl.cpp"
)
target_link_libraries(gaudi_legacy PUBLIC gaudi::gl)
target_include_directories(gaudi_legacy PUBLIC "${LIBGAUDI_ROOT}/include")
target_compile_features(gaudi_legacy PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(gaudi_legacy PRIVATE /FS)
endif()
gaudi_msvc_disable_parallel_cl_writers(gaudi_legacy)
# --- Convenience function for GL projects (mirrors add_wasm_target pattern) ---
# Usage:
#   add_gaudi_gl_project(my_project SOURCES main.cpp)
#   # Optional: CHOLMOD, TBB, OPENMP flags handled automatically
function(add_gaudi_gl_project TARGET_NAME)
    cmake_parse_arguments(GL "" "" "SOURCES" ${ARGN})

    if(NOT GL_SOURCES)
        message(FATAL_ERROR "add_gaudi_gl_project: SOURCES is required")
    endif()

    add_executable(${TARGET_NAME} ${GL_SOURCES})
    target_link_libraries(${TARGET_NAME} PRIVATE gaudi_legacy)
    target_compile_features(${TARGET_NAME} PRIVATE cxx_std_20)
    if(MSVC)
        target_compile_options(${TARGET_NAME} PRIVATE /bigobj /FS)
    endif()
    gaudi_msvc_disable_parallel_cl_writers(${TARGET_NAME})

    if(GUADI_WITH_CHOLMOD)
        target_link_libraries(${TARGET_NAME} PRIVATE gaudi::cholmod)
    endif()
endfunction()

message(STATUS "gaudi::gl configured (nanogui + OpenGL + add_gaudi_gl_project)")
