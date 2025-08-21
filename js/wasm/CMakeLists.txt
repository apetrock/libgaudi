cmake_minimum_required(VERSION 3.22)

project(
  GaudiWasm
  VERSION 0.1
  LANGUAGES CXX)

# Only build WASM targets when using Emscripten
if (NOT DEFINED EMSCRIPTEN)
    message(FATAL_ERROR "This CMakeLists.txt is intended for Emscripten builds only. Use 'emcmake cmake ..' to configure.")
endif()

# Include common configuration
include(cmake-templates/common.cmake)

# Global demo directory setup
set(DEMO_WASM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src/wasm")
set(PUBLIC_WASM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../public/wasm")
file(MAKE_DIRECTORY "${DEMO_WASM_DIR}")
file(MAKE_DIRECTORY "${PUBLIC_WASM_DIR}")

# Add subdirectories for different WASM components
add_subdirectory(modules/morton_tree)
add_subdirectory(modules/rod_constraints)
add_subdirectory(modules/strand_test)

# Create a combined target that builds everything
add_custom_target(all_wasm_targets
    DEPENDS morton_tree rod_constraints strand_test
    COMMENT "Building all WebAssembly targets"
)

# Global copy target to copy all WASM files to public directory
add_custom_target(copy_all_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:morton_tree>"
        "$<TARGET_FILE:rod_constraints>"
        "$<TARGET_FILE:strand_test>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS all_wasm_targets
    COMMENT "Copying all WASM files to public directory"
) 