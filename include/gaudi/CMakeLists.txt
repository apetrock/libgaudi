# Gaudi Module Support for WASM
# This file provides module support for the WASM build system

# Include module support functions
include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/ModuleSupport.cmake)

# Set up ARP module for all builds including WASM
# Emscripten supports C++20 modules directly

# Set up module build directory
set(ARP_MODULE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/arp_module)
set(ARP_MODULE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/arp)

# Create module build directory
file(MAKE_DIRECTORY ${ARP_MODULE_BUILD_DIR})

# Copy Ninja build file to build directory
configure_file(
    ${ARP_MODULE_SOURCE_DIR}/build.ninja
    ${ARP_MODULE_BUILD_DIR}/build.ninja
    COPYONLY
)

# Build ARP module using Ninja
add_custom_target(arp_module_build ALL
    COMMAND ${CMAKE_COMMAND} -E chdir ${ARP_MODULE_BUILD_DIR}
            ninja
    WORKING_DIRECTORY ${ARP_MODULE_BUILD_DIR}
    COMMENT "Building ARP module with Ninja"
    VERBATIM
)

# Create ARP module as a static library that includes the implementation
add_library(arp_module STATIC)

# Set up include directories for module consumers
target_include_directories(arp_module PUBLIC
    ${ARP_MODULE_BUILD_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../ext/eigen-3.4.0/eigen-3.4.0
)

# Add module-specific compile flags for consumers
target_compile_features(arp_module PUBLIC cxx_std_20)
target_compile_options(arp_module PUBLIC -fmodules-ts)

# Add the implementation source file
target_sources(arp_module PRIVATE ${ARP_MODULE_SOURCE_DIR}/arp_impl.cpp)

# Add dependency to ensure module is built before use
add_dependencies(arp_module arp_module_build)

# Export ARP module for use by other targets
set(ARP_MODULE_LIBRARY arp_module CACHE INTERNAL "ARP module library") 