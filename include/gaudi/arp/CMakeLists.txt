# ARP Module CMake Integration
# This file integrates the Ninja-built ARP module with the main CMake project

# Find Ninja executable
find_program(NINJA_EXECUTABLE ninja)

if(NOT NINJA_EXECUTABLE)
    message(FATAL_ERROR "Ninja build system not found. Please install ninja.")
endif()

# Set up module build directory
set(ARP_MODULE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/arp_module)
set(ARP_MODULE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

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
            ${NINJA_EXECUTABLE}
    WORKING_DIRECTORY ${ARP_MODULE_BUILD_DIR}
    COMMENT "Building ARP module with Ninja"
    VERBATIM
)

# Create imported target for ARP module
add_library(arp_module INTERFACE)

# Set up include directories for module consumers
target_include_directories(arp_module INTERFACE
    ${ARP_MODULE_BUILD_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/ext/eigen-3.4.0/eigen-3.4.0
)

# Add module-specific compile flags for consumers
target_compile_features(arp_module INTERFACE cxx_std_20)
target_compile_options(arp_module INTERFACE -fmodules-ts)

# Add dependency to ensure module is built before use
add_dependencies(arp_module arp_module_build)

# Export ARP module for use by other targets
set(ARP_MODULE_LIBRARY arp_module CACHE INTERNAL "ARP module library")

# Optional: Create a test target for ARP module
if(BUILD_TESTING)
    add_executable(arp_test
        ${ARP_MODULE_SOURCE_DIR}/arp_test.cpp
    )
    target_link_libraries(arp_test PRIVATE arp_module)
    target_compile_features(arp_test PRIVATE cxx_std_20)
    target_compile_options(arp_test PRIVATE -fmodules-ts)
    
    # Add test
    add_test(NAME arp_module_test COMMAND arp_test)
endif() 