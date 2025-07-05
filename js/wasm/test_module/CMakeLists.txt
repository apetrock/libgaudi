cmake_minimum_required(VERSION 3.16)
project(test_module)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Emscripten specific settings
set(CMAKE_EXECUTABLE_SUFFIX ".js")

# Test Module WASM
# Simple test module for validating WASM functionality

# Test Module
add_wasm_target(test_module
    SOURCES
        ./src/test_module.cpp
    EXPORT_NAME "TestModule"
    ES6
    OUTPUT_DIR "${PUBLIC_WASM_DIR}"
    LINK_FLAGS "-s EXPORTED_FUNCTIONS='[\"_test_function\",\"_get_test_message\",\"_add_numbers\"]' -s EXPORTED_RUNTIME_METHODS='[\"ccall\",\"cwrap\"]'"
)

# Set target properties for IDE organization
set_target_properties(test_module PROPERTIES FOLDER "wasm/test_module")

# Copy target for this module
add_custom_target(copy_test_module_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:test_module>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS test_module
    COMMENT "Copying test module WASM files to public directory"
) 