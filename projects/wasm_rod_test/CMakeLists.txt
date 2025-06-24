cmake_minimum_required(VERSION 3.16)

project(RodSimulationWasm)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Emscripten-specific settings
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".js")
    
    # Emscripten compiler flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WASM=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s ALLOW_MEMORY_GROWTH=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s EXPORT_NAME='RodSimulationModule'")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s MODULARIZE=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s ENVIRONMENT=web")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --bind")
    
    # Enable exceptions and RTTI (needed for some C++ features)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s DISABLE_EXCEPTION_CATCHING=0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions")
    
    # Memory settings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s INITIAL_MEMORY=32MB")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s STACK_SIZE=1MB")
endif()

# Include directories - point to your existing libgaudi structure
set(LIBGAUDI_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)

include_directories(
    ${LIBGAUDI_ROOT}/include
    ${LIBGAUDI_ROOT}/ext/eigen
    ${LIBGAUDI_ROOT}/ext/nanogui/include
)

# Source files
set(SOURCES
    src/rod_wasm_wrapper.cpp
    src/wasm_bindings.cpp
    # We might need to add some of the libgaudi source files if they're not header-only
    ${LIBGAUDI_ROOT}/src/SIMPLE_PARSER.cpp
    ${LIBGAUDI_ROOT}/src/TIMER.cpp
)

# Create the WASM module
add_executable(rod_simulation ${SOURCES})

# Set output directory
set_target_properties(rod_simulation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/dist
)

# Copy the generated .wasm file to the JavaScript project
if(EMSCRIPTEN)
    add_custom_command(TARGET rod_simulation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/dist/rod_simulation.js
        ${CMAKE_CURRENT_SOURCE_DIR}/../js/rod-simulation-component/public/rod_simulation.js
        
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/dist/rod_simulation.wasm
        ${CMAKE_CURRENT_SOURCE_DIR}/../js/rod-simulation-component/public/rod_simulation.wasm
        
        COMMENT "Copying WASM files to React project"
    )
endif()
