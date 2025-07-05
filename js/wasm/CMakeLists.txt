cmake_minimum_required(VERSION 3.22)

project(
  GaudiWasm
  VERSION 0.1
  LANGUAGES CXX)

# Only build WASM targets when using Emscripten
if (NOT DEFINED EMSCRIPTEN)
    message(FATAL_ERROR "This CMakeLists.txt is intended for Emscripten builds only. Use 'emcmake cmake ..' to configure.")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Let's nicely support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Include directories from main libgaudi project
set(LIBGAUDI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
set(EIGEN_INCLUDE_DIR "${LIBGAUDI_ROOT}/ext/eigen-3.4.0/eigen-3.4.0")
set(SPECTRA_INCLUDE_DIR "${LIBGAUDI_ROOT}/ext/spectra/include")
set(GAUDI_INCLUDE_DIR "${LIBGAUDI_ROOT}/include")
# Debug print directory variables
message(STATUS "LIBGAUDI_ROOT: ${LIBGAUDI_ROOT}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(
    "${EIGEN_INCLUDE_DIR}"
    "${SPECTRA_INCLUDE_DIR}"
    "${GAUDI_INCLUDE_DIR}"
)

# Common WASM logger files that can be shared across targets
set(COMMON_WASM_LOGGER_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/API/src/wasm_logger.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/API/src/wasm_logger_api.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/API/src/logger.cpp"
)

# Common Emscripten compiler flags for all WASM targets
set(COMMON_EMCC_COMPILE_FLAGS
    "-std=c++17"
    "-O2"
)

# Common Emscripten linker flags for all WASM targets
set(COMMON_EMCC_LINK_FLAGS
    "-sWASM=1"
    "-sSINGLE_FILE=1"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sENVIRONMENT=web"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sUSE_ES6_IMPORT_META=0"
    "--emit-symbol-map"
    "--bind"
    "-sINITIAL_MEMORY=268435456"
    "-sMAXIMUM_MEMORY=536870912"
    "-sSTACK_SIZE=134217728"
)

# Debug vs Release flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND COMMON_EMCC_COMPILE_FLAGS
        "-gsource-map"
        "-O1"
    )
    list(APPEND COMMON_EMCC_LINK_FLAGS
        "-gsource-map"
        "-O1"
        "-sASSERTIONS=1"
        "-sSAFE_HEAP=1"
        "-sSTACK_OVERFLOW_CHECK=2"
        "-sRETAIN_COMPILER_SETTINGS=1"
    )
else()
    list(APPEND COMMON_EMCC_COMPILE_FLAGS
        "-O2"
    )
    list(APPEND COMMON_EMCC_LINK_FLAGS
        "-O2"
        "-sASSERTIONS=0"
        "-sRETAIN_COMPILER_SETTINGS=1"
    )
endif()

# Function to create WASM targets with common settings
function(add_wasm_target TARGET_NAME)
    set(options ES6 USE_LOGGER)
    set(oneValueArgs EXPORT_NAME OUTPUT_DIR)
    set(multiValueArgs SOURCES LINK_FLAGS)
    cmake_parse_arguments(WASM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Combine sources with common component sources if requested
    set(ALL_SOURCES ${WASM_SOURCES})
    
    if(WASM_USE_LOGGER)
        list(APPEND ALL_SOURCES ${COMMON_WASM_LOGGER_SOURCES})
        message(STATUS "Adding logger sources to ${TARGET_NAME}")
    endif()
    
    # Create the executable
    add_executable(${TARGET_NAME} ${ALL_SOURCES})
    
    # Set output to .js
    set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".js")
    
    # Combine common compile flags with target-specific flags
    set(ALL_COMPILE_FLAGS ${COMMON_EMCC_COMPILE_FLAGS})
    
    # Combine common link flags with target-specific flags
    set(ALL_LINK_FLAGS ${COMMON_EMCC_LINK_FLAGS})
    
    if(WASM_EXPORT_NAME)
        list(APPEND ALL_LINK_FLAGS "-s EXPORT_NAME=\"${WASM_EXPORT_NAME}\"")
    endif()
    
    if(WASM_LINK_FLAGS)
        list(APPEND ALL_LINK_FLAGS ${WASM_LINK_FLAGS})
    endif()
    
    # Apply compile and link options separately
    target_compile_options(${TARGET_NAME} PRIVATE ${ALL_COMPILE_FLAGS})
    target_link_options(${TARGET_NAME} PRIVATE ${ALL_LINK_FLAGS})
    
    # Set output directory if specified
    if(WASM_OUTPUT_DIR)
        set_target_properties(${TARGET_NAME} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${WASM_OUTPUT_DIR}"
        )
    endif()
    
    # Print summary of what was included
    message(STATUS "Created WASM target: ${TARGET_NAME}")
    message(STATUS "  Sources: ${WASM_SOURCES}")
    message(STATUS "  Total sources: ${ALL_SOURCES}")
    if(WASM_EXPORT_NAME)
        message(STATUS "  Export name: ${WASM_EXPORT_NAME}")
    endif()
endfunction()

# Global demo directory setup
set(DEMO_WASM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src/wasm")
set(PUBLIC_WASM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../public/wasm")
file(MAKE_DIRECTORY "${DEMO_WASM_DIR}")
file(MAKE_DIRECTORY "${PUBLIC_WASM_DIR}")

# Add subdirectories for different WASM components
add_subdirectory(logger)
add_subdirectory(bar_demo)
add_subdirectory(foo_demo)
# add_subdirectory(rod_simulation)  # Directory doesn't exist
add_subdirectory(rod_constraints)
add_subdirectory(examples)
add_subdirectory(test_module)

# Create a combined target that builds everything (after subdirectories are processed)
add_custom_target(all_wasm_targets
    DEPENDS gaudi_logger_test rod_constraints_test hello_world bar_demo foo_demo test_module
    COMMENT "Building all WebAssembly targets"
)

# Global copy target to copy all WASM files to public directory
add_custom_target(copy_all_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:gaudi_logger_test>"
        "$<TARGET_FILE:rod_constraints_test>"
        "$<TARGET_FILE:hello_world>"
        "$<TARGET_FILE:bar_demo>"
        "$<TARGET_FILE:foo_demo>"
        "$<TARGET_FILE:test_module>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS all_wasm_targets
    COMMENT "Copying all WASM files to public directory"
)
