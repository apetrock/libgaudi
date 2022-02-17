if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(
  Gaudi
  VERSION 0.1
  LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

# Build Types
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}
    CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel tsan asan lsan msan ubsan"
    FORCE)

# ThreadSanitizer
set(CMAKE_C_FLAGS_TSAN
    "-fsanitize=thread -g -O1"
    CACHE STRING "Flags used by the C compiler during ThreadSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_TSAN
    "-fsanitize=thread -g -O1"
    CACHE STRING "Flags used by the C++ compiler during ThreadSanitizer builds."
    FORCE)

# AddressSanitize
set(CMAKE_C_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g3 -O1"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g3 -O1"
    CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
    FORCE)

if(CMAKE_BUILD_TYPE STREQUAL "ASAN")
  set(
    ASAN_OPTIONS
    symbolize=1
  )
  message("Set ASAN_OPTIONS to: ${ASAN_OPTIONS}")

  find_program(LLVM_SYMB_PATH llvm-symbolizer)
  set(
    ASAN_SYMBOLIZER_PATH
    ${LLVM_SYMB_PATH}
  )
  message("Set ASAN_SYMBOLIZER_PATH to: ${ASAN_SYMBOLIZER_PATH}")
endif()

# LeakSanitizer
set(CMAKE_C_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C compiler during LeakSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_LSAN
    "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
    CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
    FORCE)

# MemorySanitizer
set(CMAKE_C_FLAGS_MSAN
    "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
    CACHE STRING "Flags used by the C compiler during MemorySanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_MSAN
    "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
    CACHE STRING "Flags used by the C++ compiler during MemorySanitizer builds."
    FORCE)

# UndefinedBehaviour
set(CMAKE_C_FLAGS_UBSAN
    "-fsanitize=undefined"
    CACHE STRING "Flags used by the C compiler during UndefinedBehaviourSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_UBSAN
    "-fsanitize=undefined"
    CACHE STRING "Flags used by the C++ compiler during UndefinedBehaviourSanitizer builds."
    FORCE)



# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here

  # Let's ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Testing only available if this is the main app
  # Note this needs to be done in the main CMakeLists
  # since it calls enable_testing, which must be in the
  # main CMakeLists.
  include(CTest)

  # Docs only available if this is the main app
  find_package(Doxygen)
  if(Doxygen_FOUND)
    add_subdirectory(docs)
  else()
    message(STATUS "Doxygen not found, not building docs")
  endif()
endif()

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()

set(EIGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen-3.4.0")
set(SPECTRA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/spectra/include")

# Disable building extras we won't need (pure C++ project)
set(NANOGUI_BUILD_SHARED OFF CACHE BOOL " " FORCE)

set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)
set(NANOGUI_EIGEN_INCLUDE_DIR "${EIGEN_INCLUDE_DIR}" CACHE FILEPATH " " FORCE)

# Add the configurations from nanogui
add_subdirectory(ext/nanogui)


# For reliability of parallel build, make the NanoGUI targets dependencies
set_property(TARGET glfw glfw_objects nanogui PROPERTY FOLDER "dependencies")

include_directories("${NANOGUI_EXTRA_INCS}")
add_definitions("${NANOGUI_EXTRA_DEFS}")

#set(SPECTRA_EIGEN_INCLUDE_DIR "${EIGEN_INCLUDE_DIR}" CACHE FILEPATH " " FORCE)
#add_subdirectory(ext/spectra)

include_directories("${EIGEN_INCLUDE_DIR}" "${SPECTRA_INCLUDE_DIR}" )
add_subdirectory(src)
# The executable code is here

#add_subdirectory(projects/hello_manifold)
#add_subdirectory(projects/hello)
#add_subdirectory(projects/facemergetest)
add_subdirectory(projects/rxdiff)
add_subdirectory(projects/rmedge)
#add_subdirectory(projects/tanvec)
#add_subdirectory(projects/vortextest)
#add_subdirectory(projects/sdftest)
#add_subdirectory(projects/spectratest)
#add_subdirectory(projects/zithertest)
