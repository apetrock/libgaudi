if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(
  Gaudi
  VERSION 0.1
  LANGUAGES CXX)

cmake_minimum_required(VERSION 3.10)

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


# Formatting library
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 5.3.0)
FetchContent_MakeAvailable(fmtlib)
# Adds fmt::fmt


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
add_subdirectory(projects/hello)
add_subdirectory(projects/facemergetest)
add_subdirectory(projects/vortextest)
add_subdirectory(projects/sdftest)
add_subdirectory(projects/spectratest)
add_subdirectory(projects/zithertest)
