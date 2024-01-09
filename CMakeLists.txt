if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(
  Gaudi
  VERSION 0.1
  LANGUAGES CXX)

cmake_minimum_required(VERSION 3.22)

# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()




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
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O0"
    CACHE STRING "Flags used by the C compiler during AddressSanitizer builds."
    FORCE)
set(CMAKE_CXX_FLAGS_ASAN
    "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O0"
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
  set (CMAKE_ENABLE_EXPORTS TRUE)
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


option(WITH_OPENMP "Use OpenMP" ON)
if(WITH_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(ENV{OMP_NUM_THREADS} 8)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

set(EIGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen-3.4.0")
set(SPECTRA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/spectra/include")

option(WITH_CUDA "Build with CUDA support" ON)
if(WITH_CUDA)
    
    #find_program(NVCC_PATH nvcc /usr/local/cuda/bin)
    find_program(NVCC_PATH nvcc)
    if(NOT NVCC_PATH)
        message(STATUS "nvcc not found in /usr/bin, using default path")
    else()
        set(CMAKE_CUDA_COMPILER ${NVCC_PATH})
        message(STATUS "Using CUDA compiler: ${CMAKE_CUDA_COMPILER}")
    endif()
    
    if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
      set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
    endif()

    message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")
    
    message(STATUS "Using CUDA")  
    
    enable_language(CUDA)
endif()

option(GUADI_WITH_CHOLMOD   "Use Cholmod" ON)
if(GUADI_WITH_CHOLMOD)

  #find_package(CHOLMOD)
  message(STATUS "Using fastest solver: Cholmod")
  set(WITH_FORTRAN OFF CACHE BOOL "")
  set(WITH_PARTITION OFF CACHE BOOL "")
  set(WITH_DEMOS OFF CACHE BOOL "")
  set(BUILD_CXSPARSE OFF CACHE BOOL "")
  set(WITH_CHOLMOD ON CACHE BOOL "")
  set(WITH_SUPERNODAL ON CACHE BOOL "")
  
  #add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/SuiteSparse")
  # push binary dir
  set(SAVED_BINARY_DIR ${CMAKE_BINARY_DIR})
  set(SUITESPARSE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/SuiteSparse/build")
  set(CMAKE_BINARY_DIR ${SUITESPARSE_BINARY_DIR})
  # pop it back
  set(CMAKE_BINARY_DIR ${SAVED_BINARY_DIR})
  
  add_subdirectory("ext/SuiteSparse/")
  
  #gather all the headers
  file(GLOB_RECURSE HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/ext/SuiteSparse/*.h")

  foreach(FILE_PATH ${HEADER_LIST})
      get_filename_component(DIR_PATH ${FILE_PATH} DIRECTORY)
      list(APPEND SUITESPARSE_INCLUDE_DIRS ${DIR_PATH})
  endforeach()

  list(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)
  #message(STATUS "cholmod include dirs: ${SUITESPARSE_INCLUDE_DIRS}")
  include_directories(${SUITESPARSE_INCLUDE_DIRS})
endif()


# Disable building extras we won't need (pure C++ project)
set(SAVED_BINARY_DIR ${CMAKE_BINARY_DIR})
set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/nanogui/build")
  
set(NANOGUI_BUILD_SHARED OFF CACHE BOOL " " FORCE)

set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
set(NANOGUI_BUILD_PYTHON  OFF CACHE BOOL " " FORCE)
set(NANOGUI_INSTALL       OFF CACHE BOOL " " FORCE)
set(NANOGUI_EIGEN_INCLUDE_DIR "${EIGEN_INCLUDE_DIR}" CACHE FILEPATH " " FORCE)
# Add the configurations from nanogui
add_subdirectory(ext/nanogui)
set(CMAKE_BINARY_DIR ${SAVED_BINARY_DIR})


# For reliability of parallel build, make the NanoGUI targets dependencies
set_property(TARGET glfw glfw_objects nanogui PROPERTY FOLDER "dependencies")

include_directories("${NANOGUI_EXTRA_INCS}")
add_definitions("${NANOGUI_EXTRA_DEFS}")

#set(SPECTRA_EIGEN_INCLUDE_DIR "${EIGEN_INCLUDE_DIR}" CACHE FILEPATH " " FORCE)
#add_subdirectory(ext/spectra)

include_directories("${EIGEN_INCLUDE_DIR}" "${SPECTRA_INCLUDE_DIR}" "ext/AudioFile")
add_subdirectory(src)



add_subdirectory(projects/blank_template)
#add_subdirectory(projects/dynamic_shell_test)
add_subdirectory(projects/rod_constraints_test)
add_subdirectory(projects/rod_constraints_coupled_test)
#add_subdirectory(projects/shell_constraints_test)
#add_subdirectory(projects/shell_constraints_align_test)

#add_subdirectory(projects/rod_constraints_coupled_test)

#add_subdirectory(projects/mighty_morphin_test)
add_subdirectory(projects/rod_guided_deformation)
add_subdirectory(projects/growth_study)

#add_subdirectory(projects/reaction_diffusion_test)
#add_subdirectory(projects/vortex_test)
#add_subdirectory(projects/tangent_point_test)

#add_subdirectory(projects/shell_walk_test)

#add_subdirectory(projects/heat_geodesics_test)
#add_subdirectory(projects/fast_summation_test)
#add_subdirectory(projects/wav2fft)


