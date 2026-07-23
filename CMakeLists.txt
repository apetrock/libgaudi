if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(
  Gaudi
  VERSION 0.1
  LANGUAGES CXX)

# trace
cmake_minimum_required(VERSION 3.16)

include(FetchContent)
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()

# Build Types
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}
    CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel tsan asan lsan msan ubsan"
    FORCE)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(GAUDI_WITH_GL "Build OpenGL viewer projects (nanogui)" OFF)
option(GAUDI_WITH_VERMEER "Build WebGPU viewer projects (libvermeer)" OFF)
option(GAUDI_USE_EXPLICIT_INSTANTIATIONS "Use explicit template instantiations for faster compilation" ON)
option(WITH_OPENMP "Use OpenMP" ON)
option(WITH_CUDA "Build with CUDA support" OFF)
option(GUADI_WITH_CHOLMOD "Use Cholmod" OFF)

# ---------------------------------------------------------------------------
# Sanitizer flags (compiler-aware: MSVC, Clang, GCC)
# ---------------------------------------------------------------------------
include(cmake/gaudi_sanitizers.cmake)
include(cmake/gaudi_msvc_z7.cmake)

# ---------------------------------------------------------------------------
# Main project settings (only when this is the top-level project)
# ---------------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_ENABLE_EXPORTS TRUE)
  set(CMAKE_CXX_EXTENSIONS ON)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # MSVC C1041: multiple cl.exe writers to the same vc*.pdb (/Zi, parallel MSBuild, /MP) need /FS.
  if(MSVC)
    add_compile_options(/FS)
  endif()

  include(CTest)

  find_package(Doxygen)
  if(Doxygen_FOUND)
    add_subdirectory(docs)
  else()
    message(STATUS "Doxygen not found, not building docs")
  endif()
endif()

# ---------------------------------------------------------------------------
# OpenMP (optional)
# ---------------------------------------------------------------------------
if(WITH_OPENMP)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _GAUDI_BT_CHECK)
    if(_GAUDI_BT_CHECK MATCHES "ASAN|TSAN|MSAN|LSAN")
        message(WARNING "OpenMP is enabled but a sanitizer build is active. "
            "This combination can produce false positives or link errors on some platforms. "
            "Consider -DWITH_OPENMP=OFF for sanitizer builds.")
    endif()
    unset(_GAUDI_BT_CHECK)

    # AppleClang does not ship libomp; Homebrew's keg-only libomp needs hints
    # before find_package(OpenMP), otherwise discovery silently fails.
    if(APPLE)
        set(_gaudi_need_libomp_hints FALSE)
        if(NOT OpenMP_omp_LIBRARY OR OpenMP_omp_LIBRARY MATCHES "NOTFOUND")
            set(_gaudi_need_libomp_hints TRUE)
        endif()
        if(NOT OpenMP_CXX_FLAGS OR OpenMP_CXX_FLAGS MATCHES "NOTFOUND")
            set(_gaudi_need_libomp_hints TRUE)
        endif()
        if(_gaudi_need_libomp_hints)
            find_program(GAUDI_BREW_EXECUTABLE brew)
            if(GAUDI_BREW_EXECUTABLE)
                execute_process(
                    COMMAND "${GAUDI_BREW_EXECUTABLE}" --prefix libomp
                    OUTPUT_VARIABLE GAUDI_LIBOMP_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                    RESULT_VARIABLE GAUDI_LIBOMP_BREW_RC)
                if(GAUDI_LIBOMP_BREW_RC EQUAL 0 AND EXISTS "${GAUDI_LIBOMP_PREFIX}/include/omp.h")
                    message(STATUS "Using Homebrew libomp at ${GAUDI_LIBOMP_PREFIX}")
                    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${GAUDI_LIBOMP_PREFIX}/include"
                        CACHE STRING "OpenMP C flags" FORCE)
                    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${GAUDI_LIBOMP_PREFIX}/include"
                        CACHE STRING "OpenMP CXX flags" FORCE)
                    set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "OpenMP C libs" FORCE)
                    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP CXX libs" FORCE)
                    set(OpenMP_omp_LIBRARY "${GAUDI_LIBOMP_PREFIX}/lib/libomp.dylib"
                        CACHE FILEPATH "libomp" FORCE)
                endif()
            endif()
        endif()
        unset(_gaudi_need_libomp_hints)
    endif()

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found (CXX version ${OpenMP_CXX_VERSION})")
        # Prefer modern imported target; also keep legacy flags for older paths.
        set(GAUDI_HAS_OPENMP ON)
    elseif(OPENMP_FOUND)
        message(STATUS "OpenMP found (legacy OPENMP_FOUND)")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        set(GAUDI_HAS_OPENMP ON)
    else()
        message(WARNING "WITH_OPENMP=ON but OpenMP was not found. "
            "On macOS install libomp (brew install libomp) and reconfigure.")
        set(GAUDI_HAS_OPENMP OFF)
    endif()
else()
    set(GAUDI_HAS_OPENMP OFF)
endif()

# ---------------------------------------------------------------------------
# CUDA (optional)
# ---------------------------------------------------------------------------
if(WITH_CUDA)
    find_program(NVCC_PATH nvcc)
    if(NOT NVCC_PATH)
        message(STATUS "nvcc not found, disabling CUDA")
        set(WITH_CUDA OFF)
    else()
        set(CMAKE_CUDA_COMPILER ${NVCC_PATH})
        message(STATUS "Using CUDA compiler: ${CMAKE_CUDA_COMPILER}")
        if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
          set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
        endif()
        message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")
        enable_language(CUDA)
    endif()
endif()

# ---------------------------------------------------------------------------
# Cholmod / SuiteSparse (optional)
# ---------------------------------------------------------------------------
if(GUADI_WITH_CHOLMOD)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  message(STATUS "Using fastest solver: Cholmod")
  set(WITH_FORTRAN OFF CACHE BOOL "")
  set(WITH_PARTITION OFF CACHE BOOL "")
  set(WITH_DEMOS OFF CACHE BOOL "")
  set(BUILD_CXSPARSE OFF CACHE BOOL "")
  set(WITH_CHOLMOD ON CACHE BOOL "")
  set(WITH_SUPERNODAL ON CACHE BOOL "")

  set(SAVED_BINARY_DIR ${CMAKE_BINARY_DIR})
  set(SUITESPARSE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/SuiteSparse/build")
  set(CMAKE_BINARY_DIR ${SUITESPARSE_BINARY_DIR})
  set(CMAKE_BINARY_DIR ${SAVED_BINARY_DIR})

  add_subdirectory("ext/SuiteSparse/")

  file(GLOB_RECURSE HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/ext/SuiteSparse/*.h")
  foreach(FILE_PATH ${HEADER_LIST})
      get_filename_component(DIR_PATH ${FILE_PATH} DIRECTORY)
      list(APPEND SUITESPARSE_INCLUDE_DIRS ${DIR_PATH})
  endforeach()
  list(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)

  add_library(gaudi_cholmod INTERFACE)
  add_library(gaudi::cholmod ALIAS gaudi_cholmod)
  target_include_directories(gaudi_cholmod INTERFACE ${SUITESPARSE_INCLUDE_DIRS})
  target_link_libraries(gaudi_cholmod INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  message(STATUS "GUADI_WITH_CHOLMOD=OFF, skipping SuiteSparse/CHOLMOD")
endif()

# ---------------------------------------------------------------------------
# Core module (headers + Eigen + Spectra, always available)
# ---------------------------------------------------------------------------
include(cmake/gaudi_core.cmake)

# ---------------------------------------------------------------------------
# Explicit instantiations library (optional, for faster compilation)
# ---------------------------------------------------------------------------
if(GAUDI_USE_EXPLICIT_INSTANTIATIONS)
    add_library(gaudi_library STATIC
        "src/gaudi_arp_instantiations.cpp"
    )
    target_link_libraries(gaudi_library PUBLIC gaudi::core)
    target_compile_definitions(gaudi_library PUBLIC
        GAUDI_USE_EXPLICIT_INSTANTIATIONS=1
    )
    gaudi_msvc_disable_parallel_cl_writers(gaudi_library)
    message(STATUS "Created gaudi_library with explicit instantiations")
else()
    message(STATUS "Explicit instantiations disabled")
endif()

# ---------------------------------------------------------------------------
# GL module + viewer projects (optional)
# ---------------------------------------------------------------------------
if(GAUDI_WITH_GL)
    include(cmake/gaudi_gl.cmake)

    add_subdirectory(projects/morton_test)
    add_subdirectory(projects/blank_template)
    add_subdirectory(projects/aabb_test)
    add_subdirectory(projects/gravitys_rainbow)
    add_subdirectory(projects/rod_guided_deformation)
    add_subdirectory(projects/rod_constraints_test)
    # Uncomment as needed:
    #add_subdirectory(projects/dynamic_shell_test)
    #add_subdirectory(projects/rod_constraints_coupled_test)
    #add_subdirectory(projects/shell_constraints_test)
    #add_subdirectory(projects/shell_constraints_align_test)
    #add_subdirectory(projects/mighty_morphin_test)
    add_subdirectory(projects/growth_study)
    add_subdirectory(projects/arnoldi_test)
    add_subdirectory(projects/spectral_modes_test)
    add_subdirectory(projects/frame_curvature_demo)
    add_subdirectory(projects/darboux_cyclide_demo)
    add_subdirectory(projects/darboux_cyclide_medial_demo)
    #add_subdirectory(projects/reaction_diffusion_test)
    #add_subdirectory(projects/vortex_test)
    #add_subdirectory(projects/tangent_point_test)
    #add_subdirectory(projects/shell_walk_test)
    #add_subdirectory(projects/heat_geodesics_test)
    #add_subdirectory(projects/fast_summation_test)
    #add_subdirectory(projects/wav2fft)
else()
    message(STATUS "GAUDI_WITH_GL=OFF, skipping OpenGL projects")
endif()

# ---------------------------------------------------------------------------
# Vermeer / WebGPU module + viewer projects (optional)
# ---------------------------------------------------------------------------
if(GAUDI_WITH_VERMEER)
    include(cmake/gaudi_vermeer.cmake)
    add_subdirectory(projects/vermeer_fast_summation_demo)
    add_subdirectory(projects/vermeer_rod_guided_demo)
    add_subdirectory(projects/vermeer_darboux_cyclide_medial_demo)
    add_subdirectory(projects/vermeer_rod_constraints_demo)
    add_subdirectory(projects/vermeer_torque_loop_demo)
    add_subdirectory(projects/vermeer_dipole_tunneling_demo)
else()
    message(STATUS "GAUDI_WITH_VERMEER=OFF, skipping WebGPU projects")
endif()

# ---------------------------------------------------------------------------
# Tests (headless, no GL required)
# ---------------------------------------------------------------------------
if(BUILD_TESTING)
    include(cmake/gaudi_headless.cmake)
    add_subdirectory(tests)
endif()
