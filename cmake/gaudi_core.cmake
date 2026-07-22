# cmake/gaudi_core.cmake
# INTERFACE library providing gaudi headers + Eigen + Spectra + C++20
# This is the foundation that every build mode (headless, GL, WASM) uses.

# VS 2019+ MSBuild UseMultiToolTask can run multiple cl.exe against one vc*.pdb; C1041 may
# persist even with /FS. Disabling MTT is the reliable fix (see MSVC developer community).
if(NOT COMMAND gaudi_msvc_disable_parallel_cl_writers)
    function(gaudi_msvc_disable_parallel_cl_writers TARGET)
        if(MSVC AND CMAKE_GENERATOR MATCHES "Visual Studio")
            set_target_properties(${TARGET} PROPERTIES
                VS_GLOBAL_UseMultiToolTask "false"
            )
        endif()
    endfunction()
endif()

if(TARGET gaudi_core)
    return()
endif()

# Resolve LIBGAUDI_ROOT from this file's location
get_filename_component(LIBGAUDI_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

add_library(gaudi_core INTERFACE)
add_library(gaudi::core ALIAS gaudi_core)

target_include_directories(gaudi_core INTERFACE
    "${LIBGAUDI_ROOT}/include"
    "${LIBGAUDI_ROOT}/ext/eigen"
    "${LIBGAUDI_ROOT}/ext/spectra/include"
    "${LIBGAUDI_ROOT}/scripts/sympy/generated"
)

target_compile_features(gaudi_core INTERFACE cxx_std_20)

target_compile_definitions(gaudi_core INTERFACE
    GAUDI_REPO_ROOT="${LIBGAUDI_ROOT}"
)

# Propagate OpenMP so Eigen parallelizes matmul / reductions / some sparse ops.
if(GAUDI_HAS_OPENMP)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(gaudi_core INTERFACE OpenMP::OpenMP_CXX)
        message(STATUS "gaudi::core linking OpenMP::OpenMP_CXX")
    else()
        target_compile_options(gaudi_core INTERFACE ${OpenMP_CXX_FLAGS})
        target_link_options(gaudi_core INTERFACE ${OpenMP_CXX_FLAGS})
        if(OpenMP_omp_LIBRARY)
            target_link_libraries(gaudi_core INTERFACE ${OpenMP_omp_LIBRARY})
        endif()
        message(STATUS "gaudi::core linking OpenMP via legacy flags")
    endif()
endif()

message(STATUS "gaudi::core configured (headers + Eigen + Spectra, C++20)")
