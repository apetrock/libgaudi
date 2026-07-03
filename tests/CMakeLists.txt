# tests/CMakeLists.txt
# Headless test executable -- no GL, no nanogui, just Eigen + gaudi headers.
# Assets are resolved at compile time via GAUDI_REPO_ROOT (from gaudi::core).

add_executable(gaudi_tests gaudi_tests.cpp liblombardi_tests.cpp)
target_include_directories(gaudi_tests PRIVATE ${CMAKE_SOURCE_DIR}/ext/liblombardi/include)
target_link_libraries(gaudi_tests PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(gaudi_tests)

add_test(NAME gaudi_tests COMMAND gaudi_tests)

add_executable(vector_dirichlet_smoke vector_dirichlet_smoke.cpp)
target_link_libraries(vector_dirichlet_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(vector_dirichlet_smoke)
add_test(NAME vector_dirichlet_smoke COMMAND vector_dirichlet_smoke)

add_executable(spectral_smoke spectral_smoke.cpp)
target_link_libraries(spectral_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(spectral_smoke)
add_test(NAME spectral_smoke COMMAND spectral_smoke)

add_executable(anisotropic_cotan_smoke anisotropic_cotan_smoke.cpp)
target_link_libraries(anisotropic_cotan_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(anisotropic_cotan_smoke)
add_test(NAME anisotropic_cotan_smoke COMMAND anisotropic_cotan_smoke)

# liblombardi tests need access to ext/liblombardi/include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../ext/liblombardi/include)

message(STATUS "Configured gaudi_tests (headless, no GL)")
