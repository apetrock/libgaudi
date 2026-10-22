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

add_executable(slab_prism_smoke slab_prism_smoke.cpp)
target_include_directories(slab_prism_smoke PRIVATE ${CMAKE_SOURCE_DIR}/ext/liblombardi/include)
target_link_libraries(slab_prism_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(slab_prism_smoke)
add_test(NAME slab_prism_smoke COMMAND slab_prism_smoke)

add_executable(pierce_audit pierce_audit_main.cpp)
target_link_libraries(pierce_audit PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(pierce_audit)

add_executable(monge_fairing_smoke monge_fairing_smoke.cpp)
target_link_libraries(monge_fairing_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(monge_fairing_smoke)
add_test(NAME monge_fairing_smoke COMMAND monge_fairing_smoke)

add_executable(frame_spin_probe frame_spin_probe.cpp)
target_include_directories(frame_spin_probe PRIVATE ${CMAKE_SOURCE_DIR}/ext/liblombardi/include)
target_link_libraries(frame_spin_probe PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(frame_spin_probe)

add_executable(principal_curvature_tests principal_curvature_tests.cpp)
target_link_libraries(principal_curvature_tests PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(principal_curvature_tests)
add_test(NAME principal_curvature_tests COMMAND principal_curvature_tests)

add_executable(weingarten_tests weingarten_tests.cpp)
target_link_libraries(weingarten_tests PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(weingarten_tests)
add_test(NAME weingarten_tests COMMAND weingarten_tests)

add_executable(shell_qc_tests shell_qc_tests.cpp)
target_link_libraries(shell_qc_tests PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(shell_qc_tests)
add_test(NAME shell_qc_tests COMMAND shell_qc_tests)

add_executable(principal_curvature_probe principal_curvature_probe.cpp)
target_include_directories(principal_curvature_probe PRIVATE ${CMAKE_SOURCE_DIR}/ext/liblombardi/include)
target_link_libraries(principal_curvature_probe PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(principal_curvature_probe)

add_executable(anisotropic_cotan_smoke anisotropic_cotan_smoke.cpp)
target_link_libraries(anisotropic_cotan_smoke PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(anisotropic_cotan_smoke)
add_test(NAME anisotropic_cotan_smoke COMMAND anisotropic_cotan_smoke)

# liblombardi tests need access to ext/liblombardi/include
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../ext/liblombardi/include)

message(STATUS "Configured gaudi_tests (headless, no GL)")
