# tests/CMakeLists.txt
# Headless test executable -- no GL, no nanogui, just Eigen + gaudi headers.
# Assets are resolved at compile time via GAUDI_REPO_ROOT (from gaudi::core).

add_executable(gaudi_tests gaudi_tests.cpp)
target_link_libraries(gaudi_tests PRIVATE gaudi::core gaudi::headless)
gaudi_msvc_disable_parallel_cl_writers(gaudi_tests)

add_test(NAME gaudi_tests COMMAND gaudi_tests)

message(STATUS "Configured gaudi_tests (headless, no GL)")
