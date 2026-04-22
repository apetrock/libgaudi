# Canonical viewer sources live in `spectral_modes_main.cpp`. Some Windows setups
# keep an exclusive lock on the in-tree `main.cpp` (EBUSY), so we compile a
# generated copy under the build tree named `main.cpp` for a normal-looking TU.
set(_spectral_modes_main_cpp "${CMAKE_CURRENT_BINARY_DIR}/main.cpp")
add_custom_command(
  OUTPUT ${_spectral_modes_main_cpp}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/spectral_modes_main.cpp"
    ${_spectral_modes_main_cpp}
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/spectral_modes_main.cpp"
  COMMENT "spectral_modes_test: sync build-tree main.cpp from spectral_modes_main.cpp"
  VERBATIM)
add_gaudi_gl_project(spectral_modes_test SOURCES ${_spectral_modes_main_cpp})
