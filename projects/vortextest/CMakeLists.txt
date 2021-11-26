add_executable(vortextest main.cpp)
target_compile_features(vortextest PRIVATE cxx_std_17)
target_link_libraries(vortextest PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})