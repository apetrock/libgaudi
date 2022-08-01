add_executable(tanvec main.cpp)
target_compile_features(tanvec PRIVATE cxx_std_17)
target_link_libraries(tanvec PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})