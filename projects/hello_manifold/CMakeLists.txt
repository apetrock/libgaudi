add_executable(manifold main.cpp)
target_compile_features(manifold PRIVATE cxx_std_17)
target_link_libraries(manifold PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})