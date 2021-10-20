add_executable(hello main.cpp)
target_compile_features(hello PRIVATE cxx_std_17)
target_link_libraries(hello PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})