add_executable(sdftest main.cpp)
target_compile_features(sdftest PRIVATE cxx_std_17)
target_link_libraries(sdftest PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})