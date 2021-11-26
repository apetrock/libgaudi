add_executable(spectratest main.cpp)
target_compile_features(spectratest PRIVATE cxx_std_17)
target_link_libraries(spectratest PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})