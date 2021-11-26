add_executable(zithertest main.cpp)
target_compile_features(zithertest PRIVATE cxx_std_17)
target_link_libraries(zithertest PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})