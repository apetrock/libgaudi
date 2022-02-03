add_executable(rxdiff main.cpp)
target_compile_features(rxdiff PRIVATE cxx_std_17)
target_link_libraries(rxdiff PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})