add_executable(rmedge main.cpp)
target_compile_features(rmedge PRIVATE cxx_std_17)
target_link_libraries(rmedge PRIVATE -fsanitize=address)
target_link_libraries(rmedge PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})