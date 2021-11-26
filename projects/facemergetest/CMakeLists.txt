add_executable(facemerge main.cpp)
target_compile_features(facemerge PRIVATE cxx_std_17)
target_link_libraries(facemerge PRIVATE libGaudi fmt::fmt nanogui ${NANOGUI_EXTRA_LIBS})