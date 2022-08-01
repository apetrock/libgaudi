add_executable(bending main.cpp)
target_compile_features(bending PRIVATE cxx_std_17)
target_link_libraries(bending PRIVATE -g -ltbb)

#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(stretch PUBLIC OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(bending PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})