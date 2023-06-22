set(SOURCES
#    ../../src/GaudiGraphics/geometry_logger.cpp
    main.cpp
)

set(PROJECT_NAME shell_constraints_align_test)
add_executable(${PROJECT_NAME} ${SOURCES})

#add_executable(bending main.cpp GuadiGraphics/geometry_logger.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} PRIVATE -g -ltbb)

#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(stretch PUBLIC OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(shell_constraints_align_test PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})