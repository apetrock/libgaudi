set(SOURCES
#    ../../src/GaudiGraphics/geometry_logger.cpp
    main.cpp
)

set(PROJECT_NAME wav2fft)
add_executable(${PROJECT_NAME} ${SOURCES})

#add_executable(bending main.cpp GuadiGraphics/geometry_logger.cpp)
setup_module_target(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} PRIVATE -g -ltbb)

#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(stretch PUBLIC OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(wav2fft PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})