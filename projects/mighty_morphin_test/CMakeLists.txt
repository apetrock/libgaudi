set(SOURCES
#    ../../src/GaudiGraphics/geometry_logger.cpp
    main.cpp
)

set(PROJECT_NAME mighty_morphin_test)
add_executable(${PROJECT_NAME} ${SOURCES})

#add_executable(bending main.cpp GuadiGraphics/geometry_logger.cpp)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} PRIVATE -g -ltbb)

target_link_libraries(${PROJECT_NAME} PRIVATE libGaudi nanogui ${NANOGUI_EXTRA_LIBS})