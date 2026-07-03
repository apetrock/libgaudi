cmake_minimum_required(VERSION 3.20)
project(liblombardi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create test executable
add_executable(liblombardi_test tests/liblombardi_test.cpp)

# Include liblombardi directory (parent of test_nodes)
target_include_directories(liblombardi_test PRIVATE ${PROJECT_SOURCE_DIR}/include/liblombardi)

# Link with libgaudi core
target_link_libraries(liblombardi_test gaudi_core)
