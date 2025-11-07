cmake_minimum_required(VERSION 3.10roject(morton_test)

set(CMAKE_CXX_STANDARD17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/ext/eigen-3.4.0)

# Create executable
add_executable(morton_test main.cpp)

# Link libraries
target_link_libraries(morton_test Eigen3::Eigen)

# Set compiler flags
if(MSVC)
    target_compile_options(morton_test PRIVATE /W4)
else()
    target_compile_options(morton_test PRIVATE -Wall -Wextra -Wpedantic)
endif() 