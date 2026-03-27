# morton_test: headless project, no GL dependency
add_executable(morton_test main.cpp)
target_link_libraries(morton_test PRIVATE gaudi::core)

if(MSVC)
    target_compile_options(morton_test PRIVATE /W4)
else()
    target_compile_options(morton_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
