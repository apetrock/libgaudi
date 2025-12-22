# CMakeLists.txt for module: bvh_test
# BVH testing module that compares BVH nearest search against brute force

# Include common configuration
include(../../cmake-templates/common.cmake)

# Create the WASM target
add_wasm_target(bvh_test
    SOURCES "src/bvh_test.cpp"
    OUTPUT_FILENAME "bvh_test"
    USE_LOGGER
    EXPORT_NAME "BvhTestModule"
)


