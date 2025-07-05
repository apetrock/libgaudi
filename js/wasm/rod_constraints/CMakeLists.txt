# Rod Constraints WASM Module
# Contains the rod constraints test implementation

# Rod Constraints Test Module
add_wasm_target(rod_constraints_test
    SOURCES
        ./src/rod_constraints_wasm.cpp
    USE_LOGGER
    EXPORT_NAME "RodConstraintsModule"
    ES6
    OUTPUT_DIR "${PUBLIC_WASM_DIR}"
)

# Set target properties for IDE organization
set_target_properties(rod_constraints_test PROPERTIES FOLDER "wasm/rod_constraints")
target_include_directories(rod_constraints_test PRIVATE ${EIGEN_INCLUDE_DIR})

# Copy target for this module
add_custom_target(copy_rod_constraints_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:rod_constraints_test>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS rod_constraints_test
    COMMENT "Copying rod constraints WASM files to public directory"
)
