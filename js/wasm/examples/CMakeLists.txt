# Examples WASM Modules
# Contains simple example and test modules

# Hello World Example
add_wasm_target(hello_world
    SOURCES
        ./src/hello_world.cpp
    EXPORT_NAME "HelloWorldModule"
    ES6
    OUTPUT_DIR "${PUBLIC_WASM_DIR}"
)

# Set target properties for IDE organization
set_target_properties(hello_world PROPERTIES FOLDER "wasm/examples")

# Future examples can be added here
# add_wasm_target(geometry_utils
#     SOURCES
#         ../geometry_utils.cpp
#     EXPORT_NAME "GeometryUtilsModule"
#     ES6
#     OUTPUT_DIR "${DEMO_WASM_DIR}"
# )

# Copy target for this module
add_custom_target(copy_hello_world_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:hello_world>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS hello_world
    COMMENT "Copying hello world WASM files to public directory"
)
