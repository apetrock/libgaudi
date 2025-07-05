# CMakeLists.txt for foo_demo WASM project
add_wasm_target(foo_demo
    SOURCES src/foo_demo.cpp
    LINK_FLAGS
        "-s EXPORT_NAME=\"Foo_demoModule\""
        "-s USE_ES6_IMPORT_META=0"
)

# Copy files to public directory after build
add_custom_command(TARGET foo_demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:foo_demo>"
        "${PUBLIC_WASM_DIR}/"
    COMMENT "Copying foo_demo files to public directory"
)
