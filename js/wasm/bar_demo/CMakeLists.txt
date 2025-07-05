# CMakeLists.txt for bar_demo WASM project
add_wasm_target(bar_demo
    SOURCES src/bar_demo.cpp
    LINK_FLAGS
        "-s EXPORT_NAME=\"Bar_demoModule\""
        "-s USE_ES6_IMPORT_META=0"
)

# Copy files to public directory after build
add_custom_command(TARGET bar_demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:bar_demo>"
        "${PUBLIC_WASM_DIR}/"
    COMMENT "Copying bar_demo files to public directory"
)
