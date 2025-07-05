# Logger WASM Module
# Contains the Gaudi logger implementation and test modules

# Gaudi Logger Test Module
add_wasm_target(gaudi_logger_test
    SOURCES
        ./src/gaudi_logger_test.cpp
    USE_LOGGER
    EXPORT_NAME "GaudiLoggerModule"
    ES6
    OUTPUT_DIR "${PUBLIC_WASM_DIR}"
)

# Set target properties for IDE organization
set_target_properties(gaudi_logger_test PROPERTIES FOLDER "wasm/logger")

# Copy target for this module
add_custom_target(copy_logger_to_public
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "$<TARGET_FILE:gaudi_logger_test>"
        "${PUBLIC_WASM_DIR}/"
    DEPENDS gaudi_logger_test
    COMMENT "Copying logger WASM files to public directory"
)
