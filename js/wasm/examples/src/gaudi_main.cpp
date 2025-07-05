#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
#include <string>
#include <iostream>

// Include our new API
#include "api/gaudi_wasm_api.hpp"

/**
 * Main WebAssembly Entry Point with Gaudi Line Logger Integration
 * 
 * This demonstrates the integration between Gaudi's C++ line logger
 * and JavaScript using shared memory buffers for efficient data transfer.
 */

// Simple hello world functions for basic testing
std::string hello_world() {
    return "Hello from Gaudi WebAssembly with Line Logger API!";
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    int add_numbers(int a, int b) {
        return a + b;
    }
    
    EMSCRIPTEN_KEEPALIVE
    double multiply_floats(double x, double y) {
        return x * y;
    }
}

// Main binding
EMSCRIPTEN_BINDINGS(gaudi_main) {
    emscripten::function("hello_world", &hello_world);
    emscripten::function("add_numbers", &add_numbers);
    emscripten::function("multiply_floats", &multiply_floats);
}
