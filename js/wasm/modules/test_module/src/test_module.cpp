#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <string>

using namespace emscripten;

// Simple test function that returns a success message
extern "C" {
    EMSCRIPTEN_KEEPALIVE
    const char* test_function() {
        return "WASM test function executed successfully!";
    }
    
    EMSCRIPTEN_KEEPALIVE
    const char* get_test_message() {
        return "Hello from WASM test module!";
    }
    
    EMSCRIPTEN_KEEPALIVE
    int add_numbers(int a, int b) {
        return a + b;
    }
}

// Bind functions for easier JavaScript access
EMSCRIPTEN_BINDINGS(test_module) {
    function("testFunction", &test_function, allow_raw_pointer<ret_val>());
    function("getTestMessage", &get_test_message, allow_raw_pointer<ret_val>());
    function("addNumbers", &add_numbers);
} 