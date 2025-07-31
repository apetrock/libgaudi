#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
#include <string>
#include <iostream>
#include <vector>

// Include the actual gaudi logger
#include "../../../include/gaudi/logger.hpp"
#include "gaudi/geometry_logger.hpp"

/**
 * WASM Integration with Real Gaudi Logger + Shared Memory Buffer
 * 
 * This demonstrates C++ to JavaScript communication using the actual gaudi
 * geometry logger with a shared memory buffer for efficient line data transfer.
 */

// Shared memory buffer for line data
// Format: [x0, y0, z0, x1, y1, z1, r, g, b, a] per line (10 floats per line)
std::vector<float> lineDataBuffer;
const int FLOATS_PER_LINE = 10;

// Simple function that returns a greeting
std::string hello_world() {
    return "Hello from C++ WebAssembly with Real Gaudi Logger!";
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

// Using Embind for more advanced binding
class MathHelper {
public:
    MathHelper(double initial_value) : value(initial_value) {}
    
    double getValue() const { return value; }
    void setValue(double new_value) { value = new_value; }
    
    double calculate(double input) {
        return value * input * 2.0;
    }
    
    std::string getStatus() {
        return "MathHelper ready - value: " + std::to_string(value);
    }

private:
    double value;
};

// Gaudi Logger Integration Functions
void log_line(double x0, double y0, double z0, double x1, double y1, double z1, 
             double r, double g, double b, double a) {
    gaudi::vec3 p0(x0, y0, z0);
    gaudi::vec3 p1(x1, y1, z1);
    gaudi::vec4 color(r, g, b, a);
    
    // Log to the actual gaudi logger
    gaudi::geometry_logger::line(p0, p1, color);
    
    // Also store in our shared buffer for JavaScript access
    lineDataBuffer.insert(lineDataBuffer.end(), {
        static_cast<float>(x0), static_cast<float>(y0), static_cast<float>(z0),
        static_cast<float>(x1), static_cast<float>(y1), static_cast<float>(z1),
        static_cast<float>(r), static_cast<float>(g), static_cast<float>(b), static_cast<float>(a)
    });
}

void clear_logger() {
    // Clear both the gaudi logger and our buffer
    gg::geometry_logger::clear();
    lineDataBuffer.clear();
    std::cout << "Gaudi logger and buffer cleared\n";
}

std::string test_logger() {
    gaudi::vec3 p0(0.0, 0.0, 0.0);
    gaudi::vec3 p1(1.0, 1.0, 1.0);
    gaudi::vec4 red(1.0, 0.0, 0.0, 1.0);
    gaudi::geometry_logger::line(p0, p1, red);
    
    // Also add to buffer
    lineDataBuffer.insert(lineDataBuffer.end(), {
        0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f
    });
    
    return "Test line logged to both Gaudi logger and shared buffer";
}

// Shared Memory Access Functions
size_t get_line_count() {
    return lineDataBuffer.size() / FLOATS_PER_LINE;
}

// Get pointer to the line data buffer for JavaScript access
float* get_line_data_ptr() {
    return lineDataBuffer.empty() ? nullptr : lineDataBuffer.data();
}

size_t get_buffer_size() {
    return lineDataBuffer.size();
}

// Frame simulation - generates animated lines
void simulate_frame(double time) {
    // Clear previous frame data
    lineDataBuffer.clear();
    
    // Generate some animated lines for this frame
    const int numLines = 20;
    for (int i = 0; i < numLines; ++i) {
        double angle = (2.0 * M_PI * i / numLines) + time;
        double radius = 2.0 + sin(time * 2.0 + i * 0.5);
        
        // Start point at origin
        gaudi::vec3 p0(0.0, 0.0, 0.0);
        
        // End point rotating around
        gaudi::vec3 p1(
            radius * cos(angle),
            radius * sin(angle),
            sin(time + i * 0.3)
        );
        
        // Color cycling through spectrum
        gaudi::vec4 color(
            0.5 + 0.5 * cos(time + i * 0.2),
            0.5 + 0.5 * cos(time + i * 0.2 + 2.094),
            0.5 + 0.5 * cos(time + i * 0.2 + 4.188),
            1.0
        );
        
        // Log to gaudi logger
        gaudi::geometry_logger::line(p0, p1, color);
        
        // Add to shared buffer
        lineDataBuffer.insert(lineDataBuffer.end(), {
            static_cast<float>(p0.x()), static_cast<float>(p0.y()), static_cast<float>(p0.z()),
            static_cast<float>(p1.x()), static_cast<float>(p1.y()), static_cast<float>(p1.z()),
            static_cast<float>(color.x()), static_cast<float>(color.y()), 
            static_cast<float>(color.z()), static_cast<float>(color.w())
        });
    }
}

std::string get_logger_status() {
    return "Gaudi Logger: " + std::to_string(get_line_count()) + " lines in buffer";
}

// Bind the class and functions to JavaScript
EMSCRIPTEN_BINDINGS(hello_module) {
    emscripten::function("hello_world", &hello_world);
    emscripten::function("add_numbers", &add_numbers);
    emscripten::function("multiply_floats", &multiply_floats);
    
    // Logger functions
    emscripten::function("log_line", &log_line);
    emscripten::function("clear_logger", &clear_logger);
    emscripten::function("test_logger", &test_logger);
    emscripten::function("get_line_count", &get_line_count);
    emscripten::function("get_logger_status", &get_logger_status);
    
    // Shared memory functions
    emscripten::function("get_line_data_ptr", &get_line_data_ptr, emscripten::allow_raw_pointers());
    emscripten::function("get_buffer_size", &get_buffer_size);
    emscripten::function("simulate_frame", &simulate_frame);
    
    emscripten::class_<MathHelper>("MathHelper")
        .constructor<double>()
        .function("getValue", &MathHelper::getValue)
        .function("setValue", &MathHelper::setValue)
        .function("calculate", &MathHelper::calculate)
        .function("getStatus", &MathHelper::getStatus);
}
