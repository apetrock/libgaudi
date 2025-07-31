#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
#include <string>
#include <cmath>
#include "../../../include/gaudi/logger.hpp"
#include "gaudi/geometry_logger.hpp"

/**
 * Main WebAssembly module with real Gaudi Logger API integration
 * 
 * This demonstrates the real gaudi::logger interface working in WebAssembly
 * with efficient line data transfer to Three.js.
 */

// Simple frame counter for animation
static int frame_count = 0;

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    std::string hello_world() {
        return "Hello from libgaudi WebAssembly with real Gaudi Logger API!";
    }    EMSCRIPTEN_KEEPALIVE
    void simulate_frame() {
        // Clear previous frame's lines from the Gaudi logger
        gaudi::geometry_logger::clear();
        
        // Add animated lines for this frame using the real Gaudi logger
        float time = frame_count * 0.1f;
        
        // Rotating coordinate frame
        for (int i = 0; i < 3; i++) {
            float angle = time + i * 2.0f * M_PI / 3.0f;
            float x = cos(angle);
            float y = sin(angle);
            
            double r = (i == 0) ? 1.0 : 0.0;
            double g = (i == 1) ? 1.0 : 0.0;
            double b = (i == 2) ? 1.0 : 0.0;
            
            gaudi::vec3 p0(0, 0, 0);
            gaudi::vec3 p1(x, y, 0.5);
            gaudi::vec4 color(r, g, b, 1.0);
            gaudi::geometry_logger::line(p0, p1, color);
        }
        
        // Sine wave
        const int wave_segments = 20;
        for (int i = 0; i < wave_segments - 1; i++) {
            float t0 = float(i) / float(wave_segments - 1);
            float t1 = float(i + 1) / float(wave_segments - 1);
            
            float x0 = t0 * 4.0f - 2.0f;
            float x1 = t1 * 4.0f - 2.0f;            float y0 = sin(x0 + time) * 0.5f;
            float y1 = sin(x1 + time) * 0.5f;
            
            gaudi::vec3 p0(x0, y0, 1.0);
            gaudi::vec3 p1(x1, y1, 1.0);
            gaudi::vec4 color(0.5, 1.0, 0.5, 1.0);
            gaudi::geometry_logger::line(p0, p1, color);
        }
        
        frame_count++;
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_frame_count() {
        return frame_count;
    }    EMSCRIPTEN_KEEPALIVE
    void reset_animation() {
        frame_count = 0;
        gaudi::geometry_logger::clear();
    }
    
    // Data access functions for JavaScript
    EMSCRIPTEN_KEEPALIVE
    int get_line_count() {
        return gaudi::logger::get_lines().size() / 2;
    }
    
    EMSCRIPTEN_KEEPALIVE
    uintptr_t get_line_positions() {
        const auto &lines = gaudi::logger::get_lines();
        return reinterpret_cast<uintptr_t>(lines.data());
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_line_positions_size() {
        return gaudi::logger::get_lines().size() * 3 * sizeof(double);
    }
    
    EMSCRIPTEN_KEEPALIVE
    uintptr_t get_line_colors() {
        const auto &colors = gaudi::logger::get_line_colors();
        return reinterpret_cast<uintptr_t>(colors.data());
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_line_colors_size() {
        return gaudi::logger::get_line_colors().size() * 4 * sizeof(double);
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_point_count() {
        return gaudi::logger::get_points().size();
    }
    
    EMSCRIPTEN_KEEPALIVE
    uintptr_t get_point_positions() {
        const auto &points = gaudi::logger::get_points();
        return reinterpret_cast<uintptr_t>(points.data());
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_point_positions_size() {
        return gaudi::logger::get_points().size() * 3 * sizeof(double);
    }
    
    EMSCRIPTEN_KEEPALIVE
    uintptr_t get_point_colors() {
        const auto &colors = gaudi::logger::get_point_colors();
        return reinterpret_cast<uintptr_t>(colors.data());
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_point_colors_size() {
        return gaudi::logger::get_point_colors().size() * 4 * sizeof(double);
    }
}

// Embind bindings
EMSCRIPTEN_BINDINGS(main_module) {
    emscripten::function("hello_world", &hello_world);
    emscripten::function("simulate_frame", &simulate_frame);
    emscripten::function("get_frame_count", &get_frame_count);
    emscripten::function("reset_animation", &reset_animation);
    
    // Logger data access
    emscripten::function("get_line_count", &get_line_count);
    emscripten::function("get_line_positions", &get_line_positions);
    emscripten::function("get_line_positions_size", &get_line_positions_size);
    emscripten::function("get_line_colors", &get_line_colors);
    emscripten::function("get_line_colors_size", &get_line_colors_size);
    
    emscripten::function("get_point_count", &get_point_count);
    emscripten::function("get_point_positions", &get_point_positions);
    emscripten::function("get_point_positions_size", &get_point_positions_size);
    emscripten::function("get_point_colors", &get_point_colors);
    emscripten::function("get_point_colors_size", &get_point_colors_size);
}
