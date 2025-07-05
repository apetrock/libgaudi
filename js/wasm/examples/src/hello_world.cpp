#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Phase 1: Basic WASM "Hello World" Integration + Line Logger + Three.js Integration
 * 
 * This demonstrates C++ to JavaScript communication and includes
 * a simplified line logger for visual debugging with Three.js rendering.
 */

// Simplified vector types for WebAssembly
struct Vec3 {
    double x, y, z;
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

struct Vec4 {
    double r, g, b, a;
    Vec4(double r = 0, double g = 0, double b = 0, double a = 1) : r(r), g(g), b(b), a(a) {}
};

// Simple line logger for WebAssembly
class SimpleLineLogger {
private:
    struct LoggedLine {
        Vec3 p0, p1;
        Vec4 color;
        LoggedLine(const Vec3& start, const Vec3& end, const Vec4& col) 
            : p0(start), p1(end), color(col) {}
    };
    
    std::vector<LoggedLine> lines;
    
public:
    static SimpleLineLogger& getInstance() {
        static SimpleLineLogger instance;
        return instance;
    }
    
    void logLine(const Vec3& p0, const Vec3& p1, const Vec4& color) {
        lines.emplace_back(p0, p1, color);
        std::cout << "Logged line from (" << p0.x << "," << p0.y << "," << p0.z 
                  << ") to (" << p1.x << "," << p1.y << "," << p1.z 
                  << ") color RGBA(" << color.r << "," << color.g << "," << color.b << "," << color.a << ")\n";
    }
    
    void clear() {
        lines.clear();
        std::cout << "Line logger cleared\n";
    }
    
    size_t getLineCount() const {
        return lines.size();
    }
      std::string getStatus() const {
        return "SimpleLineLogger: " + std::to_string(lines.size()) + " lines logged";
    }
    
    // Get line data for rendering
    const LoggedLine& getLine(size_t index) const {
        return lines[index];
    }
    
    // Get all line data as flat arrays for JavaScript
    std::vector<double> getLinePositions() const {
        std::vector<double> positions;
        positions.reserve(lines.size() * 6); // 2 points * 3 coords per line
        for (const auto& line : lines) {
            positions.push_back(line.p0.x);
            positions.push_back(line.p0.y);
            positions.push_back(line.p0.z);
            positions.push_back(line.p1.x);
            positions.push_back(line.p1.y);
            positions.push_back(line.p1.z);
        }
        return positions;
    }
    
    std::vector<double> getLineColors() const {
        std::vector<double> colors;
        colors.reserve(lines.size() * 8); // 2 points * 4 color components per line
        for (const auto& line : lines) {
            // Color for start point
            colors.push_back(line.color.r);
            colors.push_back(line.color.g);
            colors.push_back(line.color.b);
            colors.push_back(line.color.a);
            // Color for end point (same)
            colors.push_back(line.color.r);
            colors.push_back(line.color.g);
            colors.push_back(line.color.b);
            colors.push_back(line.color.a);
        }
        return colors;
    }
};

// Simple function that returns a greeting
std::string hello_world() {
    return "Hello from C++ WebAssembly!";
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

// Line Logger Wrapper Functions
void log_line(double x0, double y0, double z0, double x1, double y1, double z1, 
             double r, double g, double b, double a) {
    Vec3 p0(x0, y0, z0);
    Vec3 p1(x1, y1, z1);
    Vec4 color(r, g, b, a);
    SimpleLineLogger::getInstance().logLine(p0, p1, color);
}

void clear_logger() {
    SimpleLineLogger::getInstance().clear();
}

std::string test_logger() {
    Vec3 p0(0.0, 0.0, 0.0);
    Vec3 p1(1.0, 1.0, 1.0);
    Vec4 red(1.0, 0.0, 0.0, 1.0);
    SimpleLineLogger::getInstance().logLine(p0, p1, red);
    return "Test line logged from (0,0,0) to (1,1,1) in red";
}

size_t get_line_count() {
    return SimpleLineLogger::getInstance().getLineCount();
}

std::string get_logger_status() {
    return SimpleLineLogger::getInstance().getStatus();
}

// Get line data for Three.js rendering
emscripten::val get_line_positions() {
    auto positions = SimpleLineLogger::getInstance().getLinePositions();
    return emscripten::val(emscripten::typed_memory_view(positions.size(), positions.data()));
}

emscripten::val get_line_colors() {
    auto colors = SimpleLineLogger::getInstance().getLineColors();
    return emscripten::val(emscripten::typed_memory_view(colors.size(), colors.data()));
}

// Frame simulation - generate animated lines
void simulate_frame(double time) {
    // Clear previous frame
    SimpleLineLogger::getInstance().clear();
    
    // Generate some animated lines
    double radius = 2.0;
    int numLines = 8;
    
    for (int i = 0; i < numLines; i++) {
        double angle = (2.0 * M_PI * i / numLines) + time * 0.5;
        double nextAngle = (2.0 * M_PI * (i + 1) / numLines) + time * 0.5;
        
        // Rotating lines around origin
        Vec3 p0(radius * cos(angle), radius * sin(angle), sin(time + i) * 0.5);
        Vec3 p1(0.0, 0.0, 0.0); // Center
        
        // Rainbow colors based on time and index
        double hue = (double)i / numLines + time * 0.1;
        Vec4 color(
            0.5 + 0.5 * cos(2.0 * M_PI * (hue + 0.0)),
            0.5 + 0.5 * cos(2.0 * M_PI * (hue + 0.33)),
            0.5 + 0.5 * cos(2.0 * M_PI * (hue + 0.66)),
            1.0
        );
        
        SimpleLineLogger::getInstance().logLine(p0, p1, color);
        
        // Add outer ring connections
        Vec3 p2(radius * cos(nextAngle), radius * sin(nextAngle), sin(time + i + 1) * 0.5);
        SimpleLineLogger::getInstance().logLine(p0, p2, color);
    }
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
    
    // Frame simulation and data access
    emscripten::function("get_line_positions", &get_line_positions);
    emscripten::function("get_line_colors", &get_line_colors);
    emscripten::function("simulate_frame", &simulate_frame);
    
    emscripten::class_<MathHelper>("MathHelper")
        .constructor<double>()
        .function("getValue", &MathHelper::getValue)
        .function("setValue", &MathHelper::setValue)
        .function("calculate", &MathHelper::calculate)
        .function("getStatus", &MathHelper::getStatus);
}
