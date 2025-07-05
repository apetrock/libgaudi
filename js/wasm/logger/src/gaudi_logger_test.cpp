#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include <cmath>
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>

using namespace gaudi;

class GaudiLoggerTest {
private:
    int frame_count;
    double time;
    bool is_animating;

public:
    GaudiLoggerTest() : frame_count(0), time(0.0), is_animating(false) {
        // Initialize with some test data
        add_test_data();
    }

    void step() {
        if (!is_animating) return;
        
        frame_count++;
        time += 0.016; // ~60fps
        
        // Clear previous frame
        logger::clear();
        
        // Add animated test data
        add_animated_data();
    }

    void reset() {
        frame_count = 0;
        time = 0.0;
        is_animating = false;
        logger::clear();
        add_test_data();
    }

    void start() {
        is_animating = true;
    }

    void stop() {
        is_animating = false;
    }

    int get_frame_count() const { return frame_count; }
    double get_time() const { return time; }
    bool get_is_animating() const { return is_animating; }

private:
    void add_test_data() {
        // Add some static test lines
        logger::line(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec4(1.0, 0.0, 0.0, 1.0)); // Red line
        logger::line(vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), vec4(0.0, 1.0, 0.0, 1.0)); // Green line
        logger::line(vec3(0.0, 0.0, 0.0), vec3(0.0, 0.0, 1.0), vec4(0.0, 0.0, 1.0, 1.0)); // Blue line
        
        // Add some test points
        logger::point(vec3(0.5, 0.5, 0.5), vec4(1.0, 1.0, 0.0, 1.0)); // Yellow point
        logger::point(vec3(-0.5, -0.5, -0.5), vec4(1.0, 0.0, 1.0, 1.0)); // Magenta point
    }

    void add_animated_data() {
        // Add animated lines based on time
        double t = time;
        double radius = 0.5;
        
        for (int i = 0; i < 8; i++) {
            double angle = t + i * M_PI / 4;
            double x1 = radius * cos(angle);
            double y1 = radius * sin(angle);
            double x2 = radius * cos(angle + M_PI / 8);
            double y2 = radius * sin(angle + M_PI / 8);
            
            double r = 0.5 + 0.5 * cos(angle);
            double g = 0.5 + 0.5 * sin(angle);
            double b = 0.5 + 0.5 * cos(angle + M_PI / 3);
            
            logger::line(vec3(x1, y1, 0.0), vec3(x2, y2, 0.0), vec4(r, g, b, 1.0));
        }
        
        // Add animated points
        for (int i = 0; i < 5; i++) {
            double angle = t * 2 + i * M_PI / 2.5;
            double x = 0.3 * cos(angle);
            double y = 0.3 * sin(angle);
            double z = 0.2 * sin(t * 3 + i);
            
            logger::point(vec3(x, y, z), vec4(1.0, 0.5, 0.0, 1.0));
        }
    }
};

EMSCRIPTEN_BINDINGS(gaudi_logger_test) {
    emscripten::class_<GaudiLoggerTest>("GaudiLoggerTest")
        .constructor<>()
        .function("step", &GaudiLoggerTest::step)
        .function("reset", &GaudiLoggerTest::reset)
        .function("start", &GaudiLoggerTest::start)
        .function("stop", &GaudiLoggerTest::stop)
        .function("get_frame_count", &GaudiLoggerTest::get_frame_count)
        .function("get_time", &GaudiLoggerTest::get_time)
        .function("get_is_animating", &GaudiLoggerTest::get_is_animating);
}
