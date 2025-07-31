#include "gaudi/duchamp/path_constraint_test.hpp"
#include <emscripten/bind.h>

using namespace gaudi;

class PathTest {
private:
    duchamp::block_test::ptr test_instance;
    int frame_count;

public:
PathTest() 
        : frame_count(0) {
        test_instance = duchamp::block_test::create();
    }

    void step() {
        frame_count++;
        test_instance->step(frame_count);
    }

    void reset() {
        frame_count = 0;
        test_instance = duchamp::block_test::create();
    }

    // Getters
    int get_frame_count() const { return frame_count; }
    double get_time() const { return frame_count * 0.016; } // ~60fps

private:
    // No private methods needed - all logic is in the duchamp module
};

EMSCRIPTEN_BINDINGS(path_test) {
    emscripten::class_<PathTest>("PathTest")
        .constructor<>()
        .function("step", &PathTest::step)
        .function("reset", &PathTest::reset)
        .function("get_frame_count", &PathTest::get_frame_count)
        .function("get_time", &PathTest::get_time);
}
