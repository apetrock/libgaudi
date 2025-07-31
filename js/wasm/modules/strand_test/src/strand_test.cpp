#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "gaudi/duchamp/strand_constraints_test.hpp"
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>

using namespace gaudi;

class StrandTest {
private:
    duchamp::block_test::ptr test_instance;
    int frame_count;

public:
StrandTest() 
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

EMSCRIPTEN_BINDINGS(strand_test) {
    emscripten::class_<StrandTest>("StrandTest")
        .constructor<>()
        .function("step", &StrandTest::step)
        .function("reset", &StrandTest::reset)
        .function("get_frame_count", &StrandTest::get_frame_count)
        .function("get_time", &StrandTest::get_time);
}
