// WASM Unit Test Runner
// Runs all gaudi tests and outputs results to console

#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <string>
#include <sstream>

// Include test framework
#include "gaudi/test/test.hpp"

// Include test files that don't require file system access
#include "gaudi/test/view_tests.hpp"
#include "gaudi/test/simplex_algorithm_tests.hpp"
#include "gaudi/test/morton_simplex_tests.hpp"

// BVH tests use procedurally generated meshes (no asset injection needed)
#include "gaudi/test/bvh_tests.hpp"

using namespace gaudi;

// Run all tests and return a summary string
std::string run_all_tests() {
    std::ostringstream out;
    
    out << "=== GAUDI Unit Tests ===" << std::endl;
    out << std::endl;
    
    auto result = test::Registry::run_all(false);
    
    out << std::endl;
    out << "========================" << std::endl;
    out << "Total:  " << result.total << std::endl;
    out << "Passed: " << result.passed << std::endl;
    out << "Failed: " << result.failed << std::endl;
    out << "Assertion failures: " << result.assertions_failed << std::endl;
    out << "========================" << std::endl;
    
    if (result.failed == 0) {
        out << "ALL TESTS PASSED!" << std::endl;
    } else {
        out << "SOME TESTS FAILED!" << std::endl;
    }
    
    return out.str();
}

// Run tests without mesh-dependent tests
std::string run_core_tests() {
    std::ostringstream out;
    
    out << "=== GAUDI Core Tests (no mesh data required) ===" << std::endl;
    out << std::endl;
    
    // Note: The GAUDI_TEST macro auto-registers all tests, so we run all
    // Tests that require mesh data will fail gracefully if data not set
    auto result = test::Registry::run_all(false);
    
    out << std::endl;
    out << "========================" << std::endl;
    out << "Total:  " << result.total << std::endl;
    out << "Passed: " << result.passed << std::endl;
    out << "Failed: " << result.failed << std::endl;
    out << "========================" << std::endl;
    
    return out.str();
}

// Bind functions for JavaScript access
EMSCRIPTEN_BINDINGS(unit_tests) {
    emscripten::function("runAllTests", &run_all_tests);
    emscripten::function("runCoreTests", &run_core_tests);
}
