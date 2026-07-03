# Testing

## Overview

libGaudi uses a custom test framework (not Google Test) that's integrated into the main test runner. Tests are automatically discovered and executed when running `gaudi_tests`.

## Testing Paradigm

### Test Framework

The test framework is defined in `include/gaudi/test/test.hpp` and uses:

- **`GAUDI_TEST(name)`** - Macro to register a test function
- **`GAUDI_ASSERT(expr)`** - Fatal assertion (throws on failure, stops test)
- **`GAUDI_EXPECT(expr)`** - Non-fatal assertion (records failure but continues)
- **`GAUDI_ASSERT_THROW(expr)`** - Fatal assertion that expects an exception

### Test Registration

Tests are registered using the `GAUDI_TEST(name)` macro. The macro automatically adds the test to the global test registry when the file is compiled. No manual registration needed.

**Example:**
```cpp
#include "gaudi/test/test.hpp"

GAUDI_TEST(MyTest) {
    // Test implementation here
    GAUDI_ASSERT(condition1);
    GAUDI_EXPECT(condition2);
}
```

## Test Organization

### Test Files

Tests are organized in `libgaudi/tests/` directory. Each test file typically includes the framework header and defines one or more test cases using `GAUDI_TEST()`.

**Pattern from existing tests:**
- `gaudi_tests.cpp` - Includes all other test files
- Individual test files (e.g., `bvh_tests.cpp`, `shell_dynamic_tests.cpp`) - Each contains multiple test functions

### Running Tests

#### Option 1: Run all tests
```bash
cd build
./gaudi_tests
```

#### Option 2: Run specific test via ctest
```bash
cd build
ctest -R MyTest  # Run tests matching "MyTest"
```

#### Option 3: Run all tests with ctest
```bash
cd build
ctest
```

### Adding a New Test

1. Create a new .cpp file in `libgaudi/tests/` or add test functions to an existing file

2. Use `GAUDI_TEST(name)` to register tests:
```cpp
#include "gaudi/test/test.hpp"

GAUDI_TEST(MyNewTest) {
    GraphContext ctx;
    auto node = ctx.create_node<SomeNode>();
    GAUDI_ASSERT(node != nullptr);
}
```

3. Include this file in `libgaudi/tests/gaudi_tests.cpp`:
```cpp
#include "liblombardi_tests.cpp"  // Or whatever you named your file
```

4. Rebuild and run:
```bash
cmake --build . --target gaudi_tests
./gaudi_tests
```

## Test Infrastructure

### Dependencies

All tests depend on:
- `gaudi::core` - Core library with common types
- `gaudi::headless` - Headless logger backends

### Test Dependencies

When adding tests that depend on new components, ensure the component is linked to the test target. See `libgaudi/tests/CMakeLists.txt` for examples.

## Debugging Tests

### Failures

If a test fails:
- The test framework prints the failure location (file:line)
- For fatal assertions, the test stops immediately
- The test framework collects all failures and reports them at the end

### Verbose Output

The test framework uses the console logger to output test results. You can check the build output for detailed information about which tests passed/failed.

## Known Limitations

- Tests currently use a custom framework (not Google Test)
- Some tests use plain `assert()` instead of the framework assertions
- **Future:** Custom framework is planned to be made compatible with Google Test
- The testing infrastructure requires `BUILD_TESTING=ON` when configuring cmake
