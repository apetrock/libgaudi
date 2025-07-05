# Rod Constraints Debug Improvements

This document describes the enhancements made to the rod constraints simulation for better debugging and visualization.

## Overview

The rod constraints demo has been enhanced with:
1. **Better error handling and stack traces** in WASM builds
2. **Comprehensive line logging** throughout the simulation
3. **Real-time debug visualization** in the React interface
4. **Logger API integration** for C++ to JavaScript communication

## Debug Features

### 1. Enhanced WASM Build Configuration

The `js/wasm/rod_constraints/CMakeLists.txt` now includes debug flags:

```cmake
"-s ASSERTIONS=1"           # Enable runtime assertions
"-s SAFE_HEAP=1"           # Memory safety checks
"-s DEMANGLE_SUPPORT=1"    # Better stack traces
"-s STACK_OVERFLOW_CHECK=2" # Stack overflow detection
"-s ABORT_ON_STACK_OVERFLOW=1" # Abort on stack overflow
"-s EMIT_SOURCE_MAP=1"     # Source maps for debugging
"-g4"                      # Full debug symbols
"-O0"                      # No optimization for debugging
```

### 2. Line Logging in C++ Simulation

The `include/gaudi/duchamp/rod_constraints_test.hpp` now includes comprehensive logging:

- **Coordinate axes**: Red (X), Green (Y), Blue (Z) axes
- **Rod geometry**: Gray lines showing the current rod shape
- **Tangent point gradients**: Orange lines showing gradient directions
- **Boundary gradients**: Green (positive) / Red (negative) force vectors
- **Constraint forces**: Yellow lines showing applied forces

### 3. Logger API Integration

The WASM module (`js/wasm/rod_constraints/src/rod_constraints_wasm.cpp`) exports:

```cpp
// Logger data access
get_line_count(): number
get_line_positions_ptr(): number
get_line_positions_size(): number
get_line_colors_ptr(): number
get_line_colors_size(): number
get_point_count(): number
get_point_positions_ptr(): number
get_point_positions_size(): number
get_point_colors_ptr(): number
get_point_colors_size(): number

// Manual logging
add_line(x0, y0, z0, x1, y1, z1, r, g, b, a): void
add_point(x, y, z, r, g, b, a): void
clear_logger(): void
```

### 4. React Visualization

The `js/src/components/RodConstraintsTest.tsx` component now includes:

- **Debug line visualization**: Real-time rendering of C++ logged lines
- **Debug point visualization**: Real-time rendering of C++ logged points
- **Toggle controls**: Show/hide debug lines and points
- **Debug info overlay**: Live counters for lines, points, and vertices
- **Error handling**: Better error messages and recovery

## Building the Debug Version

### Windows (Batch)
```batch
cd js/wasm
build_rod_constraints_test.bat
```

### Windows (PowerShell)
```powershell
cd js/wasm
./build_rod_constraints_test.ps1
```

### Manual Build
```bash
cd js/wasm
mkdir build && cd build
emcmake cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . --target rod_constraints_test
cmake --build . --target copy_rod_constraints_to_demo
```

## Debug Visualization Guide

### Color Coding
- **Red lines**: X-axis and negative boundary gradients
- **Green lines**: Y-axis and positive boundary gradients  
- **Blue lines**: Z-axis
- **Orange lines**: Tangent point gradients
- **Yellow lines**: Constraint forces
- **Gray lines**: Rod geometry

### Controls
- **Play/Pause**: Control simulation playback
- **Step**: Advance simulation one frame
- **Reset**: Restart simulation
- **Lines toggle**: Show/hide debug lines
- **Points toggle**: Show/hide debug points
- **Growth Rate**: Adjust rod growth parameter
- **Constraint Weight**: Adjust constraint solver weight

## Troubleshooting

### Common Issues

1. **WASM module not found**
   - Run the build script first
   - Check that files are copied to `js/wasm/demo/`

2. **C++ crashes without stack trace**
   - Ensure debug build is used (`-DCMAKE_BUILD_TYPE=Debug`)
   - Check browser console for detailed error messages

3. **No debug lines visible**
   - Verify logger API is working
   - Check that debug visualization toggles are enabled
   - Look for errors in browser console

### Debug Console Commands

In the browser console, you can manually test the logger:

```javascript
// Get the WASM instance
const instance = window.wasmInstance;

// Add a test line
instance.add_line(0, 0, 0, 1, 1, 1, 1, 0, 0, 1);

// Check line count
console.log(instance.get_line_count());

// Clear logger
instance.clear_logger();
```

## Performance Notes

- Debug builds are slower due to safety checks
- Line logging adds some overhead
- For production, use release builds with `-O2` optimization
- Debug visualization can be disabled for better performance

## Future Enhancements

- Add more detailed constraint visualization
- Implement force magnitude scaling
- Add animation controls for debug visualization
- Support for different debug visualization modes
- Export debug data for external analysis 