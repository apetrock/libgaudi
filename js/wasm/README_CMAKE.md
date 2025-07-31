# Gaudi WASM Components

This directory contains WebAssembly modules for the Gaudi library using a proper CMake build system.

## Building with CMake

The WASM components now use CMake with Emscripten for a cleaner, more maintainable build process.

### Prerequisites

- Emscripten SDK installed and activated
- CMake 3.8 or later  
- Ninja build system (recommended)

### Quick Start

```bash
# Build all targets (Release)
.\build_wasm.bat
# or
.\build_wasm.ps1

# Build in Debug mode
.\build_wasm.bat --debug
# or  
.\build_wasm.ps1 -Debug

# Build specific target
.\build_wasm.bat --target gaudi_logger_test
# or
.\build_wasm.ps1 -Target gaudi_logger_test

# Clean build
.\build_wasm.bat --clean
# or
.\build_wasm.ps1 -Clean
```

### Available Targets

The CMake structure is organized into subdirectories, each with their own CMakeLists.txt:

**logger/** - Logger functionality
- **gaudi_logger_test** - Test module for Gaudi logger functionality

**rod_simulation/** - Rod simulation modules  
- **rod_simulation_main** - Main rod simulation module

**examples/** - Simple examples and tests
- **hello_world** - Simple hello world test module

**Meta targets:**
- **all_wasm_targets** - Builds all the above targets
- **copy_to_demo** - Copies built files to ../demo/wasm/

### Manual CMake Usage

```bash
mkdir build && cd build
emcmake cmake .. -G Ninja
cmake --build . --target all_wasm_targets
cmake --build . --target copy_to_demo
```

### Build Configuration

- **Release** (default): Optimized builds with `-Os` 
- **Debug**: Debug builds with `-O2 -g -gsource-map`
- **ES6**: Enable ES6 module exports with `--es6` flag

### Output

Built files are placed in:
- `build/` directory (build artifacts)
- `../demo/wasm/` directory (copied for use)

### Files

**Root directory:**
- `CMakeLists.txt` - Main CMake configuration with shared functions
- `build_wasm.bat` - Windows batch build script
- `build_wasm.ps1` - PowerShell build script (recommended)
- `wasm_geometry_logger.cpp` - Logger implementation for WASM (shared)
- `gaudi_logger_test.cpp` - Test for the logger functionality
- `main.cpp` - Main rod simulation module
- `hello.cpp` - Simple hello world test

**Subdirectories:**
- `logger/CMakeLists.txt` - Logger module build configuration
- `rod_simulation/CMakeLists.txt` - Rod simulation build configuration  
- `examples/CMakeLists.txt` - Examples build configuration

## Project Structure

```
wasm/
├── CMakeLists.txt           # Root CMake configuration
├── build_wasm.ps1           # PowerShell build script
├── build_wasm.bat           # Batch build script
├── wasm_geometry_logger.cpp          # Shared logger implementation
├── gaudi_logger_test.cpp    # Logger test
├── main.cpp                 # Rod simulation main
├── hello.cpp                # Hello world example
├── logger/
│   └── CMakeLists.txt       # Logger module configuration
├── rod_simulation/
│   └── CMakeLists.txt       # Rod simulation configuration
└── examples/
    └── CMakeLists.txt       # Examples configuration
```

## Migration from Ad-hoc Scripts

The old PowerShell scripts (`build_gaudi_logger_test.ps1`, etc.) are now replaced by this CMake-based approach, which provides:

- Consistent build configuration across targets
- Better dependency management
- Easier maintenance and extension
- Integration with the main libgaudi CMake structure
- Support for Debug/Release builds
- Automated file copying to demo directory

## Original Phase Documentation

For historical reference, the original phase-based development documentation has been preserved below:

---

## Phase 1 Goals (Completed)

- ✅ Basic C++ to JavaScript communication
- ✅ Function calls with different parameter types
- ✅ Class instantiation and method calls
- ✅ Interactive React component integration
- ✅ Ready for Phase 2: Gaudi geometry_logger integration

## Prerequisites

You need the Emscripten SDK installed to build the WebAssembly module.

### Installing Emscripten (Windows)

1. Clone the Emscripten SDK and install following the official instructions
2. Run `emsdk_env.bat` to set up the environment before building
