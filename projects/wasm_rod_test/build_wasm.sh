#!/bin/bash

# Build script for WebAssembly rod simulation
# Make sure you have Emscripten installed and activated:
# source /path/to/emsdk/emsdk_env.sh

set -e

echo "Building Rod Simulation for WebAssembly..."

# Create build directory
mkdir -p build
cd build

# Configure with Emscripten
emcmake cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-std=c++17"

# Build
emmake make -j$(nproc)

echo "Build complete! WASM files should be in dist/ directory"
echo "Files will be automatically copied to the React project"

# List generated files
echo "Generated files:"
ls -la ../dist/
