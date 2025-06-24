@echo off
REM Build script for WebAssembly rod simulation on Windows
REM Make sure you have Emscripten installed and the environment activated

echo Building Rod Simulation for WebAssembly...

REM Create build directory
if not exist build mkdir build
cd build

REM Configure with Emscripten
emcmake cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-std=c++17"

REM Build
emmake make

echo Build complete! WASM files should be in dist/ directory
echo Files will be automatically copied to the React project

REM List generated files
echo Generated files:
dir ..\dist\

pause
