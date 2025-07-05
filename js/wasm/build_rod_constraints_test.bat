@echo off
echo Building Rod Constraints WASM Module with Debug Support...

REM Set up Emscripten environment
call emsdk_env.bat

REM Create build directory and clean if it exists
if exist "build" (
    echo Cleaning previous build...
    rmdir /s /q build
)
mkdir build
cd build

REM Configure with CMake in Debug mode
echo Configuring with CMake...
emcmake cmake .. -DCMAKE_BUILD_TYPE=Debug

REM Build the rod constraints test module
echo Building rod_constraints_test...
cmake --build . --target rod_constraints_test

REM Copy files to demo directory
echo Copying files to demo directory...
cmake --build . --target copy_rod_constraints_to_demo

echo.
echo Build complete! Files copied to js/wasm/demo/
echo.
echo Debug features enabled:
echo - Stack traces and source maps
echo - Memory safety checks
echo - Assertions enabled
echo - Line logging visualization
echo.
echo To run the demo, start the development server and navigate to the rod constraints test. 