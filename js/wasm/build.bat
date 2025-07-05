@echo off
REM Batch script to build all WebAssembly modules
REM Builds logger, rod simulation, and examples using CMake

echo Building All WebAssembly Modules...

REM Check if emcc is available
where emcc >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: emcc not found. Please install Emscripten SDK.
    echo.
    echo Installation instructions:
    echo 1. git clone https://github.com/emscripten-core/emsdk.git
    echo 2. cd emsdk
    echo 3. emsdk install latest
    echo 4. emsdk activate latest
    echo 5. emsdk_env.bat
    exit /b 1
)

for /f "tokens=*" %%i in ('emcc --version 2^>^&1 ^| findstr /n "^" ^| findstr "^1:"') do (
    set "EMCC_VERSION=%%i"
    set "EMCC_VERSION=!EMCC_VERSION:~2!"
)
echo Found Emscripten: %EMCC_VERSION%

REM Create build directory
if not exist build mkdir build
cd build

REM Configure with Emscripten
echo Configuring CMake with Emscripten...
call emcmake cmake .. -DCMAKE_BUILD_TYPE=Debug
echo CMake configuration complete.
if %ERRORLEVEL% NEQ 0 (
    echo ❌ CMake configuration failed
    cd ..
    exit /b %ERRORLEVEL%
)

REM Build all WASM targets (copy commands are part of POST_BUILD steps)
echo Building all WebAssembly targets...
call cmake --build . --target all_wasm_targets

if %ERRORLEVEL% EQU 0 (
    echo Copying files to public directory...
    call cmake --build . --target copy_all_to_public
    
    if %ERRORLEVEL% EQU 0 (
        echo ✅ All builds successful!
        echo.
        echo Generated WASM modules:
        echo - Logger API (gaudi_logger_test.js/.wasm)
        echo - Rod Simulation (rod_constraints_test.js/.wasm) 
        echo - Hello World Example (hello_world.js/.wasm)
        echo - Foo Demo (foo_demo.js/.wasm)
        echo.
        echo Files copied to: ../public/wasm/
        cd ..
        dir /b ..\public\wasm\*.js ..\public\wasm\*.wasm 2>nul
    ) else (
        echo ❌ Copy to public directory failed
        cd ..
        exit /b %ERRORLEVEL%
    )
) else (
    echo ❌ Build failed
    cd ..
    exit /b %ERRORLEVEL%
)
