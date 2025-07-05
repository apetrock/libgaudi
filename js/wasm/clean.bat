@echo off
REM Clean all WebAssembly build artifacts

echo Cleaning WebAssembly build artifacts...

REM Remove build directory
if exist build (
    echo Removing build directory...
    rmdir /s /q build
)

REM Remove individual build artifacts (in case of standalone builds)
if exist hello.js del hello.js
if exist hello.wasm del hello.wasm
if exist gaudi_logger_test.js del gaudi_logger_test.js
if exist gaudi_logger_test.wasm del gaudi_logger_test.wasm
if exist rod_simulation_main.js del rod_simulation_main.js
if exist rod_simulation_main.wasm del rod_simulation_main.wasm

REM Clean demo directory WASM files
if exist ..\demo\wasm (
    echo Cleaning demo directory...
    del /q ..\demo\wasm\*.js 2>nul
    del /q ..\demo\wasm\*.wasm 2>nul
)

REM Clean public/wasm directory (where files are served from for frontend)
if exist ..\public\wasm (
    echo Cleaning public/wasm directory...
    del /q ..\public\wasm\*.js 2>nul
    del /q ..\public\wasm\*.wasm 2>nul
)

echo ✅ Clean complete!
