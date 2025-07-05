Write-Host "Building Rod Constraints WASM Module with Debug Support..." -ForegroundColor Green

# Set up Emscripten environment
& emsdk_env.ps1

# Create build directory and clean if it exists
if (Test-Path "build") {
    Write-Host "Cleaning previous build..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force "build"
}
New-Item -ItemType Directory -Path "build"
Set-Location build

# Configure with CMake in Debug mode
Write-Host "Configuring with CMake..." -ForegroundColor Yellow
emcmake cmake .. -DCMAKE_BUILD_TYPE=Debug

# Build the rod constraints test module
Write-Host "Building rod_constraints_test..." -ForegroundColor Yellow
cmake --build . --target rod_constraints_test

# Copy files to demo directory
Write-Host "Copying files to demo directory..." -ForegroundColor Yellow
cmake --build . --target copy_rod_constraints_to_demo

Write-Host ""
Write-Host "Build complete! Files copied to js/wasm/demo/" -ForegroundColor Green
Write-Host ""
Write-Host "Debug features enabled:" -ForegroundColor Yellow
Write-Host "- Stack traces and source maps"
Write-Host "- Memory safety checks"
Write-Host "- Assertions enabled"
Write-Host "- Line logging visualization"
Write-Host ""
Write-Host "To run the demo, start the development server and navigate to the rod constraints test." -ForegroundColor Cyan 