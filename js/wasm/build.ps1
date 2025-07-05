# PowerShell script to build WebAssembly "Hello World" module
# Phase 1: Basic WASM integration

Write-Host "Building WebAssembly Hello World Module..." -ForegroundColor Green

# Check if emcc is available
try {
    $emccVersion = & emcc --version 2>&1 | Select-Object -First 1
    Write-Host "Found Emscripten: $emccVersion" -ForegroundColor Yellow
} catch {
    Write-Host "Error: emcc not found. Please install Emscripten SDK." -ForegroundColor Red
    Write-Host ""
    Write-Host "Installation instructions:" -ForegroundColor Yellow
    Write-Host "1. git clone https://github.com/emscripten-core/emsdk.git" -ForegroundColor White
    Write-Host "2. cd emsdk" -ForegroundColor White
    Write-Host "3. .\emsdk install latest" -ForegroundColor White
    Write-Host "4. .\emsdk activate latest" -ForegroundColor White
    Write-Host "5. .\emsdk_env.ps1" -ForegroundColor White
    exit 1
}

# Build parameters
$sources = "hello.cpp"
$output = "hello.js"
$flags = @(
    "-std=c++17",
    "-O2",
    "-s", "WASM=1",
    "-s", 'EXPORTED_RUNTIME_METHODS=["ccall", "cwrap"]',
    "-s", "EXPORT_ES6=1",
    "-s", "MODULARIZE=1",
    "-s", 'EXPORT_NAME="HelloModule"',
    "-s", "ALLOW_MEMORY_GROWTH=1",
    "-s", "ENVIRONMENT='web,node'",
    "--bind"
)

# Build command
$buildCommand = @("emcc", $sources, "-o", $output) + $flags

Write-Host "Build command: $($buildCommand -join ' ')" -ForegroundColor Yellow

# Execute build
try {
    & $buildCommand[0] $buildCommand[1..($buildCommand.Length-1)]
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Build successful!" -ForegroundColor Green
        Write-Host "Generated files:" -ForegroundColor Yellow
        Get-ChildItem -Path "." -Filter "hello.*" | ForEach-Object {
            Write-Host "  - $($_.Name)" -ForegroundColor White
        }
    } else {
        Write-Host "❌ Build failed with exit code $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }
} catch {
    Write-Host "❌ Build failed: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}
