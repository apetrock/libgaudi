@echo off
REM Quick setup script for the React Rod Simulation component

echo 🚀 Setting up Rod Simulation React Component with pnpm...

REM Check if pnpm is installed
pnpm --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ pnpm is not installed. Please install it first:
    echo npm install -g pnpm
    pause
    exit /b 1
)

echo ✅ pnpm found

REM 🧹 Clean previous install artifacts
if exist node_modules rmdir /s /q node_modules
if exist package-lock.json del package-lock.json
if exist pnpm-lock.yaml del pnpm-lock.yaml

echo 🧹 Cleaning WASM build artifacts
pnpm run wasm:clean

REM Install dependencies
echo 📦 Installing dependencies...
pnpm install

REM Check if installation was successful
if %errorlevel% neq 0 (
    echo ❌ Failed to install dependencies
    pause
    exit /b 1
)

echo ✅ Dependencies installed successfully

REM Type check
echo 🔍 Running type check...
pnpm run type-check

echo 🎉 Setup complete!
echo.
echo Next steps:
echo 1. Build the WASM module: cd ..\..\projects\wasm_rod_test ^&^& build_wasm.bat
echo 2. Start development server: pnpm run dev
echo 3. Build for production: pnpm run build:lib
echo.
echo Available commands:
echo   pnpm run dev       - Start development server
echo   pnpm run build:lib - Build library for distribution  
echo   pnpm run type-check - Run TypeScript type checking
echo   pnpm run lint      - Run ESLint

pause
