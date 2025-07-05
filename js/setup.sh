#!/bin/bash

# Quick setup script for the React Rod Simulation component

echo "🚀 Setting up Rod Simulation React Component with pnpm..."

# Check if pnpm is installed
if ! command -v pnpm &> /dev/null; then
    echo "❌ pnpm is not installed. Please install it first:"
    echo "npm install -g pnpm"
    exit 1
fi

echo "✅ pnpm found"

echo "🧹 Clean previous install artifacts"
if [ -d node_modules ]; then rm -rf node_modules; fi
if [ -f package-lock.json ]; then rm package-lock.json; fi
if [ -f pnpm-lock.yaml ]; then rm pnpm-lock.yaml; fi

echo "🧹 Cleaning WASM build artifacts"
pnpm run wasm:clean

# Install dependencies
echo "📦 Installing dependencies..."
pnpm install

# Check if installation was successful
if [ $? -eq 0 ]; then
    echo "✅ Dependencies installed successfully"
else
    echo "❌ Failed to install dependencies"
    exit 1
fi

# Type check
echo "🔍 Running type check..."
pnpm run type-check

echo "🎉 Setup complete!"
echo ""
echo "Next steps:"
echo "1. Build the WASM module: cd ../../projects/wasm_rod_test && ./build_wasm.sh"
echo "2. Start development server: pnpm run dev"
echo "3. Build for production: pnpm run build:lib"
echo ""
echo "Available commands:"
echo "  pnpm run dev       - Start development server"
echo "  pnpm run build:lib - Build library for distribution"
echo "  pnpm run type-check - Run TypeScript type checking"
echo "  pnpm run lint      - Run ESLint"
