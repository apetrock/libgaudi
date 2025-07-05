# WebAssembly Hello World - Phase 1

This directory contains the basic WebAssembly "Hello World" integration for the Rod Simulation Component project.

## Phase 1 Goals

- ✅ Basic C++ to JavaScript communication
- ✅ Function calls with different parameter types
- ✅ Class instantiation and method calls
- ✅ Interactive React component integration
- ⏳ Ready for Phase 2: Gaudi geometry_logger integration

## Prerequisites

You need the Emscripten SDK installed to build the WebAssembly module.

### Installing Emscripten (Windows)

1. Clone the Emscripten SDK:
   ```powershell
   git clone https://github.com/emscripten-core/emsdk.git
   cd emsdk
   ```

2. Install and activate the latest version:
   ```powershell
   .\emsdk install latest
   .\emsdk activate latest
   ```

3. Set up the environment:
   ```powershell
   .\emsdk_env.ps1
   ```

### Installing Emscripten (Linux/Mac)

1. Clone the Emscripten SDK:
   ```bash
   git clone https://github.com/emscripten-core/emsdk.git
   cd emsdk
   ```

2. Install and activate the latest version:
   ```bash
   ./emsdk install latest
   ./emsdk activate latest
   ```

3. Set up the environment:
   ```bash
   source ./emsdk_env.sh
   ```

## Building

### Windows (PowerShell)
```powershell
.\build.ps1
```

### Linux/Mac (Make)
```bash
make
```

### Debug Build
```bash
make debug
```

## Generated Files

After building, you'll have:
- `hello.js` - JavaScript module loader
- `hello.wasm` - WebAssembly binary

## Usage

The WASM module is automatically loaded by the React component when you navigate to `/wasm-hello` in the demo application.

## Functions Available

### C Functions (exported via `extern "C"`)
- `hello_world()` - Returns a greeting string
- `add_numbers(a, b)` - Adds two integers
- `multiply_floats(x, y)` - Multiplies two floating-point numbers

### C++ Classes (exported via Embind)
- `MathHelper` - A class that demonstrates object-oriented WASM usage
  - Constructor: `new MathHelper(initialValue)`
  - `getValue()` - Get current value
  - `setValue(value)` - Set new value
  - `calculate(input)` - Perform calculation: `value * input * 2`
  - `getStatus()` - Get status string

## Next Steps (Phase 2)

1. Integrate with Gaudi's `geometry_logger`
2. Stream debug lines from C++ to JavaScript
3. Render geometry in React Three Fiber
4. Performance optimization for real-time updates

## Troubleshooting

### WASM module fails to load
- Ensure Emscripten is properly installed
- Check that the build completed successfully
- Verify the generated files exist in the `wasm/` directory
- Check browser console for detailed error messages

### Build fails
- Ensure Emscripten environment is activated
- Check that all paths are correct
- Try a clean build: `make clean && make`
