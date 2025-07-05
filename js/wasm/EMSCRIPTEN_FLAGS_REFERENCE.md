# Emscripten Compiler Flags Reference

This document provides a comprehensive reference for Emscripten compiler flags, organized by category. It includes both current valid flags and deprecated flags to avoid.

## Core WASM Settings

### Basic WASM Configuration
```bash
-sWASM=1                    # Enable WebAssembly (default: 1)
-sSTANDALONE_WASM=0         # Generate standalone WASM (default: 0)
-sSINGLE_FILE=1             # Embed WASM binary in JS file (default: 0)
-sALLOW_MEMORY_GROWTH=1     # Allow dynamic memory growth (default: 0)
-sINITIAL_MEMORY=16777216   # Initial memory size in bytes (default: calculated)
-sMAXIMUM_MEMORY=2147483648 # Maximum memory size in bytes (default: 2GB)
```

### Memory Management
```bash
-sMALLOC=dlmalloc           # Malloc implementation: dlmalloc, emmalloc, mimalloc (default: dlmalloc)
-sABORTING_MALLOC=1         # Abort on malloc failure (default: 1)
-sINITIAL_HEAP=16777216     # Initial heap size (default: 16MB)
-sMEMORY_GROWTH_GEOMETRIC_STEP=0.20  # Memory growth rate (default: 0.20)
-sMEMORY_GROWTH_GEOMETRIC_CAP=100663296  # Max geometric growth (default: 96MB)
```

## Debugging and Development

### Debug Flags (Current)
```bash
-gsource-map                # Generate source maps (replaces deprecated -g4)
-sASSERTIONS=1              # Enable runtime assertions (default: 1)
-sSAFE_HEAP=1               # Check heap access (default: 0)
-sDEMANGLE_SUPPORT=1        # Enable symbol demangling (default: 0)
-sSTACK_OVERFLOW_CHECK=2    # Stack overflow detection: 0=none, 1=cookie, 2=binaryen (default: 0)
-sRETAIN_COMPILER_SETTINGS=1 # Retain compiler settings in output (default: 0)
-sEXCEPTION_STACK_TRACES=1  # Include stack traces in exceptions (default: 0)
```

### Debug Flags (Deprecated - DO NOT USE)
```bash
-g4                         # DEPRECATED: Use -gsource-map instead
-sABORT_ON_STACK_OVERFLOW=1 # DEPRECATED: No longer a valid flag
```

## Module System and Exports

### ES6 Module Support
```bash
-sMODULARIZE=1              # Wrap output in factory function (default: 0)
-sEXPORT_ES6=1              # Use ES6 module exports (default: 0)
-sEXPORT_NAME=Module        # Name of exported factory function (default: Module)
-sUSE_ES6_IMPORT_META=0     # Use ES6 import.meta (default: 0)
```

### Function Exports
```bash
-sEXPORTED_FUNCTIONS=_main,_myFunction  # Explicitly export functions
-sEXPORT_ALL=1              # Export all symbols (default: 0)
-sEXPORT_KEEPALIVE=1        # Export kept-alive symbols (default: 1)
-sEXPORTED_RUNTIME_METHODS=ccall,cwrap  # Export runtime methods
```

## Environment and Runtime

### Environment Configuration
```bash
-sENVIRONMENT=web,worker,node  # Target environments (default: web,webview,worker,node)
-sMINIMAL_RUNTIME=0         # Use minimal runtime (default: 0)
-sDYNAMIC_EXECUTION=1       # Allow eval() and dynamic code (default: 1)
-sINVOKE_RUN=1              # Automatically run main() (default: 1)
-sEXIT_RUNTIME=0            # Exit runtime after main() (default: 0)
```

### Browser Integration
```bash
-sUSE_WEBGL2=0              # DEPRECATED: Use -sMAX_WEBGL_VERSION=2
-sMIN_WEBGL_VERSION=1       # Minimum WebGL version (default: 1)
-sMAX_WEBGL_VERSION=2       # Maximum WebGL version (default: 1)
-sOFFSCREENCANVAS_SUPPORT=0 # Enable OffscreenCanvas support (default: 0)
```

## Performance and Optimization

### Optimization Levels
```bash
-O0                        # No optimization (debug builds)
-O1                        # Basic optimization
-O2                        # Full optimization (default for release)
-Os                        # Size optimization
-Oz                        # Aggressive size optimization
```

### Binaryen Optimizations
```bash
-sBINARYEN_IGNORE_IMPLICIT_TRAPS=0  # Ignore implicit traps in optimization (default: 0)
-sBINARYEN_EXTRA_PASSES=            # Extra optimization passes
-sWASM_ASYNC_COMPILATION=1          # Async WASM compilation (default: 1)
```

## Exception Handling

### Exception Configuration
```bash
-sDISABLE_EXCEPTION_CATCHING=1      # Disable exception catching (default: 1)
-sEXCEPTION_CATCHING_ALLOWED=[]     # Allow exceptions in specific functions
-sWASM_LEGACY_EXCEPTIONS=1          # Use legacy exception handling (default: 1)
-sEXPORT_EXCEPTION_HANDLING_HELPERS=0 # Export exception helpers (default: 0)
```

## Threading and Workers

### Pthread Support
```bash
-sUSE_PTHREADS=1            # Enable pthread support (default: 0)
-sPTHREAD_POOL_SIZE=4       # Pre-created thread pool size (default: 0)
-sPTHREAD_POOL_SIZE_STRICT=1 # Strict thread pool enforcement (default: 1)
-sALLOW_BLOCKING_ON_MAIN_THREAD=1  # Allow blocking on main thread (default: 1)
```

### WASM Workers
```bash
-sWASM_WORKERS=0            # Enable WASM Workers (default: 0)
-sSHARED_MEMORY=0           # Enable shared memory (default: 0)
```

## File System and I/O

### File System Configuration
```bash
-sFILESYSTEM=1              # Enable file system (default: 1)
-sFORCE_FILESYSTEM=0        # Force file system inclusion (default: 0)
-sNODERAWFS=0               # Use Node.js raw file system (default: 0)
-sCASE_INSENSITIVE_FS=0     # Case-insensitive file system (default: 0)
```

## Network and Fetch

### Fetch API
```bash
-sFETCH=0                   # Enable fetch API (default: 0)
-sFETCH_SUPPORT_INDEXEDDB=1 # Use IndexedDB for fetch (default: 1)
-sFETCH_DEBUG=0             # Enable fetch debugging (default: 0)
```

## Security and CSP

### Content Security Policy
```bash
-sDYNAMIC_EXECUTION=0       # Disable eval() for CSP compliance (default: 1)
-sTRUSTED_TYPES=0           # Enable Trusted Types support (default: 0)
```

## Legacy and Compatibility

### Legacy Support
```bash
-sLEGACY_VM_SUPPORT=0       # Enable legacy VM compatibility (default: 0)
-sPOLYFILL=1                # Include polyfills (default: 1)
-sLEGACY_RUNTIME=0          # Include legacy runtime symbols (default: 0)
```

## Common Flag Combinations

### Debug Build
```bash
-O0 -gsource-map -sASSERTIONS=1 -sSAFE_HEAP=1 -sDEMANGLE_SUPPORT=1 -sSTACK_OVERFLOW_CHECK=2
```

### Release Build
```bash
-O2 -sASSERTIONS=0 -sDEMANGLE_SUPPORT=1
```

### ES6 Module with WASM
```bash
-sMODULARIZE=1 -sEXPORT_ES6=1 -sUSE_ES6_IMPORT_META=0 -sSINGLE_FILE=1
```

### Threaded Application
```bash
-sUSE_PTHREADS=1 -sPTHREAD_POOL_SIZE=4 -sSHARED_MEMORY=1
```

## Simplified Module Loading (Recommended for Easy Integration)

### Basic Configuration (No Modularization)
```bash
-sWASM=1                    # Enable WebAssembly
-sSINGLE_FILE=1             # Embed WASM binary in JS file (single file output)
-sALLOW_MEMORY_GROWTH=1     # Allow dynamic memory growth
-sENVIRONMENT=web           # Target web environment
--bind                      # Enable Embind for C++/JS interop
```

### Usage Pattern
```javascript
// Simple loading - no factory functions needed
import '/wasm/my_module.js';  // Module automatically initializes
// Module is available as global 'Module' object
console.log(Module._myFunction());
```

### Advantages
- **Simpler loading**: No async factory functions
- **Direct access**: Functions available immediately after import
- **Less boilerplate**: No need to await module initialization
- **Compatible**: Works with standard script tags and ES6 imports

### Trade-offs
- **Global scope**: Module pollutes global namespace
- **Single instance**: Only one instance per page
- **Less isolation**: Multiple modules can conflict

## Notes

- **Flag Format**: Use `-sFLAG_NAME=value` for settings, `-FLAG` for compiler options
- **Boolean Values**: Use `1` for true, `0` for false
- **Lists**: Use comma-separated values without spaces
- **Deprecated Flags**: Always check the latest Emscripten documentation for current flags
- **Version Compatibility**: Some flags may not be available in older Emscripten versions

## Resources

- [Emscripten Compiler Settings Documentation](https://emscripten.org/docs/tools_reference/emcc.html#emcc-compiler-settings)
- [Emscripten Porting Guide](https://emscripten.org/docs/porting/index.html)
- [WebAssembly Documentation](https://webassembly.org/) 