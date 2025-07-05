# WASM Demo Workflow - Formalized Process

## Overview

This document establishes a **consistent, predictable process** for creating WASM demos. Follow this process exactly to avoid import/MIME type issues.

## Directory Structure (Fixed)

```
js/
├── public/wasm/                    # WASM files served by Vite (DO NOT CHANGE)
│   ├── hello_world.js
│   ├── hello_world.wasm
│   ├── gaudi_logger_test.js
│   ├── gaudi_logger_test.wasm
│   ├── rod_constraints_test.js
│   └── rod_constraints_test.wasm
│
├── src/components/                 # React components (DO NOT CHANGE)
│   ├── WasmHelloWorld.tsx
│   ├── GaudiLoggerTest.tsx
│   ├── RodConstraintsTest.tsx
│   └── LineLoggerIntegrationTest.tsx
│
├── wasm/                          # C++ source and build system
│   ├── CMakeLists.txt            # Root build config
│   ├── build.bat                 # Build all WASM modules
│   ├── build/                    # Build outputs (generated)
│   │   ├── examples/
│   │   ├── logger/
│   │   └── rod_constraints/
│   │
│   ├── examples/                 # Example WASM project
│   │   ├── CMakeLists.txt
│   │   ├── src/
│   │   │   └── hello_world.cpp
│   │   └── example_endpoints.ts
│   │
│   ├── logger/                   # Logger WASM project
│   │   ├── CMakeLists.txt
│   │   ├── src/
│   │   │   └── gaudi_logger_test.cpp
│   │   └── logger_endpoints.ts
│   │
│   └── rod_constraints/          # Rod constraints WASM project
│       ├── CMakeLists.txt
│       ├── src/
│       │   └── rod_constraints_wasm.cpp
│       └── rod_endpoints.ts
│
└── demo/wasm/                    # DEPRECATED - files copied to public/wasm
```

## Step-by-Step Process for New WASM Demo

### Step 1: Create C++ Source
1. Create directory: `js/wasm/{project_name}/`
2. Create `js/wasm/{project_name}/CMakeLists.txt`
3. Create `js/wasm/{project_name}/src/{project_name}.cpp`
4. Create `js/wasm/{project_name}/{project_name}_endpoints.ts`

### Step 2: Add to Build System
1. Add project to `js/wasm/CMakeLists.txt`
2. Add copy command to `js/wasm/build.bat`

### Step 3: Build WASM
```bash
cd js/wasm
./build.bat
```

### Step 4: Create React Component
1. Create `js/src/components/{ProjectName}.tsx`
2. Import using: `import('@wasmbuilds/{project_name}.js')`
3. Add to `js/src/App.tsx` routing

### Step 5: Test
1. Start dev server: `cd js && npm run dev`
2. Navigate to `http://localhost:5173/{project-name}`

## Import Paths (Fixed)

### In React Components
```typescript
// ALWAYS use this pattern:
const wasmModule = await import('@wasmbuilds/{project_name}.js');
```

### Vite Configuration (Fixed)
```typescript
// vite.config.ts - DO NOT CHANGE
resolve: {
  alias: {
    '@wasmbuilds': path.resolve(__dirname, './public/wasm'),
  },
},
```

## Build Process (Fixed)

### WASM Build Script
```batch
# js/wasm/build.bat - DO NOT CHANGE
@echo off
REM Build all WASM modules and copy to public/wasm/

echo Building WebAssembly modules...

REM Build with CMake
call emcmake cmake -B build -S .
call cmake --build build

REM Copy to public/wasm/ (Vite serves from here)
echo Copying to public/wasm/...
copy build\examples\hello_world.* public\wasm\
copy build\logger\gaudi_logger_test.* public\wasm\
copy build\rod_constraints\rod_constraints_test.* public\wasm\

echo Build complete!
```

### CMake Configuration
```cmake
# js/wasm/CMakeLists.txt - Add new projects here
add_subdirectory(examples)
add_subdirectory(logger)
add_subdirectory(rod_constraints)
# add_subdirectory(new_project)  # Add new projects here
```

## Common Issues & Solutions

### Issue: "Cannot find module '@wasmbuilds/...'"
**Solution**: Run `cd js/wasm && ./build.bat` to copy files to `public/wasm/`

### Issue: "MIME type error for .wasm"
**Solution**: Files must be in `public/wasm/` for Vite to serve with correct MIME type

### Issue: "WebAssembly.instantiate() failed"
**Solution**: Check that WASM file exists in `public/wasm/` and is not corrupted

### Issue: "HelloModule.default is not a function"
**Solution**: Ensure C++ code exports default function via Embind

## Validation Checklist

Before declaring a demo "working":

- [ ] WASM files exist in `public/wasm/`
- [ ] React component imports from `@wasmbuilds/`
- [ ] Dev server starts without errors
- [ ] Component loads without console errors
- [ ] WASM module initializes successfully
- [ ] Demo functionality works as expected

## Example: Creating "MyDemo"

### 1. Create C++ Source
```cpp
// js/wasm/mydemo/src/mydemo.cpp
#include <emscripten/bind.h>

using namespace emscripten;

int add(int a, int b) {
    return a + b;
}

EMSCRIPTEN_BINDINGS(mydemo) {
    function("add", &add);
}
```

### 2. Create TypeScript Endpoints
```typescript
// js/wasm/mydemo/mydemo_endpoints.ts
export interface MyDemoModule {
  add(a: number, b: number): number;
}
```

### 3. Create React Component
```typescript
// js/src/components/MyDemo.tsx
import { useState, useEffect } from 'react';

interface MyDemoModule {
  add(a: number, b: number): number;
}

export function MyDemo() {
  const [module, setModule] = useState<MyDemoModule | null>(null);
  const [result, setResult] = useState<number | null>(null);

  useEffect(() => {
    const loadWasm = async () => {
      const wasmModule = await import('@wasmbuilds/mydemo.js');
      const instance = await wasmModule.default();
      setModule(instance);
    };
    loadWasm();
  }, []);

  const testAdd = () => {
    if (module) {
      setResult(module.add(5, 3));
    }
  };

  return (
    <div>
      <h1>My Demo</h1>
      <button onClick={testAdd}>Test Add</button>
      {result !== null && <p>Result: {result}</p>}
    </div>
  );
}
```

### 4. Add to Build System
```cmake
# js/wasm/CMakeLists.txt
add_subdirectory(mydemo)
```

```batch
# js/wasm/build.bat
copy build\mydemo\mydemo.* public\wasm\
```

### 5. Add to Routing
```typescript
// js/src/App.tsx
import { MyDemo } from './components/MyDemo';

// In router:
<Route path="/mydemo" element={<MyDemo />} />
```

## Golden Rules

1. **NEVER change import paths** - Always use `@wasmbuilds/`
2. **NEVER move WASM files manually** - Always use build script
3. **ALWAYS put new WASM files in public/wasm/** - Vite serves from here
4. **ALWAYS follow the directory structure** - No exceptions
5. **ALWAYS run build.bat after C++ changes** - No manual copying

This process ensures consistency and prevents the "shifting ground" problem. 