# @gaudi/physics-components

React Three Fiber components for physics simulation with WebAssembly, featuring rod dynamics and constraint systems.

## Installation

```bash
npm install @gaudi/physics-components
```

## Quick Start

### Basic Rod Constraints Simulation

```tsx
import React from 'react';
import { Canvas } from '@react-three/fiber';
import { RodConstraintsTest } from '@gaudi/physics-components';

function App() {
  return (
    <div style={{ height: '100vh' }}>
      <RodConstraintsTest />
    </div>
  );
}
```

### Using Individual Components

```tsx
import React from 'react';
import { Canvas } from '@react-three/fiber';
import { GaudiLoggerRenderer, useWasmModule } from '@gaudi/physics-components';

function MyPhysicsApp() {
  const { module, instance } = useWasmModule('rod_constraints_test');
  
  return (
    <Canvas>
      <GaudiLoggerRenderer 
        wasmInstance={instance}
        wasmModule={module}
        isPlaying={true}
        showLines={true}
        showPoints={true}
      />
    </Canvas>
  );
}
```

## Components

### RodConstraintsTest
A complete rod constraints simulation with controls and visualization.

**Props:**
- `className?: string` - Additional CSS classes
- `style?: React.CSSProperties` - Inline styles

### GaudiLoggerRenderer
Renders physics simulation data as lines and points in 3D space.

**Props:**
- `wasmInstance: WasmLoggerAPI | null` - WASM instance
- `wasmModule: WasmModule | null` - WASM module with logger functions
- `isPlaying?: boolean` - Whether simulation is running
- `showLines?: boolean` - Show line segments
- `showPoints?: boolean` - Show points
- `lineWidth?: number` - Line thickness
- `pointSize?: number` - Point size
- `lineMaterial?: THREE.LineBasicMaterialParameters` - Line material props
- `pointMaterial?: THREE.PointsMaterialParameters` - Point material props

### GaudiLoggerTest
A test component for the logger system with animated visualization.

## Utilities

### useWasmModule
Hook for loading and managing WASM modules.

```tsx
import { useWasmModule } from '@gaudi/physics-components';

function MyComponent() {
  const { loading, error, module, instance } = useWasmModule(
    'rod_constraints_test',
    (module) => new module.RodConstraintsTest()
  );
  
  if (loading) return <div>Loading...</div>;
  if (error) return <div>Error: {error}</div>;
  
  return <div>Module loaded!</div>;
}
```

## WASM Modules

The library includes pre-built WASM modules for physics simulation:

- `rod_constraints_test` - Rod dynamics with constraints
- `gaudi_logger_test` - Logger system test

### Building Custom WASM Modules

1. Navigate to the `wasm/` directory
2. Run `build.bat` to build all modules
3. Copy the generated `.js` and `.wasm` files to your project's public directory

## Development

### Building the Library

```bash
npm run build:lib
```

This creates:
- `dist/index.js` - Main library bundle
- `dist/components.js` - Components only
- `dist/utils.js` - Utilities only
- TypeScript declarations in `dist/`

### Development Server

```bash
npm run dev
```

### Building WASM Modules

```bash
npm run wasm:build
```

## Dependencies

### Peer Dependencies
- React 18+
- @react-three/fiber 8+
- @react-three/drei 9+
- Three.js 0.160+

### Internal Dependencies
- Emscripten (for WASM compilation)
- Eigen (linear algebra)
- Custom physics engine (included)

## License

MIT

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request 