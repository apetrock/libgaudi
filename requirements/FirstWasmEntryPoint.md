# Rod Constraints Test - WebAssembly Integration Requirements

## Project Overview

Convert the `rod_constraints_test.hpp` C++ physics simulation to run in the browser using Emscripten and visualize the results with Three.js. This will create an interactive web-based physics simulation of dynamic rod constraints with growth mechanics and collision detection.

## Current System Analysis

### Core Components

- **Rod Physics Simulation**: Dynamic rod with constraints, bending, stretching, and collision detection
- **Growth System**: Adaptive rod growth based on SDF (Signed Distance Field) interactions
- **Constraint Solver**: Projection-based constraint solver for physics integration
- **Geometry**: 3D rod represented as connected vertices with tangent vectors
- **Collision Detection**: Rod-to-rod collision handling via dynamic collision system

### Key Features to Port

1. **Rod Dynamics**: 256-vertex loop rod with physics simulation
2. **SDF Interactions**: Sphere and multi-sphere signed distance fields
3. **Growth Mechanics**: Adaptive growth based on distance to SDF objects
4. **Visual Debugging**: Geometry logging for visualization
5. **Constraint System**: Stretch, shear, bend, twist, and collision constraints

## Technical Requirements

### Emscripten Integration

- **Build System**: CMake configuration for Emscripten compilation
- **Memory Management**: Convert smart pointers to WASM-compatible memory management
- **API Design**: C++ functions exposed to JavaScript via Embind
- **Data Transfer**: Efficient vertex/normal data exchange between C++ and JS

### Required Emscripten APIs

```cpp
// Core simulation control
EMSCRIPTEN_BINDINGS(rod_simulation) {
    class_<block_test>("RodSimulation")
        .constructor()
        .function("step", &block_test::step)
        .function("getVertices", &block_test::getVertices)      // New method needed
        .function("getNormals", &block_test::getNormals)        // New method needed
        .function("getConstraintInfo", &block_test::getConstraintInfo) // New method needed
        .function("reset", &block_test::reset)                  // New method needed
        .function("setParameters", &block_test::setParameters); // New method needed
}
```

### Debug Infrastructure Requirements

- **Interface-Only Logger**: Use the existing `gaudi::logger.hpp` interface (header-only, no implementation dependencies)
- **WASM Logger Implementation**: Create `wasm_geometry_logger.cpp` in the wasm directory implementing the `gaudi::logger` interface
- **Clean Separation**: No dependencies on `GaudiGraphics::geometry_logger` or other graphics libraries
- **WASM Line Bridge**: Expose logger data directly to JavaScript via WASM bindings
- **Line Rendering**: Real-time rendering of logged lines in React Three Fiber from the WASM logger implementation
- **Router-based Navigation**: React Router for switching between test scenarios and debug modes
- **Logger API Preservation**: Maintain compatibility with existing `gaudi::logger` interface patterns

## Line Logger Architecture Analysis

### Approach 1: React Context-based Line Logger

**Implementation**: React Context Provider with useContext hook for line accumulation and rendering.

**Benefits**:

- Native React integration - follows React patterns
- Automatic component re-rendering when lines change
- Easy provider/consumer pattern - wrap app in `<LineLoggerProvider>`
- TypeScript-friendly with context typing
- Natural fit for component-based architecture
- Built-in React dev tools support

**Drawbacks**:

- Context re-renders all consuming components when lines change
- May not be optimal for high-frequency line updates (60fps+)
- More React boilerplate code required
- Requires provider hierarchy setup
- Less optimal for WASM → JS communication patterns

**Use Cases**: Best for infrequent debug line updates, component-based line accumulation, traditional React app integration.

### Approach 2: Zustand-based Line Logger

**Implementation**: Zustand store with singleton pattern, accessed via useStore hook.

**Benefits**:

- High performance - minimal re-renders with selective subscriptions
- Global singleton access without provider setup
- Better for high-frequency updates and animation loops
- Simpler WASM integration - direct store access from anywhere
- Smaller bundle size and less overhead
- Excellent for frame-by-frame line accumulation
- Can easily persist debug state across route changes

**Drawbacks**:

- Another dependency to manage
- Less "React-like" - external state management
- Requires understanding Zustand patterns
- Might be overkill for simple debug logging

**Use Cases**: Best for real-time updates, WASM communication, performance-critical line rendering, singleton debug patterns.

### Hybrid Approach: Context + Zustand

**Implementation**: React Context wraps Zustand store, providing both patterns.

**Benefits**:

- Best of both worlds - Context API familiarity with Zustand performance
- Context provides React integration, Zustand handles state management
- Easy migration path - start with Context, optimize with Zustand
- Flexible access patterns - use Context in components, direct Zustand access for WASM

**Structure**:

```typescript
// Zustand store for performance
const useLineStore = create<LineStore>((set, get) => ({ ... }));

// Context for React integration
const LineLoggerContext = createContext<LineLoggerAPI>(null);

// Provider combines both
function LineLoggerProvider({ children }) {
  const store = useLineStore();
  return (
    <LineLoggerContext.Provider value={store}>
      {children}
    </LineLoggerContext.Provider>
  );
}
```

### Recommendation

For our use case (debug visualization with future WASM integration):

1. **Phase 0**: Start with **Zustand-only** approach for simplicity and performance
2. **Phase 1**: When WASM integration comes, Zustand provides better direct access patterns
3. **Future**: Add Context wrapper if component integration becomes complex

**Rationale**: Debug line logging is inherently singleton-oriented (like Gaudi's geometry_logger), benefits from global access, and will need high-performance updates when WASM integration happens. Zustand matches this pattern better than Context.

### React Three Fiber Visualization Requirements

- **Component Architecture**: Self-contained React component for easy integration
- **Geometry Representation**: Convert rod vertices to R3F BufferGeometry
- **Real-time Updates**: Efficient vertex buffer updates using R3F hooks
- **Visual Elements**:
  - Rod mesh with smooth curves (TubeGeometry or custom BufferGeometry)
  - SDF sphere visualizations with materials
  - Constraint force vectors (optional debug mode)
  - Growth visualization (color-coded segments using vertex colors)
- **Performance**: Non-real-time physics simulation - no need for 60fps, optimized for quality over speed
- **Interactivity**: Orbital controls, parameter adjustment via React state, play/pause
- **Exportability**: Clean component interface for integration into other React projects

### Data Structures for React Interface

```typescript
// Rod simulation data
interface RodSimulationData {
  vertices: Float32Array; // xyz positions
  normals: Float32Array; // normal vectors
  tangents: Float32Array; // tangent vectors
  growthWeights: Float32Array; // for color visualization
  frame: number;
  totalLength: number;
}

// Component props interface
interface RodSimulationProps {
  initialParams?: SimParams;
  onFrameUpdate?: (data: RodSimulationData) => void;
  debug?: boolean;
  autoPlay?: boolean;
  style?: React.CSSProperties;
  className?: string;
}

// Simulation parameters
interface SimParams {
  growthRate: number;
  constraintStrength: number;
  collisionStrength: number;
  timeStep: number;
  sdfTransitionFrame: number;
}
```

## Implementation Plan

### Phase 0: Basic Rendering Validation (FIRST MILESTONE)

1. **Simple Box Render**: Create minimal React Three Fiber scene with animated rotating box
2. **React Router Setup**: Navigation between different test scenarios and debug modes
3. **Zustand Line Logger**: Singleton debug line store with accumulate/render/flush API
4. **Line Rendering Test**: Validate debug line geometry rendering with Zustand store
5. **Component Architecture**: Verify React component structure and TypeScript integration
6. **Build System Validation**: Confirm Vite development server and build process work correctly

### Phase 1: WASM Line Logger Integration

1. **Simple WASM Test Project**: Minimal C++ project using Gaudi's geometry_logger for debug lines
2. **C++ Line Logger Bridge**: Connect existing Gaudi line logger to WASM bindings
3. **WASM→Zustand Communication**: Transfer debug line data from C++ directly to Zustand store
4. **Real-time Line Rendering**: Display WASM-generated debug lines in React Three Fiber
5. **Communication Validation**: Prove the C++ geometry_logger → WASM → Zustand → R3F pipeline works

### Phase 2: Minimal Rod Simulation Integration

1. **Thin WASM Wrapper**: Create wrapper around existing `block_test` class (no C++ changes)
2. **Basic Rod Data**: Extract vertex positions from existing simulation
3. **Simple Rod Visualization**: Basic line segments or points showing rod structure
4. **Manual Stepping**: Single-step through simulation frames manually
5. **Data Validation**: Confirm rod vertex data transfers correctly

### Phase 3: Enhanced Rod Visualization

1. **Rod Mesh Generation**: Convert rod vertices to smooth tube geometry
2. **Growth Visualization**: Color-coded segments based on growth weights
3. **SDF Object Rendering**: Visualize the sphere constraints
4. **Animation Controls**: Play/pause/step controls for simulation
5. **Parameter Interface**: Basic sliders for simulation parameters

### Phase 4: Full Integration & Polish

1. **Complete Component**: Finalize the exportable React component
2. **TypeScript Definitions**: Full type safety and documentation
3. **Performance Optimization**: Memory management and rendering efficiency
4. **Error Handling**: Graceful degradation and error boundaries
5. **NPM Package**: Publishable component with examples

## Technical Challenges & Solutions

### Challenge 1: Eigen Library Compatibility

- **Issue**: Eigen may have WASM compilation issues
- **Solution**: Use Emscripten-compatible linear algebra or selective Eigen compilation

### Challenge 2: Memory Management

- **Issue**: Interfacing C++ smart pointers with JavaScript
- **Solution**: Keep all C++ code unchanged - create thin WASM wrapper that manages the smart pointer lifecycle internally and exposes simple data extraction methods

### Challenge 3: Performance

- **Issue**: Data transfer overhead between WASM and JavaScript
- **Solution**: Minimize copying, use TypedArrays, batch updates

### Challenge 4: Complex Dependencies

- **Issue**: Multiple header dependencies may complicate compilation
- **Solution**: Create isolated compilation unit with minimal dependencies

## Proposed File Structure (Refined)

```
js/                           # Main JS/React project
├── package.json               # React component package
├── tsconfig.json             # TypeScript configuration
├── vite.config.ts            # Global build configuration
├── tailwind.config.js        # Shared styling config
├── postcss.config.js         # Shared PostCSS config
│
├── src/                      # Core reusable library
│   ├── index.ts              # Main library export
│   ├── components/           # Reusable components
│   │   ├── RodSimulation.tsx # Main rod simulation component
│   │   ├── ControlPanel.tsx  # Parameter controls
│   │   ├── RodMesh.tsx       # R3F geometry component
│   │   ├── LineRenderer.tsx  # Debug line renderer
│   │   └── BoxDemo.tsx       # Basic demo component
│   ├── hooks/                # Reusable React hooks
│   │   ├── useWasmModule.ts  # WASM loading hook
│   │   ├── useSimulation.ts  # Simulation state management
│   │   └── useLineLogger.ts  # Zustand line logger hook
│   ├── stores/               # Global state management
│   │   └── lineLoggerStore.ts # Zustand store for debug lines
│   ├── types/                # TypeScript definitions
│   │   ├── simulation.ts     # Rod simulation types
│   │   ├── lineLogger.ts     # Line logger types
│   │   └── wasm.ts           # WASM module interfaces
│   ├── utils/                # General utilities (including WASM loading)
│   │   ├── wasmLoader.ts     # WASM module loading utilities
│   │   ├── math.ts           # Math utilities
│   │   └── constants.ts      # App constants
│   └── styles/               # Shared CSS/styling
│       └── globals.css       # Global styles
│
├── Projects/                 # Demo applications and examples (flat structure)
│   ├── phase0_box/           # Phase 0 validation example
│   │   ├── BoxExample.tsx
│   │   ├── assets/
│   │   └── README.md
│   ├── phase1_logger/        # Phase 1 logger integration
│   │   ├── LoggerExample.tsx
│   │   ├── assets/
│   │   └── README.md
│   ├── phase2_rod/           # Phase 2 rod simulation
│   │   ├── RodExample.tsx
│   │   ├── assets/
│   │   └── README.md
│   └── wasm_hello/           # Basic WASM hello world
│       ├── HelloWorldDemo.tsx
│       ├── assets/
│       └── README.md
│
├── wasm/                     # Complete WASM ecosystem
│   ├── CMakeLists.txt        # Root WASM build config (orchestrates all projects)
│   ├── build.bat/.ps1        # Build scripts for all projects
│   ├── build/                # All build outputs (generated)
│   │   ├── hello_world.js/.wasm
│   │   ├── gaudi_logger_test.js/.wasm
│   │   └── rod_simulation_main.js/.wasm
│   │
│   ├── examples/             # Example WASM project
│   │   ├── CMakeLists.txt    # Example-specific build config
│   │   ├── src/              # C++ source for examples
│   │   │   ├── hello_world.cpp
│   │   │   └── basic_math.cpp
│   │   └── example_endpoints.ts  # TypeScript API for examples
│   │
│   ├── logger/               # Logger WASM project
│   │   ├── CMakeLists.txt    # Logger-specific build config
│   │   ├── src/              # C++ source for logger
│   │   │   ├── wasm_geometry_logger.cpp
│   │   │   └── logger_test.cpp
│   │   └── logger_endpoints.ts   # TypeScript API for logger
│   │
│   └── rod_simulation/       # Rod simulation WASM project
│       ├── CMakeLists.txt    # Rod simulation build config
│       ├── src/              # C++ source for rod simulation
│       │   ├── rod_wrapper.cpp
│       │   └── simulation_main.cpp
│       └── rod_endpoints.ts  # TypeScript API for rod simulation
│
└── dist/                     # Built library package
    ├── index.js              # Built library entry
    ├── index.d.ts            # TypeScript declarations
    └── components/           # Built components
```

### **Key Organizational Principles:**

1. **Separation of Concerns**:
   - `src/` = Reusable library code
   - `Projects/` = Applications and examples
   - `wasm/` = C++ build system
   - `demo/` = Generated assets

2. **Shared vs Project-Specific**:
   - Common configs at root level
   - Project-specific configs in project folders
   - Shared components in `src/components/`
   - Project-specific components in `Projects/*/`

3. **WASM Integration**:
   - WASM utilities in `src/utils/` (loading, error handling)
   - WASM projects in `wasm/{project_name}/` (C++ + TypeScript co-located)
   - WASM builds in `wasm/build/` (all generated modules)
   - Projects import WASM from `wasm/{project_name}/{project_name}_endpoints`

4. **Examples Organization**:
   - Each phase gets its own example folder
   - Clear progression from Phase 0 → Phase 2
   - Integration tests separate from examples

### **Structure Decision Points:**

**Q1: Should we keep the current demo structure or move to Projects?**
- Current: Root-level `App.tsx`, `main.tsx` for development
- Proposed: Move to `Projects/demo_app/` 
- **Decision needed**: Keep development at root or move to Projects?

**Q2: WASM file imports - relative or absolute?**
- Current: Components import `'../../demo/wasm/module.js'`
- Alternative: Use Vite alias like `'@wasm/module.js'`
- **Decision needed**: Import strategy?

**Q3: Shared configuration**
- Root `vite.config.ts` with multiple entry points?
- Or individual configs per project?
- **Decision needed**: Configuration sharing strategy?

**Q4: Component categorization**
- Which components are "library" vs "example"?
- Current mix in `src/components/` needs sorting
- **Decision needed**: What stays in core library?

**Q5: Build outputs**
- Should Projects build to their own dist folders?
- Or everything builds to root dist with different entry points?
- **Decision needed**: Build output strategy?

### **Immediate Action Items:**

1. **Decide on development workflow**: Keep current root-level development or move to Projects/demo_app?
2. **Component audit**: Which components belong in core library vs examples?
3. **Configuration strategy**: Shared configs or project-specific?
4. **Import paths**: How should projects import WASM modules and library components?

## Success Criteria

1. **Phase 0 Milestone**: Smooth animated box with React Router navigation and working line logger in browser
2. **Phase 1 Milestone**: C++ geometry_logger successfully sends debug lines to JS for rendering
3. **Phase 1.5 Milestone**: Line Logger-Three.js Integration Test - Real-time C++ line generation with JavaScript rendering pipeline
4. **Functional Simulation**: Rod physics running at stable framerate in browser
5. **Clean Component Interface**: Easy to integrate `<RodSimulation />` component
6. **Visual Quality**: Smooth, appealing rod visualization with growth animation
7. **Interactive Controls**: Real-time parameter adjustment via React state
8. **Performance**: Stable framerate (no 60fps requirement) with React concurrent features
9. **TypeScript Support**: Full type safety and IntelliSense
10. **Reusability**: Easily importable into other React projects
11. **Documentation**: Clear usage examples and API documentation

## Phase 1.5 Milestone - Line Logger Integration Test

**Objective**: Demonstrate real-time C++ to Three.js line rendering pipeline using a clean WASM logger implementation.

**Technical Requirements**:

- **Interface-Only Logger**: Use `gaudi::logger.hpp` interface with WASM-specific implementation
- **WASM Logger Implementation**: Create `wasm_geometry_logger.cpp` implementing the logger interface in the wasm directory
- **Frame-based Line Generation**: C++ generates animated lines every frame via logger interface calls (rotating axes, sine waves)
- **WASM Data Bridge**: Expose logger data directly to JavaScript via WASM bindings
- **Three.js Real-time Renderer**: BufferGeometry with dynamic position/color attributes from WASM logger data
- **Clean Architecture**: No graphics library dependencies, pure logger interface implementation

**Success Criteria**:

```typescript
// Integration test using gaudi::logger interface
interface GaudiLoggerWASMBridge {
  // Existing logger interface usage (unchanged C++ patterns)
  // gaudi::logger::line(start, end, color);

  // WASM bridge methods
  get_logged_lines_count(): number;
  get_logged_lines_positions(): number; // Float32Array view of logger data
  get_logged_lines_colors(): number; // Float32Array view of color data
  clear_logged_lines(): void; // Call gaudi::logger::clear()

  // Animation test methods
  animate_axes(time: number): void; // Generate rotating coordinate axes
  animate_sine_wave(time: number): void; // Generate animated sine wave
}
```

**Demo Features**:

- **Clean Logger Interface Usage** with `gaudi::logger::line()`, `gaudi::logger::clear()`, etc.
- **Animated coordinate axes** generated via logger interface calls
- **Sine wave animation** logged through the interface
- **Live statistics** showing logged line count from WASM logger state
- **Interactive controls** that call logger methods (clear, add test patterns)
- **React Router endpoint** `/logger-test` for testing the logger integration

This milestone validates the core infrastructure needed for rod physics visualization before implementing the full simulation.

## Phase 0 Milestone - Basic Box Demo

```tsx
import { BoxDemo } from "rod-simulation-component";

function App() {
  return (
    <div style={{ width: "100%", height: "100vh" }}>
      <BoxDemo color="#4a90e2" rotationSpeed={1.0} showControls={true} />
    </div>
  );
}
```

## Component Usage Example

```tsx
import { RodSimulation } from "rod-simulation-component";

function MyApp() {
  const [params, setParams] = useState({
    growthRate: 1.0,
    constraintStrength: 0.1,
    // ... other params
  });

  return (
    <div style={{ width: "100%", height: "100vh" }}>
      <RodSimulation
        initialParams={params}
        onFrameUpdate={(data) => console.log("Frame:", data.frame)}
        debug={false}
        autoPlay={true}
      />
    </div>
  );
}
```

## Future Extensions

- Multi-rod interactions within the same component
- Custom SDF shape designer as separate component
- Simulation recording/playback capabilities
- WebXR support for VR/AR visualization
- Real-time collaborative simulation editing
- NPM registry publication for broader ecosystem use
- Integration with popular React UI libraries (MUI, Chakra, etc.)

---

_This document will be updated as implementation progresses and requirements are refined._
