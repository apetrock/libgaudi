// Core simulation data structures
export interface RodSimulationData {
  vertices: Float32Array;      // xyz positions of rod vertices
  normals: Float32Array;       // normal vectors at each vertex
  tangents: Float32Array;      // tangent vectors along the rod
  growthWeights: Float32Array; // growth weights for visualization
  frame: number;               // current simulation frame
  totalLength: number;         // total length of the rod
  vertexCount: number;         // number of vertices
}

// Simulation parameters that can be adjusted at runtime
export interface SimulationParams {
  growthRate: number;          // Rate of rod growth (default: 1.0)
  constraintStrength: number;  // Strength of constraint forces (default: 0.1)
  collisionStrength: number;   // Strength of collision response (default: 1.0)
  timeStep: number;           // Physics timestep (default: 0.05)
  sdfTransitionFrame: number; // Frame at which SDF changes (default: 400)
  bendingStiffness: number;   // Rod bending resistance (default: 0.05)
  stretchStiffness: number;   // Rod stretch resistance (default: 0.06)
}

// Props for the main RodSimulation component
export interface RodSimulationProps {
  initialParams?: Partial<SimulationParams>;
  onFrameUpdate?: (data: RodSimulationData) => void;
  debug?: boolean;            // Show debug visualization
  autoPlay?: boolean;         // Start simulation automatically
  style?: React.CSSProperties;
  className?: string;
  width?: number;
  height?: number;
}

// Props for the control panel component
export interface ControlPanelProps {
  params: SimulationParams;
  onParamsChange: (params: Partial<SimulationParams>) => void;
  isPlaying: boolean;
  onPlayPause: () => void;
  onReset: () => void;
  onStep: () => void;
}

// WebAssembly module interface - matches our Emscripten bindings
export interface WasmModule {
  RodSimulation: new () => WasmSimulation;
}

// WebAssembly simulation instance interface - matches our C++ wrapper
export interface WasmSimulation {
  step(frame: number): void;
  reset(): void;
  getVertices(): number[];        // Returns vector<float> as JS array
  getNormals(): number[];         // Returns vector<float> as JS array  
  getTangets(): number[];         // Returns vector<float> as JS array
  getGrowthWeights(): number[];   // Returns vector<float> as JS array
  getVertexCount(): number;
  getTotalLength(): number;
  getCurrentFrame(): number;
  setGrowthRate(rate: number): void;
  setConstraintStrength(strength: number): void;
  delete(): void;                 // Emscripten cleanup method
}

// Material properties for rod visualization
export interface RodMaterialProps {
  baseColor?: string;
  growthColorLow?: string;
  growthColorHigh?: string;
  metalness?: number;
  roughness?: number;
  wireframe?: boolean;
}

// SDF visualization props
export interface SdfVisualizationProps {
  frame: number;
  opacity?: number;
  visible?: boolean;
}
