// TypeScript endpoints for BVH test module

export interface BvhTestModule {
  BvhTest: {
    new (): BvhTestInstance;
  };
}

export interface BvhTestInstance {
  // Load mesh from OBJ string
  loadMeshFromString(objContent: string): boolean;
  
  // Get mesh info
  getVertexCount(): number;
  getEdgeCount(): number;
  getTriangleCount(): number;
  
  // Build BVH structures
  buildEdgeBvh(): boolean;
  buildTriBvh(): boolean;
  
  // Run comparison tests
  testPointToEdge(tolerance: number): boolean;
  testPointToTri(tolerance: number): boolean;
  testEdgeToEdge(tolerance: number): boolean;
  
  // Get test results
  getLastBvhResult(): number;
  getLastBruteResult(): number;
  getLastTestPassed(): boolean;

  // Test harness
  setSphereObjData(objContent: string): void;
  runAllTests(): boolean;
  getLastSuiteTotal(): number;
  getLastSuiteFailed(): number;
  
  // Visualization
  visualizeResults(): void;
  clearVisualization(): void;
  
  // Logger API for GaudiLoggerRenderer
  get_line_count(): number;
  get_point_count_logger(): number;
  get_line(index: number, start: Float64Array, end: Float64Array, color: Float64Array): void;
  get_point_logger(index: number, position: Float64Array, color: Float64Array): void;
}

// Helper to load the WASM module
export async function loadBvhTestModule(): Promise<BvhTestModule> {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const factory = (await import('../../../public/wasm/bvh_test.js' as any)).default;
  return await factory();
}



