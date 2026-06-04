// TypeScript endpoints for unit_tests WASM module

export interface UnitTestsModule {
  runAllTests(): string;
  runCoreTests(): string;
  setSphereObj(objData: string): void;
}

// Helper function to load and run all tests
export async function loadAndRunTests(
  wasmLoader: () => Promise<UnitTestsModule>,
  sphereObjData?: string
): Promise<string> {
  const module = await wasmLoader();
  
  // If mesh data provided, set it for BVH tests
  if (sphereObjData) {
    module.setSphereObj(sphereObjData);
  }
  
  return module.runAllTests();
}

// Run only core tests (no mesh data required)
export async function loadAndRunCoreTests(
  wasmLoader: () => Promise<UnitTestsModule>
): Promise<string> {
  const module = await wasmLoader();
  return module.runCoreTests();
}
