// TypeScript endpoints for unit_tests WASM module

export interface UnitTestsModule {
  runAllTests(): string;
  runCoreTests(): string;
}

// Helper function to load and run all tests
export async function loadAndRunTests(
  wasmLoader: () => Promise<UnitTestsModule>
): Promise<string> {
  const module = await wasmLoader();
  return module.runAllTests();
}

// Run only core tests (no mesh data required)
export async function loadAndRunCoreTests(
  wasmLoader: () => Promise<UnitTestsModule>
): Promise<string> {
  const module = await wasmLoader();
  return module.runCoreTests();
}
