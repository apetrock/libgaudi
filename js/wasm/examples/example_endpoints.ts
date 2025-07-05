// Example WASM TypeScript endpoints
export interface ExampleWASMModule {
  hello_world(): string;
  add_numbers(a: number, b: number): number;
}

export async function loadExampleWASM(): Promise<ExampleWASMModule> {
  const module = await import('../build/hello_world.js');
  return module.default();
}
