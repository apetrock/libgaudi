/**
 * TypeScript interface for bar_demo WASM module
 * Generated automatically - do not edit manually
 */

export interface Bar_demoModule {
  add(a: number, b: number): number;
  greet(name: string): string;
}

/**
 * A simple bar_demo WASM demo
 */
export const bar_demoEndpoints = {
  moduleName: 'bar_demo',
  functions: [
    { name: 'add', returnType: 'number', parameters: [{ name: 'a', type: 'number' }, { name: 'b', type: 'number' }] },
    { name: 'greet', returnType: 'string', parameters: [{ name: 'name', type: 'string' }] }
  ]
} as const;
