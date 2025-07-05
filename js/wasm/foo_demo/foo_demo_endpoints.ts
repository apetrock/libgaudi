/**
 * TypeScript interface for foo_demo WASM module
 * Generated automatically - do not edit manually
 */

export interface Foo_demoModule {
  add(a: number, b: number): number;
  greet(name: string): string;
}

/**
 * A simple foo_demo WASM demo
 */
export const foo_demoEndpoints = {
  moduleName: 'foo_demo',
  functions: [
    { name: 'add', returnType: 'number', parameters: [{ name: 'a', type: 'number' }, { name: 'b', type: 'number' }] },
    { name: 'greet', returnType: 'string', parameters: [{ name: 'name', type: 'string' }] }
  ]
} as const;
