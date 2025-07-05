import { useState, useEffect } from 'react';
import { Bar_demoModule } from '../../wasm/bar_demo/bar_demo_endpoints';
import { loadWasmModule } from '../utils/wasmLoader';

/**
 * A simple bar_demo WASM demo
 * Generated automatically - do not edit manually
 */
export function Bar_demo() {
  const [module, setModule] = useState<Bar_demoModule | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);

  useEffect(() => {
    const loadWasm = async () => {
      try {
        setLoading(true);
        const instance = await loadWasmModule('bar_demo') as Bar_demoModule;
        setModule(instance);
        setLoading(false);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load WASM');
        setLoading(false);
      }
    };
    loadWasm();
  }, []);

  const testAdd = () => {
    if (module) {
      const result = module.add(5, 3);
      setResult(`add(5, 3) = ${result}`);
    }
  };

  const testGreet = () => {
    if (module) {
      const result = module.greet('World');
      setResult(`greet('World') = ${result}`);
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500 mx-auto mb-4"></div>
          <p>Loading bar_demo WASM module...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center text-red-600">
          <h3 className="text-xl font-bold mb-2">❌ Load Failed</h3>
          <p>{error}</p>
          <p className="text-sm mt-2">Make sure to run: cd wasm && ./build.bat</p>
        </div>
      </div>
    );
  }

  return (
    <div className="p-8 max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold mb-2">Bar_demo</h1>
        <p className="text-gray-600">A simple bar_demo WASM demo</p>
      </div>

      <div className="bg-gray-50 p-6 rounded-lg mb-6">
        <h2 className="text-xl font-semibold mb-4">Test Functions</h2>
        <div className="space-y-2">
          <button onClick={testAdd} className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600">
            Test add(5, 3)
          </button>
          <button onClick={testGreet} className="px-4 py-2 bg-green-500 text-white rounded hover:bg-green-600 ml-2">
            Test greet('World')
          </button>
        </div>
      </div>

      {result && (
        <div className="bg-green-50 border border-green-200 p-4 rounded-lg">
          <h3 className="font-semibold text-green-800 mb-2">Result:</h3>
          <p className="text-green-700 font-mono">{result}</p>
        </div>
      )}

      <div className="mt-8 p-4 bg-blue-50 rounded-lg">
        <h3 className="font-semibold mb-2">Generated Files:</h3>
        <ul className="text-sm space-y-1">
          <li>• <code>wasm/bar_demo/src/bar_demo.cpp</code> - C++ source</li>
          <li>• <code>wasm/bar_demo/bar_demo_endpoints.ts</code> - TypeScript interface</li>
          <li>• <code>src/components/Bar_demo.tsx</code> - React component</li>
        </ul>
      </div>
    </div>
  );
}
