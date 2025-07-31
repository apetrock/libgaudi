import React, { useState } from "react";
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { MortonTreeTestModule, MortonTreeTestClass } from '../../wasm/modules/morton_tree/morton_tree_endpoints';
import { loadWasmModule } from '../utils/wasmLoader';
import { GaudiLoggerRenderer, WasmLoggerAPI, WasmModule as LoggerWasmModule } from './GaudiLoggerRenderer';

interface MortonTreeTestModuleWithLogger extends MortonTreeTestModule, LoggerWasmModule {}
interface MortonTreeTestClassWithLogger extends MortonTreeTestClass, WasmLoggerAPI {}

export const MortonTreeTest = () => {
  const [module, setModule] = useState<MortonTreeTestModuleWithLogger | null>(null);
  const [instance, setInstance] = useState<MortonTreeTestClassWithLogger | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [results, setResults] = useState<string[]>([]);
  const [pointCount, setPointCount] = useState(50);

  const loadWasm = async () => {
    try {
      setLoading(true);
      setError(null);
      const wasmModule = await loadWasmModule('morton_tree_test') as MortonTreeTestModuleWithLogger;
      setModule(wasmModule);
      
      const testInstance = new wasmModule.MortonTreeTest() as MortonTreeTestClassWithLogger;
      setInstance(testInstance);
      
      // Enhanced error handling for WASM exceptions
      if (typeof wasmModule.addOnAbort === 'function') {
        wasmModule.addOnAbort((what: string) => {
          console.error('WASM Abort:', what);
          // Try to get the stack trace
          if (typeof wasmModule.getExceptionMessage === 'function') {
            console.error('Exception details:', wasmModule.getExceptionMessage());
          }
          addResult(`WASM Abort: ${what}`);
        });
      }
      
      addResult('WASM module loaded successfully');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load WASM module');
    } finally {
      setLoading(false);
    }
  };

  const addResult = (message: string) => {
    setResults(prev => [...prev, `${new Date().toLocaleTimeString()}: ${message}`]);
  };

  // Helper function to capture stack trace at call site
  const captureStackTrace = () => {
    const error = new Error();
    if (error.stack) {
      const stack = error.stack.split('\n').slice(2, 8); // Skip Error and this function
      console.log('JavaScript call stack:', stack.join('\n'));
      return stack.join(' -> ');
    }
    return 'No stack trace available';
  };

  // Point generation and visualization
  const generateRandomPoints = () => {
    if (!instance) return;
    
    try {
      addResult(`Generating ${pointCount} random points...`);
      instance.generate_random_points(pointCount);
      const newPointCount = instance.get_point_count();
      addResult(`Generated ${newPointCount} random points with Morton ordering`);
    } catch (err) {
      console.error('Detailed error in generateRandomPoints:', err);
      const errorMessage = err instanceof Error ? err.message : String(err);
      addResult(`Generate points failed: ${errorMessage}`);
      
      // Try to get more detailed error info
      if (err instanceof Error && err.stack) {
        console.error('Stack trace:', err.stack);
        addResult(`Stack trace: ${err.stack.split('\n').slice(0, 3).join(' -> ')}`);
      }
    }
  };

  // First visualization: mk_hash_tree and log_hierarchy
  const buildTreeAndLogHierarchy = () => {
    if (!instance) return;
    
    try {
      const pointCount = instance.get_point_count();
      if (pointCount === 0) {
        addResult('No points available! Please generate points first.');
        return;
      }
      
      addResult(`Building hash tree for ${pointCount} points...`);
      console.log('About to call mk_hash_tree, JavaScript call stack:');
      captureStackTrace();
      
      // Explicitly call mk_hash_tree to build the tree structure
      const success = instance.mk_hash_tree();
      if (success) {
        addResult('Hash tree built successfully, logging hierarchy...');
        // Then call log_hierarchy to visualize the tree structure
        instance.log_hierarchy();
        addResult('Successfully logged hierarchy visualization');
      } else {
        addResult('Failed to build hash tree');
      }
    } catch (err) {
      console.error('Detailed error in buildTreeAndLogHierarchy:', err);
      const callStack = captureStackTrace();
      const errorMessage = err instanceof Error ? err.message : String(err);
      addResult(`Log hierarchy failed: ${errorMessage}`);
      addResult(`Call stack: ${callStack}`);
      
      // Try to get more detailed error info
      if (err instanceof Error && err.stack) {
        console.error('Full error stack:', err.stack);
      }
    }
  };

  // Second visualization: mk_hash_tree and log_bvh
  const buildTreeAndLogBVH = () => {
    if (!instance) return;
    
    try {
      const pointCount = instance.get_point_count();
      if (pointCount === 0) {
        addResult('No points available! Please generate points first.');
        return;
      }
      
      addResult(`Building hash tree for ${pointCount} points...`);
      // Explicitly call mk_hash_tree to build the tree structure
      const success = instance.mk_hash_tree();
      if (success) {
        addResult('Hash tree built successfully, logging BVH...');
        // Then call log_bvh to visualize bounding volume hierarchy
        instance.log_bvh();
        addResult('Successfully logged BVH visualization');
      } else {
        addResult('Failed to build hash tree');
      }
    } catch (err) {
      console.error('Detailed error in buildTreeAndLogBVH:', err);
      const errorMessage = err instanceof Error ? err.message : String(err);
      addResult(`Log BVH failed: ${errorMessage}`);
      
      // Try to get more detailed error info
      if (err instanceof Error && err.stack) {
        console.error('Stack trace:', err.stack);
        addResult(`Stack trace: ${err.stack.split('\n').slice(0, 3).join(' -> ')}`);
      }
    }
  };

  const clearPoints = () => {
    if (!instance) return;
    
    try {
      instance.clear_points();
      const remainingPoints = instance.get_point_count();
      addResult(`Cleared points: ${remainingPoints} points remaining`);
    } catch (err) {
      addResult(`Clear points failed: ${err}`);
    }
  };

  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="p-6 bg-white border-b">
        <h1 className="text-3xl font-bold mb-6">Morton Tree Test</h1>
        
        <div className="mb-6">
          <button
            onClick={loadWasm}
            disabled={loading}
            className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded disabled:opacity-50"
          >
            {loading ? 'Loading...' : 'Load WASM Module'}
          </button>
          
          {error && (
            <div className="mt-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
              Error: {error}
            </div>
          )}
          
          {instance && (
            <div className="mt-4 p-4 bg-green-100 border border-green-400 text-green-700 rounded">
              WASM module loaded successfully!
            </div>
          )}
        </div>

        {instance && (
          <div className="space-y-4">
            {/* Point Generation Controls */}
            <div className="p-4 bg-gray-50 rounded">
              <h3 className="text-lg font-semibold mb-2">Step 1: Generate Points</h3>
              <div className="flex items-center gap-4 mb-2">
                <label className="text-sm font-medium">Point Count:</label>
                <input
                  type="number"
                  min="1"
                  max="1000"
                  value={pointCount}
                  onChange={(e) => setPointCount(parseInt(e.target.value) || 50)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <button
                  onClick={generateRandomPoints}
                  className="bg-green-500 hover:bg-green-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate & Visualize Points
                </button>
                <button
                  onClick={clearPoints}
                  className="bg-red-500 hover:bg-red-700 text-white font-bold py-2 px-4 rounded"
                >
                  Clear Points
                </button>
              </div>
            </div>

            {/* Tree Visualization Controls */}
            <div className="p-4 bg-gray-50 rounded">
              <h3 className="text-lg font-semibold mb-2">Step 2: Build Tree & Visualize</h3>
              <div className="grid grid-cols-2 gap-4">
                <button
                  onClick={buildTreeAndLogHierarchy}
                  className="bg-orange-500 hover:bg-orange-700 text-white font-bold py-2 px-4 rounded"
                >
                  Build Tree → Log Hierarchy
                </button>
                <button
                  onClick={buildTreeAndLogBVH}
                  className="bg-yellow-500 hover:bg-yellow-700 text-white font-bold py-2 px-4 rounded"
                >
                  Build Tree → Log BVH
                </button>
              </div>
            </div>
          </div>
        )}
      </div>

      {/* 3D Visualization using GaudiLoggerRenderer */}
      <div className="flex-1 relative min-h-[600px]">
        <Canvas
          camera={{ position: [5, 5, 5], fov: 60, near: 0.01, far: 1000 }}
          className="!absolute !inset-0"
        >
          <ambientLight intensity={0.4} />
          <pointLight position={[10, 10, 10]} />
          
          {/* All visualizations are handled by GaudiLoggerRenderer */}
          <GaudiLoggerRenderer 
            wasmInstance={instance}
            wasmModule={module}
            isPlaying={false}
            showLines={true}
            showPoints={true}
            lineWidth={2}
            pointSize={0.1}
          />
          
          <OrbitControls />
          <gridHelper args={[10, 10]} />
          <axesHelper args={[2]} />
        </Canvas>
      </div>

      {/* Test Results */}
      {results.length > 0 && (
        <div className="p-6 bg-white border-t">
          <h2 className="text-xl font-semibold mb-4">Test Results</h2>
          <div className="bg-gray-100 p-4 rounded max-h-96 overflow-y-auto">
            {results.map((result, index) => (
              <div key={index} className="text-sm font-mono mb-1 text-gray-800">
                {result}
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
};
