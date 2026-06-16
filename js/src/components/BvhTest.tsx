import React, { useState, useEffect } from "react";
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { BvhTestModule, BvhTestInstance } from '../../wasm/modules/bvh_test/bvh_test_endpoints';
import { loadWasmModule } from '../utils/wasmLoader';
import { GaudiLoggerRenderer, WasmLoggerAPI, WasmModule as LoggerWasmModule } from './GaudiLoggerRenderer';
import { ConsoleLoggerPanel } from './ConsoleLoggerPanel';
import { consoleLogger } from '../stores/consoleLoggerStore';

interface BvhTestModuleWithLogger extends BvhTestModule, LoggerWasmModule {}
interface BvhTestInstanceWithLogger extends BvhTestInstance, WasmLoggerAPI {}

// Procedural radial (UV) sphere as an OBJ string -- mirrors asawa::make_sphere
// so the demo has no asset-file dependency.
const generateSphereObj = (radius = 1.0, uSeg = 24, vSeg = 16): string => {
  const v: string[] = [];
  const f: string[] = [];
  v.push(`v 0 ${radius} 0`); // top pole -> index 1
  for (let r = 1; r < vSeg; r++) {
    const phi = (Math.PI * r) / vSeg;
    const y = radius * Math.cos(phi);
    const rr = radius * Math.sin(phi);
    for (let s = 0; s < uSeg; s++) {
      const th = (2 * Math.PI * s) / uSeg;
      v.push(`v ${rr * Math.cos(th)} ${y} ${rr * Math.sin(th)}`);
    }
  }
  v.push(`v 0 ${-radius} 0`); // bottom pole
  const top = 1;
  const bottom = v.length;
  const ring = (r: number, s: number) =>
    2 + (r - 1) * uSeg + (((s % uSeg) + uSeg) % uSeg);
  for (let s = 0; s < uSeg; s++)
    f.push(`f ${top} ${ring(1, s)} ${ring(1, s + 1)}`);
  for (let r = 1; r < vSeg - 1; r++)
    for (let s = 0; s < uSeg; s++) {
      const a = ring(r, s);
      const b = ring(r + 1, s);
      const c = ring(r + 1, s + 1);
      const d = ring(r, s + 1);
      f.push(`f ${a} ${b} ${c}`);
      f.push(`f ${a} ${c} ${d}`);
    }
  const last = vSeg - 1;
  for (let s = 0; s < uSeg; s++)
    f.push(`f ${bottom} ${ring(last, s + 1)} ${ring(last, s)}`);
  return [...v, ...f].join("\n") + "\n";
};

export const BvhTest = () => {
  const [module, setModule] = useState<BvhTestModuleWithLogger | null>(null);
  const [instance, setInstance] = useState<BvhTestInstanceWithLogger | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [results, setResults] = useState<string[]>([]);
  const [meshLoaded, setMeshLoaded] = useState(false);
  const [tolerance, setTolerance] = useState(10.0);
  const [showDebugLines, setShowDebugLines] = useState(true);
  const [showDebugPoints, setShowDebugPoints] = useState(true);

  const addResult = (message: string) => {
    setResults(prev => [...prev, `${new Date().toLocaleTimeString()}: ${message}`]);
  };

  const loadWasm = async () => {
    try {
      setLoading(true);
      setError(null);
      const wasmModule = await loadWasmModule('bvh_test') as BvhTestModuleWithLogger;
      setModule(wasmModule);
      
      // Set up console logger callbacks for WASM
      consoleLogger.setupCallbacks(wasmModule);
      consoleLogger.info("BVH Test module loaded!");
      
      // Create an instance of the BvhTest class
      const testInstance = new wasmModule.BvhTest() as BvhTestInstanceWithLogger;
      setInstance(testInstance);
      
      addResult(`WASM module loaded successfully`);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load WASM module');
    } finally {
      setLoading(false);
    }
  };

  const loadSphere = async () => {
    if (!instance) return;
    
    try {
      // Procedurally generated sphere -- no asset file needed.
      const objContent = generateSphereObj();

      const success = instance.loadMeshFromString(objContent);
      if (success) {
        setMeshLoaded(true);
        addResult(`Loaded sphere: ${instance.getVertexCount()} vertices, ${instance.getEdgeCount()} edges, ${instance.getTriangleCount()} triangles`);
        instance.visualizeResults();
      } else {
        addResult('Failed to load mesh');
      }
    } catch (err) {
      addResult(`Load sphere failed: ${err}`);
    }
  };

  const buildEdgeBvh = () => {
    if (!instance) return;
    
    try {
      const success = instance.buildEdgeBvh();
      addResult(success ? 'Edge BVH built successfully' : 'Failed to build edge BVH');
      if (success) {
        instance.visualizeResults();
      }
    } catch (err) {
      addResult(`Build edge BVH failed: ${err}`);
    }
  };

  const buildTriBvh = () => {
    if (!instance) return;
    
    try {
      const success = instance.buildTriBvh();
      addResult(success ? 'Triangle BVH built successfully' : 'Failed to build triangle BVH');
    } catch (err) {
      addResult(`Build triangle BVH failed: ${err}`);
    }
  };

  const runPointToEdgeTest = () => {
    if (!instance) return;
    
    try {
      const passed = instance.testPointToEdge(tolerance);
      const bvhResult = instance.getLastBvhResult();
      const bruteResult = instance.getLastBruteResult();
      addResult(`Point-to-Edge: BVH=${bvhResult}, Brute=${bruteResult}, ${passed ? '✓ PASSED' : '✗ FAILED'}`);
    } catch (err) {
      addResult(`Point-to-Edge test failed: ${err}`);
    }
  };

  const runEdgeToEdgeTest = () => {
    if (!instance) return;
    
    try {
      const passed = instance.testEdgeToEdge(tolerance);
      const bvhResult = instance.getLastBvhResult();
      const bruteResult = instance.getLastBruteResult();
      addResult(`Edge-to-Edge: BVH=${bvhResult}, Brute=${bruteResult}, ${passed ? '✓ PASSED' : '✗ FAILED'}`);
    } catch (err) {
      addResult(`Edge-to-Edge test failed: ${err}`);
    }
  };

  const runPointToTriTest = () => {
    if (!instance) return;
    
    try {
      const passed = instance.testPointToTri(tolerance);
      addResult(`Point-to-Triangle: ${passed ? '✓ PASSED' : '✗ FAILED (not fully implemented)'}`);
    } catch (err) {
      addResult(`Point-to-Triangle test failed: ${err}`);
    }
  };

  const runAllTests = () => {
    addResult('--- Running all BVH tests ---');
    runPointToEdgeTest();
    runEdgeToEdgeTest();
    runPointToTriTest();
    addResult('--- All tests complete ---');
  };

  const runTestSuite = () => {
    if (!instance) return;

    try {
      const passed = instance.runAllTests();
      const total = instance.getLastSuiteTotal();
      const failed = instance.getLastSuiteFailed();
      addResult(`Test Suite: total=${total}, failed=${failed}, ${passed ? '✓ PASSED' : '✗ FAILED'}`);
    } catch (err) {
      addResult(`Test Suite failed: ${err}`);
    }
  };

  const clearVisualization = () => {
    if (!instance) return;
    instance.clearVisualization();
    addResult('Visualization cleared');
  };

  return (
    <div className="h-full flex flex-col">
      {/* Controls */}
      <div className="p-6 bg-white border-b">
        <h1 className="text-3xl font-bold mb-6">BVH Test</h1>
        
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
            {/* Mesh Loading */}
            <div className="p-4 bg-gray-50 rounded">
              <h3 className="text-lg font-semibold mb-2">Mesh Loading</h3>
              <div className="flex items-center gap-4">
                <button
                  onClick={loadSphere}
                  className="bg-green-500 hover:bg-green-700 text-white font-bold py-2 px-4 rounded"
                >
                  Load Sphere OBJ
                </button>
                {meshLoaded && (
                  <span className="text-green-600">
                    ✓ Mesh loaded: {instance.getVertexCount()} verts, {instance.getEdgeCount()} edges
                  </span>
                )}
              </div>
            </div>

            {/* BVH Building */}
            {meshLoaded && (
              <div className="p-4 bg-gray-50 rounded">
                <h3 className="text-lg font-semibold mb-2">Build BVH Structures</h3>
                <div className="flex items-center gap-4">
                  <button
                    onClick={buildEdgeBvh}
                    className="bg-purple-500 hover:bg-purple-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Build Edge BVH
                  </button>
                  <button
                    onClick={buildTriBvh}
                    className="bg-purple-500 hover:bg-purple-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Build Triangle BVH
                  </button>
                </div>
              </div>
            )}

            {/* Test Controls */}
            {meshLoaded && (
              <div className="p-4 bg-gray-50 rounded">
                <h3 className="text-lg font-semibold mb-2">BVH vs Brute Force Tests</h3>
                <div className="flex items-center gap-4 mb-3">
                  <label className="text-sm font-medium">Tolerance:</label>
                  <input
                    type="number"
                    step="0.1"
                    min="0.01"
                    value={tolerance}
                    onChange={(e) => setTolerance(parseFloat(e.target.value) || 10.0)}
                    className="border rounded px-2 py-1 w-24 text-black bg-white"
                  />
                </div>
                <div className="grid grid-cols-2 lg:grid-cols-5 gap-4">
                  <button
                    onClick={runPointToEdgeTest}
                    className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Point→Edge
                  </button>
                  <button
                    onClick={runEdgeToEdgeTest}
                    className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Edge→Edge
                  </button>
                  <button
                    onClick={runPointToTriTest}
                    className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Point→Triangle
                  </button>
                  <button
                    onClick={runAllTests}
                    className="bg-orange-500 hover:bg-orange-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Run All Tests
                  </button>
                  <button
                    onClick={runTestSuite}
                    className="bg-emerald-500 hover:bg-emerald-700 text-white font-bold py-2 px-4 rounded"
                  >
                    Run Test Suite
                  </button>
                </div>
              </div>
            )}

            {/* Display Controls */}
            <div className="flex flex-wrap gap-2">
              <button
                onClick={() => setShowDebugPoints(!showDebugPoints)}
                className={`px-3 py-1 rounded text-sm ${showDebugPoints ? 'bg-blue-500 text-white' : 'bg-gray-300 text-gray-700'}`}
              >
                {showDebugPoints ? 'Hide' : 'Show'} Debug Points
              </button>
              
              <button
                onClick={() => setShowDebugLines(!showDebugLines)}
                className={`px-3 py-1 rounded text-sm ${showDebugLines ? 'bg-blue-500 text-white' : 'bg-gray-300 text-gray-700'}`}
              >
                {showDebugLines ? 'Hide' : 'Show'} Debug Lines
              </button>
              
              <button
                onClick={clearVisualization}
                className="px-3 py-1 rounded text-sm bg-red-500 text-white hover:bg-red-700"
              >
                Clear Visualization
              </button>
            </div>
          </div>
        )}
      </div>

      {/* 3D Visualization */}
      <div className="flex-1 relative min-h-[400px]">
        <Canvas
          camera={{ position: [3, 3, 3], fov: 60, near: 0.01, far: 1000 }}
          className="!absolute !inset-0"
        >
          <ambientLight intensity={0.4} />
          <pointLight position={[10, 10, 10]} />
          
          {module && instance && (
            <GaudiLoggerRenderer 
              wasmInstance={instance}
              wasmModule={module}
              showLines={showDebugLines}
              showPoints={showDebugPoints}
            />
          )}
          
          <OrbitControls />
          <gridHelper args={[10, 10]} />
          <axesHelper args={[2]} />
        </Canvas>
      </div>

      {/* Test Results */}
      <div className="p-6 bg-white border-t">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Test Results Log */}
          <div>
            <h3 className="text-lg font-semibold mb-2">Test Results</h3>
            <div className="bg-gray-900 text-green-400 p-4 rounded font-mono text-sm h-48 overflow-y-auto">
              {results.length === 0 ? (
                <div className="text-gray-500">No results yet. Load a mesh and run tests.</div>
              ) : (
                results.map((result, index) => (
                  <div key={index} className={result.includes('FAILED') ? 'text-red-400' : result.includes('PASSED') ? 'text-green-400' : 'text-gray-300'}>
                    {result}
                  </div>
                ))
              )}
            </div>
          </div>
          
          {/* Console Logger */}
          <div>
            <ConsoleLoggerPanel 
              title="Console Logger Output"
              height={192}
              showTimestamps={true}
              showLevelIcons={true}
              autoScroll={true}
            />
          </div>
        </div>
      </div>
    </div>
  );
};


