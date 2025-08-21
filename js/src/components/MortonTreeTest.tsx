import React, { useState, useRef } from "react";
import { Canvas } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import { MortonTreeTestModule, MortonTreeTestClass } from '../../wasm/modules/morton_tree/morton_tree_endpoints';
import { loadWasmModule } from '../utils/wasmLoader';
import { GaudiLoggerRenderer, WasmLoggerAPI, WasmModule as LoggerWasmModule } from './GaudiLoggerRenderer';
import { ConsoleLoggerPanel } from './ConsoleLoggerPanel';
import { consoleLogger } from '../stores/consoleLoggerStore';

interface MortonTreeTestModuleWithLogger extends MortonTreeTestModule, LoggerWasmModule {}
interface MortonTreeTestClassWithLogger extends MortonTreeTestClass, WasmLoggerAPI {
  // No additional properties needed - generate_grid_points is already in MortonTreeTestClass
}

export const MortonTreeTest = () => {
  const [module, setModule] = useState<MortonTreeTestModuleWithLogger | null>(null);
  const [instance, setInstance] = useState<MortonTreeTestClassWithLogger | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [results, setResults] = useState<string[]>([]);
  const [pointCount, setPointCount] = useState(1000);
  const [gridSize, setGridSize] = useState(5);
  const [knotSegments, setKnotSegments] = useState(100);
  const [knotType, setKnotType] = useState(2); // K=2 for cinquefoil
  const [randomKnotControlPoints, setRandomKnotControlPoints] = useState(8);
  const [randomKnotSegmentsPerChord, setRandomKnotSegmentsPerChord] = useState(10);
  const [showDebugLines, setShowDebugLines] = useState(true);
  const [showDebugPoints, setShowDebugPoints] = useState(true);
  const currentTimeRef = useRef(0);
  const [animationFrameId, setAnimationFrameId] = useState<number | null>(null);
  const isAnimatingRef = useRef(false);
  const [isAnimating, setIsAnimating] = useState(false); // Keep for UI display

  const loadWasm = async () => {
    try {
      setLoading(true);
      setError(null);
      const wasmModule = await loadWasmModule('morton_tree_test') as MortonTreeTestModuleWithLogger;
      setModule(wasmModule);
      
      // Set up console logger callbacks for WASM
      consoleLogger.setupCallbacks(wasmModule);
      
      // Test the console logger is working
      consoleLogger.info("Console logger initialized successfully!");
      consoleLogger.debug("WASM module loaded and ready");
      
      // Create an instance of the MortonTreeTest class
      const testInstance = new wasmModule.MortonTreeTest() as MortonTreeTestClassWithLogger;
      setInstance(testInstance);
      
      // Debug: Log all available functions
      console.log('Available functions on WASM instance:', Object.getOwnPropertyNames(testInstance));
      console.log('Prototype functions:', Object.getOwnPropertyNames(Object.getPrototypeOf(testInstance)));
      
      addResult(`WASM module loaded. Available functions: ${Object.getOwnPropertyNames(testInstance).join(', ')}`);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load WASM module');
    } finally {
      setLoading(false);
    }
  };

  const addResult = (message: string) => {
    setResults(prev => [...prev, `${new Date().toLocaleTimeString()}: ${message}`]);
  };

  const clearPoints = () => {
    if (!instance) return;
    
    try {
      instance.clear_points();
      const pointCount = instance.get_point_count();
      addResult(`Cleared points: ${pointCount} points remaining`);
    } catch (err) {
      addResult(`Clear points failed: ${err}`);
    }
  };

  const generateRandomPoints = () => {
    if (!instance) return;
    
    try {
      instance.generate_random_points(pointCount);
      const newPointCount = instance.get_point_count();
      addResult(`Generated ${newPointCount} random points`);
    } catch (err) {
      addResult(`Generate points failed: ${err}`);
    }
  };

  const generateGridPoints = () => {
    if (!instance) return;
    
    try {
      if (instance.generate_grid_points) {
        instance.generate_grid_points(gridSize);
        const newPointCount = instance.get_point_count();
        addResult(`Generated ${newPointCount} grid points (${gridSize}x${gridSize}x${gridSize} grid)`);
      } else {
        addResult(`Grid points function not available`);
      }
    } catch (err) {
      addResult(`Generate grid points failed: ${err}`);
    }
  };

  const generateTrefoilKnot = () => {
    if (!instance) return;
    
    try {
      if (instance.generate_trefoil_knot) {
        instance.generate_trefoil_knot(knotSegments);
        const newPointCount = instance.get_point_count();
        addResult(`Generated trefoil knot with ${newPointCount} points (${knotSegments} segments)`);
      } else {
        addResult(`Trefoil knot function not available`);
      }
    } catch (err) {
      addResult(`Generate trefoil knot failed: ${err}`);
    }
  };

  const generateKnot = () => {
    if (!instance) return;
    
    try {
      if (instance.generate_knot) {
        instance.generate_knot(knotType, knotSegments);
        const newPointCount = instance.get_point_count();
        const knotNames = ['', 'Trefoil', 'Cinquefoil', 'Septafoil', 'Nonafoil'];
        const knotName = knotType < knotNames.length ? knotNames[knotType] : `K=${knotType}`;
        addResult(`Generated ${knotName} knot with ${newPointCount} points (K=${knotType}, ${knotSegments} segments)`);
      } else {
        addResult(`General knot function not available`);
      }
    } catch (err) {
      addResult(`Generate knot failed: ${err}`);
    }
  };

  const generateRandomKnot = () => {
    if (!instance) return;
    
    try {
      if (instance.generate_random_knot) {
        instance.generate_random_knot(randomKnotControlPoints, randomKnotSegmentsPerChord);
        const newPointCount = instance.get_point_count();
        addResult(`Generated random knot with ${newPointCount} points (${randomKnotControlPoints} control points, ${randomKnotSegmentsPerChord} segments per chord)`);
      } else {
        addResult(`Random knot function not available`);
      }
    } catch (err) {
      addResult(`Generate random knot failed: ${err}`);
    }
  };


  const logZorder = () => {
    if (!instance) return;
    try {
      instance.draw_all_visualizations();
      addResult('Z-order visualization logged');
    } catch (err) {
      addResult(`Log Z-order failed: ${err}`);
    }
  };

  const logHierarchy = () => {
    if (!instance) return;
    
    try {
      instance.log_hierarchy();
      addResult('Hierarchy visualization logged');
    } catch (err) {
      addResult(`Log hierarchy failed: ${err}`);
    }
  };

  const logBVH = () => {
    if (!instance) return;
    
    try {
      instance.log_bvh();
      addResult('BVH visualization logged');
    } catch (err) {
      addResult(`Log BVH failed: ${err}`);
    }
  };

  const testLogNearest = () => {
    if (!instance) return;
    
    try {
      // Use current time as parameter
      const timeParam = currentTimeRef.current;
      instance.log_nearest(timeParam);
      addResult(`Log nearest test called with time parameter: ${timeParam}`);
    } catch (err) {
      addResult(`Log nearest test failed: ${err}`);
    }
  };

  // Animation frame callback for continuous nearest neighbor testing
  const animationFrameCallback = () => {
    if (!instance || !isAnimatingRef.current) return;
    
    try {
      instance.log_nearest(currentTimeRef.current);
      currentTimeRef.current += 1; // Increment time for smooth animation
      
      // Schedule next frame
      const nextFrameId = requestAnimationFrame(animationFrameCallback);
      setAnimationFrameId(nextFrameId);
    } catch (err) {
      addResult(`Animation frame callback failed: ${err}`);
      stopAnimation();
    }
  };

  const startAnimation = () => {
    if (!instance || isAnimatingRef.current) return;
    
    try {
      currentTimeRef.current = 0;
      isAnimatingRef.current = true;
      setIsAnimating(true);
      addResult('Starting nearest neighbor animation...');
      
      // Start the animation loop
      const frameId = requestAnimationFrame(animationFrameCallback);
      setAnimationFrameId(frameId);
    } catch (err) {
      addResult(`Failed to start animation: ${err}`);
      isAnimatingRef.current = false;
      setIsAnimating(false);
    }
  };

  const stopAnimation = () => {
    if (animationFrameId !== null) {
      cancelAnimationFrame(animationFrameId);
      setAnimationFrameId(null);
    }
    isAnimatingRef.current = false;
    setIsAnimating(false);
    addResult('Stopped nearest neighbor animation');
  };

  // Cleanup animation on unmount
  React.useEffect(() => {
    return () => {
      if (animationFrameId !== null) {
        cancelAnimationFrame(animationFrameId);
      }
    };
  }, [animationFrameId]);

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
              <h3 className="text-lg font-semibold mb-2">Point Generation</h3>
              
              {/* Random Points Row */}
              <div className="flex items-center gap-4 mb-3">
                <label className="text-sm font-medium">Random Count:</label>
                <input
                  type="number"
                  min="1"
                  max="50000"
                  value={pointCount}
                  onChange={(e) => setPointCount(parseInt(e.target.value) || 1000)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <button
                  onClick={generateRandomPoints}
                  className="bg-green-500 hover:bg-green-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate Random Points
                </button>
              </div>
              
              {/* Grid Points Row */}
              <div className="flex items-center gap-4 mb-3">
                <label className="text-sm font-medium">Grid Size:</label>
                <input
                  type="number"
                  min="2"
                  max="20"
                  value={gridSize}
                  onChange={(e) => setGridSize(parseInt(e.target.value) || 5)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <button
                  onClick={generateGridPoints}
                  className="bg-purple-500 hover:bg-purple-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate Grid Points
                </button>
              </div>
              
              {/* Trefoil Knot Row */}
              <div className="flex items-center gap-4 mb-3">
                <label className="text-sm font-medium">Knot Segments:</label>
                <input
                  type="number"
                  min="10"
                  max="1000"
                  value={knotSegments}
                  onChange={(e) => setKnotSegments(parseInt(e.target.value) || 100)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <button
                  onClick={generateTrefoilKnot}
                  className="bg-indigo-500 hover:bg-indigo-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate Trefoil Knot
                </button>
              </div>
              
              {/* General Knot Row */}
              <div className="flex items-center gap-4 mb-3">
                <label className="text-sm font-medium">K:</label>
                <input
                  type="number"
                  min="1"
                  max="10"
                  value={knotType}
                  onChange={(e) => setKnotType(parseInt(e.target.value) || 2)}
                  className="border rounded px-2 py-1 w-16 text-black bg-white"
                />
                <label className="text-sm font-medium">Segments:</label>
                <input
                  type="number"
                  min="10"
                  max="1000"
                  value={knotSegments}
                  onChange={(e) => setKnotSegments(parseInt(e.target.value) || 100)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <button
                  onClick={generateKnot}
                  className="bg-teal-500 hover:bg-teal-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate Knot
                </button>
              </div>
              
              {/* Random Knot Row */}
              <div className="flex items-center gap-4 mb-3">
                <label className="text-sm font-medium">Control Points:</label>
                <input
                  type="number"
                  min="4"
                  max="20"
                  value={randomKnotControlPoints}
                  onChange={(e) => setRandomKnotControlPoints(parseInt(e.target.value) || 8)}
                  className="border rounded px-2 py-1 w-20 text-black bg-white"
                />
                <label className="text-sm font-medium">Segments/Chord:</label>
                <input
                  type="number"
                  min="5"
                  max="50"
                  value={randomKnotSegmentsPerChord}
                  onChange={(e) => setRandomKnotSegmentsPerChord(parseInt(e.target.value) || 10)}
                  className="border rounded px-2 py-1 w-24 text-black bg-white"
                />
                <button
                  onClick={generateRandomKnot}
                  className="bg-pink-500 hover:bg-pink-700 text-white font-bold py-2 px-4 rounded"
                >
                  Generate Random Knot
                </button>
              </div>
              
              {/* Knot Type Help */}
              <div className="text-xs text-gray-600 mb-3">
                K=1: Trefoil, K=2: Cinquefoil, K=3: Septafoil, K=4: Nonafoil
              </div>
              
              {/* Random Knot Help */}
              <div className="text-xs text-gray-600 mb-3">
                Random knots use Bezier splines to connect Gaussian control points
              </div>
              
              
              {/* Clear Button Row */}
              <div className="flex items-center gap-4">
                <button
                  onClick={clearPoints}
                  className="bg-red-500 hover:bg-red-700 text-white font-bold py-2 px-4 rounded"
                >
                  Clear Points
                </button>
              </div>
            </div>

            {/* Visualization Controls */}
            <div className="p-4 bg-gray-50 rounded">
              <h3 className="text-lg font-semibold mb-2">Visualizations</h3>
              <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
                <button
                  onClick={logZorder}
                  className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded"
                >
                  Log Z Order
                </button>
                <button
                  onClick={logHierarchy}
                  className="bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded"
                >
                  Log Hierarchy
                </button>
                <button
                  onClick={logBVH}
                  className="bg-orange-500 hover:bg-orange-700 text-white font-bold py-2 px-4 rounded"
                >
                  Log BVH
                </button>
                <button
                  onClick={testLogNearest}
                  className="bg-green-500 hover:bg-green-700 text-white font-bold py-2 px-4 rounded"
                >
                  Test Log Nearest
                </button>
                <button
                  onClick={startAnimation}
                  disabled={isAnimating}
                  className="bg-purple-500 hover:bg-purple-700 text-white font-bold py-2 px-4 rounded disabled:opacity-50"
                >
                  Start Animation
                </button>
                <button
                  onClick={stopAnimation}
                  disabled={!isAnimating}
                  className="bg-red-500 hover:bg-red-700 text-white font-bold py-2 px-4 rounded disabled:opacity-50"
                >
                  Stop Animation
                </button>
              </div>
              
              {/* Time Parameter Display */}
              <div className="mt-3 text-sm text-gray-600">
                Current time parameter: {currentTimeRef.current} {isAnimating && '(Animating)'}
              </div>
            </div>

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
            </div>
          </div>
        )}
      </div>

      {/* 3D Visualization */}
      <div className="flex-1 relative min-h-[600px]">
        <Canvas
          camera={{ position: [5, 5, 5], fov: 60, near: 0.01, far: 1000 }}
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
          
          {/* Console Logger Column */}
          <div>
            <ConsoleLoggerPanel 
              title="Console Logger Output"
              height={384}
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