import React, { useState, useEffect, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { Button } from './ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { Play, Pause, RotateCcw, Plus, Trash2, Activity, Clock } from 'lucide-react';
import { loadWasmModule } from '../utils/wasmLoader';

interface GaudiLoggerTestModule {
  GaudiLoggerTest: new () => GaudiLoggerTestInstance;
  HEAPF64: Float64Array;
}

interface GaudiLoggerTestInstance {
  // Animation control
  step(): void;
  reset(): void;
  get_frame_count(): number;
  get_time(): number;
  
  // Line data access
  get_line_count(): number;
  get_line_positions_ptr(): number;
  get_line_positions_size(): number;
  get_line_colors_ptr(): number;
  get_line_colors_size(): number;
  
  // Point data access
  get_point_count(): number;
  get_point_positions_ptr(): number;
  get_point_positions_size(): number;
  get_point_colors_ptr(): number;
  get_point_colors_size(): number;
  
  // Manual testing
  add_test_line(x1: number, y1: number, z1: number, x2: number, y2: number, z2: number,
                r: number, g: number, b: number, a: number): void;
  add_test_point(x: number, y: number, z: number, r: number, g: number, b: number, a: number): void;
  clear_logger(): void;
}

interface WasmState {
  loading: boolean;
  error: string | null;
  module: GaudiLoggerTestModule | null;
  instance: GaudiLoggerTestInstance | null;
}

function GaudiLoggerGeometry({ 
  wasmModule, 
  wasmInstance, 
  isAnimating 
}: { 
  wasmModule: GaudiLoggerTestModule | null;
  wasmInstance: GaudiLoggerTestInstance | null;
  isAnimating: boolean;
}) {
  const linesRef = useRef<THREE.LineSegments>(null);
  const pointsRef = useRef<THREE.Points>(null);

  useFrame(() => {
    if (!wasmModule || !wasmInstance) return;

    // Step the animation
    if (isAnimating) {
      wasmInstance.step();
    }

    // Update Lines
    const lineCount = wasmInstance.get_line_count();
    if (linesRef.current) {
      if (lineCount > 0) {
        const positionsPtr = wasmInstance.get_line_positions_ptr();
        const colorsPtr = wasmInstance.get_line_colors_ptr();
        
        // Get line data from WASM heap
        const positions = new Float64Array(wasmModule.HEAPF64.buffer, positionsPtr, lineCount * 2 * 3);
        const colors = new Float64Array(wasmModule.HEAPF64.buffer, colorsPtr, lineCount * 2 * 4);
        
        // Convert to Three.js format
        const threePositions = new Float32Array(positions);
        const threeColors = new Float32Array(lineCount * 2 * 3);

        for (let i = 0; i < lineCount * 2; i++) {
          threeColors[i * 3 + 0] = colors[i * 4 + 0];
          threeColors[i * 3 + 1] = colors[i * 4 + 1];
          threeColors[i * 3 + 2] = colors[i * 4 + 2];
        }

        linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(threePositions, 3));
        linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(threeColors, 3));
        linesRef.current.geometry.attributes.position.needsUpdate = true;
        linesRef.current.geometry.attributes.color.needsUpdate = true;
      } else {
        linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 3));
      }
    }

    // Update Points
    const pointCount = wasmInstance.get_point_count();
    if (pointsRef.current) {
      if (pointCount > 0) {
        const positionsPtr = wasmInstance.get_point_positions_ptr();
        const colorsPtr = wasmInstance.get_point_colors_ptr();
        
        const positions = new Float64Array(wasmModule.HEAPF64.buffer, positionsPtr, pointCount * 3);
        const colors = new Float64Array(wasmModule.HEAPF64.buffer, colorsPtr, pointCount * 4);

        const threePositions = new Float32Array(positions);
        const threeColors = new Float32Array(pointCount * 3);

        for (let i = 0; i < pointCount; i++) {
          threeColors[i * 3 + 0] = colors[i * 4 + 0];
          threeColors[i * 3 + 1] = colors[i * 4 + 1];
          threeColors[i * 3 + 2] = colors[i * 4 + 2];
        }

        pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(threePositions, 3));
        pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(threeColors, 3));
        pointsRef.current.geometry.attributes.position.needsUpdate = true;
        pointsRef.current.geometry.attributes.color.needsUpdate = true;
      } else {
        pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
        pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 3));
      }
    }
  });

  return (
    <>
      <lineSegments ref={linesRef}>
        <bufferGeometry />
        <lineBasicMaterial vertexColors={true} />
      </lineSegments>
      <points ref={pointsRef}>
        <bufferGeometry />
        <pointsMaterial size={5} vertexColors={true} />
      </points>
    </>
  );
}

function ControlPanel({ 
  wasmInstance, 
  isAnimating, 
  setIsAnimating 
}: { 
  wasmInstance: GaudiLoggerTestInstance | null;
  isAnimating: boolean;
  setIsAnimating: (value: boolean) => void;
}) {
  const [frameCount, setFrameCount] = useState(0);
  const [time, setTime] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => {
      if (wasmInstance) {
        setFrameCount(wasmInstance.get_frame_count());
        setTime(wasmInstance.get_time());
      }
    }, 100);

    return () => clearInterval(interval);
  }, [wasmInstance]);

  const handleReset = () => {
    if (wasmInstance) {
      wasmInstance.reset();
      setFrameCount(0);
      setTime(0);
    }
  };

  const handleAddTestLine = () => {
    if (wasmInstance) {
      const x1 = (Math.random() - 0.5) * 4;
      const y1 = (Math.random() - 0.5) * 4;
      const z1 = (Math.random() - 0.5) * 4;
      const x2 = x1 + (Math.random() - 0.5) * 2;
      const y2 = y1 + (Math.random() - 0.5) * 2;
      const z2 = z1 + (Math.random() - 0.5) * 2;
      wasmInstance.add_test_line(x1, y1, z1, x2, y2, z2, Math.random(), Math.random(), Math.random(), 1.0);
    }
  };

  const handleClearLogger = () => {
    if (wasmInstance) {
      wasmInstance.clear_logger();
    }
  };

  return (
    <div style={{
      position: 'absolute',
      top: '20px',
      left: '20px',
      background: 'rgba(0, 0, 0, 0.8)',
      color: 'white',
      padding: '15px',
      borderRadius: '8px',
      fontFamily: 'monospace',
      fontSize: '12px',
      zIndex: 1000,
      minWidth: '250px'
    }}>
      <h3 style={{ margin: '0 0 10px 0', color: '#4a90e2' }}>Gaudi Logger Test</h3>
      
      <div style={{ marginBottom: '10px' }}>
        <strong>Frame:</strong> {frameCount}<br />
        <strong>Time:</strong> {time.toFixed(2)}s
      </div>
      
      <div style={{ display: 'flex', gap: '10px', marginBottom: '10px' }}>
        <button
          onClick={() => setIsAnimating(!isAnimating)}
          style={{
            padding: '5px 10px',
            background: isAnimating ? '#e74c3c' : '#2ecc71',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '11px'
          }}
        >
          {isAnimating ? 'Pause' : 'Play'}
        </button>
        
        <button
          onClick={handleReset}
          style={{
            padding: '5px 10px',
            background: '#3498db',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '11px'
          }}
        >
          Reset
        </button>
      </div>
      
      <div style={{ display: 'flex', gap: '10px', marginBottom: '10px' }}>
        <button
          onClick={handleAddTestLine}
          style={{
            padding: '5px 10px',
            background: '#9b59b6',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '11px'
          }}
        >
          Add Test Line
        </button>
        
        <button
          onClick={handleClearLogger}
          style={{
            padding: '5px 10px',
            background: '#e67e22',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '11px'
          }}
        >
          Clear Logger
        </button>
      </div>
      
      <div style={{ fontSize: '10px', color: '#bbb', marginTop: '10px' }}>
        <strong>Features:</strong><br />
        • Rotating coordinate frame<br />
        • Animated sine wave<br />
        • Orbiting box<br />
        • Manual line addition<br />
        • Real-time gaudi::logger interface
      </div>
    </div>
  );
}

export function GaudiLoggerTest() {
  const [wasmState, setWasmState] = useState<WasmState>({
    loading: true,
    error: null,
    module: null,
    instance: null
  });
  const [isAnimating, setIsAnimating] = useState(true);

  useEffect(() => {
    const loadWasm = async () => {
      try {
        const module = await loadWasmModule('gaudi_logger_test');
        const instance = new module.GaudiLoggerTest();
        
        setWasmState({
          loading: false,
          error: null,
          module,
          instance
        });
      } catch (error) {
        setWasmState({
          loading: false,
          error: error instanceof Error ? error.message : 'Unknown error',
          module: null,
          instance: null
        });
      }
    };

    loadWasm();
  }, []);

  if (wasmState.loading) {
    return (
      <div style={{
        width: '100vw',
        height: '100vh',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        background: '#1a1a1a',
        color: 'white',
        fontFamily: 'monospace'
      }}>
        <div style={{ textAlign: 'center' }}>
          <div style={{ fontSize: '18px', marginBottom: '10px' }}>Loading Gaudi Logger Test...</div>
          <div style={{ fontSize: '12px', color: '#888' }}>Initializing WASM module</div>
        </div>
      </div>
    );
  }

  if (wasmState.error) {
    return (
      <div style={{
        width: '100vw',
        height: '100vh',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        background: '#1a1a1a',
        color: 'white',
        fontFamily: 'monospace'
      }}>
        <div style={{ textAlign: 'center', maxWidth: '500px' }}>
          <div style={{ fontSize: '18px', marginBottom: '10px', color: '#e74c3c' }}>
            Failed to load WASM module
          </div>
          <div style={{ fontSize: '12px', color: '#888', marginBottom: '15px' }}>
            {wasmState.error}
          </div>
          <div style={{ fontSize: '11px', color: '#666' }}>
            Make sure the WASM module is built: run <code>build_gaudi_logger_test.bat</code>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div style={{ width: '100vw', height: '100vh' }}>
      <Canvas
        camera={{ position: [8, 6, 8], fov: 60 }}
        style={{ background: '#1a1a1a' }}
      >
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} />
        
        <GaudiLoggerGeometry 
          wasmModule={wasmState.module}
          wasmInstance={wasmState.instance}
          isAnimating={isAnimating}
        />
        
        <OrbitControls 
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
        />
        
        {/* Grid helper */}
        <gridHelper args={[20, 20, '#444', '#222']} />
        <axesHelper args={[3]} />
      </Canvas>
      
      <ControlPanel 
        wasmInstance={wasmState.instance}
        isAnimating={isAnimating}
        setIsAnimating={setIsAnimating}
      />
    </div>
  );
}
