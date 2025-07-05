import { useState, useEffect, useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { loadWasmModule } from '../utils/wasmLoader';

interface LoggerApiModule {
  // Animation
  step(): void;
  get_frame_count(): number;
  reset(): void;

  // Logger API
  add_point_to_gaudi_logger(x0: number, y0: number, z0: number, r: number, g: number, b: number, a: number): void;
  add_line_to_gaudi_logger(x0: number, y0: number, z0: number, x1: number, y1: number, z1: number, r: number, g: number, b: number, a: number): void;
  clear_gaudi_logger(): void;

  // Line data
  get_line_count(): number;
  get_line_positions_ptr(): number; // Pointer
  get_line_positions_size(): number; // Size in bytes
  get_line_colors_ptr(): number; // Pointer
  get_line_colors_size(): number; // Size in bytes

  // Point data
  get_point_count(): number;
  get_point_positions_ptr(): number; // Pointer
  get_point_positions_size(): number; // Size in bytes
  get_point_colors_ptr(): number; // Pointer
  get_point_colors_size(): number; // Size in bytes

  // Emscripten heap
  HEAPF64: Float64Array;
}

interface WasmState {
  loading: boolean;
  error: string | null;
  module: LoggerApiModule | null;
}

function WasmGeometry({ wasmModule }: { wasmModule: LoggerApiModule | null }) {
  const linesRef = useRef<THREE.LineSegments>(null);
  const pointsRef = useRef<THREE.Points>(null);

  useFrame(() => {
    if (!wasmModule) return;

    wasmModule.step();

    // Update Lines
    const lineCount = wasmModule.get_line_count();
    if (linesRef.current) {
      if (lineCount > 0) {
        const positionsPtr = wasmModule.get_line_positions_ptr();
        const colorsPtr = wasmModule.get_line_colors_ptr();
        const positions = new Float64Array(wasmModule.HEAPF64.buffer, positionsPtr, lineCount * 2 * 3);
        const colors = new Float64Array(wasmModule.HEAPF64.buffer, colorsPtr, lineCount * 2 * 4);
        
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
      }
    }

    // Update Points
    const pointCount = wasmModule.get_point_count();
    if (pointsRef.current) {
      if (pointCount > 0) {
        const positionsPtr = wasmModule.get_point_positions_ptr();
        const colorsPtr = wasmModule.get_point_colors_ptr();
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
      }
    }
  });

  return (
    <>
      <lineSegments ref={linesRef}>
        <bufferGeometry />
        <lineBasicMaterial vertexColors />
      </lineSegments>
      <points ref={pointsRef}>
        <bufferGeometry />
        <pointsMaterial vertexColors size={0.1} />
      </points>
    </>
  );
}


export function LoggerIntegrationTest() {
  const [wasmState, setWasmState] = useState<WasmState>({
    loading: true,
    error: null,
    module: null
  });

  const [stats, setStats] = useState({
    frameCount: 0,
    lineCount: 0,
    pointCount: 0,
  });

  useEffect(() => {
    const loadWasm = async () => {
      try {
        const wasmModule = await loadWasmModule('gaudi_logger_test');
        setWasmState({ 
          loading: false, 
          error: null, 
          module: wasmModule as LoggerApiModule 
        });
      } catch (error) {
        console.error('Failed to load WebAssembly module:', error);
        setWasmState({ 
          loading: false, 
          error: 'Failed to load WebAssembly module.', 
          module: null 
        });
      }
    };
    loadWasm();
  }, []);

  useFrame(() => {
    if (wasmState.module) {
      setStats({
        frameCount: wasmState.module.get_frame_count(),
        lineCount: wasmState.module.get_line_count(),
        pointCount: wasmState.module.get_point_count(),
      });
    }
  });

  const handleReset = () => {
    wasmState.module?.reset();
  };

  return (
    <div style={{ width: '100vw', height: '100vh' }}>
      <div style={{ position: 'absolute', top: 10, left: 10, color: 'white', zIndex: 1 }}>
        <div>Frame: {stats.frameCount}</div>
        <div>Lines: {stats.lineCount}</div>
        <div>Points: {stats.pointCount}</div>
        <button onClick={handleReset}>Reset Animation</button>
      </div>
      <Canvas>
        <color attach="background" args={['black']} />
        <OrbitControls />
        <axesHelper args={[5]} />
        {wasmState.module && <WasmGeometry wasmModule={wasmState.module} />}
      </Canvas>
      {wasmState.loading && <div style={{ position: 'absolute', top: '50%', left: '50%', color: 'white' }}>Loading WASM...</div>}
      {wasmState.error && <div style={{ position: 'absolute', top: '50%', left: '50%', color: 'red' }}>{wasmState.error}</div>}
    </div>
  );
}
