import React, { Suspense } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import { useWasmModule, useWasmSimulation } from '../hooks/useWasmModule';
import { useSimulation } from '../hooks/useSimulation';
import { RodMesh, SdfVisualization } from './RodMesh';
import { ControlPanel } from './ControlPanel';
import type { RodSimulationProps } from '../types/simulation';

/**
 * Main Rod Simulation Component
 * Self-contained React component that can be easily imported into other projects
 */
export function RodSimulation({
  initialParams,
  onFrameUpdate,
  debug = false,
  autoPlay = true,
  style,
  className,
  width,
  height
}: RodSimulationProps) {
  // Load WebAssembly module
  const { module, isLoading: wasmLoading, error: wasmError } = useWasmModule();
  
  // Create simulation instance
  const { simulation, isReady: simReady, error: simError } = useWasmSimulation(module);
  
  // Manage simulation state
  const {
    data,
    params,
    isPlaying,
    frame,
    play,
    pause,
    reset,
    step,
    updateParams
  } = useSimulation(simulation, initialParams, onFrameUpdate);

  // Auto-play on ready
  React.useEffect(() => {
    if (simReady && autoPlay && !isPlaying) {
      play();
    }
  }, [simReady, autoPlay, isPlaying, play]);

  // Error handling
  if (wasmError || simError) {
    return (
      <div 
        className={className}
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          background: '#f8f9fa',
          color: '#dc3545',
          fontFamily: 'monospace',
          width: width || '100%',
          height: height || '400px',
          ...style
        }}
      >
        <div>
          <h3>Rod Simulation Error</h3>
          <p>{wasmError || simError}</p>
        </div>
      </div>
    );
  }

  // Loading state
  if (wasmLoading || !simReady) {
    return (
      <div 
        className={className}
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          background: '#f8f9fa',
          color: '#6c757d',
          fontFamily: 'monospace',
          width: width || '100%',
          height: height || '400px',
          ...style
        }}
      >
        <div>
          <h3>Loading Rod Simulation...</h3>
          <p>Initializing WebAssembly module and physics engine</p>
        </div>
      </div>
    );
  }

  return (
    <div 
      className={className}
      style={{
        position: 'relative',
        width: width || '100%',
        height: height || '400px',
        ...style
      }}
    >
      <Canvas
        camera={{ position: [3, 3, 3], fov: 60 }}
        style={{ background: '#1a1a1a' }}
      >
        <Suspense fallback={null}>
          {/* Lighting */}
          <ambientLight intensity={0.4} />
          <directionalLight 
            position={[5, 5, 5]} 
            intensity={0.8}
            castShadow
            shadow-mapSize-width={2048}
            shadow-mapSize-height={2048}
          />
          <pointLight position={[-5, -5, -5]} intensity={0.3} />

          {/* Environment */}
          <Environment preset="city" />

          {/* Rod visualization */}
          <RodMesh
            data={data}
            baseColor="#4a90e2"
            growthColorLow="#2ecc71"
            growthColorHigh="#e74c3c"
            metalness={0.1}
            roughness={0.4}
          />

          {/* SDF visualization (debug mode) */}
          {debug && (
            <SdfVisualization 
              frame={frame}
              opacity={0.2}
              visible={true}
            />
          )}

          {/* Camera controls */}
          <OrbitControls 
            enablePan={true}
            enableZoom={true}
            enableRotate={true}
            maxDistance={10}
            minDistance={1}
          />
        </Suspense>
      </Canvas>

      {/* Control Panel */}
      <ControlPanel
        params={params}
        onParamsChange={updateParams}
        isPlaying={isPlaying}
        onPlayPause={isPlaying ? pause : play}
        onReset={reset}
        onStep={step}
      />

      {/* Debug Info */}
      {debug && data && (
        <div style={{
          position: 'absolute',
          bottom: 10,
          right: 10,
          background: 'rgba(0, 0, 0, 0.8)',
          color: 'white',
          padding: '10px',
          borderRadius: '4px',
          fontFamily: 'monospace',
          fontSize: '11px'
        }}>
          <div>Frame: {frame}</div>
          <div>Vertices: {data.vertexCount}</div>
          <div>Length: {data.totalLength.toFixed(2)}</div>
        </div>
      )}
    </div>
  );
}
