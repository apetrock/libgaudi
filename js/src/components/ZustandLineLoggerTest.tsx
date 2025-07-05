import { useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { useDebugLines } from '../stores/lineLoggerStore';
import { LineRenderer, EnhancedLineRenderer, LineRendererDebugInfo, useFlushAfterRender } from './ZustandLineRenderer';

/**
 * Animation component that generates test debug lines using Zustand store
 */
function LineAnimator({ mode = 'spiral' }: { mode?: 'spiral' | 'cube' | 'random' | 'grid' }) {
  const debugLines = useDebugLines();

  // Auto-flush after each frame
  useFlushAfterRender();
  useFrame((state) => {
    const time = state.clock.elapsedTime;
    
    // Clear old lines every few seconds
    if (Math.floor(time) % 4 === 0 && time % 1 < 0.02) {
      debugLines.clear();
    }

    // Generate different patterns based on mode
    switch (mode) {
      case 'spiral':
        generateSpiralLines(debugLines, time);
        break;
      case 'cube':
        generateCubeLines(debugLines, time);
        break;
      case 'random':
        generateRandomLines(debugLines, time);
        break;
      case 'grid':
        generateGridLines(debugLines, time);
        break;
    }
    
    // Debug: log line count occasionally
    if (Math.floor(time * 10) % 30 === 0) {
      console.log(`Line count: ${debugLines.lineCount}, Mode: ${mode}, Time: ${time.toFixed(2)}`);
    }
  });

  return null;
}

/**
 * Generate spiral pattern lines
 */
function generateSpiralLines(debugLines: any, time: number) {
  const numLines = 6;
  for (let i = 0; i < numLines; i++) {
    const angle = time + (i * Math.PI * 2) / numLines;
    const radius = 2 + Math.sin(time * 0.5) * 0.5;
    
    const start = new THREE.Vector3(
      Math.cos(angle) * radius,
      Math.sin(time * 0.3 + i) * 0.5,
      Math.sin(angle) * radius
    );
    
    const end = new THREE.Vector3(
      Math.cos(angle + 0.3) * (radius + 0.4),
      Math.sin(time * 0.3 + i + 0.3) * 0.5,
      Math.sin(angle + 0.3) * (radius + 0.4)
    );
    
    const hue = (i * 60 + time * 20) % 360;
    const color = `hsl(${hue}, 80%, 60%)`;
    debugLines.addLine(start, end, color);
  }
}

/**
 * Generate animated cube wireframe
 */
function generateCubeLines(debugLines: any, time: number) {
  const size = 1.5 + Math.sin(time * 0.7) * 0.3;
  const vertices = [
    [-size, -size, -size], [size, -size, -size], [size, size, -size], [-size, size, -size],
    [-size, -size, size], [size, -size, size], [size, size, size], [-size, size, size]
  ];
  
  const edges = [
    [0,1], [1,2], [2,3], [3,0], // bottom face
    [4,5], [5,6], [6,7], [7,4], // top face
    [0,4], [1,5], [2,6], [3,7]  // vertical edges
  ];
  
  edges.forEach(([i, j], edgeIndex) => {
    const start = new THREE.Vector3(...vertices[i]);
    const end = new THREE.Vector3(...vertices[j]);
    const hue = (edgeIndex * 30 + time * 15) % 360;
    debugLines.addLine(start, end, `hsl(${hue}, 70%, 50%)`);
  });
}

/**
 * Generate random lines
 */
function generateRandomLines(debugLines: any, time: number) {
  const numLines = 10;
  for (let i = 0; i < numLines; i++) {
    const start = new THREE.Vector3(
      (Math.random() - 0.5) * 4,
      (Math.random() - 0.5) * 4,
      (Math.random() - 0.5) * 4
    );
    const end = new THREE.Vector3(
      start.x + (Math.random() - 0.5) * 2,
      start.y + (Math.random() - 0.5) * 2,
      start.z + (Math.random() - 0.5) * 2
    );
    
    const colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#ffc107', '#e056fd'];
    const color = colors[Math.floor(time * 2 + i) % colors.length];
    debugLines.addLine(start, end, color);
  }
}

/**
 * Generate animated grid pattern
 */
function generateGridLines(debugLines: any, time: number) {
  const gridSize = 4;
  const spacing = 0.4;
  const offset = Math.sin(time * 0.8) * 0.5;
  
  // Horizontal lines
  for (let i = -gridSize; i <= gridSize; i++) {
    const y = i * spacing + offset;
    debugLines.addLine(
      new THREE.Vector3(-gridSize * spacing, y, 0),
      new THREE.Vector3(gridSize * spacing, y, 0),
      '#00ff88'
    );
  }
  
  // Vertical lines  
  for (let i = -gridSize; i <= gridSize; i++) {
    const x = i * spacing;
    debugLines.addLine(
      new THREE.Vector3(x, -gridSize * spacing + offset, 0),
      new THREE.Vector3(x, gridSize * spacing + offset, 0),
      '#ff8800'
    );
  }
}

/**
 * Control panel for Zustand line logger test
 */
function ZustandLineLoggerControls({ 
  pattern, 
  setPattern,
  useEnhanced,
  setUseEnhanced 
}: {
  pattern: string;
  setPattern: (pattern: string) => void;
  useEnhanced: boolean;
  setUseEnhanced: (use: boolean) => void;
}) {
  const debugLines = useDebugLines();
  const [fadeLines, setFadeLines] = useState(false);
  const [frameColors, setFrameColors] = useState(false);

  return (
    <div style={{
      position: 'absolute',
      top: 10,
      left: 10,
      background: 'rgba(0, 0, 0, 0.9)',
      color: 'white',
      padding: '15px',
      borderRadius: '8px',
      fontFamily: 'monospace',
      fontSize: '12px',
      zIndex: 1000,
      minWidth: '220px'
    }}>
      <h3 style={{ margin: '0 0 15px 0', color: '#4a90e2' }}>
        Zustand Line Logger
      </h3>
      
      <div style={{ marginBottom: '10px' }}>
        <label>Pattern:</label>
        <select 
          value={pattern}
          onChange={(e) => setPattern(e.target.value)}
          style={{ width: '100%', marginTop: '5px' }}
        >
          <option value="spiral">Spiral Animation</option>
          <option value="cube">Rotating Cube</option>
          <option value="random">Random Lines</option>
          <option value="grid">Animated Grid</option>
        </select>
      </div>

      <div style={{ marginBottom: '10px' }}>
        <label>
          <input
            type="checkbox"
            checked={debugLines.enabled}
            onChange={(e) => debugLines.setEnabled(e.target.checked)}
          />
          {' '}Enable Logging
        </label>
      </div>

      <div style={{ marginBottom: '10px' }}>
        <label>
          <input
            type="checkbox"
            checked={useEnhanced}
            onChange={(e) => setUseEnhanced(e.target.checked)}
          />
          {' '}Enhanced Renderer
        </label>
      </div>

      {useEnhanced && (
        <>
          <div style={{ marginBottom: '10px' }}>
            <label>
              <input
                type="checkbox"
                checked={fadeLines}
                onChange={(e) => setFadeLines(e.target.checked)}
              />
              {' '}Fade Old Lines
            </label>
          </div>

          <div style={{ marginBottom: '10px' }}>
            <label>
              <input
                type="checkbox"
                checked={frameColors}
                onChange={(e) => setFrameColors(e.target.checked)}
              />
              {' '}Frame Colors
            </label>
          </div>
        </>
      )}

      <div style={{ marginBottom: '10px' }}>
        <label>Max Lines: {debugLines.maxLines}</label>
        <input
          type="range"
          min="100"
          max="2000"
          step="100"
          value={debugLines.maxLines}
          onChange={(e) => debugLines.setMaxLines(parseInt(e.target.value))}
          style={{ width: '100%', marginTop: '5px' }}
        />
      </div>

      <button 
        onClick={() => debugLines.clear()}
        style={{
          width: '100%',
          padding: '8px',
          marginTop: '10px',
          background: '#e74c3c',
          color: 'white',
          border: 'none',
          borderRadius: '4px',
          cursor: 'pointer',
          fontSize: '12px'
        }}
      >
        Clear All Lines
      </button>

      <div style={{ 
        marginTop: '15px',
        fontSize: '11px',
        opacity: 0.8,
        lineHeight: '1.4'
      }}>
        <p style={{ margin: '0 0 5px 0' }}>
          ✅ Zustand singleton store
        </p>
        <p style={{ margin: '0 0 5px 0' }}>
          ✅ Global line accumulation
        </p>
        <p style={{ margin: '0 0 5px 0' }}>
          ✅ Frame-based flushing
        </p>
        <p style={{ margin: 0 }}>
          ✅ Ready for WASM integration
        </p>
      </div>
    </div>
  );
}

/**
 * Main Zustand line logger test component
 */
export function ZustandLineLoggerTest() {
  const [useEnhanced, setUseEnhanced] = useState(false);
  const [pattern, setPattern] = useState('spiral');

  return (
    <div style={{ width: '100%', height: '100%', position: 'relative' }}>
      <Canvas
        camera={{ position: [6, 6, 6], fov: 60 }}
        style={{ background: 'linear-gradient(135deg, #0a0a0a 0%, #1a1a3a 100%)' }}
      >
        {/* Lighting */}
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} intensity={0.4} />

        {/* Line animator - generates debug lines */}
        <LineAnimator mode={pattern as 'spiral' | 'cube' | 'random' | 'grid'} />

        {/* Line renderer - displays accumulated lines */}
        {useEnhanced ? (
          <EnhancedLineRenderer 
            fadeOldLines={true}
            maxFrameAge={8}
            showFrameColors={false}
          />
        ) : (
          <LineRenderer />
        )}

        {/* Camera controls */}
        <OrbitControls 
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
          maxDistance={15}
          minDistance={2}
        />
      </Canvas>

      {/* Controls */}
      <ZustandLineLoggerControls 
        pattern={pattern}
        setPattern={setPattern}
        useEnhanced={useEnhanced}
        setUseEnhanced={setUseEnhanced}
      />
      
      {/* Debug info */}
      <LineRendererDebugInfo />

      {/* Status indicator */}
      <div style={{
        position: 'absolute',
        bottom: 10,
        right: 10,
        background: 'rgba(74, 144, 226, 0.8)',
        color: 'white',
        padding: '8px 12px',
        borderRadius: '4px',
        fontFamily: 'monospace',
        fontSize: '11px'
      }}>
        Phase 0: Zustand Line Logger ✅
      </div>
    </div>
  );
}
