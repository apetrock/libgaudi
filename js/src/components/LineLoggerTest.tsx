import React, { useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { LineLoggerProvider, useDebugLines } from '../contexts/LineLoggerContext';
import { LineRenderer, EnhancedLineRenderer, LineRendererDebugInfo } from './LineRenderer';

/**
 * Animation component that generates test debug lines
 */
function LineAnimator() {
  const debugLines = useDebugLines();
  const [mode, setMode] = useState<'spiral' | 'cube' | 'random' | 'grid'>('spiral');

  useFrame((state) => {
    const time = state.clock.elapsedTime;
    
    // Clear old lines every few seconds
    if (Math.floor(time) % 3 === 0 && time % 1 < 0.02) {
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

    // Flush lines for this frame
    debugLines.flush();
  });

  // Expose mode setter for external control
  React.useImperativeHandle(React.createRef(), () => ({
    setMode
  }), []);

  return null;
}

/**
 * Generate spiral pattern lines
 */
function generateSpiralLines(debugLines: any, time: number) {
  const numLines = 5;
  for (let i = 0; i < numLines; i++) {
    const angle = time + (i * Math.PI * 2) / numLines;
    const radius = 2 + Math.sin(time * 0.5) * 0.5;
    
    const start = new THREE.Vector3(
      Math.cos(angle) * radius,
      Math.sin(time * 0.3 + i) * 0.5,
      Math.sin(angle) * radius
    );
    
    const end = new THREE.Vector3(
      Math.cos(angle + 0.2) * (radius + 0.3),
      Math.sin(time * 0.3 + i + 0.2) * 0.5,
      Math.sin(angle + 0.2) * (radius + 0.3)
    );
    
    const hue = (i * 60 + time * 30) % 360;
    debugLines.addLine(start, end, `hsl(${hue}, 80%, 60%)`);
  }
}

/**
 * Generate cube wireframe lines
 */
function generateCubeLines(debugLines: any, time: number) {
  const size = 1.5 + Math.sin(time) * 0.3;
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
    const hue = (edgeIndex * 30 + time * 20) % 360;
    debugLines.addLine(start, end, `hsl(${hue}, 70%, 50%)`);
  });
}

/**
 * Generate random lines
 */
function generateRandomLines(debugLines: any, time: number) {
  const numLines = 8;
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
    const color = colors[Math.floor(time + i) % colors.length];
    debugLines.addLine(start, end, color);
  }
}

/**
 * Generate grid pattern lines
 */
function generateGridLines(debugLines: any, time: number) {
  const gridSize = 3;
  const spacing = 0.5;
  const offset = Math.sin(time * 0.5) * 0.3;
  
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
 * Control panel for line logger test
 */
function LineLoggerControls() {
  const debugLines = useDebugLines();
  const [useEnhanced, setUseEnhanced] = useState(false);
  const [fadeLines, setFadeLines] = useState(false);
  const [frameColors, setFrameColors] = useState(false);

  return (
    <div style={{
      position: 'absolute',
      top: 10,
      left: 10,
      background: 'rgba(0, 0, 0, 0.8)',
      color: 'white',
      padding: '15px',
      borderRadius: '8px',
      fontFamily: 'monospace',
      fontSize: '12px',
      zIndex: 1000,
      minWidth: '200px'
    }}>
      <h3 style={{ margin: '0 0 15px 0', color: '#4a90e2' }}>
        Line Logger Test
      </h3>
      
      <div style={{ marginBottom: '10px' }}>
        <label>Pattern:</label>
        <select 
          onChange={(e) => {
            // Would need to implement pattern switching
            console.log('Pattern:', e.target.value);
          }}
          style={{ width: '100%', marginTop: '5px' }}
        >
          <option value="spiral">Spiral</option>
          <option value="cube">Cube</option>
          <option value="random">Random</option>
          <option value="grid">Grid</option>
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

      <button 
        onClick={() => debugLines.clear()}
        style={{
          width: '100%',
          padding: '5px',
          marginTop: '10px',
          background: '#e74c3c',
          color: 'white',
          border: 'none',
          borderRadius: '4px',
          cursor: 'pointer'
        }}
      >
        Clear Lines
      </button>

      <div style={{ 
        marginTop: '15px',
        fontSize: '11px',
        opacity: 0.7,
        lineHeight: '1.4'
      }}>
        <p style={{ margin: '0 0 5px 0' }}>
          Lines: {debugLines.lineCount}
        </p>
        <p style={{ margin: '0 0 5px 0' }}>
          Frame: {debugLines.frameCount}
        </p>
        <p style={{ margin: 0 }}>
          Status: {debugLines.enabled ? '✅ Active' : '❌ Disabled'}
        </p>
      </div>
    </div>
  );
}

/**
 * Main line logger test component
 */
function LineLoggerTestInner() {
  const [useEnhanced] = useState(false);

  return (
    <div style={{ width: '100%', height: '100%', position: 'relative' }}>
      <Canvas
        camera={{ position: [5, 5, 5], fov: 60 }}
        style={{ background: 'linear-gradient(135deg, #0c0c0c 0%, #1a1a2e 100%)' }}
      >
        {/* Lighting */}
        <ambientLight intensity={0.3} />
        <pointLight position={[10, 10, 10]} intensity={0.5} />

        {/* Line animator */}
        <LineAnimator />

        {/* Line renderer */}
        {useEnhanced ? (
          <EnhancedLineRenderer 
            fadeOldLines={true}
            maxFrameAge={5}
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
        />
      </Canvas>

      {/* Controls */}
      <LineLoggerControls />
      
      {/* Debug info */}
      <LineRendererDebugInfo />
    </div>
  );
}

/**
 * Line Logger Test with provider wrapper
 */
export function LineLoggerTest() {
  return (
    <LineLoggerProvider enabled={true} maxLines={500}>
      <LineLoggerTestInner />
    </LineLoggerProvider>
  );
}
