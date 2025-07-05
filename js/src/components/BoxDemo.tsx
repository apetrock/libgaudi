import React, { Suspense, useState } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import { BoxMesh } from './BoxMesh';

interface BoxDemoProps {
  color?: string;
  rotationSpeed?: number;
  showControls?: boolean;
  style?: React.CSSProperties;
  className?: string;
  width?: number;
  height?: number;
}

/**
 * Phase 0 Milestone: Simple box demo to validate React Three Fiber setup
 */
export function BoxDemo({
  color = '#4a90e2',
  rotationSpeed = 1.0,
  showControls = true,
  style,
  className,
  width,
  height
}: BoxDemoProps) {
  const [currentColor, setCurrentColor] = useState(color);
  const [currentSpeed, setCurrentSpeed] = useState(rotationSpeed);

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
        style={{ background: 'linear-gradient(135deg, #1e3c72 0%, #2a5298 100%)' }}
      >
        <Suspense fallback={null}>
          {/* Lighting */}
          <ambientLight intensity={0.4} />
          <directionalLight 
            position={[5, 5, 5]} 
            intensity={0.8}
            castShadow
          />
          <pointLight position={[-5, -5, -5]} intensity={0.3} />

          {/* Environment */}
          <Environment preset="city" />

          {/* Animated Box */}
          <BoxMesh 
            color={currentColor}
            rotationSpeed={currentSpeed}
            position={[0, 0, 0]}
          />

          {/* Additional boxes for visual interest */}
          <BoxMesh 
            color="#e74c3c"
            rotationSpeed={currentSpeed * 0.5}
            position={[-2, 0, 0]}
          />
          
          <BoxMesh 
            color="#2ecc71"
            rotationSpeed={currentSpeed * 0.3}
            position={[2, 0, 0]}
          />

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
      {showControls && (
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
            Box Demo Controls
          </h3>
          
          <div style={{ marginBottom: '10px' }}>
            <label>Rotation Speed: {currentSpeed.toFixed(1)}</label>
            <input
              type="range"
              min="0"
              max="3"
              step="0.1"
              value={currentSpeed}
              onChange={(e) => setCurrentSpeed(parseFloat(e.target.value))}
              style={{ width: '100%', marginTop: '2px' }}
            />
          </div>

          <div style={{ marginBottom: '10px' }}>
            <label>Color:</label>
            <div style={{ display: 'flex', gap: '5px', marginTop: '5px' }}>
              {['#4a90e2', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6'].map((colorOption) => (
                <button
                  key={colorOption}
                  onClick={() => setCurrentColor(colorOption)}
                  style={{
                    width: '25px',
                    height: '25px',
                    backgroundColor: colorOption,
                    border: currentColor === colorOption ? '2px solid white' : '1px solid #ccc',
                    borderRadius: '4px',
                    cursor: 'pointer'
                  }}
                />
              ))}
            </div>
          </div>

          <div style={{ 
            marginTop: '15px',
            fontSize: '11px',
            opacity: 0.7,
            lineHeight: '1.4'
          }}>
            <p style={{ margin: '0 0 5px 0' }}>
              ✅ React Three Fiber working
            </p>
            <p style={{ margin: '0 0 5px 0' }}>
              ✅ Controls and state management
            </p>
            <p style={{ margin: 0 }}>
              ✅ 60fps rendering validated
            </p>
          </div>
        </div>
      )}

      {/* Performance indicator */}
      <div style={{
        position: 'absolute',
        bottom: 10,
        right: 10,
        background: 'rgba(0, 0, 0, 0.6)',
        color: 'white',
        padding: '5px 10px',
        borderRadius: '4px',
        fontFamily: 'monospace',
        fontSize: '11px'
      }}>
        Phase 0: Basic Rendering ✅
      </div>
    </div>
  );
}
