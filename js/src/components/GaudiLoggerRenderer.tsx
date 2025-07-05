import React, { useRef } from 'react';
import { useFrame, useThree } from '@react-three/fiber';
import * as THREE from 'three';

export interface WasmLoggerAPI {
  // The WASM instance doesn't need to have logger methods
  // Logger functions are called globally
}

export interface WasmModule {
  get_line_count(): number;
  get_line_colors_data(): Float32Array;
  get_lines_data(): Float32Array;
  get_point_count(): number;
  get_point_colors_data(): Float32Array;
  get_points_data(): Float32Array;
}

interface GaudiLoggerRendererProps {
  wasmInstance: WasmLoggerAPI | null;
  wasmModule: WasmModule | null;
  isPlaying?: boolean;
  showLines?: boolean;
  showPoints?: boolean;
  lineWidth?: number;
  pointSize?: number;
  lineMaterial?: THREE.LineBasicMaterialParameters;
  pointMaterial?: THREE.PointsMaterialParameters;
}

export function GaudiLoggerRenderer({
  wasmInstance,
  wasmModule,
  isPlaying = true,
  showLines = true,
  showPoints = true,
  lineWidth = 2,
  pointSize = 0.1,
  lineMaterial = { vertexColors: true },
  pointMaterial = { vertexColors: true, size: 0.1 }
}: GaudiLoggerRendererProps) {
  const linesRef = useRef<THREE.LineSegments>(null);
  const pointsRef = useRef<THREE.Points>(null);
  const { camera } = useThree();

  useFrame(() => {
    if (!wasmInstance || !wasmModule || !isPlaying) return;

    try {
      // Update Lines - call functions through the WASM module
      if (linesRef.current && showLines) {
        const lineCount = wasmModule.get_line_count();
        if (lineCount > 0) {
          const linePositions = wasmModule.get_lines_data();
          const lineColors = wasmModule.get_line_colors_data();
          
          linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(linePositions, 3));
          linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(lineColors, 4)); // RGBA
          linesRef.current.geometry.attributes.position.needsUpdate = true;
          linesRef.current.geometry.attributes.color.needsUpdate = true;
        } else {
          // Clear geometry if no lines
          linesRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
          linesRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 4));
        }
      }

      // Update Points - call functions through the WASM module
      if (pointsRef.current && showPoints) {
        const pointCount = wasmModule.get_point_count();
        if (pointCount > 0) {
          const pointPositions = wasmModule.get_points_data();
          const pointColors = wasmModule.get_point_colors_data();
          
          pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(pointPositions, 3));
          pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(pointColors, 4)); // RGBA
          pointsRef.current.geometry.attributes.position.needsUpdate = true;
          pointsRef.current.geometry.attributes.color.needsUpdate = true;
          
          // Adjust point size based on camera distance to prevent disappearing
          const points = pointsRef.current.geometry.attributes.position.array as Float32Array;
          if (points.length > 0) {
            // Calculate average distance to camera
            let totalDistance = 0;
            for (let i = 0; i < points.length; i += 3) {
              const point = new THREE.Vector3(points[i], points[i + 1], points[i + 2]);
              totalDistance += camera.position.distanceTo(point);
            }
            const avgDistance = totalDistance / (points.length / 3);
            
            // Scale point size based on distance (closer = larger)
            const distanceScale = Math.max(0.1, Math.min(2.0, avgDistance * 0.1));
            const adjustedPointSize = pointSize * distanceScale;
            
            if (pointsRef.current.material instanceof THREE.PointsMaterial) {
              pointsRef.current.material.size = adjustedPointSize;
              pointsRef.current.material.sizeAttenuation = true;
            }
          }
        } else {
          // Clear geometry if no points
          pointsRef.current.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(0), 3));
          pointsRef.current.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(0), 4));
        }
      }
    } catch (error) {
      console.error('Error updating GaudiLoggerRenderer:', error);
    }
  });

  return (
    <>
      {showLines && (
        <lineSegments ref={linesRef}>
          <bufferGeometry />
          <lineBasicMaterial 
            {...lineMaterial} 
            linewidth={lineWidth}
            depthTest={true}
            depthWrite={true}
            transparent={true}
            alphaTest={0.1}
          />
        </lineSegments>
      )}
      
      {showPoints && (
        <points ref={pointsRef}>
          <bufferGeometry />
          <pointsMaterial 
            {...pointMaterial} 
            size={pointSize}
            sizeAttenuation={true}
            depthTest={true}
            depthWrite={true}
            transparent={true}
            alphaTest={0.1}
          />
        </points>
      )}
    </>
  );
} 