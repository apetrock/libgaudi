import { useMemo, useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';
import { useDebugLinesForRendering, useDebugLineStats, useLineLoggerStore } from '../stores/lineLoggerStore';

/**
 * React Three Fiber component that renders debug lines from Zustand store
 */
export function LineRenderer() {
  const groupRef = useRef<THREE.Group>(null);
  
  // Subscribe only to lines for optimal performance
  const lines = useDebugLinesForRendering();

  // Create line geometries and materials
  const lineObjects = useMemo(() => {
    return lines.map((line) => {
      // Create geometry for the line
      const geometry = new THREE.BufferGeometry();
      const points = [line.start, line.end];
      geometry.setFromPoints(points);

      // Create material with line color and width
      const material = new THREE.LineBasicMaterial({
        color: line.color,
        linewidth: line.width
      });

      return {
        id: line.id,
        geometry,
        material,
        frame: line.frame
      };
    });
  }, [lines]);

  // Optional: Add frame-based updates
  useFrame(() => {
    // Could add animation effects here if needed
  });

  return (
    <group ref={groupRef}>
      {lineObjects.map((lineObj) => (
        <primitive 
          key={lineObj.id} 
          object={new THREE.Line(lineObj.geometry, lineObj.material)} 
        />
      ))}
    </group>
  );
}

/**
 * Enhanced line renderer with visual effects and filtering
 */
interface EnhancedLineRendererProps {
  fadeOldLines?: boolean;
  maxFrameAge?: number;
  showFrameColors?: boolean;
  wireframe?: boolean;
}

export function EnhancedLineRenderer({ 
  fadeOldLines = false,
  maxFrameAge = 10,
  showFrameColors = false,
  wireframe = false
}: EnhancedLineRendererProps) {
  const groupRef = useRef<THREE.Group>(null);
  
  // Subscribe to both lines and frame count
  const lines = useDebugLinesForRendering();
  const { frameCount } = useDebugLineStats();

  // Create enhanced line objects with effects
  const enhancedLineObjects = useMemo(() => {
    return lines.map((line) => {
      const frameAge = frameCount - line.frame;
      
      // Skip old lines if fade is enabled
      if (fadeOldLines && frameAge > maxFrameAge) {
        return null;
      }

      // Create geometry
      const geometry = new THREE.BufferGeometry();
      const points = [line.start, line.end];
      geometry.setFromPoints(points);

      // Calculate opacity based on age
      let opacity = 1.0;
      if (fadeOldLines && frameAge > 0) {
        opacity = Math.max(0.1, 1.0 - (frameAge / maxFrameAge));
      }

      // Determine color
      let color = line.color;
      if (showFrameColors) {
        // Color-code by frame age
        const hue = (frameAge * 30) % 360;
        color = `hsl(${hue}, 70%, 60%)`;
      }      // Create material with effects
      const material = new THREE.LineBasicMaterial({
        color,
        linewidth: line.width,
        transparent: opacity < 1.0,
        opacity
      });

      return {
        id: line.id,
        geometry,
        material,
        frame: line.frame,
        frameAge
      };
    }).filter(Boolean);
  }, [lines, frameCount, fadeOldLines, maxFrameAge, showFrameColors, wireframe]);

  return (
    <group ref={groupRef}>
      {enhancedLineObjects.map((lineObj) => 
        lineObj && (
          <primitive 
            key={lineObj.id} 
            object={new THREE.Line(lineObj.geometry, lineObj.material)} 
          />
        )
      )}
    </group>
  );
}

/**
 * Debug info overlay showing line logger statistics
 */
export function LineRendererDebugInfo() {
  const stats = useDebugLineStats();
  
  return (
    <div style={{
      position: 'absolute',
      top: 10,
      right: 10,
      background: 'rgba(0, 0, 0, 0.8)',
      color: 'white',
      padding: '10px',
      borderRadius: '4px',
      fontFamily: 'monospace',
      fontSize: '12px',
      zIndex: 1000
    }}>
      <div>Lines: {stats.lineCount}</div>
      <div>Frame: {stats.frameCount}</div>
      <div>Max: {stats.maxLines}</div>
      <div>Status: {stats.enabled ? '✅' : '❌'}</div>
    </div>
  );
}

/**
 * Hook for flushing lines after rendering (to be called in useFrame)
 */
export function useFlushAfterRender() {
  // Import the store directly to access flush function
  const { flush } = useLineLoggerStore();
  
  useFrame(() => {
    // Flush after each frame to increment frame counter
    flush();
  });
}
