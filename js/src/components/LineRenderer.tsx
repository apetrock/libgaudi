import { useMemo, useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import * as THREE from 'three';
import { useDebugLinesForRendering, useDebugLineStats } from '../stores/lineLoggerStore';
import type { LineRendererConfig } from '../types/lineLogger';

/**
 * React Three Fiber component that renders accumulated debug lines from Zustand store
 */
export function LineRenderer() {
  const lines = useDebugLinesForRendering();
  const groupRef = useRef<THREE.Group>(null);

  // Create line objects for rendering
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

  return (
    <group ref={groupRef}>
      {lineObjects.map((lineObj) => (
        <primitive key={lineObj.id} object={new THREE.Line(lineObj.geometry, lineObj.material)} />
      ))}
    </group>
  );
}

/**
 * Enhanced line renderer with additional visual effects
 */
interface EnhancedLineRendererProps {
  fadeOldLines?: boolean;
  maxFrameAge?: number;
  showFrameColors?: boolean;
}

export function EnhancedLineRenderer({ 
  fadeOldLines = false,
  maxFrameAge = 10,
  showFrameColors = false
}: EnhancedLineRendererProps) {
  const lineLogger = useLineLogger();
  const groupRef = useRef<THREE.Group>(null);
    const lines = lineLogger.getLines();
  const currentFrame = lineLogger.getFrameCount();

  // Create enhanced line objects with effects
  const enhancedLineObjects = useMemo(() => {
    return lines.map((line) => {
      const frameAge = currentFrame - line.frame;
      
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
        const hue = (frameAge * 30) % 360; // Cycle through hues
        color = `hsl(${hue}, 70%, 60%)`;
      }

      // Create material with effects
      const material = new THREE.LineBasicMaterial({
        color,
        linewidth: line.width || 1,
        transparent: opacity < 1.0,
        opacity
      });

      return {
        id: line.id,
        geometry,
        material,
        frame: line.frame
      };
    }).filter(Boolean); // Remove null entries
  }, [lines, currentFrame, fadeOldLines, maxFrameAge, showFrameColors]);

  return (
    <group ref={groupRef}>
      {enhancedLineObjects.map((lineObj) => 
        lineObj && (
          <primitive key={lineObj.id} object={new THREE.Line(lineObj.geometry, lineObj.material)} />
        )
      )}
    </group>
  );
}

/**
 * Debug info overlay for line rendering
 */
export function LineRendererDebugInfo() {
  const lineLogger = useLineLogger();
  const lines = lineLogger.getLines();
  
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
      <div>Lines: {lines.length}</div>
      <div>Frame: {lineLogger.getFrameCount()}</div>
      <div>Enabled: {lineLogger.isEnabled ? '✅' : '❌'}</div>
    </div>
  );
}
