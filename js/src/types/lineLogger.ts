import * as THREE from 'three';

/**
 * Debug line data structure
 */
export interface DebugLine {
  id: string;
  start: THREE.Vector3;
  end: THREE.Vector3;
  color: string;
  width: number;
  frame: number;
  timestamp: number;
}

/**
 * Line logger API interface (matches Gaudi's geometry_logger style)
 */
export interface LineLoggerAPI {
  // Line creation methods
  addLine(start: THREE.Vector3, end: THREE.Vector3, color?: string, width?: number): void;
  addLine(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, color?: string, width?: number): void;
  
  // Frame management
  flush(): void;
  clear(): void;
  
  // Configuration
  setEnabled(enabled: boolean): void;
  setMaxLines(maxLines: number): void;
  
  // State access
  getLines(): DebugLine[];
  getLineCount(): number;
  getFrameCount(): number;
  
  // Properties
  enabled: boolean;
  maxLines: number;
}

/**
 * Line renderer configuration
 */
export interface LineRendererConfig {
  fadeOldLines?: boolean;
  maxFrameAge?: number;
  showFrameColors?: boolean;
  lineWidth?: number;
  opacity?: number;
}

/**
 * Debug line statistics
 */
export interface LineLoggerStats {
  lineCount: number;
  frameCount: number;
  enabled: boolean;
  maxLines: number;
  oldestFrame?: number;
  newestFrame?: number;
  memoryUsage?: number; // estimated bytes
}

/**
 * Line animation pattern types
 */
export type LinePattern = 'spiral' | 'cube' | 'random' | 'grid' | 'wave' | 'custom';

/**
 * Line pattern generator function type
 */
export type LinePatternGenerator = (time: number, frameCount: number, api: LineLoggerAPI) => void;
