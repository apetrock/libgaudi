import { create } from 'zustand';
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
 * Line logger store state
 */
interface LineLoggerState {
  lines: DebugLine[];
  frameCount: number;
  enabled: boolean;
  maxLines: number;
  
  // Actions
  addLine: (start: THREE.Vector3, end: THREE.Vector3, color?: string, width?: number) => void;
  addLineCoords: (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, color?: string, width?: number) => void;
  flush: () => void;
  clear: () => void;
  setEnabled: (enabled: boolean) => void;
  setMaxLines: (maxLines: number) => void;
  
  // Getters
  getLines: () => DebugLine[];
  getLineCount: () => number;
  getFrameCount: () => number;
}

/**
 * Zustand store for debug line logging
 * Singleton pattern matching Gaudi's geometry_logger
 */
export const useLineLoggerStore = create<LineLoggerState>((set, get) => {
  let lineIdCounter = 0;

  return {
    // State
    lines: [],
    frameCount: 0,
    enabled: true,
    maxLines: 1000,

    // Add line with Vector3 parameters
    addLine: (start: THREE.Vector3, end: THREE.Vector3, color = '#ffffff', width = 1) => {
      const state = get();
      if (!state.enabled) return;

      const line: DebugLine = {
        id: `line_${lineIdCounter++}`,
        start: start.clone(),
        end: end.clone(),
        color,
        width,
        frame: state.frameCount,
        timestamp: performance.now()
      };

      set((state) => {
        const newLines = [...state.lines, line];
        
        // Limit number of lines to prevent memory issues
        if (newLines.length > state.maxLines) {
          return {
            ...state,
            lines: newLines.slice(-state.maxLines)
          };
        }
        
        return {
          ...state,
          lines: newLines
        };
      });
    },

    // Add line with coordinate parameters
    addLineCoords: (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, color = '#ffffff', width = 1) => {
      const start = new THREE.Vector3(startX, startY, startZ);
      const end = new THREE.Vector3(endX, endY, endZ);
      get().addLine(start, end, color, width);
    },

    // Flush lines (increment frame count)
    flush: () => {
      set((state) => ({
        ...state,
        frameCount: state.frameCount + 1
      }));
    },

    // Clear all lines
    clear: () => {
      set((state) => ({
        ...state,
        lines: []
      }));
      lineIdCounter = 0;
    },

    // Enable/disable logging
    setEnabled: (enabled: boolean) => {
      set((state) => ({
        ...state,
        enabled
      }));
    },

    // Set maximum number of lines
    setMaxLines: (maxLines: number) => {
      set((state) => {
        const newLines = state.lines.length > maxLines 
          ? state.lines.slice(-maxLines)
          : state.lines;
        
        return {
          ...state,
          maxLines,
          lines: newLines
        };
      });
    },

    // Get all lines
    getLines: () => get().lines,

    // Get line count
    getLineCount: () => get().lines.length,

    // Get frame count
    getFrameCount: () => get().frameCount
  };
});

/**
 * Convenience hook for debug line logging
 * Provides a simplified API similar to Gaudi's geometry_logger
 */
export function useDebugLines() {
  const addLine = useLineLoggerStore((state) => state.addLine);
  const addLineCoords = useLineLoggerStore((state) => state.addLineCoords);
  const flush = useLineLoggerStore((state) => state.flush);
  const clear = useLineLoggerStore((state) => state.clear);
  const setEnabled = useLineLoggerStore((state) => state.setEnabled);
  const setMaxLines = useLineLoggerStore((state) => state.setMaxLines);
  const lines = useLineLoggerStore((state) => state.lines);
  const enabled = useLineLoggerStore((state) => state.enabled);
  const maxLines = useLineLoggerStore((state) => state.maxLines);
  const lineCount = useLineLoggerStore((state) => state.lines.length);
  const frameCount = useLineLoggerStore((state) => state.frameCount);
  
  return {
    // Simple methods that directly call store
    addLine,
    addLineCoords,
    
    // Frame management
    flush,
    clear,
    
    // Configuration
    setEnabled,
    setMaxLines,
    
    // State access
    lines,
    lineCount,
    frameCount,
    enabled,
    maxLines
  };
}

/**
 * Hook for accessing just the lines (for rendering)
 * Optimized to only re-render when lines change
 */
export function useDebugLinesForRendering() {
  return useLineLoggerStore((state) => state.lines);
}

/**
 * Hook for accessing just the stats (for debug UI)
 * Optimized to only re-render when stats change
 */
export function useDebugLineStats() {
  const lineCount = useLineLoggerStore((state) => state.lines.length);
  const frameCount = useLineLoggerStore((state) => state.frameCount);
  const enabled = useLineLoggerStore((state) => state.enabled);
  const maxLines = useLineLoggerStore((state) => state.maxLines);
  
  return {
    lineCount,
    frameCount,
    enabled,
    maxLines
  };
}
