import React, { createContext, useContext, useRef, useCallback, ReactNode } from 'react';
import * as THREE from 'three';

/**
 * Debug line data structure
 */
export interface DebugLine {
  id: string;
  start: THREE.Vector3;
  end: THREE.Vector3;
  color: string;
  width?: number;
  frame: number;
}

/**
 * Line logger interface matching Gaudi's geometry_logger style
 */
export interface LineLoggerInterface {
  addLine(start: THREE.Vector3, end: THREE.Vector3, color?: string, width?: number): void;
  addLine(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, color?: string, width?: number): void;
  getLines(): DebugLine[];
  flush(): void;
  clear(): void;
  getFrameCount(): number;
}

/**
 * Line logger context data
 */
interface LineLoggerContextData extends LineLoggerInterface {
  isEnabled: boolean;
  setEnabled: (enabled: boolean) => void;
}

const LineLoggerContext = createContext<LineLoggerContextData | null>(null);

/**
 * Props for LineLoggerProvider
 */
interface LineLoggerProviderProps {
  children: ReactNode;
  enabled?: boolean;
  maxLines?: number;
}

/**
 * React Context-based singleton line logger provider
 * Provides a centralized way to accumulate debug lines across components
 */
export function LineLoggerProvider({ 
  children, 
  enabled = true,
  maxLines = 1000 
}: LineLoggerProviderProps) {
  const linesRef = useRef<DebugLine[]>([]);
  const frameCountRef = useRef(0);
  const [isEnabled, setIsEnabled] = React.useState(enabled);
  const lineIdCounterRef = useRef(0);

  /**
   * Add a debug line (overloaded function)
   */
  const addLine = useCallback((
    startOrX: THREE.Vector3 | number,
    endOrY: THREE.Vector3 | number,
    zOrColor?: number | string,
    endXOrWidth?: number,
    endY?: number,
    endZ?: number,
    color: string = '#ffffff',
    width: number = 1
  ) => {
    if (!isEnabled) return;

    let start: THREE.Vector3;
    let end: THREE.Vector3;
    let finalColor: string;
    let finalWidth: number;

    // Handle overloaded parameters
    if (startOrX instanceof THREE.Vector3 && endOrY instanceof THREE.Vector3) {
      // Vector3 overload: addLine(start: Vector3, end: Vector3, color?, width?)
      start = startOrX.clone();
      end = endOrY.clone();
      finalColor = (zOrColor as string) || color;
      finalWidth = (endXOrWidth as number) || width;
    } else if (
      typeof startOrX === 'number' && 
      typeof endOrY === 'number' && 
      typeof zOrColor === 'number' &&
      typeof endXOrWidth === 'number' &&
      typeof endY === 'number' &&
      typeof endZ === 'number'
    ) {
      // Number overload: addLine(startX, startY, startZ, endX, endY, endZ, color?, width?)
      start = new THREE.Vector3(startOrX, endOrY, zOrColor);
      end = new THREE.Vector3(endXOrWidth, endY, endZ);
      finalColor = color;
      finalWidth = width;
    } else {
      console.error('LineLogger.addLine: Invalid parameters');
      return;
    }

    const line: DebugLine = {
      id: `line_${lineIdCounterRef.current++}`,
      start,
      end,
      color: finalColor,
      width: finalWidth,
      frame: frameCountRef.current
    };

    linesRef.current.push(line);

    // Limit number of lines to prevent memory issues
    if (linesRef.current.length > maxLines) {
      linesRef.current = linesRef.current.slice(-maxLines);
    }
  }, [isEnabled, maxLines]);

  /**
   * Get all accumulated lines
   */
  const getLines = useCallback((): DebugLine[] => {
    return [...linesRef.current];
  }, []);

  /**
   * Flush lines (increment frame count, keep lines for rendering)
   */
  const flush = useCallback((): void => {
    frameCountRef.current++;
    // Note: We don't clear lines here - they persist until manually cleared
    // This allows the renderer to display accumulated debug data
  }, []);

  /**
   * Clear all accumulated lines
   */
  const clear = useCallback((): void => {
    linesRef.current = [];
    lineIdCounterRef.current = 0;
  }, []);

  /**
   * Get current frame count
   */
  const getFrameCount = useCallback((): number => {
    return frameCountRef.current;
  }, []);

  const contextValue: LineLoggerContextData = {
    addLine: addLine as any, // Type assertion to handle overloads
    getLines,
    flush,
    clear,
    getFrameCount,
    isEnabled,
    setEnabled: setIsEnabled
  };

  return (
    <LineLoggerContext.Provider value={contextValue}>
      {children}
    </LineLoggerContext.Provider>
  );
}

/**
 * Hook to access the line logger context
 */
export function useLineLogger(): LineLoggerContextData {
  const context = useContext(LineLoggerContext);
  if (!context) {
    throw new Error('useLineLogger must be used within a LineLoggerProvider');
  }
  return context;
}

/**
 * Convenience hook for simple line logging (like Gaudi's global logger)
 */
export function useDebugLines() {
  const logger = useLineLogger();
  
  return {
    addLine: logger.addLine,
    clear: logger.clear,
    flush: logger.flush,
    lineCount: logger.getLines().length,
    frameCount: logger.getFrameCount(),
    enabled: logger.isEnabled,
    setEnabled: logger.setEnabled
  };
}
