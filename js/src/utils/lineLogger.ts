import * as THREE from 'three';

// Line data structure for debug rendering
export interface DebugLine {
  start: THREE.Vector3;
  end: THREE.Vector3;
  color: THREE.Color;
  id?: string;
  timestamp?: number;
}

// Color presets matching common debug scenarios
export const DebugColors = {
  RED: new THREE.Color(1, 0, 0),
  GREEN: new THREE.Color(0, 1, 0),
  BLUE: new THREE.Color(0, 0, 1),
  YELLOW: new THREE.Color(1, 1, 0),
  CYAN: new THREE.Color(0, 1, 1),
  MAGENTA: new THREE.Color(1, 0, 1),
  WHITE: new THREE.Color(1, 1, 1),
  GRAY: new THREE.Color(0.5, 0.5, 0.5),
  ORANGE: new THREE.Color(1, 0.5, 0),
  PURPLE: new THREE.Color(0.5, 0, 1),
} as const;

/**
 * Singleton line logger for accumulating debug visualization lines
 * Similar to Gaudi's geometry_logger but for JavaScript/React Three Fiber
 */
class LineLogger {
  private static instance: LineLogger | null = null;
  private lines: DebugLine[] = [];
  private maxLines: number = 10000;
  private frameCounter: number = 0;

  private constructor() {}

  static getInstance(): LineLogger {
    if (!LineLogger.instance) {
      LineLogger.instance = new LineLogger();
    }
    return LineLogger.instance;
  }

  /**
   * Add a line between two points with specified color
   */
  line(
    start: THREE.Vector3 | [number, number, number],
    end: THREE.Vector3 | [number, number, number],
    color: THREE.Color | string = DebugColors.WHITE,
    id?: string
  ): void {
    const startVec = Array.isArray(start) ? new THREE.Vector3(...start) : start.clone();
    const endVec = Array.isArray(end) ? new THREE.Vector3(...end) : end.clone();
    const colorObj = typeof color === 'string' ? new THREE.Color(color) : color.clone();

    const debugLine: DebugLine = {
      start: startVec,
      end: endVec,
      color: colorObj,
      id,
      timestamp: performance.now()
    };

    this.lines.push(debugLine);

    // Prevent memory bloat by limiting total lines
    if (this.lines.length > this.maxLines) {
      this.lines.shift();
    }
  }

  /**
   * Add a line from origin to a point (useful for vectors)
   */
  vector(
    end: THREE.Vector3 | [number, number, number],
    color: THREE.Color | string = DebugColors.GREEN,
    origin: THREE.Vector3 | [number, number, number] = [0, 0, 0],
    id?: string
  ): void {
    this.line(origin, end, color, id);
  }

  /**
   * Add a coordinate frame (X=red, Y=green, Z=blue axes)
   */
  frame(
    position: THREE.Vector3 | [number, number, number] = [0, 0, 0],
    scale: number = 1.0,
    id?: string
  ): void {
    const pos = Array.isArray(position) ? new THREE.Vector3(...position) : position;
    
    // X axis - Red
    this.line(
      pos,
      pos.clone().add(new THREE.Vector3(scale, 0, 0)),
      DebugColors.RED,
      id ? `${id}_x` : undefined
    );
    
    // Y axis - Green  
    this.line(
      pos,
      pos.clone().add(new THREE.Vector3(0, scale, 0)),
      DebugColors.GREEN,
      id ? `${id}_y` : undefined
    );
    
    // Z axis - Blue
    this.line(
      pos,
      pos.clone().add(new THREE.Vector3(0, 0, scale)),
      DebugColors.BLUE,
      id ? `${id}_z` : undefined
    );
  }

  /**
   * Add a wireframe box outline
   */
  box(
    center: THREE.Vector3 | [number, number, number],
    size: THREE.Vector3 | [number, number, number] | number,
    color: THREE.Color | string = DebugColors.WHITE,
    id?: string
  ): void {
    const centerVec = Array.isArray(center) ? new THREE.Vector3(...center) : center;
    const sizeVec = typeof size === 'number' 
      ? new THREE.Vector3(size, size, size)
      : Array.isArray(size) 
        ? new THREE.Vector3(...size)
        : size;

    const halfSize = sizeVec.clone().multiplyScalar(0.5);
    
    // Define the 8 corners of the box
    const corners = [
      centerVec.clone().add(new THREE.Vector3(-halfSize.x, -halfSize.y, -halfSize.z)),
      centerVec.clone().add(new THREE.Vector3( halfSize.x, -halfSize.y, -halfSize.z)),
      centerVec.clone().add(new THREE.Vector3( halfSize.x,  halfSize.y, -halfSize.z)),
      centerVec.clone().add(new THREE.Vector3(-halfSize.x,  halfSize.y, -halfSize.z)),
      centerVec.clone().add(new THREE.Vector3(-halfSize.x, -halfSize.y,  halfSize.z)),
      centerVec.clone().add(new THREE.Vector3( halfSize.x, -halfSize.y,  halfSize.z)),
      centerVec.clone().add(new THREE.Vector3( halfSize.x,  halfSize.y,  halfSize.z)),
      centerVec.clone().add(new THREE.Vector3(-halfSize.x,  halfSize.y,  halfSize.z))
    ];

    // Draw the 12 edges of the box
    const edges = [
      [0, 1], [1, 2], [2, 3], [3, 0], // Bottom face
      [4, 5], [5, 6], [6, 7], [7, 4], // Top face
      [0, 4], [1, 5], [2, 6], [3, 7]  // Vertical edges
    ];

    edges.forEach(([a, b], i) => {
      this.line(corners[a], corners[b], color, id ? `${id}_edge_${i}` : undefined);
    });
  }

  /**
   * Generate random test lines for validation
   */
  generateTestLines(count: number = 50): void {
    for (let i = 0; i < count; i++) {
      const start = new THREE.Vector3(
        (Math.random() - 0.5) * 4,
        (Math.random() - 0.5) * 4,
        (Math.random() - 0.5) * 4
      );
      const end = new THREE.Vector3(
        (Math.random() - 0.5) * 4,
        (Math.random() - 0.5) * 4,
        (Math.random() - 0.5) * 4
      );
      
      const colors = Object.values(DebugColors);
      const color = colors[Math.floor(Math.random() * colors.length)];
      
      this.line(start, end, color, `test_line_${i}`);
    }
  }

  /**
   * Get all accumulated lines
   */
  getLines(): DebugLine[] {
    return [...this.lines];
  }

  /**
   * Get lines by ID pattern
   */
  getLinesByPattern(pattern: string): DebugLine[] {
    return this.lines.filter(line => line.id?.includes(pattern));
  }

  /**
   * Clear all lines
   */
  clear(): void {
    this.lines = [];
  }

  /**
   * Clear lines by ID pattern
   */
  clearByPattern(pattern: string): void {
    this.lines = this.lines.filter(line => !line.id?.includes(pattern));
  }

  /**
   * Get statistics about current lines
   */
  getStats(): { count: number; memoryUsage: string; oldestTimestamp: number | null } {
    const count = this.lines.length;
    const memoryUsage = `${(count * 100 / 1024).toFixed(2)} KB`; // Rough estimate
    const oldestTimestamp = this.lines.length > 0 ? this.lines[0].timestamp || null : null;

    return { count, memoryUsage, oldestTimestamp };
  }

  /**
   * Set maximum number of lines to prevent memory issues
   */
  setMaxLines(max: number): void {
    this.maxLines = max;
    if (this.lines.length > max) {
      this.lines = this.lines.slice(-max);
    }
  }

  /**
   * Increment frame counter (useful for frame-based debugging)
   */
  nextFrame(): void {
    this.frameCounter++;
  }

  /**
   * Get current frame number
   */
  getFrame(): number {
    return this.frameCounter;
  }
}

// Export singleton instance
export const lineLogger = LineLogger.getInstance();

// Convenience functions that match Gaudi's geometry_logger style
export const debugLine = {
  line: (start: THREE.Vector3 | [number, number, number], end: THREE.Vector3 | [number, number, number], color?: THREE.Color | string, id?: string) => 
    lineLogger.line(start, end, color, id),
  
  vector: (end: THREE.Vector3 | [number, number, number], color?: THREE.Color | string, origin?: THREE.Vector3 | [number, number, number], id?: string) =>
    lineLogger.vector(end, color, origin, id),
    
  frame: (position?: THREE.Vector3 | [number, number, number], scale?: number, id?: string) =>
    lineLogger.frame(position, scale, id),
    
  box: (center: THREE.Vector3 | [number, number, number], size: THREE.Vector3 | [number, number, number] | number, color?: THREE.Color | string, id?: string) =>
    lineLogger.box(center, size, color, id),
    
  clear: () => lineLogger.clear(),
  
  getLines: () => lineLogger.getLines(),
  
  generateTest: (count?: number) => lineLogger.generateTestLines(count)
};
