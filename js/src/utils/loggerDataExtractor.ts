/**
 * Uniform interface for extracting logger data from WASM modules
 * This utility provides a consistent way to access line and point data
 * across different WASM modules that use the gaudi::logger interface
 * 
 * Now uses explicit API endpoints that copy data to JavaScript arrays
 * instead of opaque heap pointers for better clarity and safety
 */

export interface LoggerData {
  linePositions: Float32Array;
  lineColors: Float32Array;
  pointPositions: Float32Array;
  pointColors: Float32Array;
  lineCount: number;
  pointCount: number;
}

export interface WasmLoggerInstance {
  // Explicit data copying API - no more heap pointers
  get_line_count(): number;
  get_line_positions(): number[]; // Returns JavaScript array
  get_line_colors(): number[]; // Returns JavaScript array
  get_point_count(): number;
  get_point_positions(): number[]; // Returns JavaScript array
  get_point_colors(): number[]; // Returns JavaScript array
}

/**
 * Extract logger data from a WASM module instance using explicit API
 * @param wasmInstance The WASM module instance with logger API
 * @returns LoggerData with Float32Array ready for Three.js
 */
export function extractLoggerData(wasmInstance: WasmLoggerInstance): LoggerData {
  try {
    // Get line data using explicit API
    const lineCount = wasmInstance.get_line_count();
    const linePositionsArray = wasmInstance.get_line_positions();
    const lineColorsArray = wasmInstance.get_line_colors();
    
    // Get point data using explicit API
    const pointCount = wasmInstance.get_point_count();
    const pointPositionsArray = wasmInstance.get_point_positions();
    const pointColorsArray = wasmInstance.get_point_colors();
    
    // Convert to Float32Array for Three.js
    const linePositions = new Float32Array(linePositionsArray);
    const lineColors = new Float32Array(lineColorsArray);
    const pointPositions = new Float32Array(pointPositionsArray);
    const pointColors = new Float32Array(pointColorsArray);
    
    return {
      linePositions,
      lineColors,
      pointPositions,
      pointColors,
      lineCount,
      pointCount
    };
  } catch (error) {
    console.error('Error extracting logger data:', error);
    return {
      linePositions: new Float32Array(0),
      lineColors: new Float32Array(0),
      pointPositions: new Float32Array(0),
      pointColors: new Float32Array(0),
      lineCount: 0,
      pointCount: 0
    };
  }
} 