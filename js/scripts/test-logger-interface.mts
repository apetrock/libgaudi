#!/usr/bin/env ts-node

import { loadWasmModule } from '../src/utils/wasmLoader';
import { LineLoggerAPI } from '../src/types/lineLogger';

interface TestLoggerInstance {
  // Logger API methods
  get_line_count(): number;
  get_line_positions(): number[];
  get_line_colors(): number[];
  get_point_count(): number;
  get_point_positions(): number[];
  get_point_colors(): number[];
  
  // Test methods
  add_line(x0: number, y0: number, z0: number, x1: number, y1: number, z1: number, 
           r: number, g: number, b: number, a: number): void;
  add_point(x: number, y: number, z: number, r: number, g: number, b: number, a: number): void;
  clear_logger(): void;
}

interface TestLoggerModule {
  GaudiLoggerTest: new () => TestLoggerInstance;
}

async function testLoggerInterface(): Promise<void> {
  console.log('🧪 Testing Logger Interface...\n');
  
  try {
    // Load the WASM module
    console.log('📦 Loading WASM module...');
    const wasmModule = await loadWasmModule('gaudi_logger_test') as TestLoggerModule;
    const loggerInstance = new wasmModule.GaudiLoggerTest();
    console.log('✅ WASM module loaded successfully\n');
    
    // Test 1: Initial state
    console.log('🔍 Test 1: Initial state');
    const initialLineCount = loggerInstance.get_line_count();
    const initialPointCount = loggerInstance.get_point_count();
    console.log(`  Initial line count: ${initialLineCount}`);
    console.log(`  Initial point count: ${initialPointCount}`);
    
    // Test 2: Add known test lines
    console.log('\n🔍 Test 2: Adding known test lines');
    loggerInstance.clear_logger(); // Start fresh
    
    // Add a red line from origin to (1,0,0)
    loggerInstance.add_line(0, 0, 0, 1, 0, 0, 1.0, 0.0, 0.0, 1.0);
    console.log('  Added red line: (0,0,0) → (1,0,0)');
    
    // Add a green line from origin to (0,1,0)
    loggerInstance.add_line(0, 0, 0, 0, 1, 0, 0.0, 1.0, 0.0, 1.0);
    console.log('  Added green line: (0,0,0) → (0,1,0)');
    
    // Add a blue line from origin to (0,0,1)
    loggerInstance.add_line(0, 0, 0, 0, 0, 1, 0.0, 0.0, 1.0, 1.0);
    console.log('  Added blue line: (0,0,0) → (0,0,1)');
    
    // Test 3: Retrieve and verify line data
    console.log('\n🔍 Test 3: Retrieving line data');
    const lineCount = loggerInstance.get_line_count();
    const linePositions = loggerInstance.get_line_positions();
    const lineColors = loggerInstance.get_line_colors();
    
    console.log(`  Retrieved line count: ${lineCount}`);
    console.log(`  Retrieved positions array length: ${linePositions.length}`);
    console.log(`  Retrieved colors array length: ${lineColors.length}`);
    
    // Verify we have 3 lines (6 vertices total)
    if (lineCount !== 3) {
      throw new Error(`Expected 3 lines, got ${lineCount}`);
    }
    
    // Verify positions array has 18 elements (6 vertices × 3 coordinates)
    if (linePositions.length !== 18) {
      throw new Error(`Expected 18 position elements, got ${linePositions.length}`);
    }
    
    // Verify colors array has 24 elements (6 vertices × 4 RGBA values)
    if (lineColors.length !== 24) {
      throw new Error(`Expected 24 color elements, got ${lineColors.length}`);
    }
    
    // Test 4: Verify specific line data
    console.log('\n🔍 Test 4: Verifying specific line data');
    
    // Check first line (red line from origin to (1,0,0))
    const line1Start = [linePositions[0], linePositions[1], linePositions[2]]; // (0,0,0)
    const line1End = [linePositions[3], linePositions[4], linePositions[5]];   // (1,0,0)
    const line1Color = [lineColors[0], lineColors[1], lineColors[2], lineColors[3]]; // Red
    
    console.log(`  Line 1 start: [${line1Start.join(', ')}]`);
    console.log(`  Line 1 end: [${line1End.join(', ')}]`);
    console.log(`  Line 1 color: [${line1Color.join(', ')}]`);
    
    // Verify first line
    if (line1Start[0] !== 0 || line1Start[1] !== 0 || line1Start[2] !== 0) {
      throw new Error(`Line 1 start should be (0,0,0), got (${line1Start.join(',')})`);
    }
    if (line1End[0] !== 1 || line1End[1] !== 0 || line1End[2] !== 0) {
      throw new Error(`Line 1 end should be (1,0,0), got (${line1End.join(',')})`);
    }
    if (line1Color[0] !== 1 || line1Color[1] !== 0 || line1Color[2] !== 0 || line1Color[3] !== 1) {
      throw new Error(`Line 1 color should be red [1,0,0,1], got [${line1Color.join(',')}]`);
    }
    
    // Test 5: Add and verify points
    console.log('\n🔍 Test 5: Adding and verifying points');
    loggerInstance.clear_logger();
    
    // Add a yellow point at (0.5, 0.5, 0.5)
    loggerInstance.add_point(0.5, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0);
    console.log('  Added yellow point at (0.5, 0.5, 0.5)');
    
    // Add a magenta point at (-0.5, -0.5, -0.5)
    loggerInstance.add_point(-0.5, -0.5, -0.5, 1.0, 0.0, 1.0, 1.0);
    console.log('  Added magenta point at (-0.5, -0.5, -0.5)');
    
    const pointCount = loggerInstance.get_point_count();
    const pointPositions = loggerInstance.get_point_positions();
    const pointColors = loggerInstance.get_point_colors();
    
    console.log(`  Retrieved point count: ${pointCount}`);
    console.log(`  Retrieved point positions length: ${pointPositions.length}`);
    console.log(`  Retrieved point colors length: ${pointColors.length}`);
    
    // Verify we have 2 points
    if (pointCount !== 2) {
      throw new Error(`Expected 2 points, got ${pointCount}`);
    }
    
    // Verify positions array has 6 elements (2 points × 3 coordinates)
    if (pointPositions.length !== 6) {
      throw new Error(`Expected 6 position elements, got ${pointPositions.length}`);
    }
    
    // Verify colors array has 8 elements (2 points × 4 RGBA values)
    if (pointColors.length !== 8) {
      throw new Error(`Expected 8 color elements, got ${pointColors.length}`);
    }
    
    // Verify first point (yellow at 0.5, 0.5, 0.5)
    const point1Pos = [pointPositions[0], pointPositions[1], pointPositions[2]];
    const point1Color = [pointColors[0], pointColors[1], pointColors[2], pointColors[3]];
    
    console.log(`  Point 1 position: [${point1Pos.join(', ')}]`);
    console.log(`  Point 1 color: [${point1Color.join(', ')}]`);
    
    if (point1Pos[0] !== 0.5 || point1Pos[1] !== 0.5 || point1Pos[2] !== 0.5) {
      throw new Error(`Point 1 should be at (0.5,0.5,0.5), got (${point1Pos.join(',')})`);
    }
    if (point1Color[0] !== 1 || point1Color[1] !== 1 || point1Color[2] !== 0 || point1Color[3] !== 1) {
      throw new Error(`Point 1 color should be yellow [1,1,0,1], got [${point1Color.join(',')}]`);
    }
    
    console.log('\n✅ All logger interface tests passed!');
    console.log('\n📊 Summary:');
    console.log('  - Line data retrieval: ✅');
    console.log('  - Point data retrieval: ✅');
    console.log('  - Array conversion: ✅');
    console.log('  - Data integrity: ✅');
    
  } catch (error) {
    console.error('\n❌ Logger interface test failed:');
    console.error(error);
    process.exit(1);
  }
}

// Run the test if this file is executed directly
if (import.meta.url === `file://${process.argv[1]}`) {
  testLoggerInterface();
}

// Run the test
testLoggerInterface(); 