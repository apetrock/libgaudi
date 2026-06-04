#!/usr/bin/env tsx
/**
 * WASM Unit Test Runner
 * Loads the unit_tests WASM module and runs all gaudi unit tests
 * 
 * Usage: pnpm run wasm:unit-test
 * 
 * The test output comes from two sources:
 * 1. console_logger in test.hpp -> wasm_terminal_logger -> stdout (or JS callbacks)
 * 2. The returned summary string from runAllTests()
 */

import { fileURLToPath, pathToFileURL } from 'url';
import { dirname, join } from 'path';
import { readFileSync, existsSync } from 'fs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Path to built WASM module (SINGLE_FILE mode - no separate .wasm)
const wasmDir = join(__dirname, '..', 'public', 'wasm');
const moduleJs = join(wasmDir, 'unit_tests.js');
// Convert to file:// URL for dynamic import (required on Windows)
const moduleJsUrl = pathToFileURL(moduleJs).href;

// ANSI color codes for prettier output
const colors = {
  reset: '\x1b[0m',
  bright: '\x1b[1m',
  dim: '\x1b[2m',
  red: '\x1b[31m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  cyan: '\x1b[36m',
};

async function main() {
  console.log(`${colors.bright}${colors.cyan}🧪 GAUDI WASM Unit Test Runner${colors.reset}\n`);
  
  // Check if module is built
  if (!existsSync(moduleJs)) {
    console.error(`${colors.red}❌ unit_tests.js not found at: ${moduleJs}${colors.reset}`);
    console.error(`${colors.dim}   Build with: pnpm run wasm:build${colors.reset}`);
    process.exit(1);
  }
  
  console.log(`${colors.dim}📦 Loading module from: ${wasmDir}${colors.reset}\n`);
  
  try {
    // Dynamic import of the Emscripten-generated module
    const createModule = await import(moduleJsUrl);
    const Module = await createModule.default();
    
    // Set up logging callbacks to capture C++ console_logger output
    if (typeof Module.set_info_callback === 'function') {
      Module.set_info_callback((msg: string) => {
        console.log(`${colors.blue}[INFO]${colors.reset} ${msg}`);
      });
    }
    if (typeof Module.set_error_callback === 'function') {
      Module.set_error_callback((msg: string) => {
        console.log(`${colors.red}[ERROR]${colors.reset} ${msg}`);
      });
    }
    if (typeof Module.set_warning_callback === 'function') {
      Module.set_warning_callback((msg: string) => {
        console.log(`${colors.yellow}[WARN]${colors.reset} ${msg}`);
      });
    }
    if (typeof Module.set_debug_callback === 'function') {
      Module.set_debug_callback((msg: string) => {
        console.log(`${colors.dim}[DEBUG]${colors.reset} ${msg}`);
      });
    }
    
    // Try to load sphere.obj for BVH tests (optional)
    const sphereObjPath = join(__dirname, '..', 'public', 'assets', 'models', 'sphere.obj');
    if (existsSync(sphereObjPath)) {
      console.log(`${colors.dim}📁 Loading sphere.obj for BVH tests...${colors.reset}\n`);
      const sphereData = readFileSync(sphereObjPath, 'utf-8');
      Module.setSphereObj(sphereData);
    } else {
      // Also try the assets folder at project root
      const rootSphereObjPath = join(__dirname, '..', '..', 'assets', 'models', 'sphere.obj');
      if (existsSync(rootSphereObjPath)) {
        console.log(`${colors.dim}📁 Loading sphere.obj from assets...${colors.reset}\n`);
        const sphereData = readFileSync(rootSphereObjPath, 'utf-8');
        Module.setSphereObj(sphereData);
      } else {
        console.log(`${colors.yellow}⚠️  sphere.obj not found - BVH mesh tests may fail${colors.reset}\n`);
      }
    }
    
    // Run all tests
    console.log(`${colors.bright}Running tests...${colors.reset}\n`);
    console.log('─'.repeat(60));
    
    const result = Module.runAllTests();
    
    console.log('─'.repeat(60));
    console.log(`\n${colors.bright}Summary:${colors.reset}`);
    console.log(result);
    
    // Parse result to determine exit code
    if (result.includes('ALL TESTS PASSED')) {
      console.log(`\n${colors.green}${colors.bright}✅ All tests passed!${colors.reset}`);
      process.exit(0);
    } else {
      console.log(`\n${colors.red}${colors.bright}❌ Some tests failed!${colors.reset}`);
      process.exit(1);
    }
    
  } catch (error) {
    console.error(`${colors.red}❌ Failed to load or run tests:${colors.reset}`, error);
    process.exit(1);
  }
}

main();
