#!/usr/bin/env ts-node

import { execSync } from 'child_process';
import { existsSync, mkdirSync, readdirSync, copyFileSync } from 'fs';
import { join, resolve } from 'path';

interface BuildConfig {
  name: string;
  description: string;
  outputFiles: string[];
}

const WASM_MODULES: BuildConfig[] = [
  {
    name: 'gaudi_logger_test',
    description: 'Logger API',
    outputFiles: ['gaudi_logger_test.js', 'gaudi_logger_test.wasm']
  },
  {
    name: 'rod_constraints_test', 
    description: 'Rod Simulation',
    outputFiles: ['rod_constraints_test.js', 'rod_constraints_test.wasm']
  },
  {
    name: 'hello_world',
    description: 'Hello World Example',
    outputFiles: ['hello_world.js', 'hello_world.wasm']
  },
  {
    name: 'foo_demo',
    description: 'Foo Demo',
    outputFiles: ['foo_demo.js', 'foo_demo.wasm']
  }
];

function checkEmscripten(): string {
  try {
    const version = execSync('emcc --version', { encoding: 'utf8' });
    const versionLine = version.split('\n')[0];
    console.log(`✅ Found Emscripten: ${versionLine}`);
    return versionLine;
  } catch (error) {
    console.error('❌ Error: emcc not found. Please install Emscripten SDK.');
    console.error('\nInstallation instructions:');
    console.error('1. git clone https://github.com/emscripten-core/emsdk.git');
    console.error('2. cd emsdk');
    console.error('3. emsdk install latest');
    console.error('4. emsdk activate latest');
    console.error('5. emsdk_env.bat');
    process.exit(1);
  }
}

function runCommand(command: string, cwd?: string): void {
  console.log(`Running: ${command}`);
  try {
    execSync(command, { 
      cwd, 
      stdio: 'inherit',
      encoding: 'utf8'
    });
  } catch (error) {
    console.error(`❌ Command failed: ${command}`);
    process.exit(1);
  }
}

function copyFilesToPublic(): void {
  const buildDir = resolve(__dirname, '../wasm/build');
  const publicDir = resolve(__dirname, '../public/wasm');
  
  // Ensure public directory exists
  if (!existsSync(publicDir)) {
    mkdirSync(publicDir, { recursive: true });
  }
  
  console.log('📁 Copying files to public directory...');
  
  // Copy all .js and .wasm files from build to public
  const buildFiles = readdirSync(buildDir);
  const wasmFiles = buildFiles.filter(file => 
    file.endsWith('.js') || file.endsWith('.wasm')
  );
  
  wasmFiles.forEach(file => {
    const sourcePath = join(buildDir, file);
    const destPath = join(publicDir, file);
    copyFileSync(sourcePath, destPath);
    console.log(`  📄 ${file}`);
  });
  
  console.log(`✅ Copied ${wasmFiles.length} files to ${publicDir}`);
}

function main(): void {
  console.log('🔨 Building All WebAssembly Modules...\n');
  
  // Check Emscripten
  checkEmscripten();
  
  const wasmDir = resolve(__dirname, '../wasm');
  const buildDir = join(wasmDir, 'build');
  
  // Create build directory
  if (!existsSync(buildDir)) {
    mkdirSync(buildDir, { recursive: true });
  }
  
  // Configure with Emscripten
  console.log('⚙️  Configuring CMake with Emscripten...');
  runCommand('emcmake cmake .. -DCMAKE_BUILD_TYPE=Debug', buildDir);
  
  // Build all WASM targets
  console.log('🔨 Building all WebAssembly targets...');
  runCommand('cmake --build . --target all_wasm_targets', buildDir);
  
  // Copy files to public
  copyFilesToPublic();
  
  // Summary
  console.log('\n✅ All builds successful!');
  console.log('\nGenerated WASM modules:');
  WASM_MODULES.forEach(module => {
    console.log(`  - ${module.description} (${module.outputFiles.join(', ')})`);
  });
  console.log('\nFiles copied to: public/wasm/');
}

if (require.main === module) {
  main();
} 