#!/usr/bin/env ts-node

import { execSync } from 'child_process';
import { existsSync, mkdirSync, readdirSync, copyFileSync, writeFileSync, readFileSync } from 'fs';
import { join, resolve, normalize } from 'path';
import { fileURLToPath } from 'url';
import { cpus } from 'os';
import { MODULE_CONFIGS, getAllModules, ModuleConfig } from './module-config.js';

interface BuildConfig {
  name: string;
  description: string;
  outputFiles: string[];
}

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

function validateModuleStructure(scriptDir: string): void {
  console.log('🔍 Validating module structure...');
  
  const wasmDir = resolve(scriptDir, '../wasm');
  const modulesDir = join(wasmDir, 'modules');
  
  // Check if modules directory exists
  if (!existsSync(modulesDir)) {
    console.log('📁 Creating modules directory...');
    mkdirSync(modulesDir, { recursive: true });
  }
  
  // Validate each module exists
  const modules = getAllModules();
  for (const module of modules) {
    const moduleDir = join(modulesDir, module.name);
    if (!existsSync(moduleDir)) {
      console.error(`❌ Module directory not found: ${moduleDir}`);
      console.error(`   Please create the module directory and move your source files there.`);
      process.exit(1);
    }
    
    // Check if source files exist
    for (const source of module.sources) {
      const sourcePath = join(moduleDir, source);
      if (!existsSync(sourcePath)) {
        console.error(`❌ Source file not found: ${sourcePath}`);
        console.error(`   Expected source file for module ${module.name}`);
        process.exit(1);
      }
    }
  }
  
  console.log('✅ Module structure validation passed');
}

function generateModuleCMakeLists(scriptDir: string, module: ModuleConfig): void {
  const wasmDir = resolve(scriptDir, '../wasm');
  const moduleDir = join(wasmDir, 'modules', module.name);
  const templatePath = join(scriptDir, 'cmake-templates/module.cmake');
  const template = readFileSync(templatePath, 'utf8');
  
  // Prepare template variables
  const sources = module.sources.map(s => `"${s}"`).join(' ');
  const linkFlags = module.linkFlags ? module.linkFlags.join(' ') : '';
  const outputFilename = module.outputFiles[0].replace('.js', ''); // Remove .js extension
  
  let content = template
    .replace(/{{MODULE_NAME}}/g, module.name)
    .replace(/{{SOURCES}}/g, sources)
    .replace(/{{LINK_FLAGS}}/g, linkFlags)
    .replace(/{{OUTPUT_FILENAME}}/g, outputFilename);
  
  // Handle conditional template sections
  if (module.useLogger) {
    content = content.replace(/{{#USE_LOGGER}}(.*?){{\/USE_LOGGER}}/gs, '$1');
  } else {
    content = content.replace(/{{#USE_LOGGER}}.*?{{\/USE_LOGGER}}/gs, '');
  }
  
  if (module.exportName) {
    // Replace the entire conditional block with the export name
    content = content.replace(/{{#EXPORT_NAME}}EXPORT_NAME "{{EXPORT_NAME}}"{{\/EXPORT_NAME}}/gs, `EXPORT_NAME "${module.exportName}"`);
  } else {
    content = content.replace(/{{#EXPORT_NAME}}.*?{{\/EXPORT_NAME}}/gs, '');
  }
  
  if (module.linkFlags && module.linkFlags.length > 0) {
    content = content.replace(/{{#LINK_FLAGS}}(.*?){{\/LINK_FLAGS}}/gs, `$1 ${linkFlags}`);
  } else {
    content = content.replace(/{{#LINK_FLAGS}}.*?{{\/LINK_FLAGS}}/gs, '');
  }
  
  const outputPath = join(moduleDir, 'CMakeLists.txt');
  writeFileSync(outputPath, content);
  console.log(`📄 Generated CMakeLists.txt for ${module.name}`);
}

function generateMainCMakeLists(scriptDir: string): void {
  const wasmDir = resolve(scriptDir, '../wasm');
  const templatePath = join(scriptDir, 'cmake-templates/main.cmake');
  const template = readFileSync(templatePath, 'utf8');
  
  const modules = getAllModules();
  const moduleNames = modules.map(m => m.name);
  const targetNames = moduleNames.join(' ');
  const copyFiles = modules.map(m => `"$<TARGET_FILE:${m.name}>"`).join('\n        ');
  
  // Generate the modules section
  const modulesSection = modules.map(m => `add_subdirectory(modules/${m.name})`).join('\n');
  
  let content = template
    .replace(/{{MODULES_SECTION}}/g, modulesSection)
    .replace(/{{TARGET_NAMES}}/g, targetNames)
    .replace(/{{COPY_FILES}}/g, copyFiles);
  
  const outputPath = join(wasmDir, 'CMakeLists.txt');
  writeFileSync(outputPath, content);
  console.log('📄 Generated main CMakeLists.txt');
}

function copyFilesToPublic(scriptDir: string): void {
  const wasmDir = resolve(scriptDir, '../wasm');
  const buildDir = join(wasmDir, 'build');
  const publicDir = resolve(scriptDir, '../public/wasm');
  
  // Ensure public directory exists
  if (!existsSync(publicDir)) {
    mkdirSync(publicDir, { recursive: true });
  }
  
  console.log('📁 Copying files to public directory...');
  
  // Copy all .js files and symbol files from build to public (WASM is embedded in single file)
  const buildFiles = readdirSync(buildDir);
  const wasmFiles = buildFiles.filter(file => 
    file.endsWith('.js') || file.endsWith('.symbols')
  );
  
  let copiedCount = 0;
  wasmFiles.forEach(file => {
    const sourcePath = join(buildDir, file);
    const destPath = join(publicDir, file);
    copyFileSync(sourcePath, destPath);
    console.log(`  📄 ${file}`);
    copiedCount++;
  });
  
  console.log(`✅ Copied ${copiedCount} files to ${publicDir}`);
  
  // Validate expected files were copied
  const modules = getAllModules();
  const expectedFiles = modules.flatMap(m => m.outputFiles);
  const missingFiles = expectedFiles.filter(file => !existsSync(join(publicDir, file)));
  
  if (missingFiles.length > 0) {
    console.error('❌ Missing expected output files:');
    missingFiles.forEach(file => console.error(`   - ${file}`));
    process.exit(1);
  }
}

function copyTemplateFiles(scriptDir: string): void {
  const wasmDir = resolve(scriptDir, '../wasm');
  const templatesDir = join(wasmDir, 'cmake-templates');
  
  // Create templates directory
  if (!existsSync(templatesDir)) {
    mkdirSync(templatesDir, { recursive: true });
  }
  
  // Copy template files
  const templateFiles = ['common.cmake'];
  for (const file of templateFiles) {
    const sourcePath = join(scriptDir, 'cmake-templates', file);
    const destPath = join(templatesDir, file);
    copyFileSync(sourcePath, destPath);
    console.log(`📄 Copied template: ${file}`);
  }
}

function main(): void {
  // Check for build flags
  const isRelease = process.argv.includes('--release');
  const isClean = process.argv.includes('--clean');
  const buildType = isRelease ? 'Release' : 'Debug';
  
  // Check for parallel jobs override
  const parallelArg = process.argv.find(arg => arg.startsWith('--parallel='));
  const customParallelJobs = parallelArg ? parseInt(parallelArg.split('=')[1]) : null;
  
  const startTime = Date.now();
  console.log(`🔨 Building All WebAssembly Modules (${buildType})...\n`);
  
  // Check Emscripten
  checkEmscripten();
  
  const scriptDir = fileURLToPath(new URL('.', import.meta.url));
  const wasmDir = resolve(scriptDir, '../wasm');
  const buildDir = join(wasmDir, 'build');
  
  // Validate module structure
  validateModuleStructure(scriptDir);
  
  // Copy template files
  copyTemplateFiles(scriptDir);
  
  // Generate CMake files
  console.log('📝 Generating CMake files...');
  const modules = getAllModules();
  
  // Generate individual module CMakeLists.txt files
  for (const module of modules) {
    generateModuleCMakeLists(scriptDir, module);
  }
  
  // Generate main CMakeLists.txt
  generateMainCMakeLists(scriptDir);
  
  // Handle clean build
  if (isClean && existsSync(buildDir)) {
    console.log('🧹 Cleaning build directory...');
    runCommand('rm -rf build', wasmDir);
  }
  
  // Create build directory
  if (!existsSync(buildDir)) {
    mkdirSync(buildDir, { recursive: true });
  }
  
  // Configure with Emscripten and build cache
  console.log('⚙️  Configuring CMake with Emscripten...');
  const cmakeCache = join(buildDir, 'CMakeCache.txt');
  const needsConfigure = !existsSync(cmakeCache);
  
  if (needsConfigure) {
    runCommand(`emcmake cmake .. -DCMAKE_BUILD_TYPE=${buildType}`, buildDir);
  } else {
    console.log('📋 Using existing CMake cache (incremental build)');
  }
  
  // Build all WASM targets with parallel compilation
  console.log('🔨 Building all WebAssembly targets with parallel compilation...');
  const numCores = cpus().length;
  const parallelJobs = customParallelJobs || Math.max(1, numCores - 1); // Leave one core free
  console.log(`🚀 Using ${parallelJobs} parallel jobs`);
  runCommand(`cmake --build . --target all_wasm_targets --parallel ${parallelJobs}`, buildDir);
  
  // Copy files to public
  copyFilesToPublic(scriptDir);
  
  // Summary
  const endTime = Date.now();
  const buildTime = (endTime - startTime) / 1000;
  
  console.log('\n✅ All builds successful!');
  console.log(`⏱️  Total build time: ${buildTime.toFixed(2)}s`);
  console.log('\nGenerated WASM modules:');
  modules.forEach(module => {
    console.log(`  - ${module.description} (${module.outputFiles.join(', ')})`);
  });
  console.log('\nFiles copied to: public/wasm/');
  console.log('\n💡 Build options:');
  console.log('  --release    : Build in Release mode (faster, smaller)');
  console.log('  --clean      : Clean build directory before building');
  console.log('  --parallel N : Use N parallel jobs (default: auto-detect)');
}

// ES module equivalent of require.main === module
console.log('🔍 Debug: Checking if script is running as main module...');
console.log('import.meta.url:', import.meta.url);
console.log('process.argv[1]:', process.argv[1]);

// Normalize paths for comparison
const scriptUrl = fileURLToPath(import.meta.url);
const scriptPath = normalize(process.argv[1]);
const normalizedScriptUrl = normalize(scriptUrl);

console.log('scriptUrl:', scriptUrl);
console.log('scriptPath:', scriptPath);
console.log('normalizedScriptUrl:', normalizedScriptUrl);
console.log('Are they equal?', normalizedScriptUrl === scriptPath);

if (normalizedScriptUrl === scriptPath) {
  console.log('✅ Running as main module - starting build...');
  main();
} else {
  console.log('❌ Not running as main module');
  console.log('This script should be run directly, not imported');
} 