# WASM Project Reorganization Script
# This script reorganizes the WASM structure according to the proposed file structure

Write-Host "Starting WASM project reorganization..." -ForegroundColor Green

# Set error action to stop on any error
$ErrorActionPreference = "Stop"

try {
    # Navigate to the js directory
    Set-Location "c:\Users\jdelaney\source\libgaudi\js"
    
    Write-Host "Current directory: $(Get-Location)" -ForegroundColor Yellow
    
    # Step 1: Create new directory structure
    Write-Host "`n1. Creating new directory structure..." -ForegroundColor Cyan
    
    # Create WASM project source directories
    New-Item -Path "wasm\examples\src" -ItemType Directory -Force | Out-Null
    New-Item -Path "wasm\logger\src" -ItemType Directory -Force | Out-Null
    New-Item -Path "wasm\rod_simulation\src" -ItemType Directory -Force | Out-Null
    
    # Create Projects directories with assets
    New-Item -Path "Projects\phase0_box\assets" -ItemType Directory -Force | Out-Null
    New-Item -Path "Projects\phase1_logger\assets" -ItemType Directory -Force | Out-Null
    New-Item -Path "Projects\phase2_rod\assets" -ItemType Directory -Force | Out-Null
    New-Item -Path "Projects\wasm_hello\assets" -ItemType Directory -Force | Out-Null
    
    Write-Host "   ✓ Directory structure created"
    
    # Step 2: Move WASM C++ files to appropriate project src folders
    Write-Host "`n2. Moving WASM C++ files..." -ForegroundColor Cyan
    
    # Move logger files
    if (Test-Path "wasm\wasm_geometry_logger.cpp") {
        Move-Item "wasm\wasm_geometry_logger.cpp" "wasm\logger\src\" -Force
        Write-Host "   ✓ Moved wasm_geometry_logger.cpp to logger/src/"
    }
    
    # Move example files if they exist
    if (Test-Path "wasm\examples") {
        Get-ChildItem "wasm\examples\*.cpp" -ErrorAction SilentlyContinue | ForEach-Object {
            Move-Item $_.FullName "wasm\examples\src\" -Force
            Write-Host "   ✓ Moved $($_.Name) to examples/src/"
        }
    }
    
    # Move rod simulation files if they exist
    if (Test-Path "wasm\rod_simulation") {
        Get-ChildItem "wasm\rod_simulation\*.cpp" -ErrorAction SilentlyContinue | ForEach-Object {
            Move-Item $_.FullName "wasm\rod_simulation\src\" -Force
            Write-Host "   ✓ Moved $($_.Name) to rod_simulation/src/"
        }
    }
    
    # Step 3: Create TypeScript endpoint files
    Write-Host "`n3. Creating TypeScript endpoint files..." -ForegroundColor Cyan
    
    # Create example_endpoints.ts
    $exampleEndpoints = @"
// Example WASM TypeScript endpoints
export interface ExampleWASMModule {
  hello_world(): string;
  add_numbers(a: number, b: number): number;
}

export async function loadExampleWASM(): Promise<ExampleWASMModule> {
  const module = await import('../build/hello_world.js');
  return module.default();
}
"@
    Set-Content -Path "wasm\examples\example_endpoints.ts" -Value $exampleEndpoints
    Write-Host "   ✓ Created examples/example_endpoints.ts"
    
    # Create logger_endpoints.ts
    $loggerEndpoints = @"
// Logger WASM TypeScript endpoints
export interface LoggerWASMModule {
  get_logged_lines_count(): number;
  get_logged_lines_positions(): number; // Float32Array view
  get_logged_lines_colors(): number; // Float32Array view
  clear_logged_lines(): void;
  animate_axes(time: number): void;
  animate_sine_wave(time: number): void;
}

export async function loadLoggerWASM(): Promise<LoggerWASMModule> {
  const module = await import('../build/gaudi_logger_test.js');
  return module.default();
}
"@
    Set-Content -Path "wasm\logger\logger_endpoints.ts" -Value $loggerEndpoints
    Write-Host "   ✓ Created logger/logger_endpoints.ts"
    
    # Create rod_endpoints.ts
    $rodEndpoints = @"
// Rod Simulation WASM TypeScript endpoints
export interface RodSimulationWASMModule {
  step(): void;
  getVertices(): number; // Float32Array view
  getNormals(): number; // Float32Array view
  reset(): void;
  setParameters(params: any): void;
}

export async function loadRodSimulationWASM(): Promise<RodSimulationWASMModule> {
  const module = await import('../build/rod_simulation_main.js');
  return module.default();
}
"@
    Set-Content -Path "wasm\rod_simulation\rod_endpoints.ts" -Value $rodEndpoints
    Write-Host "   ✓ Created rod_simulation/rod_endpoints.ts"
    
    # Step 4: Update import paths in existing files
    Write-Host "`n4. Updating import paths..." -ForegroundColor Cyan
    
    # Find and update WASM import paths in TypeScript/JavaScript files
    $tsFiles = Get-ChildItem -Path "src", "Projects" -Recurse -Include "*.ts", "*.tsx", "*.js", "*.jsx" -ErrorAction SilentlyContinue
    
    foreach ($file in $tsFiles) {
        $content = Get-Content $file.FullName -Raw -ErrorAction SilentlyContinue
        if ($content) {
            $originalContent = $content
            
            # Update old demo/wasm imports to new wasm project structure
            $content = $content -replace "from\s+['\`]\.\.\/\.\.\/demo\/wasm\/([^'\`]+)['\`]", "from '../wasm/logger/logger_endpoints'"
            $content = $content -replace "import\s+['\`]\.\.\/\.\.\/demo\/wasm\/([^'\`]+)['\`]", "import '../wasm/logger/logger_endpoints'"
            
            if ($content -ne $originalContent) {
                Set-Content -Path $file.FullName -Value $content
                Write-Host "   ✓ Updated imports in $($file.Name)"
            }
        }
    }
    
    # Step 5: Clean up old demo folder if it exists
    Write-Host "`n5. Cleaning up old structure..." -ForegroundColor Cyan
    
    if (Test-Path "demo") {
        Remove-Item "demo" -Recurse -Force -ErrorAction SilentlyContinue
        Write-Host "   ✓ Removed old demo folder"
    }
    
    # Step 6: Display final structure
    Write-Host "`n6. Final structure verification..." -ForegroundColor Cyan
    Write-Host "WASM project structure:" -ForegroundColor Yellow
    
    if (Test-Path "wasm") {
        tree wasm /F | Write-Host
    }
    
    Write-Host "`n✅ Reorganization completed successfully!" -ForegroundColor Green
    Write-Host "`nNext steps:" -ForegroundColor Yellow
    Write-Host "1. Update CMakeLists.txt files for each WASM project"
    Write-Host "2. Test the build process: wasm\build.bat"
    Write-Host "3. Update component imports to use new endpoint files"
    Write-Host "4. Test the React app: npm run dev"
    
} catch {
    Write-Host "`n❌ Error during reorganization: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "Stack trace: $($_.ScriptStackTrace)" -ForegroundColor Red
    exit 1
}
