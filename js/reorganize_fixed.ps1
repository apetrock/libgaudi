# WASM Project Reorganization Script
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
    if (Test-Path "wasm\wasm_logger.cpp") {
        Move-Item "wasm\wasm_logger.cpp" "wasm\logger\src\" -Force
        Write-Host "   ✓ Moved wasm_logger.cpp to logger/src/"
    } else {
        Write-Host "   - wasm_logger.cpp not found (may already be moved)"
    }
    
    # Move example files if they exist
    if (Test-Path "wasm\examples") {
        Get-ChildItem "wasm\examples\*.cpp" -ErrorAction SilentlyContinue | ForEach-Object {
            Move-Item $_.FullName "wasm\examples\src\" -Force
            Write-Host "   ✓ Moved $($_.Name) to examples/src/"
        }
    }
    
    # Move any loose cpp files in wasm root to examples
    Get-ChildItem "wasm\*.cpp" -ErrorAction SilentlyContinue | ForEach-Object {
        Move-Item $_.FullName "wasm\examples\src\" -Force
        Write-Host "   ✓ Moved $($_.Name) to examples/src/"
    }
    
    # Move rod simulation files if they exist
    if (Test-Path "wasm\rod_simulation") {
        Get-ChildItem "wasm\rod_simulation\*.cpp" -ErrorAction SilentlyContinue | ForEach-Object {
            Move-Item $_.FullName "wasm\rod_simulation\src\" -Force
            Write-Host "   ✓ Moved $($_.Name) to rod_simulation/src/"
        }
    }
    
    Write-Host "`n3. Creating TypeScript endpoint files..." -ForegroundColor Cyan
    
    # Create example_endpoints.ts with here-string that works in PowerShell
    @"
// Example WASM TypeScript endpoints
export interface ExampleWASMModule {
  hello_world(): string;
  add_numbers(a: number, b: number): number;
}

export async function loadExampleWASM(): Promise<ExampleWASMModule> {
  const module = await import('../build/hello_world.js');
  return module.default();
}
"@ | Out-File -FilePath "wasm\examples\example_endpoints.ts" -Encoding UTF8
    Write-Host "   ✓ Created examples/example_endpoints.ts"
    
    # Create logger_endpoints.ts
    @"
// Logger WASM TypeScript endpoints
export interface LoggerWASMModule {
  get_logged_lines_count(): number;
  get_logged_lines_positions(): number;
  get_logged_lines_colors(): number;
  clear_logged_lines(): void;
  animate_axes(time: number): void;
  animate_sine_wave(time: number): void;
}

export async function loadLoggerWASM(): Promise<LoggerWASMModule> {
  const module = await import('../build/gaudi_logger_test.js');
  return module.default();
}
"@ | Out-File -FilePath "wasm\logger\logger_endpoints.ts" -Encoding UTF8
    Write-Host "   ✓ Created logger/logger_endpoints.ts"
    
    # Create rod_endpoints.ts
    @"
// Rod Simulation WASM TypeScript endpoints
export interface RodSimulationWASMModule {
  step(): void;
  getVertices(): number;
  getNormals(): number;
  reset(): void;
  setParameters(params: any): void;
}

export async function loadRodSimulationWASM(): Promise<RodSimulationWASMModule> {
  const module = await import('../build/rod_simulation_main.js');
  return module.default();
}
"@ | Out-File -FilePath "wasm\rod_simulation\rod_endpoints.ts" -Encoding UTF8
    Write-Host "   ✓ Created rod_simulation/rod_endpoints.ts"
    
    # Step 4: Clean up old demo folder if it exists
    Write-Host "`n4. Cleaning up old structure..." -ForegroundColor Cyan
    
    if (Test-Path "demo") {
        Remove-Item "demo" -Recurse -Force -ErrorAction SilentlyContinue
        Write-Host "   ✓ Removed old demo folder"
    } else {
        Write-Host "   - No old demo folder to remove"
    }
    
    Write-Host "`n✅ Reorganization completed successfully!" -ForegroundColor Green
    Write-Host "`nFinal WASM structure:" -ForegroundColor Yellow
    
    if (Test-Path "wasm") {
        Get-ChildItem "wasm" -Recurse | ForEach-Object {
            $indent = "  " * (($_.FullName -split "\\").Count - (Split-Path $PWD -Leaf).Length - 2)
            Write-Host "$indent$($_.Name)" -ForegroundColor Gray
        }
    }
    
    Write-Host "`nNext steps:" -ForegroundColor Yellow
    Write-Host "1. Update CMakeLists.txt files for each WASM project"
    Write-Host "2. Update component imports to use new endpoint files"
    Write-Host "3. Test the build process"
    
} catch {
    Write-Host "`n❌ Error during reorganization: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}
