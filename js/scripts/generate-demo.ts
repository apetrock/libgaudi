#!/usr/bin/env node

import * as fs from 'fs';
import * as path from 'path';

interface DemoConfig {
  name: string;
  description: string;
  functions: Array<{
    name: string;
    returnType: string;
    parameters: Array<{ name: string; type: string }>;
    description: string;
  }>;
}

class DemoGenerator {
  private projectRoot: string;
  private wasmDir: string;
  private srcDir: string;

  constructor() {
    this.projectRoot = path.resolve(__dirname, '..');
    this.wasmDir = path.join(this.projectRoot, 'wasm');
    this.srcDir = path.join(this.projectRoot, 'src');
  }

  async generateDemo(config: DemoConfig): Promise<void> {
    console.log(`🎯 Generating WASM demo: ${config.name}`);
    
    try {
      // Step 1: Create WASM project directory
      await this.createWasmProject(config);
      
      // Step 2: Create React component
      await this.createReactComponent(config);
      
      // Step 3: Update build system
      await this.updateBuildSystem(config);
      
      // Step 4: Update routing
      await this.updateRouting(config);
      
      console.log(`✅ Demo '${config.name}' generated successfully!`);
      console.log('\n📋 Next steps:');
      console.log(`1. cd wasm && ./build.bat`);
      console.log(`2. cd .. && npm run dev`);
      console.log(`3. Navigate to http://localhost:5173/${config.name.toLowerCase()}`);
      
    } catch (error) {
      console.error('❌ Failed to generate demo:', error);
      process.exit(1);
    }
  }

  private async createWasmProject(config: DemoConfig): Promise<void> {
    const projectDir = path.join(this.wasmDir, config.name);
    const srcDir = path.join(projectDir, 'src');
    
    // Create directories
    await this.ensureDir(projectDir);
    await this.ensureDir(srcDir);
    
    // Generate CMakeLists.txt
    const cmakeContent = this.generateCMakeLists(config);
    await fs.promises.writeFile(path.join(projectDir, 'CMakeLists.txt'), cmakeContent);
    
    // Generate C++ source
    const cppContent = this.generateCppSource(config);
    await fs.promises.writeFile(path.join(srcDir, `${config.name}.cpp`), cppContent);
    
    // Generate TypeScript endpoints
    const tsContent = this.generateTypeScriptEndpoints(config);
    await fs.promises.writeFile(path.join(projectDir, `${config.name}_endpoints.ts`), tsContent);
    
    console.log(`📁 Created WASM project: ${projectDir}`);
  }

  private async createReactComponent(config: DemoConfig): Promise<void> {
    const componentPath = path.join(this.srcDir, 'components', `${this.pascalCase(config.name)}.tsx`);
    
    const componentContent = this.generateReactComponent(config);
    await fs.promises.writeFile(componentPath, componentContent);
    
    console.log(`⚛️  Created React component: ${componentPath}`);
  }

  private async updateBuildSystem(config: DemoConfig): Promise<void> {
    // Update main CMakeLists.txt
    const mainCmakePath = path.join(this.wasmDir, 'CMakeLists.txt');
    let cmakeContent = await fs.promises.readFile(mainCmakePath, 'utf-8');
    
    // Add subdirectory if not already present
    const subdirLine = `add_subdirectory(${config.name})`;
    if (!cmakeContent.includes(subdirLine)) {
      // Find the add_subdirectory section and add our line
      const lines = cmakeContent.split('\n');
      const insertIndex = lines.findIndex(line => line.includes('add_subdirectory('));
      if (insertIndex !== -1) {
        lines.splice(insertIndex + 1, 0, subdirLine);
        cmakeContent = lines.join('\n');
        await fs.promises.writeFile(mainCmakePath, cmakeContent);
        console.log(`🔧 Updated main CMakeLists.txt`);
      }
    }
  }

  private async updateRouting(config: DemoConfig): Promise<void> {
    const appPath = path.join(this.srcDir, 'App.tsx');
    let appContent = await fs.promises.readFile(appPath, 'utf-8');
    
    // Add import if not present
    const importLine = `import { ${this.pascalCase(config.name)} } from './components/${this.pascalCase(config.name)}';`;
    if (!appContent.includes(importLine)) {
      // Find the last import statement
      const lines = appContent.split('\n');
      const lastImportIndex = lines.findLastIndex(line => line.startsWith('import '));
      if (lastImportIndex !== -1) {
        lines.splice(lastImportIndex + 1, 0, importLine);
        appContent = lines.join('\n');
      }
    }
    
    // Add route if not present
    const routeLine = `        <Route path="/${config.name.toLowerCase()}" element={<${this.pascalCase(config.name)} />} />`;
    if (!appContent.includes(routeLine)) {
      // Find the Routes section and add our route
      const lines = appContent.split('\n');
      const routesIndex = lines.findIndex(line => line.includes('<Routes>'));
      if (routesIndex !== -1) {
        // Find the closing Routes tag
        const closingRoutesIndex = lines.findIndex((line, index) => 
          index > routesIndex && line.includes('</Routes>')
        );
        if (closingRoutesIndex !== -1) {
          lines.splice(closingRoutesIndex, 0, routeLine);
          appContent = lines.join('\n');
        }
      }
    }
    
    await fs.promises.writeFile(appPath, appContent);
    console.log(`🛣️  Updated App.tsx routing`);
  }

  private generateCMakeLists(config: DemoConfig): string {
    return `# CMakeLists.txt for ${config.name} WASM project

# Use the common WASM target function
add_wasm_target(${config.name}
    SOURCES src/${config.name}.cpp
    OUTPUT_DIR "\${DEMO_WASM_DIR}"
    LINK_FLAGS
        "-s EXPORT_NAME=\\"${this.pascalCase(config.name)}Module\\""
        "-s USE_ES6_IMPORT_META=0"
)

# Copy files to demo directory after build
add_custom_command(TARGET ${config.name} POST_BUILD
    COMMAND \${CMAKE_COMMAND} -E copy_if_different
        "\${CMAKE_CURRENT_BINARY_DIR}/${config.name}.js"
        "\${CMAKE_CURRENT_BINARY_DIR}/${config.name}.wasm"
        "\${DEMO_WASM_DIR}/"
    COMMENT "Copying ${config.name} files to demo directory"
)
`;
  }

  private generateCppSource(config: DemoConfig): string {
    const functionBindings = config.functions.map(func => {
      const params = func.parameters.map(p => `${p.type} ${p.name}`).join(', ');
      const paramNames = func.parameters.map(p => p.name).join(', ');
      return `    function("${func.name}", &${func.name});`;
    }).join('\n');

    const functionDeclarations = config.functions.map(func => {
      const params = func.parameters.map(p => `${p.type} ${p.name}`).join(', ');
      return `${func.returnType} ${func.name}(${params}) {
    // TODO: Implement ${func.description}
    return ${func.returnType === 'int' ? '0' : func.returnType === 'double' ? '0.0' : '""'};
}`;
    }).join('\n\n');

    return `#include <emscripten/bind.h>
#include <string>

using namespace emscripten;

${functionDeclarations}

EMSCRIPTEN_BINDINGS(${config.name}) {
${functionBindings}
}
`;
  }

  private generateTypeScriptEndpoints(config: DemoConfig): string {
    const functionSignatures = config.functions.map(func => {
      const params = func.parameters.map(p => `${p.name}: ${this.tsType(func.parameters.find(p => p.name === p.name)?.type || 'any')}`).join(', ');
      return `  ${func.name}(${params}): ${this.tsType(func.returnType)};`;
    }).join('\n');

    return `/**
 * TypeScript interface for ${config.name} WASM module
 * Generated automatically - do not edit manually
 */

export interface ${this.pascalCase(config.name)}Module {
${functionSignatures}
}

/**
 * ${config.description}
 */
export const ${config.name}Endpoints = {
  moduleName: '${config.name}',
  functions: ${JSON.stringify(config.functions, null, 2)}
} as const;
`;
  }

  private generateReactComponent(config: DemoConfig): string {
    const functionTests = config.functions.map(func => {
      const args = func.parameters.map(p => {
        switch (p.type) {
          case 'int': return '42';
          case 'double': return '3.14';
          case 'std::string': return '"test"';
          default: return 'null';
        }
      }).join(', ');
      
      return `  const test${this.pascalCase(func.name)} = () => {
    if (module) {
      const result = module.${func.name}(${args});
      setResult(\`${func.name}(${args}) = \${result}\`);
    }
  };`;
    }).join('\n\n');

    const testButtons = config.functions.map(func => 
      `        <button onClick={test${this.pascalCase(func.name)}} className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600">
          Test ${func.name}
        </button>`
    ).join('\n');

    return `import { useState, useEffect } from 'react';
import { ${this.pascalCase(config.name)}Module } from '../../wasm/${config.name}/${config.name}_endpoints';

/**
 * ${config.description}
 * Generated automatically - do not edit manually
 */
export function ${this.pascalCase(config.name)}() {
  const [module, setModule] = useState<${this.pascalCase(config.name)}Module | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);

  useEffect(() => {
    const loadWasm = async () => {
      try {
        setLoading(true);
        const wasmModule = await import('@wasmbuilds/${config.name}.js');
        const instance = await wasmModule.default();
        setModule(instance);
        setLoading(false);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load WASM');
        setLoading(false);
      }
    };
    loadWasm();
  }, []);

${functionTests}

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500 mx-auto mb-4"></div>
          <p>Loading ${config.name} WASM module...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-center text-red-600">
          <h3 className="text-xl font-bold mb-2">❌ Load Failed</h3>
          <p>{error}</p>
          <p className="text-sm mt-2">Make sure to run: cd wasm && ./build.bat</p>
        </div>
      </div>
    );
  }

  return (
    <div className="p-8 max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold mb-2">${config.name}</h1>
        <p className="text-gray-600">${config.description}</p>
      </div>

      <div className="bg-gray-50 p-6 rounded-lg mb-6">
        <h2 className="text-xl font-semibold mb-4">Test Functions</h2>
        <div className="space-y-2">
${testButtons}
        </div>
      </div>

      {result && (
        <div className="bg-green-50 border border-green-200 p-4 rounded-lg">
          <h3 className="font-semibold text-green-800 mb-2">Result:</h3>
          <p className="text-green-700 font-mono">{result}</p>
        </div>
      )}

      <div className="mt-8 p-4 bg-blue-50 rounded-lg">
        <h3 className="font-semibold mb-2">Generated Files:</h3>
        <ul className="text-sm space-y-1">
          <li>• <code>wasm/${config.name}/src/${config.name}.cpp</code> - C++ source</li>
          <li>• <code>wasm/${config.name}/${config.name}_endpoints.ts</code> - TypeScript interface</li>
          <li>• <code>src/components/${this.pascalCase(config.name)}.tsx</code> - React component</li>
        </ul>
      </div>
    </div>
  );
}
`;
  }

  private tsType(cppType: string): string {
    switch (cppType) {
      case 'int': return 'number';
      case 'double': return 'number';
      case 'std::string': return 'string';
      case 'bool': return 'boolean';
      default: return 'any';
    }
  }

  private pascalCase(str: string): string {
    return str.charAt(0).toUpperCase() + str.slice(1);
  }

  private async ensureDir(dir: string): Promise<void> {
    if (!fs.existsSync(dir)) {
      await fs.promises.mkdir(dir, { recursive: true });
    }
  }
}

// CLI interface
async function main() {
  const args = process.argv.slice(2);
  
  if (args.length === 0) {
    console.log('Usage: pnpm generate-demo <demo-name>');
    console.log('');
    console.log('Example:');
    console.log('  pnpm generate-demo calculator');
    console.log('');
    console.log('This will generate:');
    console.log('  • wasm/calculator/src/calculator.cpp');
    console.log('  • wasm/calculator/calculator_endpoints.ts');
    console.log('  • src/components/Calculator.tsx');
    console.log('  • Updated build system and routing');
    process.exit(1);
  }

  const demoName = args[0];
  
  // Default demo configuration
  const config: DemoConfig = {
    name: demoName,
    description: `A simple ${demoName} WASM demo with basic functions`,
    functions: [
      {
        name: 'add',
        returnType: 'int',
        parameters: [
          { name: 'a', type: 'int' },
          { name: 'b', type: 'int' }
        ],
        description: 'Add two integers'
      },
      {
        name: 'multiply',
        returnType: 'double',
        parameters: [
          { name: 'x', type: 'double' },
          { name: 'y', type: 'double' }
        ],
        description: 'Multiply two doubles'
      },
      {
        name: 'greet',
        returnType: 'std::string',
        parameters: [
          { name: 'name', type: 'std::string' }
        ],
        description: 'Generate a greeting message'
      }
    ]
  };

  const generator = new DemoGenerator();
  await generator.generateDemo(config);
}

if (require.main === module) {
  main().catch(console.error);
} 