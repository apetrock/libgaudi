#!/usr/bin/env ts-node

const fs = require('fs');
const path = require('path');

interface DemoConfig {
  name: string;
  description: string;
}

function pascalCase(str: string) {
  return str.charAt(0).toUpperCase() + str.slice(1);
}

function ensureDir(dir: string) {
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir, { recursive: true });
  }
}

function main() {
  const args = process.argv.slice(2);
  if (args.length === 0) {
    console.log('Usage: npx ts-node scripts/gen-demo.ts <demo_name>');
    process.exit(1);
  }
  const demoName = args[0];
  const config: DemoConfig = {
    name: demoName,
    description: `A simple ${demoName} WASM demo`,
  };

  console.log(`🎯 Generating WASM demo: ${config.name}`);

  // 1. Create WASM project directory and files
  const wasmDir = path.join(__dirname, '../wasm', config.name);
  const wasmSrcDir = path.join(wasmDir, 'src');
  ensureDir(wasmSrcDir);

  // CMakeLists.txt
  const cmakeContent = `# CMakeLists.txt for ${config.name} WASM project
add_wasm_target(${config.name}
    SOURCES src/${config.name}.cpp
    LINK_FLAGS
        "-s EXPORT_NAME=\\"${pascalCase(config.name)}Module\\""
        "-s USE_ES6_IMPORT_META=0"
)

# Copy files to demo directory after build
add_custom_command(TARGET ${config.name} POST_BUILD
    COMMAND \${CMAKE_COMMAND} -E copy_if_different
        "\$<TARGET_FILE:${config.name}>"
        "\${DEMO_WASM_DIR}/"
    COMMENT "Copying ${config.name} files to demo directory"
)
`;
  fs.writeFileSync(path.join(wasmDir, 'CMakeLists.txt'), cmakeContent);
  console.log(`📁 Created: ${path.join(wasmDir, 'CMakeLists.txt')}`);

  // C++ source
  const cppContent = `#include <emscripten/bind.h>
#include <string>

using namespace emscripten;

int add(int a, int b) {
    return a + b;
}

std::string greet(std::string name) {
    return "Hello, " + name + "!";
}

EMSCRIPTEN_BINDINGS(${config.name}) {
    function("add", &add);
    function("greet", &greet);
}
`;
  fs.writeFileSync(path.join(wasmSrcDir, `${config.name}.cpp`), cppContent);
  console.log(`📁 Created: ${path.join(wasmSrcDir, `${config.name}.cpp`)}`);

  // TypeScript endpoints
  const tsContent = `/**
 * TypeScript interface for ${config.name} WASM module
 * Generated automatically - do not edit manually
 */

export interface ${pascalCase(config.name)}Module {
  add(a: number, b: number): number;
  greet(name: string): string;
}

/**
 * ${config.description}
 */
export const ${config.name}Endpoints = {
  moduleName: '${config.name}',
  functions: [
    { name: 'add', returnType: 'number', parameters: [{ name: 'a', type: 'number' }, { name: 'b', type: 'number' }] },
    { name: 'greet', returnType: 'string', parameters: [{ name: 'name', type: 'string' }] }
  ]
} as const;
`;
  fs.writeFileSync(path.join(wasmDir, `${config.name}_endpoints.ts`), tsContent);
  console.log(`📁 Created: ${path.join(wasmDir, `${config.name}_endpoints.ts`)}`);

  // 2. Create React component
  const componentsDir = path.join(__dirname, '../src/components');
  ensureDir(componentsDir);
  const componentName = pascalCase(config.name);
  const componentPath = path.join(componentsDir, `${componentName}.tsx`);
  const reactContent = `import { useState, useEffect } from 'react';
import { ${componentName}Module } from '../../wasm/${config.name}/${config.name}_endpoints';
import { loadWasmModule } from '../utils/wasmLoader';

/**
 * ${config.description}
 * Generated automatically - do not edit manually
 */
export function ${componentName}() {
  const [module, setModule] = useState<${componentName}Module | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);

  useEffect(() => {
    const loadWasm = async () => {
      try {
        setLoading(true);
        const instance = await loadWasmModule('${config.name}') as ${componentName}Module;
        setModule(instance);
        setLoading(false);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load WASM');
        setLoading(false);
      }
    };
    loadWasm();
  }, []);

  const testAdd = () => {
    if (module) {
      const result = module.add(5, 3);
      setResult(\`add(5, 3) = \${result}\`);
    }
  };

  const testGreet = () => {
    if (module) {
      const result = module.greet('World');
      setResult(\`greet('World') = \${result}\`);
    }
  };

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
        <h1 className="text-3xl font-bold mb-2">${componentName}</h1>
        <p className="text-gray-600">${config.description}</p>
      </div>

      <div className="bg-gray-50 p-6 rounded-lg mb-6">
        <h2 className="text-xl font-semibold mb-4">Test Functions</h2>
        <div className="space-y-2">
          <button onClick={testAdd} className="px-4 py-2 bg-blue-500 text-white rounded hover:bg-blue-600">
            Test add(5, 3)
          </button>
          <button onClick={testGreet} className="px-4 py-2 bg-green-500 text-white rounded hover:bg-green-600 ml-2">
            Test greet('World')
          </button>
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
          <li>• <code>src/components/${componentName}.tsx</code> - React component</li>
        </ul>
      </div>
    </div>
  );
}
`;
  fs.writeFileSync(componentPath, reactContent);
  console.log(`📁 Created: ${componentPath}`);

  // 3. Update main CMakeLists.txt
  const mainCmakePath = path.join(__dirname, '../wasm/CMakeLists.txt');
  let cmakeMain = fs.readFileSync(mainCmakePath, 'utf-8');
  const subdirLine = `add_subdirectory(${config.name})`;
  if (!cmakeMain.includes(subdirLine)) {
    // Find the add_subdirectory section and add our line
    const lines = cmakeMain.split('\n');
    const insertIndex = lines.findIndex((line: string) => line.includes('add_subdirectory('));
    if (insertIndex !== -1) {
      lines.splice(insertIndex + 1, 0, subdirLine);
      cmakeMain = lines.join('\n');
      fs.writeFileSync(mainCmakePath, cmakeMain);
      console.log(`🔧 Updated: ${mainCmakePath}`);
    }
  }

  // 4. Print routing instructions
  console.log(`\n✅ Demo '${config.name}' generated successfully!`);
  console.log(`\n📋 Next steps:`);
  console.log(`1. Add this route to your App.tsx:`);
  console.log(`   <Route path="/${config.name}" element={<${componentName} />} />`);
  console.log(`2. Build the WASM: cd wasm && ./build.bat`);
  console.log(`3. Start dev server: cd .. && npm run dev`);
  console.log(`4. Navigate to: http://localhost:5173/${config.name}`);
  console.log(`\n🎉 Your new demo is ready!`);
}

main(); 