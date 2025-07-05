/**
 * WASM Module Loader Utility
 * 
 * Provides a unified interface for loading WebAssembly modules
 * with proper error handling and type safety.
 */

import { useState, useEffect } from 'react';

export interface WasmModule {
  [key: string]: any;
}

// Cache for loaded modules to prevent reloading
const moduleCache = new Map<string, Promise<WasmModule>>();

/**
 * React hook for loading WASM modules with state management
 * @param moduleFilename - The filename of the WASM module to load
 * @param instanceFactory - Optional function to create an instance from the module
 * @returns Object containing loading state, error, module, and instance
 */
export function useWasmModule<T extends WasmModule, I>(
  moduleFilename: string,
  instanceFactory?: (module: T) => I
) {
  const [state, setState] = useState<{
    loading: boolean;
    error: string | null;
    module: T | null;
    instance: I | null;
  }>({
    loading: true,
    error: null,
    module: null,
    instance: null
  });

  useEffect(() => {
    let isMounted = true;

    async function loadWasm() {
      try {
        setState(prev => ({ ...prev, loading: true, error: null }));
        
        const wasmModule = await loadWasmModule(moduleFilename) as T;
        
        if (!isMounted) return;

        let instance: I | null = null;
        if (instanceFactory) {
          instance = instanceFactory(wasmModule);
        }

        setState({
          loading: false,
          error: null,
          module: wasmModule,
          instance
        });
      } catch (error) {
        if (!isMounted) return;
        
        console.error(`Failed to load WASM module '${moduleFilename}':`, error);
        setState({
          loading: false,
          error: error instanceof Error ? error.message : 'Unknown error',
          module: null,
          instance: null
        });
      }
    }
    
    loadWasm();

    return () => {
      isMounted = false;
    };
  }, [moduleFilename, instanceFactory]);

  return state;
}

/**
 * Loads a WASM module by filename with caching and error recovery
 * @param moduleFilename - The filename of the WASM module to load (e.g., "gaudi_logger_test.js")
 * @param forceReload - Force reload even if cached (default: false)
 * @returns Promise that resolves to the loaded WASM module
 */
export async function loadWasmModule(moduleFilename: string, forceReload: boolean = false): Promise<WasmModule> {
  try {
    // Ensure the filename has .js extension
    const normalizedFilename = moduleFilename.endsWith('.js') ? moduleFilename : `${moduleFilename}.js`;
    
    // Check cache first (unless force reload is requested)
    if (!forceReload && moduleCache.has(normalizedFilename)) {
      return await moduleCache.get(normalizedFilename)!;
    }

    // Dynamic import using the provided filename
    // Files in public directory are served at root path, so use /wasm/ instead of /public/wasm/
    const module = await import(/* @vite-ignore */ `/wasm/${normalizedFilename}`);

    // Initialize the module if it has an initialize function
    let wasmModule: WasmModule;
    if (module.default && typeof module.default === 'function') {
      wasmModule = await module.default();
    } else if (module.default) {
      wasmModule = module.default;
    } else {
      wasmModule = module;
    }

    // Cache the successful load
    const modulePromise = Promise.resolve(wasmModule);
    moduleCache.set(normalizedFilename, modulePromise);
    
    return wasmModule;
    
  } catch (error) {
    // Remove from cache if it failed
    const normalizedFilename = moduleFilename.endsWith('.js') ? moduleFilename : `${moduleFilename}.js`;
    moduleCache.delete(normalizedFilename);
    
    // If it's a module loading error, try to clear cache and retry once
    if (!forceReload && error instanceof Error && 
        (error.message.includes('fetch') || error.message.includes('import'))) {
      return loadWasmModule(moduleFilename, true);
    }
    
    throw new Error(`Failed to load WASM module '${moduleFilename}': ${error instanceof Error ? error.message : 'Unknown error'}`);
  }
}

/**
 * Clears the module cache for a specific module or all modules
 * @param moduleFilename - Optional specific module to clear, or undefined to clear all
 */
export function clearModuleCache(moduleFilename?: string): void {
  if (moduleFilename) {
    const normalizedFilename = moduleFilename.endsWith('.js') ? moduleFilename : `${moduleFilename}.js`;
    moduleCache.delete(normalizedFilename);
  } else {
    moduleCache.clear();
  }
}

/**
 * Checks if a WASM module is available
 * @param moduleFilename - The filename of the WASM module to check
 * @returns Promise that resolves to true if the module is available
 */
export async function isWasmModuleAvailable(moduleFilename: string): Promise<boolean> {
  try {
    await loadWasmModule(moduleFilename);
    return true;
  } catch {
    return false;
  }
}

/**
 * Gets a list of available WASM modules
 * @returns Array of available module names
 */
export function getAvailableWasmModules(): string[] {
  return ['hello_world', 'foo_demo', 'bar_demo', 'gaudi_logger_test', 'rod_constraints_test', 'test_module'];
} 