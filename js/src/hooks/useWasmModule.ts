import { useState, useEffect, useRef } from 'react';
import { loadWasmModule, clearModuleCache, type WasmModule } from '../utils/wasmLoader';
import type { WasmSimulation } from '../types/simulation';

interface UseWasmModuleReturn {
  module: WasmModule | null;
  isLoading: boolean;
  error: string | null;
  reload: () => Promise<void>;
  clearCache: () => void;
}

/**
 * Hook to load the WebAssembly module for rod simulation
 * Handles loading state and error management with caching
 */
export function useWasmModule(wasmPath: string = '/rod_simulation.js'): UseWasmModuleReturn {
  const [module, setModule] = useState<WasmModule | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const loadingRef = useRef<Promise<WasmModule> | null>(null);

  const loadModule = async (forceReload: boolean = false): Promise<WasmModule> => {
    try {
      // Check if we're already loading
      if (loadingRef.current && !forceReload) {
        return loadingRef.current;
      }

      if (!forceReload) {
        setIsLoading(true);
        setError(null);
      }

      // Extract filename from path
      const filename = wasmPath.split('/').pop() || 'rod_simulation.js';
      
      // Load the module using our improved loader
      const loadPromise = loadWasmModule(filename, forceReload);
      loadingRef.current = loadPromise;
      
      const loadedModule = await loadPromise;
      setModule(loadedModule);
      setIsLoading(false);
      
      return loadedModule;
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Unknown error loading WASM module');
      setIsLoading(false);
      throw err;
    }
  };

  const reload = async (): Promise<void> => {
    try {
      await loadModule(true);
    } catch (err) {
      // Error is already handled in loadModule
    }
  };

  const clearCache = (): void => {
    clearModuleCache();
    setModule(null);
    setError(null);
  };

  useEffect(() => {
    let cancelled = false;

    const initModule = async () => {
      try {
        await loadModule();
      } catch (err) {
        if (!cancelled) {
          // Error is already handled in loadModule
        }
      }
    };

    initModule();

    return () => {
      cancelled = true;
      loadingRef.current = null;
    };
  }, [wasmPath]);

  return { module, isLoading, error, reload, clearCache };
}

interface UseWasmSimulationReturn {
  simulation: WasmSimulation | null;
  isReady: boolean;
  error: string | null;
}

/**
 * Hook to create and manage a WASM simulation instance
 */
export function useWasmSimulation(module: WasmModule | null): UseWasmSimulationReturn {
  const [simulation, setSimulation] = useState<WasmSimulation | null>(null);
  const [error, setError] = useState<string | null>(null);
  const simulationRef = useRef<WasmSimulation | null>(null);

  useEffect(() => {
    if (!module) {
      setSimulation(null);
      return;
    }

    try {
      // Clean up previous simulation
      if (simulationRef.current) {
        simulationRef.current.delete();
      }

      // Create new simulation instance
      const newSimulation = new module.RodSimulation();
      simulationRef.current = newSimulation;
      setSimulation(newSimulation);
      setError(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to create simulation instance');
    }

    // Cleanup function
    return () => {
      if (simulationRef.current) {
        simulationRef.current.delete();
        simulationRef.current = null;
      }
    };
  }, [module]);

  return {
    simulation,
    isReady: simulation !== null,
    error
  };
}
