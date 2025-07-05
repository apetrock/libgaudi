import React, { useState } from 'react';
import { loadWasmModule, clearModuleCache, isWasmModuleAvailable } from '../utils/wasmLoader';
import { Button } from './ui/button';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from './ui/card';
import { Badge } from './ui/badge';
import { RefreshCw, Trash2, Play, AlertTriangle, CheckCircle } from 'lucide-react';

interface ModuleStatus {
  name: string;
  available: boolean;
  loading: boolean;
  error: string | null;
}

export function WasmDebugPanel() {
  const [modules, setModules] = useState<ModuleStatus[]>([
    { name: 'hello_world.js', available: false, loading: false, error: null },
    { name: 'foo_demo.js', available: false, loading: false, error: null },
    { name: 'bar_demo.js', available: false, loading: false, error: null },
    { name: 'gaudi_logger_test.js', available: false, loading: false, error: null },
  ]);

  const checkModule = async (moduleName: string) => {
    const moduleIndex = modules.findIndex(m => m.name === moduleName);
    if (moduleIndex === -1) return;

    setModules(prev => prev.map((m, i) => 
      i === moduleIndex ? { ...m, loading: true, error: null } : m
    ));

    try {
      const available = await isWasmModuleAvailable(moduleName);
      setModules(prev => prev.map((m, i) => 
        i === moduleIndex ? { ...m, available, loading: false } : m
      ));
    } catch (error) {
      setModules(prev => prev.map((m, i) => 
        i === moduleIndex ? { 
          ...m, 
          available: false, 
          loading: false, 
          error: error instanceof Error ? error.message : 'Unknown error' 
        } : m
      ));
    }
  };

  const reloadModule = async (moduleName: string) => {
    const moduleIndex = modules.findIndex(m => m.name === moduleName);
    if (moduleIndex === -1) return;

    setModules(prev => prev.map((m, i) => 
      i === moduleIndex ? { ...m, loading: true, error: null } : m
    ));

    try {
      await loadWasmModule(moduleName, true); // Force reload
      setModules(prev => prev.map((m, i) => 
        i === moduleIndex ? { ...m, available: true, loading: false } : m
      ));
    } catch (error) {
      setModules(prev => prev.map((m, i) => 
        i === moduleIndex ? { 
          ...m, 
          available: false, 
          loading: false, 
          error: error instanceof Error ? error.message : 'Unknown error' 
        } : m
      ));
    }
  };

  const clearAllCache = () => {
    clearModuleCache();
    setModules(prev => prev.map(m => ({ ...m, available: false, error: null })));
  };

  const checkAllModules = () => {
    modules.forEach(m => checkModule(m.name));
  };

  return (
    <div className="p-6 max-w-4xl mx-auto">
      <Card>
        <CardHeader>
          <CardTitle className="flex items-center gap-2">
            <AlertTriangle className="h-5 w-5 text-yellow-500" />
            WASM Module Debug Panel
          </CardTitle>
          <CardDescription>
            Debug and troubleshoot WASM module loading issues. Use this to identify why modules load once then stop loading.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Control Buttons */}
          <div className="flex gap-2 flex-wrap">
            <Button onClick={checkAllModules} variant="outline" size="sm">
              <Play className="h-4 w-4 mr-2" />
              Check All Modules
            </Button>
            <Button onClick={clearAllCache} variant="outline" size="sm">
              <Trash2 className="h-4 w-4 mr-2" />
              Clear All Cache
            </Button>
          </div>

          {/* Module Status List */}
          <div className="space-y-3">
            {modules.map((module) => (
              <div key={module.name} className="flex items-center justify-between p-3 border rounded-lg">
                <div className="flex items-center gap-3">
                  <div className="flex items-center gap-2">
                    {module.loading ? (
                      <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-500" />
                    ) : module.available ? (
                      <CheckCircle className="h-4 w-4 text-green-500" />
                    ) : (
                      <AlertTriangle className="h-4 w-4 text-red-500" />
                    )}
                    <span className="font-mono text-sm">{module.name}</span>
                  </div>
                  <Badge variant={module.available ? "default" : "destructive"}>
                    {module.available ? "Available" : "Unavailable"}
                  </Badge>
                </div>
                
                <div className="flex items-center gap-2">
                  {module.error && (
                    <span className="text-xs text-red-500 max-w-48 truncate" title={module.error}>
                      {module.error}
                    </span>
                  )}
                  <Button 
                    onClick={() => checkModule(module.name)} 
                    variant="outline" 
                    size="sm"
                    disabled={module.loading}
                  >
                    <Play className="h-3 w-3" />
                  </Button>
                  <Button 
                    onClick={() => reloadModule(module.name)} 
                    variant="outline" 
                    size="sm"
                    disabled={module.loading}
                  >
                    <RefreshCw className="h-3 w-3" />
                  </Button>
                </div>
              </div>
            ))}
          </div>

          {/* Debug Information */}
          <div className="mt-6 p-4 bg-gray-50 rounded-lg">
            <h4 className="font-semibold mb-2">Debug Information</h4>
            <ul className="text-sm space-y-1 text-gray-600">
              <li>• <strong>Check All:</strong> Tests if modules are available without loading them</li>
              <li>• <strong>Reload:</strong> Forces a fresh load of the module (bypasses cache)</li>
              <li>• <strong>Clear Cache:</strong> Removes all cached modules from memory</li>
              <li>• <strong>Common Issues:</strong> Memory leaks, module state persistence, import cache conflicts</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </div>
  );
} 