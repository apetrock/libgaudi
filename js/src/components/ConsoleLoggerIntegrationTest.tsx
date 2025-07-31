import React, { useEffect, useState } from 'react';
import { ConsoleLoggerRenderer, ConsoleLoggerStats } from './ConsoleLoggerRenderer';
import { consoleLogger, useConsoleLogger, LogLevel } from '../stores/consoleLoggerStore';

/**
 * Test component demonstrating console logger integration with WASM
 */
export function ConsoleLoggerIntegrationTest() {
  const [wasmModule, setWasmModule] = useState<any>(null);
  const [isSetup, setIsSetup] = useState(false);
  const logger = useConsoleLogger();

  // Load WASM module (placeholder for actual module loading)
  useEffect(() => {
    // This would typically be where you load your WASM module
    // For testing, we'll simulate having a module
    const mockWasmModule = {
      // Mock the WASM terminal logger functions
      log_info: (message: string) => console.log('[WASM INFO]', message),
      log_warning: (message: string) => console.warn('[WASM WARNING]', message),
      log_error: (message: string) => console.error('[WASM ERROR]', message),
      log_debug: (message: string) => console.log('[WASM DEBUG]', message),
      
      // Mock callback setters
      set_info_callback: (callback: (msg: string) => void) => {
        console.log('Info callback set');
        // In real implementation, this would register the callback in WASM
      },
      set_warning_callback: (callback: (msg: string) => void) => {
        console.log('Warning callback set');
      },
      set_error_callback: (callback: (msg: string) => void) => {
        console.log('Error callback set');
      },
      set_debug_callback: (callback: (msg: string) => void) => {
        console.log('Debug callback set');
      },
      
      // Mock message getters (for batch sync)
      get_info_messages: () => ({ size: () => 0, get: () => '' }),
      get_warning_messages: () => ({ size: () => 0, get: () => '' }),
      get_error_messages: () => ({ size: () => 0, get: () => '' }),
      get_debug_messages: () => ({ size: () => 0, get: () => '' }),
    };
    
    setWasmModule(mockWasmModule);
  }, []);

  // Setup WASM callbacks when module is ready
  const setupCallbacks = () => {
    if (wasmModule && !isSetup) {
      try {
        consoleLogger.setupCallbacks(wasmModule);
        setIsSetup(true);
        
        // Add a test log to confirm setup
        consoleLogger.info('Console logger callbacks setup complete');
      } catch (error) {
        consoleLogger.error(`Failed to setup callbacks: ${error}`);
      }
    }
  };

  // Test functions
  const testJSLogging = () => {
    consoleLogger.debug('This is a debug message from JavaScript');
    consoleLogger.info('This is an info message from JavaScript');
    consoleLogger.warning('This is a warning message from JavaScript');
    consoleLogger.error('This is an error message from JavaScript');
  };

  const testWasmLogging = () => {
    if (wasmModule) {
      // These would trigger the WASM callbacks in real implementation
      wasmModule.log_debug('This is a debug message from WASM');
      wasmModule.log_info('This is an info message from WASM');
      wasmModule.log_warning('This is a warning message from WASM');
      wasmModule.log_error('This is an error message from WASM');
      
      // For testing, manually add to our store since we don't have real WASM
      logger.addLog(LogLevel.DEBUG, 'This is a debug message from WASM (simulated)');
      logger.addLog(LogLevel.INFO, 'This is an info message from WASM (simulated)');
      logger.addLog(LogLevel.WARNING, 'This is a warning message from WASM (simulated)');
      logger.addLog(LogLevel.ERROR, 'This is an error message from WASM (simulated)');
    }
  };

  const syncWithWasm = () => {
    if (wasmModule) {
      consoleLogger.sync(wasmModule);
      consoleLogger.info('Synced with WASM logger');
    }
  };

  const clearLogs = () => {
    consoleLogger.clear();
  };

  return (
    <div className="console-logger-integration-test p-4 space-y-4">
      <h2 className="text-xl font-bold">Console Logger Integration Test</h2>
      
      {/* Status */}
      <div className="bg-gray-100 p-3 rounded">
        <div className="text-sm space-y-1">
          <div>WASM Module: {wasmModule ? '✅ Loaded' : '❌ Not loaded'}</div>
          <div>Callbacks Setup: {isSetup ? '✅ Ready' : '❌ Not setup'}</div>
          <div>Total Logs: {logger.getLogCount()}</div>
        </div>
      </div>

      {/* Controls */}
      <div className="flex flex-wrap gap-2">
        <button
          onClick={setupCallbacks}
          disabled={!wasmModule || isSetup}
          className="px-3 py-2 bg-blue-500 text-white rounded disabled:bg-gray-400"
        >
          Setup WASM Callbacks
        </button>
        
        <button
          onClick={testJSLogging}
          className="px-3 py-2 bg-green-500 text-white rounded"
        >
          Test JS Logging
        </button>
        
        <button
          onClick={testWasmLogging}
          disabled={!wasmModule}
          className="px-3 py-2 bg-purple-500 text-white rounded disabled:bg-gray-400"
        >
          Test WASM Logging
        </button>
        
        <button
          onClick={syncWithWasm}
          disabled={!wasmModule}
          className="px-3 py-2 bg-orange-500 text-white rounded disabled:bg-gray-400"
        >
          Sync with WASM
        </button>
        
        <button
          onClick={clearLogs}
          className="px-3 py-2 bg-red-500 text-white rounded"
        >
          Clear Logs
        </button>
      </div>

      {/* Stats */}
      <ConsoleLoggerStats />

      {/* Logger renderer */}
      <ConsoleLoggerRenderer 
        height={400}
        showTimestamps={true}
        showFrameNumbers={true}
        showLevelIcons={true}
        autoScroll={true}
        filterLevel={LogLevel.DEBUG}
      />

      {/* Usage instructions */}
      <div className="bg-blue-50 p-4 rounded">
        <h3 className="font-semibold mb-2">Usage Instructions:</h3>
        <ol className="text-sm space-y-1 list-decimal list-inside">
          <li>Click "Setup WASM Callbacks" to register callbacks with the WASM module</li>
          <li>Use "Test JS Logging" to add logs directly from JavaScript</li>
          <li>Use "Test WASM Logging" to simulate logs coming from WASM</li>
          <li>Use "Sync with WASM" to pull any accumulated logs from WASM</li>
          <li>Logs will appear in real-time as they're generated</li>
        </ol>
      </div>

      {/* Integration code example */}
      <div className="bg-gray-100 p-4 rounded">
        <h3 className="font-semibold mb-2">Integration Code Example:</h3>
        <pre className="text-xs bg-black text-green-400 p-2 rounded overflow-x-auto">
{`// In your component:
import { consoleLogger, useConsoleLogger } from '../stores/consoleLoggerStore';

// Setup callbacks when WASM module loads
useEffect(() => {
  if (wasmModule) {
    consoleLogger.setupCallbacks(wasmModule);
  }
}, [wasmModule]);

// Log from JavaScript
consoleLogger.info('Hello from JS');

// WASM will automatically call callbacks when it logs
// wasmModule.log_info('Hello from WASM'); // -> triggers callback

// Render logs
<ConsoleLoggerRenderer />
`}
        </pre>
      </div>
    </div>
  );
}
