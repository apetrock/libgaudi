import { create } from 'zustand';
import { useEffect } from 'react';

/**
 * Console log levels matching C++ enum
 */
export enum LogLevel {
  DEBUG = 0,
  INFO = 1,
  WARNING = 2,
  ERROR = 3
}

/**
 * Console log entry structure matching C++ LogEntry
 */
export interface ConsoleLogEntry {
  level: LogLevel;
  message: string;
  timestamp: number;
  frame: number;
}

/**
 * Console logger configuration
 */
export interface ConsoleLoggerConfig {
  enabled: boolean;
  maxLogs: number;
  autoScroll: boolean;
  showTimestamps: boolean;
  showFrameNumbers: boolean;
  filterLevel: LogLevel;
}

/**
 * Console logger callback types
 */
export type ConsoleLogCallback = (entry: ConsoleLogEntry) => void;
export type ConsoleLogBatchCallback = (entries: ConsoleLogEntry[]) => void;

/**
 * Console logger store state
 */
interface ConsoleLoggerState {
  logs: ConsoleLogEntry[];
  config: ConsoleLoggerConfig;
  callbacks: ConsoleLogCallback[];
  batchCallbacks: ConsoleLogBatchCallback[];
  
  // Actions
  addLog: (level: LogLevel, message: string, timestamp?: number, frame?: number) => void;
  clearLogs: () => void;
  updateConfig: (config: Partial<ConsoleLoggerConfig>) => void;
  
  // Callback management
  addCallback: (callback: ConsoleLogCallback) => () => void;
  addBatchCallback: (callback: ConsoleLogBatchCallback) => () => void;
  
  // WASM integration
  syncWithWasm: (wasmModule: any) => void;
  setupWasmCallbacks: (wasmModule: any) => void;
  fetchWasmLogs: (wasmModule: any) => void;
  
  // Getters
  getLogs: () => ConsoleLogEntry[];
  getFilteredLogs: () => ConsoleLogEntry[];
  getLogCount: () => number;
}

/**
 * Default configuration
 */
const defaultConfig: ConsoleLoggerConfig = {
  enabled: true,
  maxLogs: 1000,
  autoScroll: true,
  showTimestamps: true,
  showFrameNumbers: true,
  filterLevel: LogLevel.DEBUG
};

/**
 * Zustand store for console logging
 */
export const useConsoleLoggerStore = create<ConsoleLoggerState>((set, get) => ({
  // State
  logs: [],
  config: defaultConfig,
  callbacks: [],
  batchCallbacks: [],
  
  // Actions
  addLog: (level, message, timestamp, frame) => {
    const state = get();
    if (!state.config.enabled) return;
    
    const entry: ConsoleLogEntry = {
      level,
      message,
      timestamp: timestamp ?? Date.now() / 1000,
      frame: frame ?? 0
    };
    
    set(state => {
      const newLogs = [...state.logs, entry];
      
      // Enforce max logs limit
      if (newLogs.length > state.config.maxLogs) {
        newLogs.splice(0, newLogs.length - state.config.maxLogs);
      }
      
      // Call callbacks
      state.callbacks.forEach(callback => callback(entry));
      
      return { logs: newLogs };
    });
  },
  
  clearLogs: () => {
    set({ logs: [] });
  },
  
  updateConfig: (newConfig) => {
    set(state => ({
      config: { ...state.config, ...newConfig }
    }));
  },
  
  // Callback management
  addCallback: (callback) => {
    set(state => ({
      callbacks: [...state.callbacks, callback]
    }));
    
    // Return unsubscribe function
    return () => {
      set(state => ({
        callbacks: state.callbacks.filter(cb => cb !== callback)
      }));
    };
  },
  
  addBatchCallback: (callback) => {
    set(state => ({
      batchCallbacks: [...state.batchCallbacks, callback]
    }));
    
    // Return unsubscribe function
    return () => {
      set(state => ({
        batchCallbacks: state.batchCallbacks.filter(cb => cb !== callback)
      }));
    };
  },
  
  // WASM integration
  syncWithWasm: (wasmModule) => {
    if (!wasmModule) return;
    
    try {
      const state = get();
      
      // Get logs from WASM using our terminal logger API
      const infoMessages = wasmModule.get_info_messages();
      const warningMessages = wasmModule.get_warning_messages();
      const errorMessages = wasmModule.get_error_messages();
      const debugMessages = wasmModule.get_debug_messages();
      
      // Convert WASM logs to our format
      const newEntries: ConsoleLogEntry[] = [];
      
      // Process each log level
      for (let i = 0; i < infoMessages.size(); i++) {
        newEntries.push({
          level: LogLevel.INFO,
          message: infoMessages.get(i),
          timestamp: Date.now() / 1000,
          frame: 0 // TODO: Get frame from WASM if available
        });
      }
      
      for (let i = 0; i < warningMessages.size(); i++) {
        newEntries.push({
          level: LogLevel.WARNING,
          message: warningMessages.get(i),
          timestamp: Date.now() / 1000,
          frame: 0
        });
      }
      
      for (let i = 0; i < errorMessages.size(); i++) {
        newEntries.push({
          level: LogLevel.ERROR,
          message: errorMessages.get(i),
          timestamp: Date.now() / 1000,
          frame: 0
        });
      }
      
      for (let i = 0; i < debugMessages.size(); i++) {
        newEntries.push({
          level: LogLevel.DEBUG,
          message: debugMessages.get(i),
          timestamp: Date.now() / 1000,
          frame: 0
        });
      }
      
      // Only add new logs (avoid duplicates based on current count)
      const currentLogCount = state.logs.length;
      const totalWasmLogs = infoMessages.size() + warningMessages.size() + 
                           errorMessages.size() + debugMessages.size();
      
      if (totalWasmLogs > currentLogCount) {
        const newLogs = newEntries.slice(currentLogCount);
        
        if (newLogs.length > 0) {
          set(state => {
            const allLogs = [...state.logs, ...newLogs];
            
            // Enforce max logs limit
            if (allLogs.length > state.config.maxLogs) {
              allLogs.splice(0, allLogs.length - state.config.maxLogs);
            }
            
            // Call batch callbacks
            state.batchCallbacks.forEach(callback => callback(newLogs));
            
            return { logs: allLogs };
          });
        }
      }
    } catch (error) {
      console.error('Failed to sync with WASM terminal logger:', error);
    }
  },

  // Setup WASM callbacks for real-time logging
  setupWasmCallbacks: (wasmModule) => {
    if (!wasmModule) return;
    
    try {
      const store = get();
      
      // Set up real-time callbacks
      wasmModule.set_info_callback((message: string) => {
        store.addLog(LogLevel.INFO, message);
      });
      
      wasmModule.set_warning_callback((message: string) => {
        store.addLog(LogLevel.WARNING, message);
      });
      
      wasmModule.set_error_callback((message: string) => {
        store.addLog(LogLevel.ERROR, message);
      });
      
      wasmModule.set_debug_callback((message: string) => {
        store.addLog(LogLevel.DEBUG, message);
      });
    } catch (error) {
      console.error('Failed to setup WASM callbacks:', error);
    }
  },
  
  fetchWasmLogs: (wasmModule) => {
    get().syncWithWasm(wasmModule);
  },
  
  // Getters
  getLogs: () => get().logs,
  
  getFilteredLogs: () => {
    const state = get();
    return state.logs.filter(log => log.level >= state.config.filterLevel);
  },
  
  getLogCount: () => get().logs.length
}));

/**
 * Hook for console logger with automatic cleanup
 */
export function useConsoleLogger() {
  const store = useConsoleLoggerStore();
  
  // Add default console callback if not already added
  useEffect(() => {
    const unsubscribe = store.addCallback((entry) => {
      const levelNames = ['DEBUG', 'INFO', 'WARNING', 'ERROR'];
      const levelName = levelNames[entry.level] || 'UNKNOWN';
      const timestamp = new Date(entry.timestamp * 1000).toISOString();
      console.log(`[${levelName}] ${timestamp} Frame:${entry.frame} ${entry.message}`);
    });
    
    return unsubscribe;
  }, [store.addCallback]); // Only depend on the specific function, not the entire store
  
  return store;
}

/**
 * Convenience functions for logging
 */
export const consoleLogger = {
  debug: (message: string) => useConsoleLoggerStore.getState().addLog(LogLevel.DEBUG, message),
  info: (message: string) => useConsoleLoggerStore.getState().addLog(LogLevel.INFO, message),
  warning: (message: string) => useConsoleLoggerStore.getState().addLog(LogLevel.WARNING, message),
  error: (message: string) => useConsoleLoggerStore.getState().addLog(LogLevel.ERROR, message),
  clear: () => useConsoleLoggerStore.getState().clearLogs(),
  sync: (wasmModule: any) => useConsoleLoggerStore.getState().syncWithWasm(wasmModule),
  setupCallbacks: (wasmModule: any) => useConsoleLoggerStore.getState().setupWasmCallbacks(wasmModule)
};
