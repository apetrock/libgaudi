import React, { useEffect, useRef, useMemo } from 'react';
import { useConsoleLogger, LogLevel, ConsoleLogEntry } from '../stores/consoleLoggerStore';

/**
 * Props for the ConsoleLoggerRenderer component
 */
interface ConsoleLoggerRendererProps {
  className?: string;
  height?: number;
  maxHeight?: number;
  showTimestamps?: boolean;
  showFrameNumbers?: boolean;
  showLevelIcons?: boolean;
  autoScroll?: boolean;
  filterLevel?: LogLevel;
}

/**
 * Get display color for each log level
 */
const getLogLevelColor = (level: LogLevel): string => {
  switch (level) {
    case LogLevel.DEBUG:
      return 'text-gray-500';
    case LogLevel.INFO:
      return 'text-blue-600';
    case LogLevel.WARNING:
      return 'text-yellow-600';
    case LogLevel.ERROR:
      return 'text-red-600';
    default:
      return 'text-gray-700';
  }
};

/**
 * Get display icon for each log level
 */
const getLogLevelIcon = (level: LogLevel): string => {
  switch (level) {
    case LogLevel.DEBUG:
      return '🐛';
    case LogLevel.INFO:
      return 'ℹ️';
    case LogLevel.WARNING:
      return '⚠️';
    case LogLevel.ERROR:
      return '❌';
    default:
      return '•';
  }
};

/**
 * Get display name for each log level
 */
const getLogLevelName = (level: LogLevel): string => {
  switch (level) {
    case LogLevel.DEBUG:
      return 'DEBUG';
    case LogLevel.INFO:
      return 'INFO';
    case LogLevel.WARNING:
      return 'WARNING';
    case LogLevel.ERROR:
      return 'ERROR';
    default:
      return 'UNKNOWN';
  }
};

/**
 * Format timestamp for display
 */
const formatTimestamp = (timestamp: number): string => {
  const date = new Date(timestamp * 1000);
  return date.toLocaleTimeString('en-US', { 
    hour12: false, 
    hour: '2-digit', 
    minute: '2-digit', 
    second: '2-digit',
    fractionalSecondDigits: 3
  });
};

/**
 * React component for rendering console logs from the WASM logger
 */
export function ConsoleLoggerRenderer({
  className = '',
  height = 300,
  maxHeight,
  showTimestamps = true,
  showFrameNumbers = true,
  showLevelIcons = true,
  autoScroll = true,
  filterLevel = LogLevel.DEBUG
}: ConsoleLoggerRendererProps) {
  const logger = useConsoleLogger();
  const scrollRef = useRef<HTMLDivElement>(null);
  
  // Filter logs based on level - memoized to prevent unnecessary re-renders
  const filteredLogs = useMemo(() => {
    return logger.logs.filter(log => log.level >= filterLevel);
  }, [logger.logs, filterLevel]);
  
  // Auto-scroll to bottom when new logs arrive
  useEffect(() => {
    if (autoScroll && scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
    }
  }, [filteredLogs, autoScroll]);
  
  return (
    <div className={`console-logger-renderer ${className}`}>
      {/* Header */}
      <div className="flex justify-between items-center p-2 bg-gray-100 border-b">
        <span className="font-semibold text-sm">
          Console Logger ({filteredLogs.length} logs)
        </span>
        <div className="flex gap-2">
          <button
            onClick={logger.clearLogs}
            className="px-2 py-1 text-xs bg-red-500 text-white rounded hover:bg-red-600"
          >
            Clear
          </button>
        </div>
      </div>
      
      {/* Logs container */}
      <div
        ref={scrollRef}
        className="overflow-y-auto bg-black text-green-400 font-mono text-xs"
        style={{ 
          height: maxHeight ? 'auto' : height,
          maxHeight: maxHeight || 'none'
        }}
      >
        {filteredLogs.length === 0 ? (
          <div className="p-4 text-gray-500 text-center">
            No logs to display
          </div>
        ) : (
          filteredLogs.map((log, index) => (
            <LogEntry
              key={`${log.timestamp}-${index}`}
              log={log}
              showTimestamp={showTimestamps}
              showFrameNumber={showFrameNumbers}
              showLevelIcon={showLevelIcons}
            />
          ))
        )}
      </div>
    </div>
  );
}

/**
 * Individual log entry component
 */
interface LogEntryProps {
  log: ConsoleLogEntry;
  showTimestamp: boolean;
  showFrameNumber: boolean;
  showLevelIcon: boolean;
}

function LogEntry({ log, showTimestamp, showFrameNumber, showLevelIcon }: LogEntryProps) {
  const levelColor = getLogLevelColor(log.level);
  const levelIcon = getLogLevelIcon(log.level);
  const levelName = getLogLevelName(log.level);
  
  return (
    <div className="px-2 py-1 border-b border-gray-800 hover:bg-gray-900">
      <div className="flex items-start gap-2">
        {/* Level indicator */}
        <span className={`${levelColor} min-w-0 flex-shrink-0`}>
          {showLevelIcon && <span className="mr-1">{levelIcon}</span>}
          <span className="font-semibold">[{levelName}]</span>
        </span>
        
        {/* Timestamp */}
        {showTimestamp && (
          <span className="text-gray-500 text-xs min-w-0 flex-shrink-0">
            {formatTimestamp(log.timestamp)}
          </span>
        )}
        
        {/* Frame number */}
        {showFrameNumber && (
          <span className="text-purple-400 text-xs min-w-0 flex-shrink-0">
            Frame:{log.frame}
          </span>
        )}
        
        {/* Message */}
        <span className="text-green-400 break-words flex-1 min-w-0">
          {log.message}
        </span>
      </div>
    </div>
  );
}

/**
 * Compact console logger for smaller spaces
 */
export function CompactConsoleLogger({ filterLevel = LogLevel.INFO }: { filterLevel?: LogLevel }) {
  const logger = useConsoleLogger();
  const filteredLogs = logger.logs.filter(log => log.level >= filterLevel);
  const recentLogs = filteredLogs.slice(-5); // Show only last 5 logs
  
  return (
    <div className="compact-console-logger bg-gray-900 text-green-400 p-2 rounded">
      <div className="text-xs font-semibold text-gray-400 mb-1">
        Console ({filteredLogs.length})
      </div>
      {recentLogs.map((log, index) => (
        <div key={`${log.timestamp}-${index}`} className="text-xs truncate">
          <span className={getLogLevelColor(log.level)}>
            [{getLogLevelName(log.level)}]
          </span>
          <span className="ml-1">{log.message}</span>
        </div>
      ))}
    </div>
  );
}

/**
 * Console logger stats component
 */
export function ConsoleLoggerStats() {
  const logger = useConsoleLogger();
  
  const debugCount = logger.logs.filter(log => log.level === LogLevel.DEBUG).length;
  const infoCount = logger.logs.filter(log => log.level === LogLevel.INFO).length;
  const warningCount = logger.logs.filter(log => log.level === LogLevel.WARNING).length;
  const errorCount = logger.logs.filter(log => log.level === LogLevel.ERROR).length;
  
  return (
    <div className="console-logger-stats flex gap-4 text-xs">
      <span className="text-gray-500">Debug: {debugCount}</span>
      <span className="text-blue-600">Info: {infoCount}</span>
      <span className="text-yellow-600">Warning: {warningCount}</span>
      <span className="text-red-600">Error: {errorCount}</span>
    </div>
  );
}
