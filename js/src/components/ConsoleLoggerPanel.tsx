import React from 'react';
import { ConsoleLoggerRenderer } from './ConsoleLoggerRenderer';
import { consoleLogger, LogLevel } from '../stores/consoleLoggerStore';

interface ConsoleLoggerPanelProps {
  title?: string;
  height?: number;
  className?: string;
  showTimestamps?: boolean;
  showLevelIcons?: boolean;
  autoScroll?: boolean;
  filterLevel?: LogLevel;
}

export const ConsoleLoggerPanel: React.FC<ConsoleLoggerPanelProps> = ({
  title = "Console Logger Output",
  height = 400,
  className = "",
  showTimestamps = true,
  showLevelIcons = true,
  autoScroll = true,
  filterLevel = LogLevel.DEBUG
}) => {
  const clearLogs = () => {
    consoleLogger.clear();
  };

  return (
    <div className={`flex flex-col ${className}`}>
      
      {/* Console Logger Renderer */}
      <div className="border border-gray-300 rounded overflow-hidden">
        <ConsoleLoggerRenderer
          height={height}
          showTimestamps={showTimestamps}
          showLevelIcons={showLevelIcons}
          autoScroll={autoScroll}
          filterLevel={filterLevel}
        />
      </div>
      
      {/* Footer with stats or controls if needed */}
      <div className="mt-2 text-xs text-gray-500">
        Console output from WASM module (simplified)
      </div>
    </div>
  );
};
