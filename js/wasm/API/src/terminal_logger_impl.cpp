#include "../include/wasm_geometry_logger.h"
#include "../include/wasm_terminal_logger.h"
#include "gaudi/common.h"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/logger.hpp"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace gaudi {
namespace logger {

// Static storage for accumulated logs
static std::vector<LogEntry> accumulated_logs;
static int current_frame = 0;

void log_info(const std::string &message) {
  wasm_terminal_logger::log_info(message);
}

void log_warning(const std::string &message) {
  wasm_terminal_logger::log_warning(message);
}

void log_error(const std::string &message) {
  wasm_terminal_logger::log_error(message);
}

void log_debug(const std::string &message) {
  wasm_terminal_logger::log_debug(message);
}

// Accumulated logging functions implementation
void add_log(LogLevel level, const std::string& message) {
    LogEntry entry;
    entry.level = level;
    entry.message = message;
    entry.timestamp = 0.0; // Could use emscripten_get_now() if needed
    entry.frame = current_frame;
    accumulated_logs.push_back(entry);
}

void clear_logs() {
    accumulated_logs.clear();
}

const std::vector<LogEntry>& get_logs() {
    return accumulated_logs;
}

int get_frame_count() {
    return current_frame;
}

void increment_frame() {
    current_frame++;
}

// LogStream class is now defined in the header file

// Global logger stream objects (replace std::cout, std::cerr, etc.)
LogStream info(LogStream::INFO);
LogStream warning(LogStream::WARNING);
LogStream error(LogStream::ERROR);
LogStream debug(LogStream::DEBUG);

} // namespace logger
} // namespace gaudi
