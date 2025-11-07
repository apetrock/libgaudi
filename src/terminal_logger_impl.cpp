#include "gaudi/logger.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>

namespace gaudi {
namespace logger {

// Static storage for accumulated logs (similar to geometry_logger)
static std::vector<LogEntry> s_logs;
static int s_frame_count = 0;

// Get current timestamp in seconds
static double get_timestamp() {
    auto now = std::chrono::steady_clock::now();
    auto duration = now.time_since_epoch();
    return std::chrono::duration<double>(duration).count();
}

void log_info(const std::string& message) {
    std::cout << "[INFO] " << message << std::endl;
}

void log_warning(const std::string& message) {
    std::cout << "[WARNING] " << message << std::endl;
}

void log_error(const std::string& message) {
    std::cerr << "[ERROR] " << message << std::endl;
}

void log_debug(const std::string& message) {
    std::cout << "[DEBUG] " << message << std::endl;
}

// Accumulated logging functions
void add_log(LogLevel level, const std::string& message) {
    LogEntry entry;
    entry.level = level;
    entry.message = message;
    entry.timestamp = get_timestamp();
    entry.frame = s_frame_count;
    s_logs.push_back(entry);
}

void clear_logs() {
    s_logs.clear();
}

const std::vector<LogEntry>& get_logs() {
    return s_logs;
}

int get_frame_count() {
    return s_frame_count;
}

void increment_frame() {
    s_frame_count++;
}

// LogStream class is now defined in the header file

// Global logger stream objects (replace std::cout, std::cerr, etc.)
LogStream info(LogStream::INFO);
LogStream warning(LogStream::WARNING);
LogStream error(LogStream::ERROR);
LogStream debug(LogStream::DEBUG);

} // namespace logger
} // namespace gaudi 