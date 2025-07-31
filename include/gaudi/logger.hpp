#ifndef __LIBGAUDI_LOGGER__
#define __LIBGAUDI_LOGGER__

#include "common.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace gaudi {
namespace logger {

// Log levels enumeration
enum class LogLevel {
    DEBUG = 0,
    INFO = 1,
    WARNING = 2,
    ERROR = 3
};

// Log entry structure for accumulated logging
struct LogEntry {
    LogLevel level;
    std::string message;
    double timestamp;
    int frame;
};

// System logging functions (for general application logging)
void log_info(const std::string& message);
void log_warning(const std::string& message);
void log_error(const std::string& message);
void log_debug(const std::string& message);

// Accumulated logging functions (similar to geometry_logger)
void add_log(LogLevel level, const std::string& message);
void clear_logs();
const std::vector<LogEntry>& get_logs();
int get_frame_count();
void increment_frame();

// Stream-like logger class implementation
class LogStream {
public:
    enum Level { INFO, WARNING, ERROR, DEBUG };
    
    LogStream(Level level) : m_level(level) {}
    
    // Template for any type - convert to string and log
    template<typename T>
    LogStream& operator<<(const T& value) {
        m_stream << value;
        return *this;
    }
    
    // Handle std::endl
    LogStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        if (manip == static_cast<std::ostream& (*)(std::ostream&)>(std::endl)) {
            flush();
        }
        return *this;
    }
    
    void flush() {
        std::string message = m_stream.str();
        if (!message.empty()) {
            LogLevel level_enum;
            switch (m_level) {
                case INFO: level_enum = LogLevel::INFO; break;
                case WARNING: level_enum = LogLevel::WARNING; break;
                case ERROR: level_enum = LogLevel::ERROR; break;
                case DEBUG: level_enum = LogLevel::DEBUG; break;
            }
            
            // Send to both immediate and accumulated logging
            switch (m_level) {
                case INFO: log_info(message); break;
                case WARNING: log_warning(message); break;
                case ERROR: log_error(message); break;
                case DEBUG: log_debug(message); break;
            }
            add_log(level_enum, message);
            
            m_stream.str("");
        }
    }
    
    ~LogStream() { flush(); }
    
private:
    Level m_level;
    std::stringstream m_stream;
};

// Global logger objects (replace std::cout, std::cerr, etc.)
extern LogStream info;
extern LogStream warning;
extern LogStream error;
extern LogStream debug;

} // namespace logger
} // namespace gaudi

#endif