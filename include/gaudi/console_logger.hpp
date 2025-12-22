#ifndef __GAUDI_CONSOLE_LOGGER__
#define __GAUDI_CONSOLE_LOGGER__

#include "gaudi/common.h"
#include <string>
#include <sstream>

#ifdef __EMSCRIPTEN__
#include "../../js/wasm/API/include/wasm_terminal_logger.h"
#else
#include <iostream>
#endif

namespace gaudi
{
  namespace console_logger
  {
    // Internal logging functions for the stream implementation
    namespace internal
    {
      inline void log_debug(const std::string& message)
      {
#ifdef __EMSCRIPTEN__
        wasm_terminal_logger::log_debug(message);
#else
        std::cout << "[DEBUG] " << message << std::endl;
#endif
      }

      inline void log_info(const std::string& message)
      {
#ifdef __EMSCRIPTEN__
        wasm_terminal_logger::log_info(message);
#else
        std::cout << "[INFO] " << message << std::endl;
#endif
      }

      inline void log_warning(const std::string& message)
      {
#ifdef __EMSCRIPTEN__
        wasm_terminal_logger::log_warning(message);
#else
        std::cout << "[WARNING] " << message << std::endl;
#endif
      }

      inline void log_error(const std::string& message)
      {
#ifdef __EMSCRIPTEN__
        wasm_terminal_logger::log_error(message);
#else
        std::cout << "[ERROR] " << message << std::endl;
#endif
      }

      // Clear all logs
      inline void clear()
      {
#ifdef __EMSCRIPTEN__
        wasm_terminal_logger::clear();
#endif
      }
    }

    // Stream-style logging (similar to existing logger)
    class LogStream
    {
    private:
      std::string level;
      std::string buffer;

    public:
      LogStream(const std::string& log_level) : level(log_level) {}

      template<typename T>
      LogStream& operator<<(const T& value)
      {
        std::ostringstream oss;
        oss << value;
        buffer += oss.str();
        return *this;
      }

      // Handle std::endl and other manipulators
      LogStream& operator<<(std::ostream& (*manipulator)(std::ostream&))
      {
        if (manipulator == static_cast<std::ostream& (*)(std::ostream&)>(std::endl))
        {
          // Flush the buffer
          if (level == "DEBUG") internal::log_debug(buffer);
          else if (level == "INFO") internal::log_info(buffer);
          else if (level == "WARNING") internal::log_warning(buffer);
          else if (level == "ERROR") internal::log_error(buffer);
          
          buffer.clear();
        }
        return *this;
      }

      ~LogStream()
      {
        // Flush any remaining buffer content
        if (!buffer.empty())
        {
          if (level == "DEBUG") internal::log_debug(buffer);
          else if (level == "INFO") internal::log_info(buffer);
          else if (level == "WARNING") internal::log_warning(buffer);
          else if (level == "ERROR") internal::log_error(buffer);
        }
      }
    };

    // Stream-style logging objects
    static LogStream debug_stream("DEBUG");
    static LogStream info_stream("INFO");
    static LogStream warning_stream("WARNING");
    static LogStream error_stream("ERROR");

    // Clean streaming interface
    static LogStream& debug = debug_stream;
    static LogStream& info = info_stream;
    static LogStream& warning = warning_stream;
    static LogStream& error = error_stream;

  } // namespace console_logger
} // namespace gaudi

#endif // __GAUDI_CONSOLE_LOGGER__
