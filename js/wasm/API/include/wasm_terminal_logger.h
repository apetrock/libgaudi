/*
 *  wasm_terminal_logger.h
 *  WASM Terminal Logger Interface
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include "gaudi/common.h"
#include <emscripten/val.h>

#ifndef __WASM_TERMINAL_LOGGER__
#define __WASM_TERMINAL_LOGGER__

namespace gaudi {

// Set logging callbacks from JavaScript using emscripten::val
void set_info_callback(const emscripten::val& callback);
void set_warning_callback(const emscripten::val& callback);
void set_error_callback(const emscripten::val& callback);
void set_debug_callback(const emscripten::val& callback);

class wasm_terminal_logger {

public:
  static wasm_terminal_logger &get_instance();
  
  // Terminal logging functions
  static void log_info(const std::string& message);
  static void log_warning(const std::string& message);
  static void log_error(const std::string& message);
  static void log_debug(const std::string& message);
  
  // Utility functions
  static void clear();
  static void render();

  bool &initialized() { return instance_flag; }
  bool initialized() const { return instance_flag; }
  
  // Storage for log messages
  std::vector<std::string> _info_messages;
  std::vector<std::string> _warning_messages;
  std::vector<std::string> _error_messages;
  std::vector<std::string> _debug_messages;

private:
  wasm_terminal_logger() {
  }

  wasm_terminal_logger(const wasm_terminal_logger &);
  wasm_terminal_logger &operator=(const wasm_terminal_logger &);

  static wasm_terminal_logger *global_instance;
  static bool instance_flag;
};

// Global terminal logger API functions for WASM modules
emscripten::val get_info_messages();
emscripten::val get_warning_messages();
emscripten::val get_error_messages();
emscripten::val get_debug_messages();
int get_info_count();
int get_warning_count();
int get_error_count();
int get_debug_count();

} // namespace gaudi

#endif 