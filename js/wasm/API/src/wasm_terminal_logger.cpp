/*
 *  wasm_terminal_logger.cpp
 *  WASM Terminal Logger Implementation
 *
 *  Created by John Delaney on 12/29/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "../include/wasm_terminal_logger.h"
#include <emscripten/val.h>
#include <emscripten/bind.h>
#include <cstddef>

namespace gaudi {

//////////////////////////
// Callback storage and functions
//////////////////////////

// Callback storage using emscripten::val instead of function pointers
static emscripten::val info_callback = emscripten::val::null();
static emscripten::val warning_callback = emscripten::val::null();
static emscripten::val error_callback = emscripten::val::null();
static emscripten::val debug_callback = emscripten::val::null();

// Set logging callbacks from JavaScript
void set_info_callback(const emscripten::val& callback) {
    info_callback = callback;
}

void set_warning_callback(const emscripten::val& callback) {
    warning_callback = callback;
}

void set_error_callback(const emscripten::val& callback) {
    error_callback = callback;
}

void set_debug_callback(const emscripten::val& callback) {
    debug_callback = callback;
}

//////////////////////////
// Terminal Logger Implementation
//////////////////////////

bool wasm_terminal_logger::instance_flag = false;
wasm_terminal_logger *wasm_terminal_logger::global_instance = nullptr;

wasm_terminal_logger &wasm_terminal_logger::get_instance() {
  static wasm_terminal_logger logger;
  if (!logger.initialized()) {
    logger.initialized() = true;
  }
  return logger;
}

void wasm_terminal_logger::log_info(const std::string& message) {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  logger._info_messages.push_back(message);
  
  if (!info_callback.isNull()) {
    info_callback(emscripten::val(message));
  } else {
    std::cout << "[INFO] " << message << std::endl;
  }
}

void wasm_terminal_logger::log_warning(const std::string& message) {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  logger._warning_messages.push_back(message);
  
  if (!warning_callback.isNull()) {
    warning_callback(emscripten::val(message));
  } else {
    std::cout << "[WARNING] " << message << std::endl;
  }
}

void wasm_terminal_logger::log_error(const std::string& message) {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  logger._error_messages.push_back(message);
  
  if (!error_callback.isNull()) {
    error_callback(emscripten::val(message));
  } else {
    std::cout << "[ERROR] " << message << std::endl;
  }
}

void wasm_terminal_logger::log_debug(const std::string& message) {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  logger._debug_messages.push_back(message);
  
  if (!debug_callback.isNull()) {
    debug_callback(emscripten::val(message));
  } else {
    std::cout << "[DEBUG] " << message << std::endl;
  }
}

void wasm_terminal_logger::clear() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  logger._info_messages.clear();
  logger._warning_messages.clear();
  logger._error_messages.clear();
  logger._debug_messages.clear();
}

void wasm_terminal_logger::render() {
  // For terminal logging, render might just flush to console
  // or prepare data for JavaScript consumption
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  // Implementation depends on how you want to expose logs to JavaScript
}

// Global API functions for WASM modules
emscripten::val get_info_messages() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return emscripten::val::array(logger._info_messages);
}

emscripten::val get_warning_messages() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return emscripten::val::array(logger._warning_messages);
}

emscripten::val get_error_messages() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return emscripten::val::array(logger._error_messages);
}

emscripten::val get_debug_messages() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return emscripten::val::array(logger._debug_messages);
}

int get_info_count() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return logger._info_messages.size();
}

int get_warning_count() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return logger._warning_messages.size();
}

int get_error_count() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return logger._error_messages.size();
}

int get_debug_count() {
  wasm_terminal_logger &logger = wasm_terminal_logger::get_instance();
  return logger._debug_messages.size();
}

// EMSCRIPTEN_BINDINGS for terminal logger
EMSCRIPTEN_BINDINGS(terminal_logger) {
  // Data access functions
  emscripten::function("get_info_messages", &get_info_messages);
  emscripten::function("get_warning_messages", &get_warning_messages);
  emscripten::function("get_error_messages", &get_error_messages);
  emscripten::function("get_debug_messages", &get_debug_messages);
  emscripten::function("get_info_count", &get_info_count);
  emscripten::function("get_warning_count", &get_warning_count);
  emscripten::function("get_error_count", &get_error_count);
  emscripten::function("get_debug_count", &get_debug_count);
  
  // Logging functions
  emscripten::function("log_info", &wasm_terminal_logger::log_info);
  emscripten::function("log_warning", &wasm_terminal_logger::log_warning);
  emscripten::function("log_error", &wasm_terminal_logger::log_error);
  emscripten::function("log_debug", &wasm_terminal_logger::log_debug);
  
  // Control functions
  emscripten::function("clear_logs", &wasm_terminal_logger::clear);
  emscripten::function("render_logs", &wasm_terminal_logger::render);
  
  // Callback functions
  emscripten::function("set_info_callback", &set_info_callback);
  emscripten::function("set_warning_callback", &set_warning_callback);
  emscripten::function("set_error_callback", &set_error_callback);
  emscripten::function("set_debug_callback", &set_debug_callback);
}

} // namespace gaudi 