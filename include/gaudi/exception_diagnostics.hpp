#pragma once

// Install a terminate handler that prints the active exception and a native
// backtrace. Sanitizers do not report std::length_error / other thrown
// exceptions — this fills that gap for asan/debug runs.

#include <cstdlib>
#include <exception>
#include <execinfo.h>
#include <iostream>
#include <typeinfo>

namespace gaudi {

inline void print_native_backtrace(std::ostream &os, int max_frames = 64) {
  void *frames[128];
  const int n = backtrace(frames, max_frames);
  char **symbols = backtrace_symbols(frames, n);
  os << "backtrace (" << n << " frames):\n";
  if (!symbols) {
    os << "  (backtrace_symbols failed)\n";
    return;
  }
  for (int i = 0; i < n; ++i)
    os << "  " << symbols[i] << "\n";
  free(symbols);
}

inline void exception_terminate_handler() {
  std::cerr << "\n=== gaudi: std::terminate ===\n";
  if (auto eptr = std::current_exception()) {
    try {
      std::rethrow_exception(eptr);
    } catch (const std::exception &e) {
      std::cerr << "exception type: " << typeid(e).name() << "\n"
                << "what(): " << e.what() << "\n";
    } catch (...) {
      std::cerr << "exception type: (non-std::exception)\n";
    }
  } else {
    std::cerr << "(no active exception)\n";
  }
  print_native_backtrace(std::cerr);
  std::cerr << "=== end terminate ===\n" << std::flush;
  std::abort();
}

inline void install_exception_diagnostics() {
  static bool installed = false;
  if (installed)
    return;
  installed = true;
  std::set_terminate(exception_terminate_handler);
}

} // namespace gaudi
