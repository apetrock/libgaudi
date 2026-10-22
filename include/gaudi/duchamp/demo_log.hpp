#ifndef GAUDI_DUCHAMP_DEMO_LOG_HPP
#define GAUDI_DUCHAMP_DEMO_LOG_HPP

#include <filesystem>
#include <fstream>
#include <iostream>

namespace gaudi {
namespace duchamp {

inline std::ofstream &demo_log_file() {
  static std::ofstream out;
  return out;
}

/// Truncate `path` once. Later `demo_log()` writes go there until the process ends.
inline void open_demo_log(const std::filesystem::path &path) {
  std::error_code ec;
  if (path.has_parent_path())
    std::filesystem::create_directories(path.parent_path(), ec);
  std::ofstream &out = demo_log_file();
  if (out.is_open())
    out.close();
  out.open(path, std::ios::out | std::ios::trunc);
  out << std::unitbuf;
}

/// The braid-demo file when `open_demo_log` succeeded, otherwise the console.
inline std::ostream &demo_log() {
  std::ofstream &out = demo_log_file();
  if (out.is_open())
    return out;
  return std::cout;
}

} // namespace duchamp
} // namespace gaudi

#endif
