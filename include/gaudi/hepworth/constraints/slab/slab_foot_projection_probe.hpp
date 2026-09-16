#ifndef GAUDI_HEPWORTH_SLAB_FOOT_PROJECTION_PROBE_HPP
#define GAUDI_HEPWORTH_SLAB_FOOT_PROJECTION_PROBE_HPP

#include "slab_pierce_stencil.hpp"
#include "slab_prism_operators.hpp"

#include "gaudi/hepworth/block/projection_diagnostics.hpp"
#include "gaudi/hepworth/block/slab_prism_log.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace gaudi {
namespace hepworth {
namespace block {

/// Optional override; otherwise GAUDI_PROBE_FOOT env (comma-separated feet ok).
inline index_t &slab_foot_probe_override() {
  static index_t foot = -1;
  return foot;
}

inline void set_slab_foot_probe(index_t foot) { slab_foot_probe_override() = foot; }

inline std::vector<index_t> slab_foot_probe_ids() {
  if (slab_foot_probe_override() >= 0)
    return {slab_foot_probe_override()};
  const char *env = std::getenv("GAUDI_PROBE_FOOT");
  if (env == nullptr || env[0] == '\0' || env[0] == '0')
    return {};
  std::vector<index_t> out;
  std::stringstream ss(env);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty())
      continue;
    out.push_back(static_cast<index_t>(std::stoi(tok)));
  }
  return out;
}

inline bool slab_foot_probe_active(index_t foot) {
  if (!projection_diag_enabled() && !slab_prism_log_enabled())
    return false;
  for (index_t id : slab_foot_probe_ids()) {
    if (id == foot)
      return true;
  }
  return false;
}

inline std::string slab_foot_probe_dump_path(index_t foot) {
  return "dump/probe_foot_" + std::to_string(foot) + ".jsonl";
}

inline void append_probe_line(index_t foot, const std::string &line) {
  const char *env = std::getenv("GAUDI_PROBE_FOOT_DUMP");
  if (env != nullptr && env[0] == '0')
    return;
  const std::string path = slab_foot_probe_dump_path(foot);
  std::error_code ec;
  std::filesystem::create_directories(
      std::filesystem::path(path).parent_path(), ec);
  std::ofstream out(path, std::ios::app);
  if (!out)
    return;
  out << line << "\n";
}

inline std::string vec3_json(const vec3 &v) {
  std::ostringstream os;
  os << "[" << v.x() << "," << v.y() << "," << v.z() << "]";
  return os.str();
}

inline void log_load_vert_world_state(index_t vert, int frame,
                                      const char *phase,
                                      const std::vector<vec3> &x) {
  if (vert < 0 || static_cast<size_t>(vert) >= x.size())
    return;
  if (!slab_foot_probe_active(vert))
    return;
  const vec3 &xi = x[static_cast<size_t>(vert)];
  const std::string msg =
      std::string("load_vert state frame=") + std::to_string(frame) + " phase=" +
      phase + " vert=" + std::to_string(vert) + " x=" + vec3_json(xi) +
      " |x|=" + std::to_string(xi.norm());
  if (projection_diag_enabled())
    projection_diag_log(msg);
  else
    slab_prism_log(msg);
  append_probe_line(vert, "{\"frame\":" + std::to_string(frame) + ",\"phase\":\"" +
                               std::string(phase) + "\",\"x\":" + vec3_json(xi) +
                               ",\"norm\":" + std::to_string(xi.norm()) + "}");
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_FOOT_PROJECTION_PROBE_HPP
