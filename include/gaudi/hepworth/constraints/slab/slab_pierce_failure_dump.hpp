#ifndef GAUDI_HEPWORTH_SLAB_PIERCE_FAILURE_DUMP_HPP
#define GAUDI_HEPWORTH_SLAB_PIERCE_FAILURE_DUMP_HPP

#include "slab_pierce_stencil.hpp"

#include <cstdlib>
#include <fstream>
#include <string>

namespace gaudi {
namespace hepworth {
namespace block {

struct pierce_failure_record {
  index_t foot = -1;
  pierce_failure_reason reason = pierce_failure_reason::ok;
  vec3 x = vec3::Zero();
  vec3 n_in = vec3::Zero();
  real area = 0.0;
  size_t valence = 0;
  index_t hit_face = -1;
  real ray_t = -1.0;
  ray_miss_detail miss_detail = ray_miss_detail::none;
};

inline std::string pierce_failure_dump_path() {
  const char *env = std::getenv("GAUDI_PIERCE_FAILURE_DUMP");
  if (env != nullptr && env[0] != '\0')
    return std::string(env);
  return "dump/pierce_failures.jsonl";
}

inline void append_pierce_failure_record(const pierce_failure_record &rec,
                                         const std::string &path) {
  std::ofstream out(path, std::ios::app);
  if (!out)
    return;
  out << "{\"foot\":" << rec.foot << ",\"reason\":\""
      << pierce_failure_reason_name(rec.reason) << "\",\"x\":["
      << rec.x.x() << "," << rec.x.y() << "," << rec.x.z() << "],\"n_in\":["
      << rec.n_in.x() << "," << rec.n_in.y() << "," << rec.n_in.z()
      << "],\"area\":" << rec.area << ",\"valence\":" << rec.valence
      << ",\"hit_face\":" << rec.hit_face << ",\"ray_t\":" << rec.ray_t
      << ",\"miss_detail\":\""
      << ray_miss_detail_name(rec.miss_detail) << "\""
      << "}\n";
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // GAUDI_HEPWORTH_SLAB_PIERCE_FAILURE_DUMP_HPP
