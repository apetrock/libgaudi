#pragma once

#include <cstdlib>

#include "lewitt/pipeline_graph.hpp"

namespace gaudi {
namespace vermeer {

using render_graph_type = lewitt::render_graph_type;
using ffmpeg_record_config = lewitt::ffmpeg_record_config;

struct vermeer_config {
  render_graph_type renderer_type = render_graph_type::ssao_deferred;
  ffmpeg_record_config record{};
};

inline void apply_record_env(vermeer_config &config) {
  if (std::getenv("VERMEER_RECORD")) {
    config.record.enabled = true;
  }
  if (const char *path = std::getenv("VERMEER_RECORD_PATH")) {
    config.record.output_path = path;
  }
}

} // namespace vermeer
} // namespace gaudi
