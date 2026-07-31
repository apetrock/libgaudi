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

  // Auto PCA/OBB camera framing (opt-out with hotkey C).
  // Capture every period; min-jerk ease between consecutive targets.
  bool fit_bbox_auto = true;
  float fit_period_s = 8.0f; // longer tween (~240 frames @ 30fps sim)
  // Retarget hysteresis (see bbox_framing_config).
  float fit_anisotropy_min = 1.35f;
  float fit_keep_alignment = 0.82f;
  float fit_min_score_improvement = 0.12f;
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
