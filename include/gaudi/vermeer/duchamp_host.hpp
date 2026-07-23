#pragma once

// Vermeer Duchamp run harness: Space = pause/resume, '.' = one step while paused.

#include <memory>
#include <utility>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/exception_diagnostics.hpp"
#include "gaudi/vermeer/duchamp_playback.hpp"
#include "gaudi/vermeer/duchamp_project.hpp"
#include "gaudi/vermeer/vermeer_config.hpp"
#include "lewitt/application.h"
#include "lewitt/performance.hpp"

namespace gaudi {
namespace vermeer {

class duchamp_host {
public:
  explicit duchamp_host(duchamp::demo_trait::ptr demo, vermeer_config config = {})
      : _demo(std::move(demo)), _config(config) {}

  int run(uint32_t initial_width = 1280, uint32_t initial_height = 720) {
    // Apply VERMEER_RECORD* here so demos that construct duchamp_host
    // directly (not via vermeer()) still pick up gaudi.py defaults.
    apply_record_env(_config);

    lewitt::app_runner app;
    auto playback = std::make_shared<duchamp_playback>();
    playback->set_window_title(_demo ? _demo->name() : "Vermeer");

    app.set_project_renderer([this, playback](lewitt::gpu_context &ctx) {
      return std::make_unique<duchamp_project>(_demo, _config, playback);
    });

    if (!app.onInit(initial_width, initial_height)) {
      return 1;
    }

    while (app.isRunning()) {
      LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_host::frame");
      app.onFrame(static_cast<uint>(playback->sim_frame.load()));
    }

    app.onFinish();
    return 0;
  }

private:
  duchamp::demo_trait::ptr _demo;
  vermeer_config _config;
};

inline int vermeer(duchamp::demo_trait::ptr demo, vermeer_config config = {},
                   uint32_t initial_width = 1280, uint32_t initial_height = 720) {
  gaudi::install_exception_diagnostics();
  apply_record_env(config);
  duchamp_host host(std::move(demo), config);
  return host.run(initial_width, initial_height);
}

} // namespace vermeer
} // namespace gaudi
