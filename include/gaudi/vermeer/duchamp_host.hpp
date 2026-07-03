#pragma once

#include <utility>

#include "gaudi/duchamp/demo_trait.hpp"
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
    lewitt::app_runner app;
    _frame = 0;

    app.set_project_renderer([this](lewitt::gpu_context &ctx) {
      return std::make_unique<duchamp_project>(_demo, _config);
    });

    if (!app.onInit(initial_width, initial_height)) {
      return 1;
    }

    while (app.isRunning()) {
      LEWITT_PERF_SCOPE_PATH("gaudi::vermeer::duchamp_host::frame");
      app.onFrame(static_cast<uint>(_frame++));
    }

    app.onFinish();
    return 0;
  }

private:
  duchamp::demo_trait::ptr _demo;
  vermeer_config _config;
  int _frame = 0;
};

inline int vermeer(duchamp::demo_trait::ptr demo, vermeer_config config = {},
                   uint32_t initial_width = 1280, uint32_t initial_height = 720) {
  apply_record_env(config);
  duchamp_host host(std::move(demo), config);
  return host.run(initial_width, initial_height);
}

} // namespace vermeer
} // namespace gaudi
