#pragma once

#include <memory>

#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vermeer/medial_axis_project.hpp"
#include "lewitt/application.h"
#include "lewitt/gpu_session.hpp"

namespace gaudi {
namespace vermeer {

class medial_duchamp_host {
public:
  explicit medial_duchamp_host(duchamp::demo_trait::ptr demo) : _demo(std::move(demo)) {}

  int run(int width = 1280, int height = 720) {
    lewitt::app_runner app;
    _frame = 0;

    app.set_project_renderer(
        [this](lewitt::gpu_context &ctx) {
          return std::make_unique<medial_axis_project>(_demo);
        });

    if (!app.onInit(width, height))
      return 1;

    while (app.isRunning()) {
      app.onFrame(static_cast<uint>(_frame++));
    }

    app.onFinish();
    return 0;
  }

private:
  duchamp::demo_trait::ptr _demo;
  int _frame = 0;
};

} // namespace vermeer
} // namespace gaudi
