#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <vector>
#if defined(WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#ifdef ERROR
#undef ERROR
#endif
#endif

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/duchamp/darboux_cyclide_medial.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"

namespace {

constexpr int kSingleFitVertex = 0;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

gaudi::duchamp::cyclide_medial_params demo_params() {
  gaudi::duchamp::cyclide_medial_params params;
  params.l0_scale = 0.1;
  params.fit_p = 3.0;
  params.max_iters = 12;
  params.tol = 1e-8;
  params.max_travel_scale = 12.0;
  return params;
}

class Scene : public gg::Scene {
public:
  static ScenePtr create(const gaudi::duchamp::cyclide_medial_params &params) {
    return std::make_shared<Scene>(params);
  }

  Scene(const gaudi::duchamp::cyclide_medial_params &params)
      : gg::Scene(), __params(params) {
    __frame =
        gaudi::test::make_torus_frame(gaudi::vec3(0.7, -0.45, 0.3),
                                      gaudi::vec3(0.25, 0.65, 0.9));
    gaudi::test::TorusMesh torus =
        gaudi::test::make_offset_torus_shell(36, 20, __major_radius,
                                             __minor_radius, __frame);
    __M = torus.shell;

    gaudi::asawa::shell::shell &M = *__M;
    __single_fit_vertex =
        std::clamp(kSingleFitVertex, 0, static_cast<int>(M.vert_count()) - 1);
    __candidates =
        gaudi::duchamp::compute_single_fit_cyclide_medial_candidates(
            M, __single_fit_vertex, __params, &__stats, &__single_fit_Q);
    std::cerr << "single-fit cyclide vertex=" << __single_fit_vertex
              << std::endl;
    gaudi::duchamp::print_cyclide_medial_stats(std::cerr, __stats);

    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    refresh_mesh_buffer();
    mSceneObjects.push_back(_objs[0]);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  void refresh_mesh_buffer() {
    gg::fillBuffer_ref(*__M, _objs[0], gg::colorRGB(0.72, 0.74, 0.78, 1.0));
  }

  void draw_medial_rays() {
    const gaudi::vec4 ray_color(1.0, 0.92, 0.15, 1.0);
    const gaudi::vec4 trace_color(1.0, 0.35, 0.05, 1.0);
    const gaudi::vec4 axis_color(0.35, 0.35, 0.35, 1.0);
    const gaudi::vec4 fit_color(0.15, 1.0, 0.25, 1.0);
    const gaudi::vec4 level_color(0.1, 0.9, 1.0, 1.0);

    gg::geometry_logger::line(
        __frame.center - 1.7 * __major_radius * __frame.z_axis,
        __frame.center + 1.7 * __major_radius * __frame.z_axis, axis_color);

    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);

    const gaudi::vec3 fit_p = x[__single_fit_vertex];
    const gaudi::real marker = 0.08 * gaudi::asawa::shell::avg_length(M, x);
    gg::geometry_logger::line(fit_p - marker * __frame.x_axis,
                              fit_p + marker * __frame.x_axis, fit_color);
    gg::geometry_logger::line(fit_p - marker * __frame.y_axis,
                              fit_p + marker * __frame.y_axis, fit_color);
    gg::geometry_logger::line(fit_p - marker * __frame.z_axis,
                              fit_p + marker * __frame.z_axis, fit_color);
    draw_fit_level_slice(fit_p, level_color);

    for (const auto &cand : __candidates) {
      const gaudi::vec3 p = x[cand.vertex];
      for (size_t i = 1; i < cand.trace_world.size(); ++i) {
        gg::geometry_logger::line(cand.trace_world[i - 1], cand.trace_world[i],
                                  trace_color);
      }
      if (cand.accepted) {
        gg::geometry_logger::line(p, cand.center_world, ray_color);
      }
    }
  }

  void draw_fit_level_slice(const gaudi::vec3 &origin,
                            const gaudi::vec4 &color) {
    gaudi::vec3 radial = origin - __frame.center;
    radial -= radial.dot(__frame.z_axis) * __frame.z_axis;
    if (radial.norm() < 1e-10) {
      radial = __frame.x_axis;
    } else {
      radial.normalize();
    }
    const gaudi::vec3 binormal = __frame.z_axis;
    const gaudi::real extent = 1.4 * __minor_radius;
    constexpr int N = 48;
    const gaudi::real h = 2.0 * extent / gaudi::real(N);

    auto sample = [&](int i, int j) {
      const gaudi::real u = -extent + h * gaudi::real(i);
      const gaudi::real v = -extent + h * gaudi::real(j);
      const gaudi::vec3 p = u * radial + v * binormal;
      return gaudi::albers::eval_darboux(__single_fit_Q, p);
    };
    auto world = [&](gaudi::real u, gaudi::real v) {
      return origin + u * radial + v * binormal;
    };
    auto interp = [](gaudi::real a, gaudi::real b) {
      const gaudi::real denom = a - b;
      if (std::abs(denom) < 1e-12) {
        return gaudi::real(0.5);
      }
      return std::clamp(a / denom, gaudi::real(0.0), gaudi::real(1.0));
    };

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        const gaudi::real f00 = sample(i, j);
        const gaudi::real f10 = sample(i + 1, j);
        const gaudi::real f11 = sample(i + 1, j + 1);
        const gaudi::real f01 = sample(i, j + 1);

        std::vector<gaudi::vec3> hits;
        if ((f00 <= 0.0) != (f10 <= 0.0)) {
          const gaudi::real t = interp(f00, f10);
          hits.push_back(world(-extent + h * (gaudi::real(i) + t),
                               -extent + h * gaudi::real(j)));
        }
        if ((f10 <= 0.0) != (f11 <= 0.0)) {
          const gaudi::real t = interp(f10, f11);
          hits.push_back(world(-extent + h * gaudi::real(i + 1),
                               -extent + h * (gaudi::real(j) + t)));
        }
        if ((f11 <= 0.0) != (f01 <= 0.0)) {
          const gaudi::real t = interp(f11, f01);
          hits.push_back(world(-extent + h * (gaudi::real(i + 1) - t),
                               -extent + h * gaudi::real(j + 1)));
        }
        if ((f01 <= 0.0) != (f00 <= 0.0)) {
          const gaudi::real t = interp(f01, f00);
          hits.push_back(world(-extent + h * gaudi::real(i),
                               -extent + h * (gaudi::real(j + 1) - t)));
        }

        for (size_t k = 1; k < hits.size(); k += 2) {
          gg::geometry_logger::line(hits[k - 1], hits[k], color);
        }
      }
    }
  }

  virtual void onAnimate(int /*frame*/) {}

  virtual void onDraw(gg::Viewer &viewer) {
    draw_medial_rays();
    gg::geometry_logger::render();
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

private:
  gaudi::asawa::shell::shell::ptr __M;
  gaudi::test::TorusFrame __frame;
  gaudi::duchamp::cyclide_medial_params __params;
  gaudi::duchamp::cyclide_medial_stats __stats;
  std::vector<gaudi::duchamp::cyclide_medial_candidate> __candidates;
  gaudi::albers::vec14 __single_fit_Q = gaudi::albers::vec14::Zero();
  int __single_fit_vertex = 0;
  gaudi::real __major_radius = 1.25;
  gaudi::real __minor_radius = 0.35;
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height,
                       const gaudi::duchamp::cyclide_medial_params &params) {
    return std::make_shared<App>(width, height, params);
  }

  App(int width, int height,
      const gaudi::duchamp::cyclide_medial_params &params)
      : gg::SimpleApp(width, height, 4.0, false, "darboux_cyclide_medial_") {
    scene = Scene::create(params);
    this->setScene(scene);
    this->initUI();
  }

  void initUI() { performLayout(); }

  ~App() override = default;

  ScenePtr scene;
};

} // namespace

int main() {
  try {
    const auto params = demo_params();
    nanogui::init();
    AppPtr app = App::create(1280, 740, params);
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();
  } catch (const std::runtime_error &e) {
    std::cerr << "Caught a fatal error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}
