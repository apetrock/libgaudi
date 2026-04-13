#include "GaudiMath/typedefs.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
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

#include "gaudi/duchamp/spectral_modes_demo.hpp"

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {
public:
  static ScenePtr create(
      gaudi::duchamp::laplace_modes_band band =
          gaudi::duchamp::laplace_modes_band::largest_magnitude,
      int num_eigenmodes = 32, double mid_band_shift = -1.0,
      double mid_slider_t = -2.0) {
    return std::make_shared<Scene>(band, num_eigenmodes, mid_band_shift,
                                   mid_slider_t);
  }

  explicit Scene(gaudi::duchamp::laplace_modes_band band =
                     gaudi::duchamp::laplace_modes_band::largest_magnitude,
                 int num_eigenmodes = 32, double mid_band_shift = -1.0,
                 double mid_slider_t = -2.0)
      : gg::Scene() {
    initScene(band, num_eigenmodes, mid_band_shift, mid_slider_t);
  }

  void initScene(gaudi::duchamp::laplace_modes_band band, int num_eigenmodes,
                 double mid_band_shift, double mid_slider_t) {
    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);

    __demo = gaudi::duchamp::spectral_modes_demo::create(
        num_eigenmodes, 1e-6, band, mid_band_shift, mid_slider_t);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);

    _mesh_dirty = true;
  }

  void mark_mesh_dirty() { _mesh_dirty = true; }

  void prev_mode() {
    if (!__demo || __demo->mode_count() <= 1)
      return;
    int i = __demo->mode_index() - 1;
    if (i < 0)
      i = __demo->mode_count() - 1;
    __demo->set_mode_index(i);
    mark_mesh_dirty();
  }

  void next_mode() {
    if (!__demo || __demo->mode_count() <= 1)
      return;
    int i = __demo->mode_index() + 1;
    if (i >= __demo->mode_count())
      i = 0;
    __demo->set_mode_index(i);
    mark_mesh_dirty();
  }

  int current_mode() const { return __demo ? __demo->mode_index() : 0; }
  int mode_count() const { return __demo ? __demo->mode_count() : 1; }

  void log_active_mode_stats() const {
    if (__demo)
      __demo->log_current_mode_quality();
  }

  /// Same pattern as `growth_study/main.cpp`: vec4 from demo, then colorRGB with
  /// explicit alpha = 1.0 (albedo only; no vertex alpha).
  void refresh_mesh_buffer() {
    std::vector<gg::colorRGB> colors;
    auto mesh_colors = __demo->get_mesh_colors();
    colors.reserve(mesh_colors.size());
    for (const auto &col : mesh_colors) {
      colors.push_back(gg::colorRGB(col[0], col[1], col[2], 1.0));
    }
    gg::fillBuffer_ref(*__demo->__M, _objs[0], colors);
  }

  /// Matches `arnoldi_test` / `aabb_test`: upload mesh while animating (key W on).
  virtual void onAnimate(int /*frame*/) {
    refresh_mesh_buffer();
    gg::geometry_logger::render();
  }

  /// One upload when mode changes or first frame — no extra nanogui / per-frame
  /// uploads that could fight the deferred pass.
  virtual void onDraw(gg::Viewer &viewer) {
    if (_mesh_dirty) {
      refresh_mesh_buffer();
      _mesh_dirty = false;
    }
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

  gaudi::duchamp::spectral_modes_demo::ptr __demo;

private:
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
  bool _mesh_dirty = true;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(
      int width, int height,
      gaudi::duchamp::laplace_modes_band band =
          gaudi::duchamp::laplace_modes_band::largest_magnitude,
      int num_eigenmodes = 32, double mid_band_shift = -1.0,
      double mid_slider_t = -2.0) {
    return std::make_shared<App>(width, height, band, num_eigenmodes,
                                 mid_band_shift, mid_slider_t);
  }

  App(int width, int height,
      gaudi::duchamp::laplace_modes_band band =
          gaudi::duchamp::laplace_modes_band::largest_magnitude,
      int num_eigenmodes = 32, double mid_band_shift = -1.0,
      double mid_slider_t = -2.0)
      : gg::SimpleApp(width, height, 4.0, false, "spectral_modes_") {
    this->setScene(scene = Scene::create(band, num_eigenmodes, mid_band_shift,
                                          mid_slider_t));
    this->initUI();
  }

  /// Same as `arnoldi_test` / `aabb_test`: layout only, no overlay widgets.
  void initUI() { performLayout(); }

  bool keyboardEvent(int key, int scancode, int action,
                     int modifiers) override {
    if (nanogui::Screen::keyboardEvent(key, scancode, action, modifiers))
      return true;

    if (action == GLFW_PRESS && scene) {
      if (key == GLFW_KEY_COMMA || key == GLFW_KEY_LEFT_BRACKET) {
        scene->prev_mode();
        scene->log_active_mode_stats();
        cerr << "[spectral_modes] mode " << scene->current_mode() << " / "
             << (scene->mode_count() - 1) << endl;
        drawContents();
        return true;
      }
      if (key == GLFW_KEY_PERIOD || key == GLFW_KEY_RIGHT_BRACKET) {
        scene->next_mode();
        scene->log_active_mode_stats();
        cerr << "[spectral_modes] mode " << scene->current_mode() << " / "
             << (scene->mode_count() - 1) << endl;
        drawContents();
        return true;
      }
    }

    return SimpleApp::keyboardEvent(key, scancode, action, modifiers);
  }

  ~App() override = default;

  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    using gaudi::duchamp::laplace_modes_band;
    bool opt_low = false;
    bool opt_high = false;
    bool opt_mid = false;
    bool shift_fraction_set = false;
    int num_modes = 32;
    double mid_shift = -1.0;
    double shift_fraction = 0.0;
    for (int i = 1; i < argc; ++i) {
      const std::string a(argv[i]);
      if (a == "--low-freq" || a == "--smooth")
        opt_low = true;
      else if (a == "--high-freq" || a == "--largest-magnitude")
        opt_high = true;
      else if (a == "--mid" || a == "--mid-freq")
        opt_mid = true;
      else if (a == "--shift" && i + 1 < argc) {
        try {
          mid_shift = std::stod(argv[++i]);
        } catch (const std::exception &) {
          cerr << "[spectral_modes] ignoring bad --shift value\n";
        }
      } else if ((a == "--shift-fraction" || a == "--shift_fraction") &&
                 i + 1 < argc) {
        try {
          shift_fraction = std::stod(argv[++i]);
          shift_fraction_set = true;
          if (!std::isfinite(shift_fraction) || shift_fraction < 0.0) {
            cerr << "[spectral_modes] --shift-fraction must be finite and >= 0; "
                    "using 1\n";
            shift_fraction = 1.0;
          }
        } catch (const std::exception &) {
          cerr << "[spectral_modes] ignoring bad --shift-fraction value\n";
        }
      } else if (a == "--modes" && i + 1 < argc) {
        try {
          int m = std::stoi(argv[++i]);
          num_modes = std::max(2, std::min(m, 4096));
        } catch (const std::exception &) {
          cerr << "[spectral_modes] ignoring bad --modes value\n";
        }
      }
    }

    laplace_modes_band band = laplace_modes_band::largest_magnitude;
    double mid_slider_t = -2.0;

    if (opt_low)
      band = laplace_modes_band::low_frequency;
    else if (opt_high)
      band = laplace_modes_band::largest_magnitude;
    else if (shift_fraction_set) {
      const double t = std::clamp(shift_fraction, 0.0, 1.0);
      if (t <= 0.0)
        band = laplace_modes_band::low_frequency;
      else if (t >= 1.0)
        band = laplace_modes_band::largest_magnitude;
      else {
        band = laplace_modes_band::shift_invert_middle;
        mid_slider_t = t;
      }
    } else if (opt_mid) {
      band = laplace_modes_band::shift_invert_middle;
      mid_slider_t = -1.0;
    }
    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";
    cout << "[spectral_modes] default: largest-|lambda| on (-L_sym+eps*I) (high-freq).\n";
    cout << "[spectral_modes] --low-freq / --smooth = smallest algebraic (smooth; "
            "mode 1 skips constant).\n";
    cout << "[spectral_modes] --shift-fraction F in [0,1]: F=0 smallest algebraic "
            "(low), F=1 largest |lambda| (high), 0<F<1 shift-invert with sigma "
            "interpolated low->high. --mid alone = shift-invert at upper sigma target. "
            "Optional --shift SIGMA (explicit; mid band only; ignores slider).\n";
    cout << "[spectral_modes] optional --modes N (default 32).\n";
    cout << "[spectral_modes] , or [ = prev mode; . or ] = next mode; "
            "W = toggle animation (like other viewers); S = single step\n";

    nanogui::init();

    AppPtr app = App::create(1280, 740, band, num_modes, mid_shift, mid_slider_t);
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();
  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());
    std::cerr << error_msg << endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), nullptr, MB_ICONERROR | MB_OK);
#endif
    return -1;
  }

  return 0;
}
