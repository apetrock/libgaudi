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

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/duchamp/spectral_modes_demo.hpp"
#include "gaudi/duchamp/spectral_projection_integrator.hpp"

#include <Eigen/Core>

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

namespace {
// ---- integrator tunables (quick prototyping knobs) --------------------------
constexpr double kNormalStepDt = 0.01;
constexpr double kNormalStepAmplitude = 4.0;

// Keep |f|_2 = 1 every frame so visible step size stays stable (basis truncation
// otherwise contracts the displacement field geometrically).
constexpr bool kRenormalizeProjectedCoeffs = true;

// Lock the tracked mode column-wise: pass σ = λ̃_tracked (first-order
// Rayleigh–Schrödinger) into the next shift-invert solve.
constexpr bool kRetargetSigmaToTrackedMode = true;

// CFL backtracking: halve dt while the predicted perturbation would cross modes
// (r_tracked ≥ kCflAlpha), up to kMaxDtBacktracks times.
constexpr bool kAdaptiveDtCfl = true;
constexpr double kCflAlpha = 0.5;
constexpr int kMaxDtBacktracks = 5;

// Hungarian column alignment is slower but optimal for large N; greedy is fine
// for small N and when modes don't cluster tightly.
constexpr bool kUseHungarianAlignment = false;

// Per-frame remeshing via `asawa::shell::dynamic` (subdivide / collapse / flip).
// When topology changes (vertex count differs), we rebuild the basis from
// scratch — the previous `_prev_evecs` snapshot has mismatched row count so
// warm-start / prediction is invalid for that frame.
constexpr bool kAnimateRemesh = true;

// Edge-length thresholds in multiples of the mesh's initial average edge length `l0`:
//   * `kRemeshCollapseScale * l0` — collapse edges shorter than this
//   * `kRemeshStretchScale  * l0` — subdivide edges longer than this
//   * `kRemeshBridgeScale   * l0` — proximity for self-contact bridge/merge
// Wider collapse/stretch band ⇒ less remeshing churn; tighter band ⇒ more aggressive.
constexpr double kRemeshCollapseScale = 2.0;
constexpr double kRemeshStretchScale = 1.75;
constexpr double kRemeshBridgeScale = 0.5;
// -----------------------------------------------------------------------------

std::vector<gg::colorRGB> integrator_colors_magenta_cyan(const Eigen::VectorXd &phi) {
  const Eigen::Matrix<double, 4, 1> col_a(1.0, 0.0, 1.0, 1.0);
  const Eigen::Matrix<double, 4, 1> col_b(0.0, 1.0, 1.0, 1.0);
  const auto mapped =
      gaudi::duchamp::scalar_field_to_mesh_colors(phi, col_a, col_b);
  std::vector<gg::colorRGB> colors;
  colors.reserve(mapped.size());
  for (const auto &c : mapped)
    colors.emplace_back(c[0], c[1], c[2], 1.0);
  return colors;
}
} // namespace

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

    gaudi::duchamp::spectral_projection_config pcfg;
    pcfg.requested_modes = num_eigenmodes;
    pcfg.epsilon_shift = 1e-6;
    pcfg.band = band;
    pcfg.mid_band_shift = mid_band_shift;
    pcfg.mid_slider_t = mid_slider_t;
    pcfg.renormalize_projected_coeffs = kRenormalizeProjectedCoeffs;
    pcfg.retarget_sigma_to_tracked_mode = kRetargetSigmaToTrackedMode;
    pcfg.adaptive_dt_cfl = kAdaptiveDtCfl;
    pcfg.cfl_alpha = kCflAlpha;
    pcfg.max_dt_backtracks = kMaxDtBacktracks;
    pcfg.use_hungarian_alignment = kUseHungarianAlignment;
    _integrator = gaudi::duchamp::spectral_projection_integrator(pcfg);
    _integrator.set_mesh(__demo->__M);

    {
      using gaudi::asawa::get_vec_data;
      using gaudi::asawa::shell::dynamic;
      std::vector<gaudi::vec3> &x = get_vec_data(*__demo->__M, 0);
      const gaudi::real l0 =
          gaudi::asawa::shell::avg_length(*__demo->__M, x);
      // Thresholds as multiples of l0 — tunables hoisted to top of file.
      __dynamic_surf = dynamic::create(
          __demo->__M, 
          static_cast<gaudi::real>(kRemeshCollapseScale) * l0,
          static_cast<gaudi::real>(kRemeshStretchScale) * l0,
          static_cast<gaudi::real>(kRemeshBridgeScale) * l0);
    }

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
    if (_integrator_ready)
      _integrator.set_coeffs_one_hot(__demo->mode_index());
    mark_mesh_dirty();
  }

  void next_mode() {
    if (!__demo || __demo->mode_count() <= 1)
      return;
    int i = __demo->mode_index() + 1;
    if (i >= __demo->mode_count())
      i = 0;
    __demo->set_mode_index(i);
    if (_integrator_ready)
      _integrator.set_coeffs_one_hot(__demo->mode_index());
    mark_mesh_dirty();
  }

  void basis_refresh() {
    if (!__demo)
      return;
    _integrator.set_mesh(__demo->__M);
    if (!_integrator.rebuild_operator()) {
      cerr << "[spectral_modes] B: rebuild_operator failed\n";
      return;
    }
    if (!_integrator.compute_modes(nullptr, true)) {
      cerr << "[spectral_modes] B: compute_modes failed\n";
      return;
    }
    _integrator.set_coeffs_one_hot(__demo->mode_index());
    _integrator_ready = true;
    _viz_integrator = true;
    _needs_spectral_refresh = false;
    mark_mesh_dirty();
  }

  void integrator_projection_step() {
    if (!__demo)
      return;
    if (!_integrator_ready) {
      cerr << "[spectral_modes] N: press B first to build the projection basis.\n";
      return;
    }
    _integrator.set_mesh(__demo->__M);
    if (!_integrator.advance_projection_step(kNormalStepDt, kNormalStepAmplitude, true)) {
      cerr << "[spectral_modes] N: advance_projection_step failed\n";
      return;
    }
    _viz_integrator = true;
    _needs_spectral_refresh = false;
    mark_mesh_dirty();
  }

  /// Per-vertex arrows `x_i -> x_i + n_i * f_i * amplitude`, where `f = evecs * coeffs`
  /// is the signed vertex field the integrator is about to apply (via `abs(f)`).
  /// Positive f → cyan (outward push), negative → magenta (inward pull). Matches the
  /// mesh coloring so you can read which direction the mode wants to deform.
  /// Lines are re-submitted each time so they stay visible across `render()`/`clear()` frames.
  void log_mode_direction_field(double amplitude) {
    if (!__demo || !_integrator_ready || !_integrator.spectrum_ok())
      return;
    auto &M = *__demo->__M;
    const std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const int nv = static_cast<int>(x.size());
    const Eigen::VectorXd f = _integrator.vertex_field();
    if (static_cast<int>(f.size()) != nv)
      return;
    const Vec4d col_pos(0.0, 1.0, 1.0, 1.0);
    const Vec4d col_neg(1.0, 0.0, 1.0, 1.0);
    for (int vi = 0; vi < nv; ++vi) {
      const gaudi::vec3 n =
          gaudi::asawa::shell::vert_normal(M, gaudi::asawa::shell::vert_id(vi), x)
              .normalized();
      const double fi = f[vi];
      const Vec3d p0 = x[static_cast<size_t>(vi)];
      const Vec3d p1 = p0 + n * (fi * amplitude);
      gg::geometry_logger::line(p0, p1, fi >= 0.0 ? col_pos : col_neg);
    }
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
    const int nv = __demo->__M->vert_count();
    std::vector<gg::colorRGB> colors(static_cast<size_t>(nv),
                                     gg::colorRGB(0.65, 0.68, 0.72, 1.0));
    if (_viz_integrator && _integrator_ready && _integrator.spectrum_ok() &&
        !_needs_spectral_refresh) {
      const Eigen::VectorXd phi = _integrator.vertex_field();
      colors = integrator_colors_magenta_cyan(phi);
    } else if (!_needs_spectral_refresh) {
      const auto mesh_colors = __demo->get_mesh_colors();
      const int n = std::min(static_cast<int>(mesh_colors.size()), nv);
      for (int i = 0; i < n; ++i) {
        const auto &c = mesh_colors[static_cast<size_t>(i)];
        colors[static_cast<size_t>(i)] = gg::colorRGB(c[0], c[1], c[2], 1.0);
      }
    }
    gg::fillBuffer_ref(*__demo->__M, _objs[0], colors);
  }

  /// Re-submit any debug-line overlays (mode direction-field arrows) each frame.
  /// `geometry_logger::clear()` drains the CPU-side buffer at the end of every draw,
  /// so these must be pushed fresh before `render()`.
  void push_debug_overlays() {
    if (_viz_integrator && _integrator_ready && _integrator.spectrum_ok() &&
        !_needs_spectral_refresh) {
      log_mode_direction_field(kNormalStepAmplitude);
    }
  }

  /// If `kAnimateRemesh`, runs one subdivide/collapse/flip pass and ensures dense
  /// packing. Returns true if the vertex count changed (basis must be rebuilt
  /// from scratch; warm-start / prediction are invalid for the post-remesh frame).
  bool remesh_if_enabled() {
    if (!kAnimateRemesh || !__dynamic_surf || !__demo)
      return false;
    const int nv_before = __demo->__M->vert_count();
    __dynamic_surf->step();
    const int nv_after = __demo->__M->vert_count();
    return nv_after != nv_before;
  }

  /// W toggles animation. Every tick:
  ///   (optional) remesh → if topology changed, rebuild basis from scratch →
  ///   `advance_projection_step` (displace, rebuild op, predict, backtrack dt,
  ///    retarget σ, solve, align). Rebuilding the basis every frame is the point —
  ///   the mode tracks the geometry as it moves. Press **B** to (re)build once;
  ///   **N** steps once without animation.
  virtual void onAnimate(int /*frame*/) {
    if (!_integrator_ready || !_integrator.spectrum_ok())
      return;

    const bool topo_changed = remesh_if_enabled();
    _integrator.set_mesh(__demo->__M);

    if (topo_changed) {
      // Previous `_prev_evecs` has the old row count — can't warm-start or
      // predict this frame. Solve the new operator from scratch, re-seat coeffs
      // on the tracked mode index, and skip the projection step this tick.
      if (!_integrator.rebuild_operator() ||
          !_integrator.compute_modes(nullptr, true)) {
        cerr << "[spectral_modes] W: post-remesh basis rebuild failed; "
                "halting animation (press B to retry).\n";
        _integrator_ready = false;
        _viz_integrator = false;
        _needs_spectral_refresh = true;
        return;
      }
      _integrator.set_coeffs_one_hot(__demo->mode_index());
      _viz_integrator = true;
      _needs_spectral_refresh = false;
      mark_mesh_dirty();
      return;
    }

    if (!_integrator.advance_projection_step(kNormalStepDt,
                                             kNormalStepAmplitude, true)) {
      cerr << "[spectral_modes] W: advance_projection_step failed this "
              "tick; retaining prior basis and retrying next frame.\n";
    } else {
      const auto &p = _integrator.last_prediction();
      if (p.valid) {
        const int k = p.tracked_index;
        const auto &prev_evals = _integrator.previous_evals();
        const auto &evals = _integrator.evals();
        const double lam_old =
            (k >= 0 && k < prev_evals.size()) ? prev_evals[k] : 0.0;
        const double lam_new =
            (k >= 0 && k < evals.size()) ? evals[k] : 0.0;
        cerr << "[spectral_modes] track mode " << k
             << ": lam " << lam_old << " -> lam_tilde=" << p.tracked_lambda_tilde
             << " (actual " << lam_new << "), gap=" << p.tracked_gap
             << ", r=" << p.tracked_robustness
             << (p.tracked_robustness >= 0.5 ? "  [!crossing risk]" : "")
             << "\n";
      }
    }
    _viz_integrator = true;
    _needs_spectral_refresh = false;
    mark_mesh_dirty();
  }

  virtual void onDraw(gg::Viewer &viewer) {
    if (_mesh_dirty) {
      refresh_mesh_buffer();
      _mesh_dirty = false;
    }
    push_debug_overlays();
    gg::geometry_logger::render();
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

  gaudi::duchamp::spectral_modes_demo::ptr __demo;
  gaudi::asawa::shell::dynamic::ptr __dynamic_surf;
  gaudi::duchamp::spectral_projection_integrator _integrator;
  bool _integrator_ready = false;
  bool _viz_integrator = false;
  /// Set after each `dynamic::step`: demo eigen-colors and cached integrator modes are invalid.
  bool _needs_spectral_refresh = false;

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
      if (key == GLFW_KEY_B) {
        scene->basis_refresh();
        drawContents();
        return true;
      }
      if (key == GLFW_KEY_N) {
        scene->integrator_projection_step();
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
    int num_modes = 8;
    double lambda_t = 0.0;
    bool lambda_set = false;
    for (int i = 1; i < argc; ++i) {
      const std::string a(argv[i]);
      if (a == "--lambda" && i + 1 < argc) {
        try {
          const double v = std::stod(argv[++i]);
          if (!std::isfinite(v)) {
            cerr << "[spectral_modes] --lambda must be finite in [0,1]; using 0\n";
          } else {
            lambda_t = std::clamp(v, 0.0, 1.0);
            lambda_set = true;
          }
        } catch (const std::exception &) {
          cerr << "[spectral_modes] ignoring bad --lambda value\n";
        }
      } else if (a == "--N" && i + 1 < argc) {
        try {
          const int m = std::stoi(argv[++i]);
          num_modes = std::max(2, std::min(m, 4096));
        } catch (const std::exception &) {
          cerr << "[spectral_modes] ignoring bad --N value\n";
        }
      }
    }

    // --lambda F in [0,1]: F=0 snaps to smallest-algebraic solver (lowest K modes);
    // F=1 snaps to largest-magnitude (highest K modes); interior F routes to
    // shift-invert with sigma = F * |A|_inf (see spectral_modes_demo auto-sigma).
    laplace_modes_band band = laplace_modes_band::low_frequency;
    double mid_slider_t = -2.0;
    const double mid_shift = -1.0;
    if (lambda_t <= 0.0) {
      band = laplace_modes_band::low_frequency;
    } else if (lambda_t >= 1.0) {
      band = laplace_modes_band::largest_magnitude;
    } else {
      band = laplace_modes_band::shift_invert_middle;
      mid_slider_t = lambda_t;
    }

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";
    cout << "[spectral_modes] --lambda F in [0,1] (default 0): F=0 lowest K modes "
            "(smallest algebraic), F=1 highest K (largest |lambda|), interior F = "
            "shift-invert at sigma = F * |A|_inf.\n";
    cout << "[spectral_modes] --N K (default 8): number of eigenmodes to compute.\n";
    if (lambda_set)
      cout << "[spectral_modes] lambda=" << lambda_t << " N=" << num_modes << "\n";
    else
      cout << "[spectral_modes] (no --lambda given; using default lambda=0, N="
           << num_modes << ")\n";
    cout << "[spectral_modes] , or [ = prev mode; . or ] = next mode; "
            "W = toggle animation; S = single step\n";
    cout << "[spectral_modes] B = projection basis refresh (rebuild + eigenmodes + "
            "colors from integrator); N key = one projected normal step (dt="
         << kNormalStepDt << ", amplitude=" << kNormalStepAmplitude << ")\n";
    cout << "[spectral_modes] Each animation frame runs `shell::dynamic::step()` "
            "(remesh + pack); press B after the mesh moves to refresh modes.\n";

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
