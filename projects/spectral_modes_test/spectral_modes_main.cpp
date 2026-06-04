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

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/laplace_spectrum.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/laplacian_anisotropic.hpp"
#include "gaudi/kusama/vector_dirichlet_guided.hpp"
#include "gaudi/duchamp/spectral_modes_demo.hpp" // scalar_field_to_mesh_colors, laplace_modes_band

#include <Eigen/Core>
#include <nanogui/nanogui.h>

using std::cerr;
using std::cout;
using std::endl;

namespace {

// -- viewer config -----------------------------------------------------------
constexpr bool kAnisotropyEnabled = false;
constexpr gaudi::kusama::anisotropic_laplacian_kind kAnisotropyKindDefault =
    gaudi::kusama::anisotropic_laplacian_kind::conductance;
constexpr double kAnisotropySigmaU = 1.0;
constexpr double kAnisotropySigmaV = 0.01;
// Match `frame_curvature_demo`'s VD lambda — coherent direction field.
constexpr double kAnisotropyVDLambda = 3.0;
constexpr gaudi::asawa::shell::face_curvature_stencil kAnisotropyStencil =
    gaudi::asawa::shell::face_curvature_stencil::one_ring;

// Always show the curvature-guidance overlay when aniso is on; useful even
// when iso to inspect the field — flip to keep behavior parallel.
constexpr bool kDrawGuidanceOverlay = true;

// Debug line scale: matches `frame_curvature_demo` (avg_length * 3, one-sided).
constexpr double kDebugLineScale = 3.0;
// ---------------------------------------------------------------------------

} // namespace

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {
public:
  static ScenePtr create(gaudi::duchamp::laplace_modes_band band, int num_modes,
                         double mid_band_shift, double mid_slider_t) {
    return std::make_shared<Scene>(band, num_modes, mid_band_shift,
                                   mid_slider_t);
  }

  Scene(gaudi::duchamp::laplace_modes_band band, int num_modes,
        double mid_band_shift, double mid_slider_t)
      : gg::Scene(), _band(band), _requested_modes(num_modes),
        _mid_band_shift(mid_band_shift), _mid_slider_t(mid_slider_t),
        _aniso_kind(kAnisotropyKindDefault) {
    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);

    load_mesh();
    rebuild_and_solve();
  }

  // ---- one-time mesh setup ------------------------------------------------
  void load_mesh() {
    __M = gaudi::asawa::shell::load_bunny();
    gaudi::asawa::shell::triangulate(*__M);
    if (!__M->verts_are_dense_packed())
      gaudi::asawa::shell::pack(*__M);
    if (!__M->verts_are_dense_packed())
      throw std::runtime_error("spectral_modes: pack failed to dense-pack mesh");
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(*__M, 0);
    gaudi::asawa::center(x, 2.0);
  }

  // ---- build operator + solve eigenmodes ----------------------------------
  void rebuild_and_solve() {
    if (!__M)
      return;
    auto &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const int nv = static_cast<int>(x.size());

    Eigen::SparseMatrix<double> L;
    if (kAnisotropyEnabled) {
      if (!compute_guidance_field(M, x))
        return;
      L = gaudi::kusama::build_anisotropic_cotan_laplacian(
          M, x, _g, _slot_map, _nE,
          static_cast<gaudi::real>(kAnisotropySigmaU),
          static_cast<gaudi::real>(kAnisotropySigmaV), _aniso_kind);
      cerr << "[spectral_modes] aniso L: nnz=" << L.nonZeros()
           << " kind=" << aniso_kind_string() << "\n";
    } else {
      gaudi::kusama::laplacian lap(__M, x);
      L = lap.stiffness();
    }

    Eigen::SparseMatrix<double> Ls = gaudi::kusama::symmetrize_sparse(L);
    Eigen::SparseMatrix<double> Lpos = Ls;
    Lpos *= -1.0;
    Lpos.makeCompressed();
    _A = gaudi::kusama::regularize_stiffness(Lpos, 1e-6);
    if (static_cast<int>(_A.rows()) != nv)
      cerr << "[spectral_modes] WARNING: operator size " << _A.rows()
           << " != nv " << nv << "\n";

    bool ok = false;
    if (_band == gaudi::duchamp::laplace_modes_band::low_frequency) {
      ok = gaudi::kusama::laplace_eigs_low_frequency(_A, _requested_modes, 0,
                                                       _evals, _evecs);
    } else if (_band == gaudi::duchamp::laplace_modes_band::largest_magnitude) {
      ok = gaudi::kusama::laplace_eigs_largest_magnitude(
          _A, _requested_modes, 0, _evals, _evecs);
    } else {
      double sigma = _mid_band_shift;
      if (sigma < 0.0) {
        const double bound = gaudi::kusama::sparse_sym_max_row_sum_abs(_A);
        const double t = (_mid_slider_t >= 0.0 && _mid_slider_t <= 1.0)
                             ? _mid_slider_t
                             : 0.35;
        sigma = t * std::max(bound, 1e-12);
      }
      ok = gaudi::kusama::laplace_eigs_shift_invert_nearest(
          _A, sigma, _requested_modes, 0, _evals, _evecs);
    }
    _spectrum_ok = ok;
    if (!ok) {
      cerr << "[spectral_modes] eigen solve failed\n";
      return;
    }
    _num_modes = static_cast<int>(_evals.size());
    _mode_idx = (_num_modes > 1 &&
                 _band == gaudi::duchamp::laplace_modes_band::low_frequency)
                    ? 1
                    : 0;
    cerr << "[spectral_modes] computed " << _num_modes
         << " eigenmodes (nnz=" << _A.nonZeros() << ")\n";
    const int ncheck = std::min(3, _num_modes);
    for (int j = 0; j < ncheck; ++j)
      gaudi::kusama::log_laplace_eigen_stats(_A, j, _evals[j], _evecs.col(j));

    if (kAnisotropyEnabled)
      log_aniso_vs_iso_frobenius_ratio();

    rebuild_overlay_lines();
    _mesh_dirty = true;
  }

  // ---- guidance solve (anisotropy and/or overlay) -------------------------
  bool compute_guidance_field(gaudi::asawa::shell::shell &M,
                              std::vector<gaudi::vec3> &x) {
    _g.resize(0);
    _slot_map.clear();
    _edge_inc.clear();
    _nE = 0;
    _nE = gaudi::kusama::build_compact_edge_dof_map(M, _slot_map);
    if (_nE <= 0)
      return false;
    gaudi::kusama::build_edge_incident_faces(M, _slot_map, _nE, _edge_inc);
    try {
      _g = gaudi::kusama::solve_curvature_guided_vector_dirichlet(
          M, x, static_cast<gaudi::real>(kAnisotropyVDLambda),
          kAnisotropyStencil, /*apply_sign_coherence=*/true);
    } catch (const std::exception &e) {
      cerr << "[spectral_modes] VD solve failed: " << e.what() << "\n";
      _g.resize(0);
      _nE = 0;
      return false;
    }
    if (_g.size() != 2 * _nE) {
      _g.resize(0);
      _nE = 0;
      return false;
    }
    cerr << "[spectral_modes] guidance: nE=" << _nE << " |g|=" << _g.norm()
         << "\n";
    return true;
  }

  // ---- precompute guidance overlay lines once -----------------------------
  void rebuild_overlay_lines() {
    _overlay_p0.clear();
    _overlay_p1.clear();
    if (!kDrawGuidanceOverlay || !__M)
      return;
    auto &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);

    // For the iso path we still want to *show* the curvature field, so solve
    // it once here if it isn't already cached.
    if (_g.size() != 2 * _nE || _nE <= 0) {
      if (!compute_guidance_field(M, x))
        return;
    }
    if (_g.size() != 2 * _nE || _nE <= 0)
      return;

    const gaudi::real scale = static_cast<gaudi::real>(kDebugLineScale) *
                              gaudi::asawa::shell::avg_length(M, x);
    _overlay_p0.reserve(static_cast<size_t>(_nE));
    _overlay_p1.reserve(static_cast<size_t>(_nE));
    for (gaudi::asawa::shell::CornerId c : M.get_edge_range()) {
      const int slot = static_cast<int>(c) / 2;
      const int e = _slot_map[static_cast<size_t>(slot)];
      if (e < 0 || e >= _nE)
        continue;
      const gaudi::vec3 mid =
          0.5 * (x[M.vert(c)] + x[M.vert(M.next(c))]);
      const gaudi::vec3 n_e = gaudi::kusama::edge_average_normal(
          M, x, _edge_inc[static_cast<size_t>(e)]);
      gaudi::vec3 dir = gaudi::kusama::edge_guidance_vector_3d(
          M, c, x, n_e, static_cast<gaudi::real>(_g(e)),
          static_cast<gaudi::real>(_g(e + _nE)));
      if (dir.norm() > 1e-12)
        dir.normalize();
      _overlay_p0.push_back(mid);
      _overlay_p1.push_back(mid + scale * dir);
    }
  }

  void push_overlay_lines() {
    if (!kDrawGuidanceOverlay || _overlay_p0.empty())
      return;
    const gaudi::vec4 col(0.2, 0.85, 0.55, 1.0);
    for (size_t i = 0; i < _overlay_p0.size(); ++i)
      gg::geometry_logger::line(_overlay_p0[i], _overlay_p1[i], col);
  }

  // ---- diagnostics --------------------------------------------------------
  void log_aniso_vs_iso_frobenius_ratio() {
    if (!__M)
      return;
    auto &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    gaudi::kusama::laplacian lap(__M, x);
    Eigen::SparseMatrix<double> L_iso =
        gaudi::kusama::symmetrize_sparse(lap.stiffness());
    Eigen::SparseMatrix<double> L_ani =
        gaudi::kusama::build_anisotropic_cotan_laplacian(
            M, x, _g, _slot_map, _nE,
            static_cast<gaudi::real>(kAnisotropySigmaU),
            static_cast<gaudi::real>(kAnisotropySigmaV), _aniso_kind);
    L_ani = gaudi::kusama::symmetrize_sparse(L_ani);
    if (L_iso.rows() != L_ani.rows() || L_iso.cols() != L_ani.cols())
      return;
    Eigen::SparseMatrix<double> D = L_ani - L_iso;
    const double r = D.norm() / std::max(1e-30, L_iso.norm());
    cerr << "[spectral_modes] aniso sanity: ||L_aniso - L_iso||_F / "
            "||L_iso||_F = "
         << r << " (≈0 ⇒ guidance is effectively isotropic)\n";
  }

  // ---- mode navigation ----------------------------------------------------
  void prev_mode() {
    if (!_spectrum_ok || _num_modes <= 1)
      return;
    _mode_idx = (_mode_idx + _num_modes - 1) % _num_modes;
    _mesh_dirty = true;
  }
  void next_mode() {
    if (!_spectrum_ok || _num_modes <= 1)
      return;
    _mode_idx = (_mode_idx + 1) % _num_modes;
    _mesh_dirty = true;
  }
  int current_mode() const { return _mode_idx; }
  int mode_count() const { return _spectrum_ok ? _num_modes : 1; }

  void log_active_mode_stats() const {
    if (_spectrum_ok && _num_modes > 0)
      gaudi::kusama::log_laplace_eigen_stats(_A, _mode_idx, _evals[_mode_idx],
                                               _evecs.col(_mode_idx));
  }

  void toggle_anisotropic_kind() {
    if (!kAnisotropyEnabled)
      return;
    using gaudi::kusama::anisotropic_laplacian_kind;
    _aniso_kind = (_aniso_kind == anisotropic_laplacian_kind::conductance)
                      ? anisotropic_laplacian_kind::fem_d
                      : anisotropic_laplacian_kind::conductance;
    cerr << "[spectral_modes] anisotropy kind = " << aniso_kind_string()
         << "\n";
    rebuild_and_solve();
  }
  const char *aniso_kind_string() const {
    return _aniso_kind ==
                   gaudi::kusama::anisotropic_laplacian_kind::conductance
               ? "conductance"
               : "fem_d";
  }

  // ---- coloring -----------------------------------------------------------
  std::vector<gg::colorRGB> colors_for_active_mode() const {
    const int nv = __M ? __M->vert_count() : 0;
    std::vector<gg::colorRGB> colors(static_cast<size_t>(nv),
                                     gg::colorRGB(0.65, 0.68, 0.72, 1.0));
    if (!_spectrum_ok || _evecs.size() == 0 || nv <= 0)
      return colors;
    const int n = std::min(nv, static_cast<int>(_evecs.rows()));
    if (n <= 0)
      return colors;
    Eigen::VectorXd phi = _evecs.col(_mode_idx).head(n);
    const gaudi::vec4 col_a(1.0, 0.0, 1.0, 1.0);
    const gaudi::vec4 col_b(0.0, 1.0, 1.0, 1.0);
    auto mapped =
        gaudi::duchamp::scalar_field_to_mesh_colors(phi, col_a, col_b);
    for (int i = 0; i < static_cast<int>(mapped.size()) && i < nv; ++i) {
      const auto &c = mapped[static_cast<size_t>(i)];
      colors[static_cast<size_t>(i)] = gg::colorRGB(c[0], c[1], c[2], 1.0);
    }
    return colors;
  }

  void refresh_mesh_buffer() {
    gg::fillBuffer_ref(*__M, _objs[0], colors_for_active_mode());
  }

  // ---- gg::Scene overrides ------------------------------------------------
  virtual void onAnimate(int /*frame*/) override {
    // Static viewer: nothing to step. Overlay lines are already cached.
  }

  virtual void onDraw(gg::Viewer &viewer) override {
    if (_mesh_dirty) {
      refresh_mesh_buffer();
      _mesh_dirty = false;
    }
    push_overlay_lines();
    gg::geometry_logger::render();
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    gg::geometry_logger::clear();
  }

  gaudi::asawa::shell::shell::ptr __M;

private:
  // mode-band config
  gaudi::duchamp::laplace_modes_band _band;
  int _requested_modes = 8;
  double _mid_band_shift = -1.0;
  double _mid_slider_t = -2.0;

  // operator + spectrum
  Eigen::SparseMatrix<double> _A;
  Eigen::VectorXd _evals;
  Eigen::MatrixXd _evecs;
  bool _spectrum_ok = false;
  int _num_modes = 0;
  int _mode_idx = 0;

  // anisotropy
  gaudi::kusama::anisotropic_laplacian_kind _aniso_kind;

  // guidance field cache (used by aniso operator and the overlay)
  Eigen::VectorXd _g;
  std::vector<int> _slot_map;
  std::vector<std::vector<gaudi::asawa::shell::FaceId>> _edge_inc;
  int _nE = 0;

  // precomputed overlay segments (avoids re-solve per frame)
  std::vector<gaudi::vec3> _overlay_p0;
  std::vector<gaudi::vec3> _overlay_p1;

  // gg state
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
  bool _mesh_dirty = true;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height,
                       gaudi::duchamp::laplace_modes_band band, int num_modes,
                       double mid_band_shift, double mid_slider_t) {
    return std::make_shared<App>(width, height, band, num_modes, mid_band_shift,
                                 mid_slider_t);
  }

  App(int width, int height, gaudi::duchamp::laplace_modes_band band,
      int num_modes, double mid_band_shift, double mid_slider_t)
      : gg::SimpleApp(width, height, 4.0, false, "spectral_modes_") {
    setScene(scene = Scene::create(band, num_modes, mid_band_shift,
                                   mid_slider_t));
    set_rotate_ball(false);
    initUI();
  }

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
             << (scene->mode_count() - 1) << "\n";
        drawContents();
        return true;
      }
      if (key == GLFW_KEY_PERIOD || key == GLFW_KEY_RIGHT_BRACKET) {
        scene->next_mode();
        scene->log_active_mode_stats();
        cerr << "[spectral_modes] mode " << scene->current_mode() << " / "
             << (scene->mode_count() - 1) << "\n";
        drawContents();
        return true;
      }
      if (key == GLFW_KEY_B) {
        scene->rebuild_and_solve();
        drawContents();
        return true;
      }
      if (key == GLFW_KEY_K && kAnisotropyEnabled) {
        scene->toggle_anisotropic_kind();
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
          if (std::isfinite(v)) {
            lambda_t = std::clamp(v, 0.0, 1.0);
            lambda_set = true;
          }
        } catch (const std::exception &) {
        }
      } else if (a == "--N" && i + 1 < argc) {
        try {
          num_modes = std::clamp(std::stoi(argv[++i]), 2, 4096);
        } catch (const std::exception &) {
        }
      }
    }

    laplace_modes_band band;
    double mid_slider_t = -2.0;
    const double mid_shift = -1.0;
    if (lambda_t <= 0.0)
      band = laplace_modes_band::low_frequency;
    else if (lambda_t >= 1.0)
      band = laplace_modes_band::largest_magnitude;
    else {
      band = laplace_modes_band::shift_invert_middle;
      mid_slider_t = lambda_t;
    }

    cout << "[spectral_modes] static viewer: load mesh -> build operator -> "
            "solve K eigenmodes -> color by selected mode.\n";
    cout << "[spectral_modes] --lambda F in [0,1] (default 0): F=0 lowest K "
            "modes, F=1 highest K, interior F = shift-invert at "
            "sigma=F*|A|_inf.\n";
    cout << "[spectral_modes] --N K (default 8): number of eigenmodes.\n";
    if (lambda_set)
      cout << "[spectral_modes] lambda=" << lambda_t << " N=" << num_modes
           << "\n";
    cout << "[spectral_modes] keys: , [ = prev mode | . ] = next mode | "
            "B = rebuild + re-solve\n";
    if (kAnisotropyEnabled)
      cout << "[spectral_modes] aniso on: K = toggle conductance <-> fem_d\n";

    nanogui::init();
    AppPtr app =
        App::create(1280, 740, band, num_modes, mid_shift, mid_slider_t);
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();
  } catch (const std::runtime_error &e) {
    const std::string msg =
        std::string("Caught a fatal error: ") + e.what();
    cerr << msg << endl;
#if defined(WIN32)
    MessageBoxA(nullptr, msg.c_str(), nullptr, MB_ICONERROR | MB_OK);
#endif
    return -1;
  }
  return 0;
}
