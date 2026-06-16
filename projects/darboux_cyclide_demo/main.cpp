#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
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

#include "gaudi/albers/darboux_cyclide.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/test/darboux_cyclide_torus_fixture.hpp"
#include "gaudi/calder/least_squares_fit.hpp"
#include "gaudi/geometry_logger.hpp"

namespace {

constexpr int kSingleFitVertex = 0;
constexpr bool kAutoIsolateWorstFit = false;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

struct FitIsolation {
  int vertex = kSingleFitVertex;
  gaudi::real signed_alignment = 1.0;
  gaudi::real abs_alignment = 1.0;
  gaudi::real gradient_norm = 0.0;
  gaudi::real value_at_origin = 0.0;
};

void sync_gaudi_debug_to_gg() {
  const auto &lines = gaudi::geometry_logger::get_lines();
  const auto &line_cols = gaudi::geometry_logger::get_line_colors();
  for (size_t i = 0; i + 1 < lines.size(); i += 2) {
    gg::geometry_logger::line(lines[i], lines[i + 1], line_cols[i]);
  }

  const auto &points = gaudi::geometry_logger::get_points();
  const auto &point_cols = gaudi::geometry_logger::get_point_colors();
  for (size_t i = 0; i < points.size(); ++i) {
    gg::geometry_logger::point(points[i], point_cols[i]);
  }
}

class Scene : public gg::Scene {
public:
  static ScenePtr create(int stride, gaudi::real axis_tilt) {
    return std::make_shared<Scene>(stride, axis_tilt);
  }

  Scene(int stride, gaudi::real axis_tilt)
      : gg::Scene(), __stride(std::max(1, stride)) {
        const gaudi::vec3 axis(0.25 + axis_tilt, 0.65, 0.9);
    __frame =
        gaudi::test::make_torus_frame(gaudi::vec3(0.7, -0.45, 0.3), axis);
        // Canonical torus at the origin, major circle in XY, tube axis +Z.
    // Use --axis-tilt to perturb the tube axis for harder cases.
    //gaudi::vec3 axis = gaudi::vec3::UnitZ();
    //if (std::abs(axis_tilt) > 1e-12) {
    //  axis = gaudi::vec3(axis_tilt, 0.0, 1.0).normalized();
    //}
    //__frame = gaudi::test::make_torus_frame(gaudi::vec3::Zero(), axis);
    gaudi::test::TorusMesh torus =
        gaudi::test::make_offset_torus_shell(36, 20, __major_radius,
                                             __minor_radius, __frame);
    __M = torus.shell;

    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    __vertex_normals = gaudi::asawa::shell::vertex_normals(M, x);
    __l0 = 0.01 * gaudi::asawa::shell::avg_length(M, x);
    // Per-vertex fit: p_pov = vertex positions, n_pov = vertex normals.
    // generic_fit binds face-area * face-normal samples on the tree internally.
    const auto fits =
        gaudi::calder::darboux_cyclide_normal_constrained(
            M, x, __vertex_normals, __l0, 3.0);
    isolate_fit(M, fits);
    evaluate_single_fit_on_vertices(M);

    _objs.resize(1);
    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    refresh_mesh_buffer();
    mSceneObjects.push_back(_objs[0]);
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  void isolate_fit(gaudi::asawa::shell::shell &M,
                   const std::vector<gaudi::albers::vec14> &fits) {
    const std::vector<gaudi::vec3> &x = gaudi::asawa::const_get_vec_data(M, 0);
    if (fits.empty() || x.empty()) {
      return;
    }

    __single_fit_vertex =
        std::clamp(kSingleFitVertex, 0, static_cast<int>(fits.size()) - 1);
    if constexpr (kAutoIsolateWorstFit) {
      FitIsolation worst;
      worst.abs_alignment = std::numeric_limits<gaudi::real>::infinity();
      for (auto vi : M.get_vert_range()) {
        const int v = static_cast<int>(vi);
        if (v < 0 || v >= static_cast<int>(fits.size()) ||
            v >= static_cast<int>(__vertex_normals.size())) {
          continue;
        }

        const gaudi::vec3 g =
            gaudi::albers::darboux_grad(fits[static_cast<size_t>(v)],
                                        gaudi::vec3::Zero());
        const gaudi::real g_norm = g.norm();
        if (!g.allFinite() || g_norm < 1e-12) {
          continue;
        }

        const gaudi::vec3 n = __vertex_normals[static_cast<size_t>(v)].normalized();
        const gaudi::real signed_alignment = g.normalized().dot(n);
        const gaudi::real abs_alignment = std::abs(signed_alignment);
        if (abs_alignment < worst.abs_alignment) {
          worst.vertex = v;
          worst.signed_alignment = signed_alignment;
          worst.abs_alignment = abs_alignment;
          worst.gradient_norm = g_norm;
          worst.value_at_origin = gaudi::albers::eval_darboux(
              fits[static_cast<size_t>(v)], gaudi::vec3::Zero());
        }
      }
      __single_fit_vertex = worst.vertex;
      std::cerr << "isolated Darboux fit vertex=" << worst.vertex
                << " signed_align=" << worst.signed_alignment
                << " abs_align=" << worst.abs_alignment
                << " grad_norm=" << worst.gradient_norm
                << " D0=" << worst.value_at_origin << " l0=" << __l0
                << std::endl;
    }

    __single_fit_origin = x[static_cast<size_t>(__single_fit_vertex)];
    __single_fit_Q = fits[static_cast<size_t>(__single_fit_vertex)];
  }

  void evaluate_single_fit_on_vertices(gaudi::asawa::shell::shell &M) {
    const std::vector<gaudi::vec3> &x = gaudi::asawa::const_get_vec_data(M, 0);
    __cyclide_normals.assign(x.size(), gaudi::vec3::Zero());

    gaudi::real max_abs_value = 0.0;
    gaudi::real rms_value = 0.0;
    gaudi::real alignment_sum = 0.0;
    int value_samples = 0;
    int alignment_samples = 0;
    for (auto vi : M.get_vert_range()) {
      const int v = static_cast<int>(vi);
      if (v < 0 || v >= static_cast<int>(x.size()) ||
          v >= static_cast<int>(__vertex_normals.size())) {
        continue;
      }

      const gaudi::vec3 local = x[static_cast<size_t>(v)] - __single_fit_origin;
      const gaudi::real value = gaudi::albers::eval_darboux(__single_fit_Q, local);
      if (std::isfinite(value)) {
        max_abs_value = std::max(max_abs_value, std::abs(value));
        rms_value += value * value;
        ++value_samples;
      }

      const gaudi::vec3 g = gaudi::albers::darboux_grad(__single_fit_Q, local);
      __cyclide_normals[static_cast<size_t>(v)] = g;
      if (g.allFinite() && g.norm() > 1e-12) {
        alignment_sum +=
            std::abs(g.normalized().dot(__vertex_normals[static_cast<size_t>(v)].normalized()));
        ++alignment_samples;
      }
    }

    if (value_samples > 0) {
      rms_value = std::sqrt(rms_value / gaudi::real(value_samples));
    }
    const gaudi::real mean_alignment =
        alignment_samples > 0 ? alignment_sum / gaudi::real(alignment_samples)
                              : 0.0;
    std::cerr << "single-fit Darboux evaluated on all vertices:"
              << " fit_vertex=" << __single_fit_vertex
              << " max_abs_D=" << max_abs_value << " rms_D=" << rms_value
              << " mean_abs_normal_alignment=" << mean_alignment
              << " samples=" << value_samples << std::endl;
  }

  void refresh_mesh_buffer() {
    gg::fillBuffer_ref(*__M, _objs[0], gg::colorRGB(0.72, 0.74, 0.78, 1.0));
  }

  void draw_normal_debug() {
    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const gaudi::real scale = 4.0 * gaudi::asawa::shell::avg_length(M, x);

    const gaudi::vec4 mesh_color(0.1, 0.85, 1.0, 1.0);
    const gaudi::vec4 cyclide_color(1.0, 0.15, 0.9, 1.0);
    const gaudi::vec4 axis_color(1.0, 0.9, 0.1, 1.0);
    const gaudi::vec4 fit_color(0.15, 1.0, 0.25, 1.0);

    gg::geometry_logger::line(
        __frame.center - 1.7 * __major_radius * __frame.z_axis,
        __frame.center + 1.7 * __major_radius * __frame.z_axis, axis_color);

    auto verts = M.get_vert_range();
    const gaudi::real marker = 0.08 * scale;
    gg::geometry_logger::line(__single_fit_origin - marker * __frame.x_axis,
                              __single_fit_origin + marker * __frame.x_axis,
                              fit_color);
    gg::geometry_logger::line(__single_fit_origin - marker * __frame.y_axis,
                              __single_fit_origin + marker * __frame.y_axis,
                              fit_color);
    gg::geometry_logger::line(__single_fit_origin - marker * __frame.z_axis,
                              __single_fit_origin + marker * __frame.z_axis,
                              fit_color);

    for (int k = 0; k < static_cast<int>(verts.size()); ++k) {
      const int vi = static_cast<int>(verts[static_cast<size_t>(k)]);
      const gaudi::vec3 p = x[vi];
      const gaudi::vec3 mesh_n = __vertex_normals[vi].normalized();
      gaudi::vec3 cyclide_n = __cyclide_normals[vi];
      if (cyclide_n.norm() < 1e-10) {
        continue;
      }
      cyclide_n.normalize();
      if (cyclide_n.dot(mesh_n) < 0.0) {
        cyclide_n *= -1.0;
      }

      gg::geometry_logger::line(p, p + scale * mesh_n, mesh_color);
      gg::geometry_logger::line(p, p + scale * cyclide_n, cyclide_color);
    }
  }

  virtual void onAnimate(int /*frame*/) {}

  virtual void onDraw(gg::Viewer &viewer) {
    draw_normal_debug();
    sync_gaudi_debug_to_gg();
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
  std::vector<gaudi::vec3> __vertex_normals;
  std::vector<gaudi::vec3> __cyclide_normals;
  gaudi::albers::vec14 __single_fit_Q = gaudi::albers::vec14::Zero();
  gaudi::vec3 __single_fit_origin = gaudi::vec3::Zero();
  gaudi::real __l0 = 0.0;
  gaudi::real __major_radius = 1.25;
  gaudi::real __minor_radius = 0.35;
  int __single_fit_vertex = kSingleFitVertex;
  int __stride = 4;
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height, int stride, gaudi::real axis_tilt) {
    return std::make_shared<App>(width, height, stride, axis_tilt);
  }

  App(int width, int height, int stride, gaudi::real axis_tilt)
      : gg::SimpleApp(width, height, 4.0, false, "darboux_cyclide_") {
    scene = Scene::create(stride, axis_tilt);
    this->setScene(scene);
    this->initUI();
  }

  void initUI() { performLayout(); }

  ~App() override = default;

  ScenePtr scene;
};

} // namespace

int main(int argc, char *argv[]) {
  try {
    int stride = 4;
    gaudi::real axis_tilt = 0.0;
    for (int i = 1; i < argc; ++i) {
      std::string arg(argv[i]);
      if (arg == "--stride" && i + 1 < argc) {
        stride = std::max(1, std::atoi(argv[++i]));
      } else if (arg == "--axis-tilt" && i + 1 < argc) {
        axis_tilt = std::atof(argv[++i]);
      }
    }

    nanogui::init();
    AppPtr app = App::create(1280, 740, stride, axis_tilt);
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
