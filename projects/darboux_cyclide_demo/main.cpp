#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
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

namespace {

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

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
        gaudi::calder::darboux_cyclide(M, x, __vertex_normals, __l0, 3.0);
    __cyclide_normals.resize(fits.size(), gaudi::vec3::Zero());
    for (size_t i = 0; i < fits.size(); ++i) {
      gaudi::vec3 g = gaudi::albers::darboux_grad(fits[i], gaudi::vec3::Zero());
      __cyclide_normals[i] = -fits[i][9] * g;
    }

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

  void draw_normal_debug() {
    gaudi::asawa::shell::shell &M = *__M;
    std::vector<gaudi::vec3> &x = gaudi::asawa::get_vec_data(M, 0);
    const gaudi::real scale = 4.0 * gaudi::asawa::shell::avg_length(M, x);

    const gaudi::vec4 mesh_color(0.1, 0.85, 1.0, 1.0);
    const gaudi::vec4 cyclide_color(1.0, 0.15, 0.9, 1.0);
    const gaudi::vec4 axis_color(1.0, 0.9, 0.1, 1.0);

    gg::geometry_logger::line(
        __frame.center - 1.7 * __major_radius * __frame.z_axis,
        __frame.center + 1.7 * __major_radius * __frame.z_axis, axis_color);

    auto verts = M.get_vert_range();
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
  gaudi::real __l0 = 0.0;
  gaudi::real __major_radius = 1.25;
  gaudi::real __minor_radius = 0.35;
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
