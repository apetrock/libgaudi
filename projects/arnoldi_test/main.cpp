#include "GaudiMath/typedefs.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
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

#include <complex>
#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "gaudi/duchamp/arnoldi_test.hpp"
#include "gaudi/geometry_logger.hpp"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

void sync_gaudi_debug_to_gg() {
  const auto &lines = gaudi::geometry_logger::get_lines();
  const auto &line_cols = gaudi::geometry_logger::get_line_colors();
  for (size_t i = 0; i + 1 < lines.size(); i += 2) {
    gg::geometry_logger::line(lines[i], lines[i + 1], line_cols[i]);
  }

  const auto &points = gaudi::geometry_logger::get_points();
  const auto &point_cols = gaudi::geometry_logger::get_point_colors();
  for (size_t i = 0; i < points.size(); i++) {
    gg::geometry_logger::point(points[i], point_cols[i]);
  }
}

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  void initScene() {
    _objs.resize(1);

    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);

    __surf = gaudi::duchamp::arnoldi_test::create();
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
    colors = {
        gg::colorRGB(0.5, 0.5, 0.5, 1.0),
        gg::colorRGB(0.0, 1.0, 1.0, 1.0),
    };
  }

  virtual void onAnimate(int frame) {
    gaudi::geometry_logger::clear();

    __surf->step(frame);

    sync_gaudi_debug_to_gg();

    gg::fillBuffer_ref(*__surf->__M, _objs[0], colors[0]);

    gg::geometry_logger::render();
  }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    gg::geometry_logger::clear();
  }

private:
  gaudi::duchamp::arnoldi_test::ptr __surf;
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
  vector<gg::colorRGB> colors;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(int width, int height, std::string file) {
    return std::make_shared<App>(width, height, file);
  }

  typedef double Real;

  App(int width, int height, std::string file)
      : gg::SimpleApp(width, height, 4.0, false, "arnoldi_test_") {
    this->setScene(scene = Scene::create());
    this->initUI();
  }

  void initUI() {
    using namespace nanogui;
    int w = 256;
    performLayout();
  }

  ~App() {}

  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc << " arguments:" << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    AppPtr app = App::create(1280, 740, std::string(argv[0]));

    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    nanogui::shutdown();

  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());
    std::cerr << error_msg << endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#endif
    return -1;
  }

  return 0;
}
