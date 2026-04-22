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
// #include "gaudi/asawa/asawa.h"

#include "gaudi/duchamp/growth_study.hpp"
   
#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  void initScene() {
    std::cerr << "[growth_study] Scene::initScene start" << std::endl;
    //_experiment = duchamp::mean_shift_experiment<growth>::create();

    _objs.resize(1);

    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);
    std::cerr << "[growth_study] primary buffer initialized" << std::endl;

    /* //a second buffer object... for curves or...
    _objs[1] = gg::BufferObject::create();
    _objs[1]->init();
    mSceneObjects.push_back(_objs[1]);
    */

    __surf = gaudi::duchamp::growth_study::create();
     __surf->set_nan_probe(true, 1); // stderr: first NaN per field + tri/edge stats
    std::cerr << "[growth_study] surface created" << std::endl;
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
    colors = {
        gg::colorRGB(0.0, 0.8, 0.4, 1.0),
        gg::colorRGB(0.0, 1.0, 1.0, 1.0),
    };
  }
  virtual void onAnimate(int frame) {
    if (frame < 3 || frame % 60 == 0) {
      std::cerr << "[growth_study] onAnimate frame " << frame << std::endl;
    }

    __surf->step(frame);
    std::vector<gg::colorRGB> colors;
    auto mesh_colors = __surf->get_mesh_colors();
    for (auto col : mesh_colors) {
      colors.push_back(gg::colorRGB(col[0], col[1], col[2], 1.0));
    }

    gg::fillBuffer_ref(*__surf->__M, _objs[0], colors);

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
  gaudi::duchamp::growth_study::ptr __surf;
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
      : gg::SimpleApp(width, height, 4.0, false, "growth_study_") {
    std::cerr << "[growth_study] App ctor after SimpleApp" << std::endl;
    this->setScene(scene = Scene::create());
    std::cerr << "[growth_study] scene attached" << std::endl;
    this->initUI();
    std::cerr << "[growth_study] UI initialized" << std::endl;
  }

  void initUI() {
    using namespace nanogui;
    int w = 256;
    performLayout();
    // window->center();
  }

  ~App() {}

  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    std::cerr << "[growth_study] main start" << std::endl;
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    std::cerr << "[growth_study] nanogui::init" << std::endl;
    nanogui::init();

    std::cerr << "[growth_study] creating app" << std::endl;
    AppPtr app = App::create(1280, 740, std::string(argv[0]));
    std::cerr << "[growth_study] app created" << std::endl;

    // app->setScene(Scene::create());
    std::cerr << "[growth_study] drawAll" << std::endl;
    app->drawAll();
    std::cerr << "[growth_study] setVisible" << std::endl;
    app->setVisible(true);
    std::cerr << "[growth_study] entering mainloop" << std::endl;
    std::cerr << "[growth_study] focus the window, then press W to toggle "
                 "animation (simulation). S = single step. Q = quit.\n"
              << std::flush;
    nanogui::mainloop();
    std::cerr << "[growth_study] mainloop returned" << std::endl;
    // delete app;
    std::cerr << "[growth_study] nanogui::shutdown" << std::endl;
    nanogui::shutdown();
    std::cerr << "[growth_study] shutdown complete" << std::endl;

  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());
    std::cerr << error_msg << std::endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#endif
    return -1;
  } catch (const std::exception &e) {
    std::string error_msg =
        std::string("Caught a fatal std::exception: ") + std::string(e.what());
    std::cerr << error_msg << std::endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#endif
    return -2;
  } catch (...) {
    std::string error_msg = "Caught an unknown fatal error.";
    std::cerr << error_msg << std::endl;
#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#endif
    return -3;
  }

  return 0;
}
