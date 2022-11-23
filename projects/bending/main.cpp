#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <vector>
#if defined(WIN32)
#include <windows.h>
#endif

#include <stdio.h> /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <complex>
#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"
#include "manifold/asawa/asawa.h"

#include "manifold/bins.hpp"
#include "manifold/duchamp/growth_experiments.hpp"
#include "manifold/duchamp/surface_experiments.hpp"

#include "manifold/vec_addendum.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

typedef duchamp::growth_space<double> growth;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  void initScene() {
    //_experiment = duchamp::mean_shift_experiment<growth>::create();
    _experiment = duchamp::covariant_flow_experiment<growth>::create();

    _experiment->init();

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  virtual void onAnimate(int frame) {
    _experiment->step(frame);
    auto colors = _experiment->getColor();
    auto surf = _experiment->get_surf();

    gg::fillBuffer(surf, _obj, colors);
    std::cout << "rendering debug" << std::endl;
    gg::geometry_logger::render();

    if (frame % 10 == 0) {
      _experiment->dump_obj("growth", frame);
      _experiment->dump_gaudi("growth", frame);
    }
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
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  duchamp::simple_experiment_base<growth>::ptr _experiment;
};

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(std::string file) { return std::make_shared<App>(file); }

  typedef double Real;

  App(std::string file) : gg::SimpleApp() {
    this->setScene(scene = Scene::create());
    this->initUI();
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
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    AppPtr app = App::create(std::string(argv[0]));

    // app->setScene(Scene::create());
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    // delete app;
    nanogui::shutdown();

  } catch (const std::runtime_error &e) {
    std::string error_msg =
        std::string("Caught a fatal error: ") + std::string(e.what());

#if defined(WIN32)
    MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
#else
    std::cerr << error_msg << endl;
#endif

    return -1;
  }

  return 0;
}
