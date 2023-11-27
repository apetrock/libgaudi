#include "GaudiMath/typedefs.hpp"
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
    //_experiment = duchamp::mean_shift_experiment<growth>::create();

    _objs.resize(1);

    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);

    /* //a second buffer object... for curves or...
    _objs[1] = gg::BufferObject::create();
    _objs[1]->init();
    mSceneObjects.push_back(_objs[1]);
    */

    __surf = gaudi::duchamp::growth_study::create();
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
    colors = {
        gg::colorRGB(0.0, 0.8, 0.4, 1.0),
        gg::colorRGB(0.0, 1.0, 1.0, 1.0),
    };
  }
  virtual void onAnimate(int frame) {

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
  static AppPtr create(int width, int height, std::string file) {
    return std::make_shared<App>(width, height, file);
  }

  typedef double Real;

  App(int width, int height, std::string file)
      : gg::SimpleApp(width, height, 4.0, true, "growth_study_") {
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

    AppPtr app = App::create(1280, 740, std::string(argv[0]));

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
