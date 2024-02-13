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

#include "gaudi/duchamp/gravitys_rainbow.hpp"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

Vec2d random_vec2_with_angle() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0f, 2.0f * M_PI);
  double angle = dis(gen);
  return Vec2d(std::cos(angle), std::sin(angle));
}

double random_normal(double mean, double std_dev) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dis(mean, std_dev);
  return dis(gen);
}

Vec2d rotate_on_disk(const Vec2d &v, double angle) {
  double radius = v.norm();
  double theta = std::atan2(v.y(), v.x());
  theta += angle;
  double x = radius * std::cos(theta);
  double y = radius * std::sin(theta);
  return Vec2d(x, y);
}

Vec3d compose_color(const Vec2d &v) {
  double angle = std::atan2(v.y(), v.x());
  double red = std::cos(angle);
  double green = std::cos(angle + M_PI * 2.0f / 3.0f);
  double blue = std::cos(angle + M_PI * 4.0f / 3.0f);
  return Vec3d(red, green, blue) * 0.5f + Vec3d(0.5f, 0.5f, 0.5f);
}

std::array<Vec3d, 3> get_rand_colors() {
  double gd = (3.0 - sqrt(5.0)) * M_PI;
  double phi = (1.0 + sqrt(5.0)) / 2.0;
  Vec2d c0 = random_vec2_with_angle();
  c0.normalize();
  gaudi::vec2 c1 = rotate_on_disk(c0, random_normal(1.0, 0.75) * M_PI);
  c1.normalize();
  gaudi::vec2 c2 = rotate_on_disk(c0, random_normal(1.0, 0.75) * M_PI);
  c2.normalize();

  Vec3d c03 = compose_color(c0);
  Vec3d c13 = compose_color(c1);
  Vec3d c23 = compose_color(c2);

  double D = random_normal(0.1, 0.05);
  if (random_normal(0.0, 0.5) > 0) {
    c03 *= 1.0 + D;
    c13 *= 1.0 - D;
  } else {
    c03 *= 1.0 - D;
    c13 *= 1.0 + D;
  }
  return {c03, c13, c23};
}

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  void initScene() {
    //_experiment = duchamp::mean_shift_experiment<growth>::create();
    __surf = gaudi::duchamp::rainbow::create();

    _objs.resize(__surf->__R.size());
    for (int i = 0; i < __surf->__R.size(); i++) {
      _objs[i] = gg::BufferObject::create();
      _objs[i]->init();
      mSceneObjects.push_back(_objs[i]);
    }
    std::array<Vec3d, 3> colors_xyz = get_rand_colors();
    // use std::transform to convert colors_xyz to std::vector<gg::colorRGB>
    std::transform(
        colors_xyz.begin(), colors_xyz.end(), this->colors.begin(),
        [](const Vec3d &v) { return gg::colorRGB(v.x(), v.y(), v.z(), 1.0); });

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  virtual void onAnimate(int frame) {

    __surf->step(frame);

    for (int i = 0; i < __surf->__R.size(); i++) {
      gg::fillBuffer_ref(*__surf->__R[i], _objs[i], colors[i]);
    }
    //   std::cout << "rendering debug" << std::endl;
    //   asawa::test();

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
  // gaudi::duchamp::fast_summation_test::ptr __surf;
  gaudi::duchamp::rainbow::ptr __surf;
  std::array<gg::colorRGB, 3> colors;
  std::vector<gg::DrawablePtr> mSceneObjects;
  std::vector<gg::BufferObjectPtr> _objs;
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
  static AppPtr create(int w, int h, std::string file) {
    return std::make_shared<App>(w, h, file);
  }

  typedef double Real;

  App(int w, int h, std::string file)
      : gg::SimpleApp(w, h, 6.0, true, "gravitys_rainbow_") {

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

    AppPtr app = App::create(1920, 1080, std::string(argv[0]));

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
