#include "GaudiMath/typedefs.hpp"
#include <algorithm>
#include <array>
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

#include "gaudi/duchamp/rod_guided_deformation.hpp"
#include "gaudi/duchamp/rod_guided_deformation_with_morph.hpp"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

// use #ifdef to defin screen size
// 4k
//  #define _SW 3840
//  #define _SH 2160

// HD
#define _SW 1920
#define _SH 1080

// 720p
// #define _SW 1280
// #define _SH 720
// #define _SW 1280
// #define _SH 720

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

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() { initScene(); }

  std::array<Vec3d, 2> get_rand_colors() {
    double gd = (3.0 - sqrt(5.0)) * M_PI;
    double phi = (1.0 + sqrt(5.0)) / 2.0;
    Vec2d c0 = random_vec2_with_angle();
    c0.normalize();
    gaudi::vec2 c1 = rotate_on_disk(c0, random_normal(1.0, 0.75) * M_PI);
    c1.normalize();

    Vec3d c03 = compose_color(c0);
    Vec3d c13 = compose_color(c1);

    double D = random_normal(0.1, 0.05);
    if (random_normal(0.0, 0.5) > 0) {
      c03 *= 1.0 + D;
      c13 *= 1.0 - D;
    } else {
      c03 *= 1.0 - D;
      c13 *= 1.0 + D;
    }
    return {c03, c13};
  }

  void initScene() {
    //_experiment = duchamp::mean_shift_experiment<growth>::create();

    _objs.resize(2);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);

    _objs[0] = gg::BufferObject::create();
    _objs[0]->init();
    mSceneObjects.push_back(_objs[0]);

    _objs[1] = gg::BufferObject::create();
    _objs[1]->init();
    mSceneObjects.push_back(_objs[1]);

    __surf = gaudi::duchamp::rod_guided_deformation_with_morph::create();

    for (int i = 0; i < 12; i++) {
      std::array<Vec3d, 2> colors = get_rand_colors();
      std::cout << " {vec3(" << colors[0].transpose() << "),"
                << "  vec3(" << colors[1].transpose() << ")}" << std::endl;
    }
    std::array<Vec3d, 2> colors_array = get_rand_colors();
    Vec3d c03 = colors_array[0];
    Vec3d c13 = colors_array[1];
    colors = {
        gg::colorRGB(c03.x(), c03.y(), c03.z(), 1.0),
        gg::colorRGB(c13.x(), c13.y(), c13.z(), 1.0),
    };
  }
  virtual void onAnimate(int frame) {
    
    __surf->step(frame);
    std::array<Vec3d, 2> colors_array = __surf->get_colors(frame);
    Vec3d c03 = colors_array[0];
    Vec3d c13 = colors_array[1];
    colors = {
        gg::colorRGB(c03.x(), c03.y(), c03.z(), 1.0),
        gg::colorRGB(c13.x(), c13.y(), c13.z(), 1.0),
    };

    gg::fillBuffer_ref(*__surf->__M, _objs[0], colors[0]);
    gg::fillBuffer_ref(*__surf->__R, _objs[1], colors[1]);

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
  // gaudi::duchamp::rod_guided_deformation::ptr __surf;
  gaudi::duchamp::rod_guided_deformation_with_morph::ptr __surf;

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
      : gg::SimpleApp(width, height, 2.0, true, "florp_drive_") {
    this->setScene(scene = Scene::create());
    this->set_rotate_ball(false);
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

    AppPtr app = App::create(_SW, _SH, std::string(argv[0]));

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
