//#include "nanoguiincludes.h"

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

#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/bins.hpp"
#include "manifold/fdt.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/moving_mesh.hpp"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

class MainScene;
using MainScenePtr = std::shared_ptr<MainScene>;

class MainScene : public gg::Scene {

public:
  static MainScenePtr create() { return std::make_shared<MainScene>(); }

  MainScene() : gg::Scene() {
    using namespace nanogui;

    m2::obj_loader<space3> load;
    m2::modify<space3> mod;

    _meshGraph = &load("assets/messer.obj");
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();

    mod.centerGeometry(*_meshGraph);
    std::cout << "creating buffer" << std::endl;
    _obj = gg::BufferObject::create();
    _obj->initBuffer();
    gg::fillBuffer(_meshGraph, _obj);
    mSceneObjects.push_back(_obj);

    std::cout << "creating Points" << std::endl;
    auto randNormalVec = [](float mean, float std) {
      auto randomFunc =
          [distribution_ = std::normal_distribution<float>(mean, std),
           random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
            return Vec4(distribution_(random_engine_),
                        distribution_(random_engine_),
                        distribution_(random_engine_), 1.0);
            ;
          };
      return randomFunc;
    };

    int N = 2 << 16;
    std::vector<GaudiMath::Vec4> positions;
    std::generate_n(std::back_inserter(positions), N, randNormalVec(0, 0.5));
    std::vector<GaudiMath::Vec4> colors;
    std::transform(positions.begin(), positions.end(),
                   std::back_inserter(colors), [](const Vec4 pos) {
                     float d = pos(Eigen::seq(0, 2), 0).norm();
                     float pi = M_PI;
                     return Vec4(0.5 + 0.5 * cos(2.0 * pi * d + 0.000 * pi),
                                 0.5 + 0.5 * cos(2.0 * pi * d + 0.666* pi),
                                 0.5 + 0.5 * cos(2.0 * pi * d + 1.333 * pi), 1.0);
                   });

    _points = gg::PointBuffer::create();
    _points->initBuffer();
    gg::insertPoints<space3>(positions, colors, _points);
    mSceneObjects.push_back(_points);
  }

  virtual void onAnimate() {

    m2::modify<space3> mod;
    gg::fillBuffer(_meshGraph, _obj);
  }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::BufferObjectPtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
  }

private:
  std::vector<gg::BufferObjectPtr> mDebug;
  std::vector<gg::BufferObjectPtr> mSceneObjects;

  gg::BufferObjectPtr _obj;
  gg::PointBufferPtr _points;
  gg::LineBufferPtr _lines;

  m2::control<space3> *_meshGraph;
};

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    gg::SimpleAppPtr app = gg::SimpleApp::create(std::string(argv[0]));

    app->setScene(MainScene::create());
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
    // delete app;
    nanogui::shutdown();

  }

  catch (const std::runtime_error &e) {
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
