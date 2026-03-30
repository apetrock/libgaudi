#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <random>
#include <vector>
#include <iostream>
#include <string>

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/calder/tree_code.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/geometry_logger.hpp"

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

using std::cerr;
using std::cout;
using std::endl;

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

namespace gaudi {
namespace duchamp {
using namespace asawa;

class aabb_build {
public:
  typedef std::shared_ptr<aabb_build> ptr;
  static ptr create() { return std::make_shared<aabb_build>(); }

  aabb_build() { load_shell(); }

  void load_shell() {
    __M = shell::load_bunny();
    shell::triangulate(*__M);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);
  }

  void step(int frame) {
    gaudi::geometry_logger::clear();

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<index_t> face_vert_ids = M.get_face_vert_ids();

    vec3 tri_center = (1.0 / 3.0) *
        (x[face_vert_ids[0]] + x[face_vert_ids[1]] + x[face_vert_ids[2]]);
    std::vector<vec3> queries = {tri_center};

    calder::visualize_shell_bvh(M, queries, 0.5);
  }

  shell::shell::ptr __M;
};

} // namespace duchamp
} // namespace gaudi

using namespace GaudiMath;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {
public:
  static ScenePtr create() { return std::make_shared<Scene>(); }
  Scene() : gg::Scene() { initScene(); }

  void initScene() {
    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    __surf = gaudi::duchamp::aabb_build::create();
    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);

    _color = gg::colorRGB(0.5, 0.5, 0.5, 1.0);
  }

  virtual void onAnimate(int frame) {
    __surf->step(frame);

    sync_gaudi_debug_to_gg();

    gg::fillBuffer_ref(*__surf->__M, _obj, _color);
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
  gaudi::duchamp::aabb_build::ptr __surf;
  std::vector<gg::DrawablePtr> mSceneObjects;
  gg::BufferObjectPtr _obj;
  gg::colorRGB _color;
};

class App;
using AppPtr = std::shared_ptr<App>;

class App : public gg::SimpleApp {
public:
  static AppPtr create(std::string file) {
    return std::make_shared<App>(file);
  }
  App(std::string file) : gg::SimpleApp(1280, 720, 4.0, true, "aabb_test_") {
    this->setScene(scene = Scene::create());
    this->set_rotate_ball(false);
    this->initUI();
  }
  void initUI() {
    using namespace nanogui;
    performLayout();
  }
  ~App() {}
  ScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    nanogui::init();
    AppPtr app = App::create(std::string(argv[0]));
    app->drawAll();
    app->setVisible(true);
    nanogui::mainloop();
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
