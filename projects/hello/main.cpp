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
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/bins.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/moving_mesh.hpp"

#include "manifold/vec_addendum.h"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

    using std::cerr;
using std::cout;
using std::endl;

template <typename SPACE>
void debugVorticity(m2::surf<SPACE> *mesh, gg::DebugBufferPtr debug) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::cout << "  rendering vorticity" << std::endl;

  for (auto f : mesh->get_faces()) {
    coordinate_type v = 0.1 * f->data;
    coordinate_type cen = f->calc_center();

    coordinate_type l0 = cen - 1.0 * v;
    coordinate_type l1 = cen + 1.0 * v;
    T mag = m2::va::norm(v);
    T pi = M_PI;
    gg::Vec4 col(0.5 + 0.5 * cos(2.0 * pi * mag + 0.000 * pi),
             0.5 + 0.5 * cos(2.0 * pi * mag + 0.666 * pi),
             0.5 + 0.5 * cos(2.0 * pi * mag + 1.333 * pi), mag);
    //gg::Vec4 col(0.5, 0.5, 0.5, 1.0);
    debug->pushLine(gg::Vec4(l0[0], l0[1], l0[2], 1.0),
                    gg::Vec4(l1[0], l1[1], l1[2], 1.0), col);
  }

  debug->renderLines();
};

class VortexScene;
using VortexScenePtr = std::shared_ptr<VortexScene>;

class VortexScene : public gg::Scene {

public:
  static VortexScenePtr create() { return std::make_shared<VortexScene>(); }

  VortexScene() : gg::Scene() {
    using namespace nanogui;

    m2::obj_loader<space3> load;
    m2::subdivide<space3> sub;
    m2::make<space3> mk;
    m2::convex_hull<space3> ch;
    m2::add_handle<space3> ah;

    m2::construct<space3> bevel;
    m2::modify<space3> mod;

    //_meshGraph = &load("assets/bunny.obj");
    //_meshGraph = &load("assets/messer.obj");
    // std::cout << "--make cube" << std::endl;
    //_meshGraph =  mk.cube(0.05,1.0,1.0);
    _meshGraph = mk.cube(1.0, 1.0, 1.0);

    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    //_meshGraph = &sub.subdivide_control(*_meshGraph);

    std::cout << "--center" << std::endl;
    std::cout << "--update_all" << std::endl;
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();
    std::cout << "--build sheet" << std::endl;
    _vortex = new m2::vortex_sheet<space3>(_meshGraph);

    mod.centerGeometry(*_meshGraph);
    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);
  }

  virtual void onAnimate() {
    std::cout << "====== " << std::endl;

    std::cout << " vertices: " << _meshGraph->get_vertices().size()
              << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;

    m2::modify<space3> mod;
    double max = 0.025;
    double min = max / 10.0;
    double dt = 0.1;

    _vortex->minLength = max;
    _vortex->minCollapseLength = min;
    _vortex->edgeJoinThresh = 0.25 * min;
    _vortex->regLength = max;

    _vortex->integrateBaroclinity(dt, 0.004);
    //_vortex->addCurveVorticity(0.2, 1);
    _meshGraph->verify();
    _vortex->updateCirculation();
    _meshGraph->verify();

    _vortex->integrateVelocityRK2(dt);
    _meshGraph->verify();
    _vortex->remesh();
    _meshGraph->verify();

    _vortex->updateVorticity();
    m2::remesh<space3> rem;
    rem.triangulate(_meshGraph);
    debugVorticity(_meshGraph, _debugLines);

    _meshGraph->update_all();
    _meshGraph->pack();
    mod.centerGeometry(*_meshGraph);

    gg::fillBuffer(_meshGraph, _obj);
  }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
    _debugLines->clear();
  }

private:
  gg::DebugBufferPtr _debugLines = NULL;
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<space3> *_meshGraph;
  m2::vortex_sheet<space3> *_vortex;
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

    app->setScene(VortexScene::create());
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
