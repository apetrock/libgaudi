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
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/triangle_operations.hpp"

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

    initDebugLines();
    createGeometry();
    // testTree();
    // testDistance();
    test();
    // testNormals();
    fillWireFrame<space3>(_meshGraph, _debugLines);
    _debugLines->renderLines();
  }

  void initDebugLines() {
    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);
  }

  virtual void onDraw(gg::Viewer &viewer) {

    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    _debugLines->clear();
  }

  void createGeometry() {
    std::cout << "creating buffer" << std::endl;
    m2::obj_loader<space3> load;
    m2::affine<space3> mod;
    m2::make<space3> mk;
    m2::subdivide<space3> sub;
    std::cout << "  loading assets" << std::endl;
    //_meshGraph = &load("assets/messer.obj");
    _meshGraph = mk.cube(1.0, 1.0, 1.0);

    //_meshGraph = &sub.subdivide_control(*_meshGraph);
    //_meshGraph = &sub.subdivide_control(*_meshGraph);
    //_meshGraph = &sub.subdivide_control(*_meshGraph);
    //_meshGraph = &sub.subdivide_control(*_meshGraph);
    m2::remesh<space3> rem;
    rem.triangulate(_meshGraph);

    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();

    std::cout << "  init buffer" << std::endl;

    mod.centerGeometry(*_meshGraph);
    _obj = gg::BufferObject::create();
    _obj->init();
    gg::fillBuffer(_meshGraph, _obj);
    mSceneObjects.push_back(_obj);
  }

  template <typename SPACE>
  void test_merge(m2::surf<SPACE> *surf, int face0, int face1) {
    M2_TYPEDEFS;

    std::cout << "merging faces: " << face0 << " " << face1 << std::endl;

    std::vector<face_ptr> faces = surf->get_faces();
    faces[0]->print();
    std::vector<triangle_pair> pairs;
    pairs.push_back(triangle_pair(m2::ci::get_tris<SPACE>(faces[face0])[0],
                                  m2::ci::get_tris<SPACE>(faces[face1])[0]));

    m2::face_merger<SPACE> merger(surf, 1.0);
    merger.merge_collected(pairs, faces);
  }

  void testMergeBridgedFaces0() { 
    test_merge<space3>(_meshGraph, 0, 1);
    //test_merge<space3>(_meshGraph, 0, 1);
  }

  void testMergeBridgedFaces1() { test_merge<space3>(_meshGraph, 0, 3); }
  void testMergeBridgedFaces2() { test_merge<space3>(_meshGraph, 2, 3); }
  void testMergeBridgedFaces3() { test_merge<space3>(_meshGraph, 4, 6); }
  void testMergeSharedEdgeFaces() { test_merge<space3>(_meshGraph, 1, 3); }

  void testMergeSharedVertex() {
    m2::subdivide<space3> sub;
    _meshGraph = &sub.subdivide_control(*_meshGraph);

    m2::remesh<space3> rem;
    rem.triangulate(_meshGraph);

    test_merge<space3>(_meshGraph, 1, 55);

    /*
    _meshGraph = &sub.subdivide_control(*_meshGraph);
    rem.triangulate(_meshGraph);
    */
  }

  void test() {
    //testMergeBridgedFaces0();
    //testMergeBridgedFaces1();
    //testMergeBridgedFaces2();
    //testMergeBridgedFaces3();
    //testMergeSharedEdgeFaces();
    testMergeSharedVertex();

    gg::fillBuffer(_meshGraph, _obj);
  }

  template <typename SPACE>
  void fillWireFrame(m2::surf<space3> *surf, gg::DebugBufferPtr debugLines) {
    M2_TYPEDEFS;
    return;
    edge_ptr e0 = surf->get_edges()[10];
    int f0 = e0->v1()->face()->position_in_set();
    int f1 = e0->v1()->vnext()->face()->position_in_set();

    std::cout << "f0: " << f0 << std::endl;
    std::cout << "f1: " << f1 << std::endl;
    f0 = 1;
    f1 = 2;
    for (auto e : surf->get_edges()) {
      if (!e)
        continue;

      std::cout << e->v1()->edge()->position_in_set() << " "
                << e->v1()->vertex()->size() << " "
                << e->v1()->face()->position_in_set() << " "
                << e->v2()->face()->position_in_set() << std::endl;
      // test 0-6 is good test;

      Vec4 col = Vec4(1.0, 0.0, 0.0, 1.0);
      if (e->v1()->face()->position_in_set() == f0 ||
          e->v2()->face()->position_in_set() == f0) {
        col = Vec4(0.0, 1.0, 0.0, 1.0);
      }

      if (e->v1()->face()->position_in_set() == f1 ||
          e->v2()->face()->position_in_set() == f1) {
        col = Vec4(0.0, 0.0, 1.0, 1.0);
      }

      if (e->position_in_set() == 6 || e->position_in_set() == 7 ||
          e->position_in_set() == 14 || e->position_in_set() == 16) {
        col = Vec4(0.0, 1.0, 1.0, 1.0);
      }

      coordinate_type c0 = m2::ci::get_coordinate<SPACE>(e->v1()->vertex());
      coordinate_type c1 = m2::ci::get_coordinate<SPACE>(e->v2()->vertex());

      _debugLines->pushLine(Vec4(c0[0], c0[1], c0[2], 1.0),
                            Vec4(c1[0], c1[1], c1[2], 0.0), col);
    }
  }

private:
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj;
  gg::PointBufferPtr _points;
  gg::DebugBufferPtr _debugLines;

  m2::surf<space3> *_meshGraph;
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

    gg::SimpleAppPtr app = gg::SimpleApp::create();

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
