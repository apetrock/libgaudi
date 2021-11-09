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

#include "manifold/harmonic_integrators.hpp"
#include "manifold/tree_code.hpp"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

template <typename SPACE>
m2::aabb_tree<SPACE, typename SPACE::triangle_type>
build_tree(m2::control<SPACE> *mesh) {
  using namespace nanogui;
  M2_TYPEDEFS;
  std::cout << "  building tree" << std::endl;

  std::vector<face_ptr> faces = mesh->get_faces();
  std::vector<triangle_type> tris;

  std::for_each(faces.begin(), faces.end(), [&tris](const face_ptr &f) {
    std::vector<triangle_type> ftris = f->get_tris();
    tris.insert(tris.end(), ftris.begin(), ftris.end());
  });

  m2::aabb_tree<SPACE, triangle_type> tree(tris);

  using tree_type = m2::aabb_tree<SPACE, triangle_type>;
  std::cout << "tree.nodes.size(): " << tree.nodes.size() << std::endl;

  return tree;
};

class MainScene;
using MainScenePtr = std::shared_ptr<MainScene>;

class MainScene : public gg::Scene {

public:
  static MainScenePtr create() { return std::make_shared<MainScene>(); }

  MainScene() : gg::Scene() {
    using namespace nanogui;

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();

    mSceneObjects.push_back(_debugLines);

    createGeometry();
    //testTree();
    //testDistance();
    testAvg();
    //testNormals();
    _debugLines->renderLines();
  }

  space3::vec3 rainbow(double d) {
    double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
    double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
    double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
    return space3::vec3(r, g, b);
  }

  void createGeometry() {
    std::cout << "creating buffer" << std::endl;
    m2::obj_loader<space3> load;
    m2::modify<space3> mod;

    std::cout << "  loading assets" << std::endl;
    _meshGraph = &load("assets/messer.obj");
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();

    std::cout << "  init buffer" << std::endl;

    mod.centerGeometry(*_meshGraph);
    //_obj = gg::BufferObject::create();
    //_obj->init();
    // gg::fillBuffer(_meshGraph, _obj);
    // mSceneObjects.push_back(_obj);
  }

  std::vector<space3::vec3> createPoints(int N) {
    auto randNormalVec = [](float mean, float std) {
      auto randomFunc =
          [distribution_ = std::normal_distribution<double>(mean, std),
           random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
            return space3::vec3(distribution_(random_engine_),
                                distribution_(random_engine_),
                                distribution_(random_engine_));
            ;
          };
      return randomFunc;
    };

    std::vector<space3::vec3> points;
    std::generate_n(std::back_inserter(points), N, randNormalVec(0, 0.5));
    return points;
  }

  void testTree() {
    auto tree = build_tree(_meshGraph);
    for (int i = 0; i < tree.nodes.size(); i++) {
      auto &node = tree.nodes[i];
      if (node.size == 0)
        continue;

      Vec4 col(0.5, 0.5, 0.5, 1.0);
      auto half = node.bbox.half;
      auto cen = node.bbox.center;
      if (node.isLeaf())
        _debugLines->pushBox(Vec4(cen[0], cen[1], cen[2], 1.0),
                             Vec4(half[0], half[1], half[2], 0.0), col);
    }
  }

  void testDistance() {
    std::cout << "creating Points" << std::endl;
    int N = 2 << 15;
    std::vector<space3::vec3> positions = createPoints(N);
    m2::mesh_calculator<space3> calc;
    std::vector<double> dist = calc.calcDistanceFromMesh(_meshGraph, positions);

    std::vector<space3::vec3> colors;
    for (int i = 0; i < positions.size(); i++) {
      double d = dist[i];
      double pi = M_PI;
      double r = d > 0.01 ? 4.0 * d : 0;
      double b = d < 0.01 ? 50.0 * d : 0;

      // std::cout << "d: " << d << std::endl;
      colors.push_back(space3::vec3(r, 0.0, b));
    }

    _points = gg::PointBuffer::create();
    _points->init();
    gg::insertPoints<space3>(positions, colors, _points);
    mSceneObjects.push_back(_points);
  }

  void testAvg() {
    std::cout << "creating Points" << std::endl;
    int N = 2 << 14;
    std::vector<space3::vec3> positions = createPoints(N);

    std::cout << "creating vecs" << std::endl;
    std::vector<space3::triangle_type> triangles;

    auto &faces = _meshGraph->get_faces();
    for (int i = 0; i < faces.size(); i++) {
      std::vector<space3::triangle_type> tris = faces[i]->get_tris();
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    int Nf = triangles.size();
    std::vector<space3::vec3> faceVectors = createPoints(Nf);
    std::for_each(faceVectors.begin(), faceVectors.end(),
                  [](auto &fv) { fv.normalize(); });

    double area = std::accumulate(
        triangles.begin(), triangles.end(), 0.0,
        [](double a, auto &triangle) { return a + triangle.area(); });
    area /= double(triangles.size());
    double reg = 2.0 * sqrt(area);

#if 0
    for (int i = 0; i < triangles.size(); i++) {

      auto p = triangles[i].center();
      auto a = faceVectors[i];
      auto pa = p + 1.0 * a;
      auto c = rainbow(m2::va::norm(a));
      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
    return;
#endif
    std::cout << "computing harmonic avg" << std::endl;
    m2::mesh_calculator<space3> calc;
    std::vector<space3::vec3> avgs =
        calc.harmonicAvg(_meshGraph, faceVectors, positions, reg);

    std::vector<space3::vec3> colors;
    for (int i = 0; i < avgs.size(); i++) {
      auto p = positions[i];
      auto a = avgs[i];
      auto pa = p + 0.25 * a;
      auto c = rainbow(m2::va::norm(a));
      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

  void testNormals() {
    std::cout << "creating Points" << std::endl;
    int N = 2 << 12;
    std::vector<space3::vec3> positions = createPoints(N);

    std::cout << "creating vecs" << std::endl;
    std::vector<space3::triangle_type> triangles;
    auto &faces = _meshGraph->get_faces();
    std::vector<space3::vec3> normals;

    for (int i = 0; i < faces.size(); i++) {
      std::vector<space3::triangle_type> tris = faces[i]->get_tris();
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    double area = 0;
    for (int i = 0; i < triangles.size(); i++) {
      normals.push_back(triangles[i].normal());
      area += triangles[i].area();
    }
    area /= double(triangles.size());
    double reg = sqrt(area);
    std::cout << "calc'd reg: " << reg << std::endl;
#if 0
    for (int i = 0; i < triangles.size(); i++) {

      auto p = triangles[i].center();
      auto a = normals[i];
      auto pa = p + 0.1 * a;
      auto c = rainbow(m2::va::norm(a));
      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
    return;
#endif

    std::cout << "computing harmonic avg" << std::endl;
    m2::mesh_calculator<space3> calc;
    std::vector<space3::vec3> avgs =
        calc.harmonicAvg(_meshGraph, normals, positions, reg);

    std::vector<space3::vec3> colors;
    double mx =
        std::accumulate(avgs.begin(), avgs.end(), 0.0, [](double a, auto &c) {
          return max(a, m2::va::norm(c));
        });

    for (int i = 0; i < avgs.size(); i++) {
      auto p = positions[i];
      auto a = avgs[i];
      a = a / mx;

      auto pa = p + 0.25 * a;
      // std::cout << a.transpose() << std::endl;
      auto c = rainbow(m2::va::norm(a));
      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

  virtual void onAnimate() {

    m2::modify<space3> mod;
    gg::fillBuffer(_meshGraph, _obj);
  }

  virtual void onDraw(gg::Viewer &viewer) {
    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });
  }

private:
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj;
  gg::PointBufferPtr _points;
  gg::DebugBufferPtr _debugLines;

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
