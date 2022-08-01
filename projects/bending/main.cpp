//#include "nanoguiincludes.h"

#include "manifold/coordinate_interface.hpp"
#include "manifold/m2.hpp"
#include <exception>
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
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/conj_grad.hpp"
#include "manifold/diffuse.hpp"
#include "manifold/laplacian.hpp"

#include "manifold/bins.hpp"
#include "manifold/conj_grad.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/position_optimization.hpp"

#include "manifold/harmonic_integrators.hpp"
#include "manifold/vec_addendum.h"

#include "manifold/triangle_operations.hpp"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

template <typename T> class stretch_space {
public:
  typedef T double_type;
  typedef T real;
  typedef std::complex<T> complex;
  // always use homogeneous coordinates, provides decent error checking
  typedef Eigen::Matrix<T, 3, 1> coordinate_type;

  typedef line<T, coordinate_type> line_type;
  typedef triangle<T, coordinate_type> triangle_type;
  typedef swept_point<T, coordinate_type> swept_point_type;
  typedef swept_triangle<T, coordinate_type> swept_triangle_type;
  typedef bounding_box<T, coordinate_type> box_type;
  typedef Eigen::Matrix<T, 2, 2> mat2;
  typedef Eigen::Matrix<T, 3, 3> mat3;
  typedef Eigen::Matrix<T, 4, 4> mat4;
  typedef Eigen::Matrix<T, 4, 3> mat43;

  typedef unsigned short ushort;
  typedef unsigned int uint;
  typedef unsigned long ulong;

  typedef double double_t;
  typedef float float_t;

  typedef Eigen::Quaternion<T> quat;
  typedef Eigen::Matrix<T, 2, 1> vec2;
  typedef Eigen::Matrix<T, 3, 1> vec3;
  typedef Eigen::Matrix<T, 4, 1> vec4;

  typedef Eigen::Matrix<uint, 2, 1> uint2;
  typedef Eigen::Matrix<uint, 2, 1> uint4;

  typedef Eigen::Matrix<T, 2, 1> int2;
  typedef Eigen::Matrix<T, 4, 1> int3;
  typedef Eigen::Matrix<T, 4, 1> int4;

  enum class face_vertex_index { BARY = 0, MAXINDEX = 1 };

  enum class edge_index { PHI0 = 0, PHI1 = 1, MAXINDEX = 2 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    SMOOTH = 2,
    MAXINDEX = 3
  };

  enum class face_index { NORMAL = 0, CENTER = 1, AREA = 2, MAXINDEX = 3 };

  static storage_type get_type(face_vertex_index idx) {
    switch (idx) {
    case face_vertex_index::BARY:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    };
  }
  static storage_type get_type(edge_index idx) {
    switch (idx) {
    case edge_index::PHI0:
      return storage_type::REAL;
    case edge_index::PHI1:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }

  static storage_type get_type(face_index idx) {
    switch (idx) {
    case face_index::NORMAL:
      return storage_type::VEC3;
    case face_index::CENTER:
      return storage_type::VEC3;
    case face_index::AREA:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }

  static storage_type get_type(vertex_index idx) {
    switch (idx) {
    case vertex_index::COORDINATE:
      return storage_type::VEC3;
    case vertex_index::COLOR:
      return storage_type::VEC3;
    case vertex_index::SMOOTH:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }
};

typedef stretch_space<double> stretch;

class Scene;
using ScenePtr = std::shared_ptr<Scene>;

class Scene : public gg::Scene {

public:
  static ScenePtr create() { return std::make_shared<Scene>(); }

  Scene() : gg::Scene() {
    initScene();
    // initUI();
  }

  void initScene() {

    m2::obj_loader<stretch> load;
    m2::subdivide<stretch> sub;
    m2::make<stretch> mk;
    m2::convex_hull<stretch> ch;

    m2::construct<stretch> bevel;
    m2::affine<stretch> mod;
    std::string start_frame = "";
    // std::string start_frame = "stretch.46.gaudi";
    if (!start_frame.empty()) {
      FILE *file;
      file = fopen(start_frame.c_str(), "rb");
      if (file != NULL) {
        this->load_gaudi(start_frame);
      }
    } else {

      _meshGraph = &load("assets/bunny.obj");
      //_meshGraph = &load("assets/messer.obj");
      // std::cout << "--make cube" << std::endl;
      //_meshGraph =  mk.cube(0.05,1.0,1.0);
      //_meshGraph = mk.cube(1.0, 1.0, 1.0);

      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      m2::remesh<stretch> rem;
      rem.triangulate(_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      std::cout << "--center" << std::endl;
      std::cout << "--update_all" << std::endl;
      _meshGraph->update_all();
      std::cout << "--pack" << std::endl;
      _meshGraph->pack();

      mod.centerGeometry(*_meshGraph);
    }

    int N = 0;
    _integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.5, 3.0, 0.35);
    //_integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.1, 3.0, 0.5);
    _integrator->add_default_vertex_policy<typename stretch::real>(
        stretch::vertex_index::SMOOTH);
    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);
  }

  template <typename SPACE>
  vector<m2::colorRGB> getColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);

    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    std::vector<m2::colorRGB> vert_colors(smooth.size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = K[i];
      typename SPACE::real N = 0;
      typename SPACE::real s = smooth[i];

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorC(0.56, 0.50, 0.60);

      typename SPACE::coordinate_type mx = m2::va::mix(s, colorS, colorC);
      vert_colors[i] = m2::colorRGB(mx[0], mx[1], mx[2], 1.0);
      i++;
    }
    return vert_colors;
  }

#if 1
  space3::vec3 rainbow(double d) {
    double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
    double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
    double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
    return space3::vec3(r, g, b);
  }

  space3::vec3 grey(double d) { return space3::vec3(0.5, 0.5, 0.5); }
  space3::vec3 red(double d) { return space3::vec3(1.0, 0.0, 0.0); }
  space3::vec3 green(double d) { return space3::vec3(0.0, 1.0, 0.0); }
  space3::vec3 blue(double d) { return space3::vec3(0.0, 0.0, 1.0); }

  template <typename SPACE>
  void print_vecs(const std::vector<typename SPACE::vec3> &p0,
                  const std::vector<typename SPACE::vec3> &p1, double D = 0.1,
                  int col = 0) {
    double mx =
        std::accumulate(p1.begin(), p1.end(), 0.0, [](double a, auto &c) {
          return max(a, m2::va::norm(c));
        });

    for (int i = 0; i < p0.size(); i++) {

      const auto &p = p0[i];
      const auto &a = p1[i];

      auto pp0 = p - D * a / mx;
      auto pp1 = p + D * a / mx;

      // std::cout << a.transpose() << std::endl;
      auto c = grey(m2::va::norm(a));
      if (col == 1)
        c = red(m2::va::norm(a));
      if (col == 2)
        c = green(m2::va::norm(a));
      if (col == 3)
        c = blue(m2::va::norm(a));

      _debugLines->pushLine(Vec4(pp0[0], pp0[1], pp0[2], 1.0),
                            Vec4(pp1[0], pp1[1], pp1[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

#endif

  template <typename SPACE>
  void
  init_edge_constraints(m2::surf<SPACE> *surf,
                        typename m2::constraint_set<SPACE>::ptr constraints) {
    M2_TYPEDEFS;

    std::vector<edge_ptr> edges = _meshGraph->get_edges();
    std::vector<size_t> indices;
    std::vector<stretch::real> lengths;

    for (auto e : edges) {
      vertex_ptr v0 = e->v1()->vertex();
      vertex_ptr v1 = e->v2()->vertex();
      vertex_ptr v2 = e->v1()->prev()->vertex();
      vertex_ptr v3 = e->v2()->prev()->vertex();

      size_t i0 = v0->position_in_set();
      size_t i1 = v1->position_in_set();
      size_t i2 = v2->position_in_set();
      size_t i3 = v3->position_in_set();

      indices.push_back(i0);
      indices.push_back(i1);
      indices.push_back(i2);
      indices.push_back(i3);

      stretch::coordinate_type c0 = m2::ci::get_coordinate<stretch>(v0);
      stretch::coordinate_type c1 = m2::ci::get_coordinate<stretch>(v1);
      stretch::coordinate_type c2 = m2::ci::get_coordinate<stretch>(v2);
      stretch::coordinate_type c3 = m2::ci::get_coordinate<stretch>(v3);
      stretch::real l01 = m2::va::norm(stretch::coordinate_type(c0 - c1));
      stretch::real l23 = m2::va::norm(stretch::coordinate_type(c3 - c2));
      lengths.push_back(l01);
      lengths.push_back(l23);
    }

    int i = 0;
    int N = 0.5 * lengths.size();

    for (i = 0; i < N; i++) {
      size_t i0 = indices[4 * i + 0];
      size_t i1 = indices[4 * i + 1];
      size_t i2 = indices[4 * i + 2];
      size_t i3 = indices[4 * i + 3];

      real l01 = lengths[2 * i + 0];
      real l23 = lengths[2 * i + 1];
#if 1
      typename m2::edge_stretch<SPACE>::ptr stretch01 =
          m2::edge_stretch<SPACE>::create(i0, i1, l01);
      constraints->add_constraint(stretch01);
#endif
#if 0
      typename m2::edge_stretch<SPACE>::ptr stretch23 =
          m2::edge_stretch<SPACE>::create(i2, i3, l23);
      constraints->add_constraint(stretch23);
#endif

#if 0
      typename edge_length<SPACE>::ptr length01 =
          edge_length<SPACE>::create(i0, i1);
      constraints->add_constraint(length01);
#endif

#if 0
      typename m2::cross_ratio<SPACE>::ptr cross =
          m2::cross_ratio<SPACE>::create(i0, i1, i2, i3, 1.1);
      constraints->add_constraint(cross);
#endif
    }
  }

  double smoothstep(double t, double l, double h, double t0, double x0) {
    double x = (l * t + t0);
    double fx = 3.0 * x * x - 2.0 * x * x * x;
    std::cout << " x/t: " << t << " " << x << std::endl;
    if (x <= 0)
      return x0;
    if (x > 0 && x < 1)
      return x0 + h * fx;
    if (x >= 1)
      return x0 + h;
  }

  virtual void onAnimate(int frame) {

    _meshGraph->update_all();
    _meshGraph->reset_flags();
    _meshGraph->pack();
    if (frame == 1)
      _integrator->integrate();

    std::cout << "frame: " << frame << std::endl;
    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<stretch>(_meshGraph)
              << std::endl;

    double dt = _params.dt;
    auto colors = getColor(_meshGraph);
#if 1

    std::vector<stretch::vec3> positions =
        m2::ci::get_coordinates<stretch>(_meshGraph);

    // build constraints to capture current config
    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    using edge_ptr = typename m2::surf<stretch>::edge *;
    using vertex_ptr = typename m2::surf<stretch>::vertex *;

    m2::optimizer<stretch> opt;

    m2::constraint_set<stretch>::ptr constraints =
        m2::constraint_set<stretch>::create(_meshGraph);

    // init_edge_constraints<stretch>(_meshGraph, constraints);
    constraints->add_constraint(m2::bend<stretch>::create(_meshGraph));

#if 1
    // perturb config
    std::vector<stretch::vec3> normals =
        m2::ci::get_vertex_normals<stretch>(_meshGraph);
    positions[0] += 0.1 * normals[0];
#endif

    constraints->set_positions(positions);
    opt.update(constraints);
    positions = constraints->get_positions();
    m2::ci::set_coordinates<stretch>(positions, _meshGraph);

#if 0
    std::vector<stretch::mat3> blocks = opt.block_diag;
    int i = 0;
    std::vector<stretch::coordinate_type> v0s;
    std::vector<stretch::coordinate_type> v1s;
    std::vector<stretch::coordinate_type> v2s;
    for (auto &block : blocks) {
      //std::cout << block << std::endl; 
      //std::cout << std::endl;
      stretch::coordinate_type p0 = positions[i++];
      stretch::coordinate_type v0 = block.block(0, 0, 3, 1);
      stretch::coordinate_type v1 = block.block(0, 1, 3, 1);
      stretch::coordinate_type v2 = block.block(0, 2, 3, 1);
      v0s.push_back(p0 + v0);
      v1s.push_back(p0 + v1);
      v2s.push_back(p0 + v2);
    }
    print_vecs<stretch>(positions, v0s, 0.1, 1);
    print_vecs<stretch>(positions, v1s, 0.1, 2);
    print_vecs<stretch>(positions, v2s, 0.1, 3);
    std::cout << normals[0] << std::endl;
#endif
#if 1
    // std::cout << "print vecs" << std::endl;

    // print_vecs<stretch>(positions, normals);
    std::cout << "rendering debug" << std::endl;
    _debugLines->renderLines();
#endif

    m2::ci::set_coordinates<stretch>(positions, _meshGraph);
    // this->dump_gaudi(frame);
    _meshGraph->print();

    gg::fillBuffer(_meshGraph, _obj, colors);
#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "stretch." << frame << ".obj";
    m2::write_obj<stretch>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    m2::flattened_surf<stretch> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    m2::flattened_surf<stretch> fsurf(_meshGraph);
    fsurf.write("stretch." + std::to_string(frame) + ".gaudi");
  }

  virtual void onDraw(gg::Viewer &viewer) {

    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    _debugLines->clear();
  }

  struct {
    double dt = 0.01;
  } _params;

private:
  double _max = 0.0;
  double _min = 0.0;

  gg::DebugBufferPtr _debugLines = NULL;
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<stretch> *_meshGraph;
  m2::surf_integrator<stretch> *_integrator;
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
