//#include "nanoguiincludes.h"

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

#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/conj_grad.hpp"
#include "manifold/laplacian.hpp"
#include "manifold/diffuse.hpp"

#include "manifold/bins.hpp"
#include "manifold/conj_grad.hpp"
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

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

template <typename T> class rxdiff_space {
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
  typedef Eigen::Matrix<T, 3, 3> mat3;
  typedef Eigen::Matrix<T, 4, 4> mat4;
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

  enum class edge_index { MAXINDEX = 0 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    RXA = 2,
    RXB = 3,
    MAXINDEX = 4
  };

  enum class face_index {
    NORMAL = 0,
    CENTER = 1,
    AREA = 2,

    MAXINDEX = 3
  };
  enum class face_vertex_index { BARY = 0, MAXINDEX = 1 };
};

typedef rxdiff_space<double> rxd3;

class RxScene;
using RxScenePtr = std::shared_ptr<RxScene>;

class RxScene : public gg::Scene {

public:
  static RxScenePtr create() { return std::make_shared<RxScene>(); }

  RxScene() : gg::Scene() {
    initScene();
    // initUI();
  }

  void initScene() {

    m2::obj_loader<rxd3> load;
    m2::subdivide<rxd3> sub;
    m2::make<rxd3> mk;
    m2::convex_hull<rxd3> ch;

    m2::construct<rxd3> bevel;
    m2::affine<rxd3> mod;

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

    m2::remesh<rxd3> rem;
    rem.triangulate(_meshGraph);
    //_meshGraph = &sub.subdivide_control(*_meshGraph);

    std::cout << "--center" << std::endl;
    std::cout << "--update_all" << std::endl;
    _meshGraph->update_all();
    std::cout << "--pack" << std::endl;
    _meshGraph->pack();

    mod.centerGeometry(*_meshGraph);

    int N = 0;
    _integrator = new m2::surf_integrator<rxd3>(_meshGraph, 0.1, 3.0, 0.25);

    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXA);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXB);
    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    initRxMesh(_meshGraph);

    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);
  }

#if 1
  space3::vec3 rainbow(double d) {
    double r = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.000));
    double g = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.333));
    double b = 0.5 + 0.5 * cos(2.0 * M_PI * (d + 0.666));
    return space3::vec3(r, g, b);
  }

  space3::vec3 grey(double d) { return space3::vec3(0.5, 0.5, 0.5); }

  template <typename SPACE>
  void print_vecs(const std::vector<typename SPACE::vec3> &p0,
                  const std::vector<typename SPACE::vec3> &p1) {
    double mx =
        std::accumulate(p1.begin(), p1.end(), 0.0, [](double a, auto &c) {
          return max(a, m2::va::norm(c));
        });

    for (int i = 0; i < p0.size(); i++) {

      const auto &p = p0[i];
      const auto &a = p1[i];

      auto pa = p + 0.025 * a / mx;
      // std::cout << a.transpose() << std::endl;
      auto c = grey(m2::va::norm(a));
      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

  template <typename SPACE>
  std::vector<typename SPACE::vec3>
  calcNormals(const std::vector<typename SPACE::vec3> &positions) {
    M2_TYPEDEFS;

    std::cout << "creating vecs" << std::endl;
    std::vector<m2::face_triangle<SPACE>> triangles;
    auto &faces = _meshGraph->get_faces();
    std::vector<typename SPACE::vec3> normals;
    std::vector<typename SPACE::real> rx(faces.size());
    std::cout << "calc_normal 1" << std::endl;

    for (int i = 0; i < faces.size(); i++) {
      if (!_meshGraph->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;

      std::vector<m2::face_triangle<SPACE>> tris =
          m2::ci::get_tris<SPACE>(faces[i]);
      real rxa = 0;
      real rxb = 0;

      m2::for_each_face<SPACE>(faces[i], [&rxa, &rxb](face_vertex_ptr fv) {
        int j = fv->next()->vertex()->position_in_set();
        real l = fv->template get<real>(SPACE::face_vertex_index::BARY);

        real ra = fv->vertex()->template get<real>(SPACE::vertex_index::RXA);
        real rb = fv->vertex()->template get<real>(SPACE::vertex_index::RXB);
        rxa += l * ra;
        rxb += l * rb;
      });
      // rxa = rxa - 0.25;
      // rxb = rxb - 0.25;
      // rx[i] = rxa * rxb * rxb;
      rx[i] = 2.5 * rxb - rxa;

      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }
    std::cout << "calc_normal 2" << std::endl;

    for (int i = 0; i < triangles.size(); i++) {
      // normals.push_back(rx[i] * triangles[i].normal());
      normals.push_back(rx[i] * triangles[i].normal());
    }

    double reg = _max;
    std::cout << "calc_normal 3" << std::endl;

    std::cout << "computing harmonic avg" << std::endl;
    m2::mesh_calculator<SPACE> calc;
    // std::vector<typename SPACE::vec3> avgs =
    //    calc.harmonicNormal(_meshGraph, normals, positions, reg);
    std::vector<typename SPACE::vec3> avgs =
        calc.harmonicAvg(_meshGraph, normals, positions, reg);
    std::cout << "calc_normal 4" << std::endl;

    return avgs;
  }
#endif



  template <typename SPACE> void setRxColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto rxA = getRx(surf, SPACE::vertex_index::RXA);
    auto rxB = getRx(surf, SPACE::vertex_index::RXB);
    int i = 0;
    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = K[i];
      typename SPACE::real N = 0;
      typename SPACE::real avgA = rxA[i];
      typename SPACE::real avgB = rxB[i];

      typename SPACE::coordinate_type colorA(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorB(0.0, 1.0, 0.0);
      typename SPACE::coordinate_type color = avgA * colorA + avgB * colorB;
      // typename SPACE::coordinate_type color = (1.0 - k) * colorA + k *
      // colorB;

      v->template set<typename SPACE::coordinate_type>(
          SPACE::vertex_index::COLOR, color);
      i++;
    }
  }

  std::vector<m2::colorRGB> getRxColors(m2::surf<rxd3> *surf) {
    auto vertices = surf->get_vertices();
    std::vector<m2::colorRGB> colors(vertices.size());
    for (int i = 0; i < vertices.size(); i++) {
      auto v = vertices[i];
      rxd3::coordinate_type color =
          v->template get<typename rxd3::coordinate_type>(
              rxd3::vertex_index::COLOR);
      colors[i] = m2::colorRGB(color[0], color[1], color[2], 1.0);
    }

    return colors;
  }

  std::vector<rxd3::real> getRx(m2::surf<rxd3> *surf, rxd3::vertex_index id) {
    auto vertices = surf->get_vertices();
    std::vector<rxd3::real> rxs(surf->get_vertices().size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      rxs[i] = v->template get<typename rxd3::real>(id);
      i++;
    }
    return rxs;
  }

  void setRx(m2::surf<rxd3> *surf, std::vector<rxd3::real> rx,
             rxd3::vertex_index id) {
    int i = 0;
    for (auto v : surf->get_vertices()) {
      v->template set<typename rxd3::real>(id, rx[i]);
      i++;
    }
  }

  template <typename SPACE>
  void curvature_to_rx(m2::surf<SPACE> *surf,
                       std::vector<typename SPACE::real> &rxA,
                       std::vector<typename SPACE::real> &rxB) {
    using namespace m2;
    M2_TYPEDEFS;
    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();

    const auto [min, max] = std::minmax_element(K.begin(), K.end());
    std::for_each(K.begin(), K.end(),
                  [min, max](auto &e) { e = (e - *min) / (*max - *min); });

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0.1, 1);
    int i = 0;

    auto vertices = surf->get_vertices();
    coordinate_array coords;
    for (auto v : vertices) {
      if (dist(e2) > 0.95 && K[i] < 0.1 && rxA[i] > 0.95) {
        coords.push_back(ci::get_coordinate<SPACE>(v));
      }
      i++;
    }

    i = 0;
    int N = coords.size();
    for (auto cj : coords) {
      i = 0;
      for (auto v : vertices) {
        coordinate_type ci = ci::get_coordinate<SPACE>(v);
        coordinate_type dc = cj - ci;
        real dist = va::norm(dc);
        // std::cout << ci.transpose() << " " << cj.transpose() << std::endl;
        real C = 1.0 * _max;
        T d2 = dist * dist * dist;
        T l2 = C * C * C;

        T kappa = exp(-0.5 * d2 / l2);
        rxB[i] += (kappa);
        rxA[i] = (1.0 - kappa) * rxA[i];
        i++;
      }
    }
    /*
    double mx = 0.0;
    for (int i = 0; i < rxB.size(); i++) {
      mx = std::max(mx, rxA[i]);
      mx = std::max(mx, rxB[i]);
    }

    for (int i = 0; i < rxB.size(); i++) {
      rxB[i] /= mx;
      rxA[i] += (1.0 - rxB[i]);
    }
    */
    /*
    mx = 0.0;
    for (int i = 0; i < rxB.size(); i++) {
      mx = std::max(mx, rxA[i]);
      mx = std::max(mx, rxB[i]);
    }

    for (auto k : K) {
      if(k< 1e-6) continue;
      //std::cout << k << std::endl;
      real k2 = k * k + 0.1;
      //k2 = std::min(k2, real(10.0));
      //k2 = std::max(k2, real(0.1));

      rxB[i] += rxA[i] *  C / k2 * dist(e2); //if low curvature and lots of A
    add some B rxA[i] += rxB[i] *  C / k2 * dist(e2); //if high curvature and
    lots of B add some A i++;
    }
    */
  }

  template <typename SPACE> void initRxMesh(m2::surf<SPACE> *surf) {
    using namespace m2;
    M2_TYPEDEFS;

    auto vertices = surf->get_vertices();

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    std::vector<typename SPACE::real> rxA(vertices.size());
    std::vector<typename SPACE::real> rxB(vertices.size());

    int i = 0;
    coordinate_array coords;
    for (auto v : vertices) {
      rxA[i] = 1.0;
      rxB[i] = 0.0;
      i++;
    }

    setRx(surf, rxA, SPACE::vertex_index::RXA);
    setRx(surf, rxB, SPACE::vertex_index::RXB);
  }

  void reset() { this->initRxMesh(_meshGraph); }

  virtual void onAnimate(int frame) {

    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<rxd3>(_meshGraph) << std::endl;

    _meshGraph->update_all();
    _meshGraph->reset_flags();
    _meshGraph->pack();

    _integrator->integrate();

    // gg::fillBuffer(_meshGraph, _obj);
    std::cout << "get rx " << std::endl;
    std::vector<rxd3::real> rxA = getRx(_meshGraph, rxd3::vertex_index::RXA);
    std::vector<rxd3::real> rxB = getRx(_meshGraph, rxd3::vertex_index::RXB);
    std::vector<rxd3::real> fA(rxA);
    std::vector<rxd3::real> fB(rxB);

    std::cout << "curve to rx " << std::endl;

    curvature_to_rx(_meshGraph, rxA, rxB);

    double An = 0.0, Bn = 0.0;
    double Du = rxParams.Du;
    double Dv = 0.4875 * Du;
    // double Dv = rxParams.Dv;

    double F = rxParams.F;
    double k = rxParams.k;
    double dt = rxParams.dt;

    std::cout << "Du: " << Du << std::endl;
    std::cout << "Dv: " << Dv << std::endl;
    std::cout << " F: " << F << std::endl;
    std::cout << " k: " << k << std::endl;
    std::cout << " dt: " << dt << std::endl;
#if 1
    for (int i = 0; i < rxA.size(); i++) {
      double u = rxA[i];
      double v = rxB[i];

      double rx = u * v * v;
      double up = (-rx + F * (1.0 - u));
      double vp = (rx - (F + k) * v);

      fA[i] = up;
      fB[i] = vp;
    }

    m2::diffuse<rxd3> diff;
    rxA = diff.second_order(_meshGraph, rxA, fA, dt, Du);
    rxB = diff.second_order(_meshGraph, rxB, fB, dt, Dv);

    for (int i = 0; i < rxA.size(); i++) {
      rxA[i] = std::clamp(rxA[i], 0.0, 1.0);
      rxB[i] = std::clamp(rxB[i], 0.0, 1.0);
      An += rxA[i];
      Bn += rxB[i];
    }

    std::cout << "rx stage norm: " << An / double(rxA.size()) << " "
              << Bn / double(rxB.size()) << std::endl;
    std::cout << "area: " << m2::ci::area<rxd3>(_meshGraph) << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<rxd3>(_meshGraph) << std::endl;

    setRx(_meshGraph, rxA, rxd3::vertex_index::RXA);
    setRx(_meshGraph, rxB, rxd3::vertex_index::RXB);

    // mod.centerGeometry(*_meshGraph);
#endif
    setRxColor(_meshGraph);
    auto colors = getRxColors(_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);
#if 1

    std::vector<rxd3::vec3> positions =
        m2::ci::get_coordinates<rxd3>(_meshGraph);
    std::vector<rxd3::vec3> normals = calcNormals<rxd3>(positions);
    std::cout << "normal length: " << normals.size() << std::endl;
    std::cout << "rx length: " << rxA.size() << std::endl;

    std::cout << "print vecs" << std::endl;
    print_vecs<rxd3>(positions, normals);
    std::cout << "rendering debug" << std::endl;
    _debugLines->renderLines();

    double Nn = 0, NRxn = 0;
    for (int i = 0; i < normals.size(); i++) {
      Nn += m2::va::norm<rxd3::real>(normals[i]);
      NRxn += m2::va::norm<rxd3::real>(normals[i]);
      positions[i] += 0.005 * normals[i];
    }

    std::cout << "pos stage norm: " << Nn / double(rxA.size()) << " "
              << NRxn / double(rxB.size()) << std::endl;

    m2::ci::set_coordinates<rxd3>(positions, _meshGraph);

#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "rxdiff." << frame << ".obj";
    m2::write_obj<rxd3>(*_meshGraph, ss.str());
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
    double dt = 2.0;
    double Du = 1.0e-4, Dv = 0.025;
    // double k = 0.061, F = 0.070; //pretty good one
    double k = 0.061, F = 0.068;
    // double k = 0.06, F = 0.082;
    // double k = 0.06132758, F = 0.037;
    // double k = 0.059, F = 0.03;
    // double k = 0.056, F = 0.098; //interesting?
    // double k = 0.0613, F = 0.06; // uskatish
    // double, k = 0.0550, F = 0.1020 //interesting?
    // double k = 0.0628, F = 0.0567;
    // double k = 0.061, F = 0.42;
    // double k = 0.045, F = 0.01; // waves

  } rxParams;

private:
  double _max = 0.0;
  double _min = 0.0;

  gg::DebugBufferPtr _debugLines = NULL;
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<rxd3> *_meshGraph;
  m2::surf_integrator<rxd3> *_integrator;
};

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

class RxApp;
using RxAppPtr = std::shared_ptr<RxApp>;

class RxApp : public gg::SimpleApp {
public:
  static RxAppPtr create(std::string file) {
    return std::make_shared<RxApp>(file);
  }

  typedef double Real;

  RxApp(std::string file) : gg::SimpleApp() {
    this->setScene(scene = RxScene::create());
    this->initUI();
  }

  void initUI() {
    using namespace nanogui;
    int w = 256;
    Window *window = new Window(this, "Button demo");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout());
    float Kd = 4.0;
    float Kf = 10.0;

    Slider *slider0 = new Slider(window);
    slider0->setValue(pow(this->scene->rxParams.Du, 1.0 / Kd));
    slider0->setFixedWidth(w);

    Slider *slider2 = new Slider(window);
    slider2->setValue(this->scene->rxParams.F * Kf);
    slider2->setFixedWidth(w);

    Slider *slider3 = new Slider(window);
    slider3->setValue(this->scene->rxParams.k * Kf);
    slider3->setFixedWidth(w);

    Slider *slider4 = new Slider(window);
    slider4->setValue(this->scene->rxParams.dt / 2.0);
    slider4->setFixedWidth(w);

    slider0->setCallback(
        [this, Kd](float value) { this->scene->rxParams.Du = pow(value, Kd); });

    slider2->setCallback(
        [this, Kf](float value) { this->scene->rxParams.F = value / Kf; });

    slider3->setCallback(
        [this, Kf](float value) { this->scene->rxParams.k = value / Kf; });

    slider4->setCallback(
        [this, Kf](float value) { this->scene->rxParams.dt = 2.0 * value; });
    Button *b = window->add<Button>("reset");
    b->setFixedSize(Vector2i(w, 22));
    b->setCallback([this] { this->scene->reset(); });

    performLayout();
    // window->center();
  }

  ~RxApp() {}

  RxScenePtr scene;
};

int main(int argc, char *argv[]) {
  try {
    cout << "You have entered " << argc << " arguments:"
         << "\n";

    for (int i = 0; i < argc; ++i)
      cout << argv[i] << "\n";

    nanogui::init();

    RxAppPtr app = RxApp::create(std::string(argv[0]));

    // app->setScene(RxScene::create());
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
