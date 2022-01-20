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

#include "manifold/harmonic_integrators.hpp"
#include "manifold/moving_mesh.hpp"
#include "manifold/vec_addendum.h"

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

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

    _meshGraph = &sub.subdivide_control(*_meshGraph);
    _meshGraph = &sub.subdivide_control(*_meshGraph);
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
    std::cout << "--build sheet" << std::endl;
    _vortex = new m2::vortex_sheet<rxd3>(_meshGraph);
    initRxMesh(_meshGraph);

    int N = 0;
    double l = m2::ci::mean_length<rxd3>(_meshGraph);
    double max = l;
    double min = max / 8.0;

    std::cout << "avg length: " << l << std::endl;

    _vortex->minLength = min;
    _vortex->minCollapseLength = min;
    _vortex->edgeJoinThresh = 0.25 * min;
    _vortex->regLength = max;

    mod.centerGeometry(*_meshGraph);
    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    _debugLines = gg::DebugBuffer::create();
    _debugLines->init();
    mSceneObjects.push_back(_debugLines);
  }

  template <typename SPACE> void cacheBary(m2::surf<rxd3> *surf) {
    M2_TYPEDEFS;
    for (auto f : surf->get_faces()) {
      typename SPACE::coordinate_type c = m2::ci::center<SPACE>(f);
      typename SPACE::coordinate_type l = m2::ci::point_to_bary<SPACE>(f, c);
      int i = 0;
      m2::for_each_face<SPACE>(f, [l, i](face_vertex_ptr fv) mutable {
        fv->template set<typename SPACE::real>(SPACE::face_vertex_index::BARY,
                                               l[i]);
        i++;
      });
    }
  }

  template <typename SPACE> void setRxColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto rxA = getRx(surf, SPACE::vertex_index::RXA);
    auto rxB = getRx(surf, SPACE::vertex_index::RXB);
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real N = 0;
      typename SPACE::real avgA = rxA[i];
      typename SPACE::real avgB = rxB[i];

      typename SPACE::coordinate_type colorA(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorB(0.0, 1.0, 0.0);
      typename SPACE::coordinate_type color = avgA * colorA + avgB * colorB;

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

  template <typename SPACE, typename TYPE>
  std::vector<TYPE> matLaplace(m2::surf<SPACE> *surf, std::vector<TYPE> U) {
    M2_TYPEDEFS;
    std::vector<TYPE> L(U.size());
    int i = 0;
    real aAvg = 0;

    for (auto v : surf->get_vertices()) {
      real ui = U[i];
      real u = 0.0;

      m2::for_each_vertex<SPACE>(
          v, [&aAvg, &u, &ui, &U](face_vertex_ptr fv) mutable {
            face_vertex_ptr fvp = fv->vprev()->next();
            face_vertex_ptr fvn = fv->vnext()->next();

            int j = fv->next()->vertex()->position_in_set();
            real uj = U[j];

            real cotp = m2::ci::cotan<SPACE>(fvp);
            real cotn = m2::ci::cotan<SPACE>(fvn);
            real area = m2::ci::area<SPACE>(fv->face());
            real l = fv->template get<real>(SPACE::face_vertex_index::BARY);
            real denom = 1.0 / (2.0 * l * area);

            aAvg += denom;
            real K = (cotp + cotn) * denom;
            u += K * (uj - ui);
          });
      L[i] = u;

      i++;
    }

    aAvg /= real(U.size());

    for (auto &l : L) {
      l /= aAvg;
    }

    return L;
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  getLaplaceScalar(m2::surf<SPACE> *surf, std::vector<typename SPACE::real> U) {
    return matLaplace<SPACE, typename SPACE::real>(surf, U);
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  diffMult(m2::surf<SPACE> *surf, std::vector<typename SPACE::real> X,
           typename SPACE::real dt, typename SPACE::real C) {
    std::vector<typename SPACE::real> Y = X;
    std::vector<typename SPACE::real> MX =
        matLaplace<SPACE, typename SPACE::real>(surf, X);

    int i = 0;
    for (int i = 0; i < Y.size(); i++) {
      Y[i] -= dt * C * MX[i];
    }
    return Y;
  }

  template <typename SPACE>
  typename SPACE::real dot(const std::vector<typename SPACE::real> &A,
                           const std::vector<typename SPACE::real> &B) {
    M2_TYPEDEFS;
    real v = 0.0;
    for (int i = 0; i < A.size(); i++) {
      v += A[i] * B[i];
    }
    return v;
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  addScaledVec(const std::vector<typename SPACE::real> &X0,
               const std::vector<typename SPACE::real> &V,
               typename SPACE::real C) {
    M2_TYPEDEFS;
    std::vector<real> X(X0);
    for (int i = 0; i < X.size(); i++) {
      X[i] = X0[i] + C * V[i];
    }
    return X;
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  diffuse(m2::surf<SPACE> *surf, std::vector<typename SPACE::real> B,
          typename SPACE::real dt, typename SPACE::real C) {
    M2_TYPEDEFS;
    /*
     x =          initial guess for solution of Ax=b
     r = b - Ax = residual, to be made small
     p = r      = initial "search direction"

     do while ( new_r , new_r ) not small
        v = Ap                        ... matrix-vector multiply
        a = ( r , r ) / ( p , v )     ... dot product
        x = x + a*p                   ... updated approximate solution
        r_new = r - a*v               ... update the residual
        g = ( r_new , r_new ) / ( r , r )
        p = r_new + g*p               ... update search direction
        r = r_new
     enddo
    */

    // x =          initial guess for solution of Ax=b
    std::vector<real> X = B;

    // std::vector<real> MU = diffMult(surf, U0, dt, C);
    std::vector<real> MX = diffMult(surf, X, dt, C);

    // r = b - Ax = residual, to be made small
    std::vector<real> R = addScaledVec<SPACE>(B, MX, -1);
    // p = r = initial "search direction"
    std::vector<real> P = R;
    real sig1 = dot<SPACE>(R, R);
    int ii = 0;
    while (ii < 10 && sig1 > 1e-3) {

      // v = Ap                        ... matrix-vector multiply
      std::vector<real> V = diffMult(surf, P, dt, C);

      // a = ( r , r ) / ( p , v )     ... dot product
      real a = sig1 / dot<SPACE>(P, V);
      // std::cout << sig1 << std::endl;
      // x = x + a*p                   ... updated approximate solution
      X = addScaledVec<SPACE>(X, P, a);
      // r_new = r - a * v... update the residual
      std::vector<real> Rn = addScaledVec<SPACE>(R, V, -a);
      // g = ( r_new , r_new ) / ( r , r )
      real sig0 = sig1;
      sig1 = dot<SPACE>(Rn, Rn);
      real g = sig1 / sig0;

      // p = r_new + g *p... update search direction
      P = addScaledVec<SPACE>(Rn, P, g);

      // r = r_new
      R = Rn;
      ii++;
    }

    return X;
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

  void initRxMesh(m2::surf<rxd3> *surf) {
    auto vertices = surf->get_vertices();

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    std::vector<typename rxd3::real> rxA(vertices.size());
    std::vector<typename rxd3::real> rxB(vertices.size());
    int i = 0;
    for (auto v : vertices) {
      // std::cout << dist(e2);
      /*
      if (dist(e2) > 0.85)
        rxA[i] = 10.0;
      else
        rxA[i] = 0.0;
      rxB[i] = 0.5;
      */
      /*
      if (dist(e2) > 0.5)
        rxA[i] = 1.0;
      else
        rxA[i] = 0.0;

      if (dist(e2) > 0.5)
        rxB[i] = 1.0;
      else
        rxB[i] = 0.0;
      */
#if 0
      typename rxd3::real t = dist(e2);
      t = t * t;
      typename rxd3::real tt = max(t, 1.0 - t);
      rxA[i] = 1.0 - tt;
      rxB[i] = tt;
#else
      if (dist(e2) > 0.99) {
        typename rxd3::real t = dist(e2);
        rxA[i] = 0.0;
        rxB[i] = 1.0;
      } else {
        rxA[i] = 1.0;
        rxB[i] = 0.0;
      }
#endif
      // rxA[i] = dist(e2);
      // rxB[i] = dist(e2);
      i++;
    }

    setRx(surf, rxA, rxd3::vertex_index::RXA);
    setRx(surf, rxB, rxd3::vertex_index::RXB);
  }

  void reset() { this->initRxMesh(_meshGraph); }

  virtual void onAnimate() {
    /*
    std::cout << "====== " << std::endl;

    std::cout << " vertices: " << _meshGraph->get_vertices().size()
              << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    */
    m2::affine<rxd3> mod;

    cacheBary<rxd3>(_meshGraph);
    std::vector<rxd3::real> rxA = getRx(_meshGraph, rxd3::vertex_index::RXA);
    std::vector<rxd3::real> rxB = getRx(_meshGraph, rxd3::vertex_index::RXB);
    // std::vector<rxd3::real> rxAd = getLaplaceScalar(_meshGraph, rxA);
    // std::vector<rxd3::real> rxBd = getLaplaceScalar(_meshGraph, rxB);

    double An = 0.0, Bn = 0.0;
    double Du = rxParams.Du;
    double Dv = 0.5 * Du;
    // double Dv = rxParams.Dv;

    double F = rxParams.F;
    double k = rxParams.k;
    double dt = rxParams.dt;

    std::cout << "Du: " << Du << std::endl;
    std::cout << "Dv: " << Dv << std::endl;
    std::cout << " F: " << F << std::endl;
    std::cout << " k: " << k << std::endl;
    std::cout << " dt: " << dt << std::endl;

    for (int i = 0; i < rxA.size(); i++) {

      double u = rxA[i];
      double v = rxB[i];

      // double Lu = rxAd[i];
      // double Lv = rxBd[i];
      double rx = u * v * v;
      // double up  = u + dt * Du * Lu;
      // double vp  = v + dt * Dv * Lv;
      double up = u + dt * (-rx + F * (1.0 - u));
      double vp = v + dt * (+rx - (F + k) * v);
      // double up = u + dt * (-rx + F * (1.0 - u));
      // double vp = v + dt * ( rx - (F + k) * v);
      rxA[i] = std::clamp(up, 0.0, 1.0);
      rxB[i] = std::clamp(vp, 0.0, 1.0);

      An += rxA[i];
      Bn += rxB[i];
    }

    rxA = diffuse(_meshGraph, rxA, dt, Du);
    rxB = diffuse(_meshGraph, rxB, dt, Dv);

    std::cout << "stage norm: " << An / double(rxA.size()) << " "
              << Bn / double(rxB.size()) << std::endl;
    std::cout << "area: " << m2::ci::area<rxd3>(_meshGraph) << std::endl;

    setRx(_meshGraph, rxA, rxd3::vertex_index::RXA);
    setRx(_meshGraph, rxB, rxd3::vertex_index::RXB);

    setRxColor(_meshGraph);
    auto colors = getRxColors(_meshGraph);
    // mod.centerGeometry(*_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);
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
    double dt = 0.1;
    double Du = 0.2, Dv = 0.1;
    //double k = 0.059, F = 0.037;
    //double k = 0.0542, F = 0.0223;
    //double k = 0.0628, F = 0.0567;
    double k = 0.045, F = 0.01; //waves

  } rxParams;

private:
  gg::DebugBufferPtr _debugLines = NULL;
  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<rxd3> *_meshGraph;
  m2::vortex_sheet<rxd3> *_vortex;
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
    float Kd = 1.0;
    float Kf = 10.0;

    Slider *slider0 = new Slider(window);
    slider0->setValue(this->scene->rxParams.Du * Kd);
    slider0->setFixedWidth(w);

    Slider *slider2 = new Slider(window);
    slider2->setValue(this->scene->rxParams.F * Kf);
    slider2->setFixedWidth(w);

    Slider *slider3 = new Slider(window);
    slider3->setValue(this->scene->rxParams.k * Kf);
    slider3->setFixedWidth(w);

    Slider *slider4 = new Slider(window);
    slider4->setValue(this->scene->rxParams.dt);
    slider4->setFixedWidth(w);

    slider0->setCallback(
        [this, Kd](float value) { this->scene->rxParams.Du = value / Kd; });

    slider2->setCallback(
        [this, Kf](float value) { this->scene->rxParams.F = value / Kf; });

    slider3->setCallback(
        [this, Kf](float value) { this->scene->rxParams.k = value / Kf; });

    slider4->setCallback(
        [this, Kf](float value) { this->scene->rxParams.dt = 10.0 * value; });
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
