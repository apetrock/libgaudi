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

  enum class edge_index { MAXINDEX = 0 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    RXA = 2,
    RXB = 3,
    RXC = 4,
    RXD = 5,
    RXE = 6,
    RXF = 7,
    SMOOTH = 8,
    MAXINDEX = 9
  };

  enum class face_index {
    NORMAL = 0,
    CENTER = 1,
    AREA = 2,

    MAXINDEX = 3
  };

  static storage_type get_type(face_vertex_index idx) {
    switch (idx) {
    case face_vertex_index::BARY:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    };
  }
  static storage_type get_type(edge_index idx) { return storage_type::SIZE; }

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

    case vertex_index::RXA:
      return storage_type::REAL;

    case vertex_index::RXB:
      return storage_type::REAL;

    case vertex_index::RXC:
      return storage_type::REAL;

    case vertex_index::RXD:
      return storage_type::REAL;

    case vertex_index::RXE:
      return storage_type::REAL;

    case vertex_index::RXF:
      return storage_type::REAL;

    case vertex_index::SMOOTH:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }
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
    ///std::string start_frame = "rxd.100.gaudi";
    std::string start_frame = "";
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

      m2::remesh<rxd3> rem;
      rem.triangulate(_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      std::cout << "--center" << std::endl;
      std::cout << "--update_all" << std::endl;
      _meshGraph->update_all();
      std::cout << "--pack" << std::endl;
      _meshGraph->pack();

      mod.centerGeometry(*_meshGraph);
      initRxMesh(_meshGraph);
    }

    int N = 0;
    _integrator = new m2::surf_integrator<rxd3>(_meshGraph, 0.125, 3.0, 0.5);
    //_integrator = new m2::surf_integrator<rxd3>(_meshGraph, 1.0, 3.0, 0.5);

    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXA);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXB);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXC);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXD);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXE);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::RXF);
    _integrator->add_default_vertex_policy<typename rxd3::real>(
        rxd3::vertex_index::SMOOTH);

    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    //

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
  space3::vec3 red(double d) { return space3::vec3(1.0, 0.0, 0.0); }
  space3::vec3 green(double d) { return space3::vec3(0.0, 1.0, 0.0); }
  space3::vec3 blue(double d) { return space3::vec3(0.0, 0.0, 1.0); }

  template <typename SPACE>
  void print_vecs(const std::vector<typename SPACE::vec3> &p0,
                  const std::vector<typename SPACE::vec3> &p1, int col = 0) {
    double mx =
        std::accumulate(p1.begin(), p1.end(), 0.0, [](double a, auto &c) {
          return max(a, m2::va::norm(c));
        });

    for (int i = 0; i < p0.size(); i++) {

      const auto &p = p0[i];
      const auto &a = p1[i];

      auto pa = p + 0.1 * a / mx;
      // std::cout << a.transpose() << std::endl;
      auto c = grey(m2::va::norm(a));
      if (col == 1)
        c = red(m2::va::norm(a));
      if (col == 2)
        c = green(m2::va::norm(a));
      if (col == 3)
        c = blue(m2::va::norm(a));

      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

  template <typename SPACE> std::vector<typename SPACE::vec3> calcNormals() {
    M2_TYPEDEFS;

    std::vector<rxd3::vec3> positions =
        m2::ci::get_coordinates<rxd3>(_meshGraph);

    auto rxA = getRx(_meshGraph, SPACE::vertex_index::RXA);
    auto rxB = getRx(_meshGraph, SPACE::vertex_index::RXB);
    auto rxC = getRx(_meshGraph, SPACE::vertex_index::RXC);

    auto rxD = getRx(_meshGraph, SPACE::vertex_index::RXD);
    auto rxE = getRx(_meshGraph, SPACE::vertex_index::RXE);
    auto rxF = getRx(_meshGraph, SPACE::vertex_index::RXF);

    auto rxV = std::vector<typename SPACE::real>(rxA.size());
    int i = 0;
    for (auto rxa : rxA) {

      rxV[i] =                                                  //
          2.0 * rxA[i] - 1.0 * rxB[i] +                         //
          2.0 * sqrt(2.0) * rxC[i] - 1.0 * sqrt(2.0) * rxD[i] + //
          2.0 * sqrt(3.0) * rxE[i] - 1.0 * sqrt(3.0) * rxF[i];
      i++;
    }

    std::cout << "creating vecs" << std::endl;
    std::vector<m2::face_triangle<SPACE>> triangles;
    auto &faces = _meshGraph->get_faces();

    for (int i = 0; i < faces.size(); i++) {
      if (!_meshGraph->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<m2::face_triangle<SPACE>> tris =
          m2::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }
    std::vector<typename SPACE::real> rxFa =
        m2::ci::verts_to_faces<SPACE, typename SPACE::real>(rxV, _meshGraph);

    std::vector<typename SPACE::vec3> normF;
    std::vector<typename SPACE::vec3> normV(positions.size());
    std::cout << "calc_normal 2" << std::endl;
    for (int i = 0; i < triangles.size(); i++) {
      normF.push_back(rxFa[i] * triangles[i].normal());
    }

    double reg = 1.0 * _max;
    std::cout << "calc_normal 3" << std::endl;

    std::cout << "computing harmonic avg" << std::endl;
    m2::mesh_calculator<SPACE> calc;

    normV = calc.harmonicAvg(_meshGraph, normF, positions, reg);
    return normV;
    std::vector<typename SPACE::real> divV =
        calc.template divergence<typename SPACE::real, typename SPACE::vec3>(
            _meshGraph, normF, positions, reg);

    std::vector<typename SPACE::real> divF =
        m2::ci::verts_to_faces<SPACE, typename SPACE::real>(divV, _meshGraph);

    std::vector<typename SPACE::vec3> grads =
        calc.template gradient<typename SPACE::vec3, typename SPACE::real>(
            _meshGraph, divF, positions, reg);

    for (int i = 0; i < normV.size(); i++) {
      normV[i] = normV[i] - grads[i];
    }
    return normV;
  }
#endif

  template <typename SPACE> void setRxColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto rxA = getRx(surf, SPACE::vertex_index::RXA);
    auto rxB = getRx(surf, SPACE::vertex_index::RXB);
    auto rxC = getRx(surf, SPACE::vertex_index::RXC);
    auto rxD = getRx(surf, SPACE::vertex_index::RXD);
    auto rxE = getRx(surf, SPACE::vertex_index::RXE);
    auto rxF = getRx(surf, SPACE::vertex_index::RXF);
    auto smooth = getRx(surf, SPACE::vertex_index::SMOOTH);

    int i = 0;
    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = K[i];
      typename SPACE::real N = 0;
      typename SPACE::real avgA = rxA[i];
      typename SPACE::real avgB = rxB[i];
      typename SPACE::real avgC = rxC[i];
      typename SPACE::real avgD = rxD[i];
      typename SPACE::real avgE = rxE[i];
      typename SPACE::real avgF = rxF[i];
      typename SPACE::real s = smooth[i];

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorA(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorB(1.0, 0.0, 1.0);
      typename SPACE::coordinate_type colorC(0.0, 0.0, 1.0);
      typename SPACE::coordinate_type colorD(0.0, 1.0, 1.0);
      typename SPACE::coordinate_type colorE(0.0, 1.0, 0.0);
      typename SPACE::coordinate_type colorF(1.0, 1.0, 0.0);

      typename SPACE::coordinate_type color = //
          0.5 * avgA * colorA +               //
          0.1 * avgB * colorD +               //
          0.5 * avgC * colorB +               //
          0.1 * avgD * colorE +               //
          0.5 * avgE * colorC +               //
          0.1 * avgF * colorF;

      // typename SPACE::coordinate_type color = (1.0 - k) * colorA + k *
      // colorB;
      typename SPACE::coordinate_type mx = m2::va::mix(s, colorS, color);

      v->template set<typename SPACE::coordinate_type>(
          SPACE::vertex_index::COLOR, mx);
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
                       std::vector<typename SPACE::real> &rxB,
                       std::vector<typename SPACE::real> &rxC,
                       std::vector<typename SPACE::real> &rxD,
                       std::vector<typename SPACE::real> &rxE,
                       std::vector<typename SPACE::real> &rxF) {
    using namespace m2;
    M2_TYPEDEFS;
    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();

    const auto [min, max] = std::minmax_element(K.begin(), K.end());
    std::for_each(K.begin(), K.end(),
                  [min, max](auto &e) { e = (e - *min) / (*max - *min); });

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    std::uniform_int_distribution<int> idist(0, 2);

    int i = 0;

    auto vertices = surf->get_vertices();
    coordinate_array coords;

    for (auto v : vertices) {
      if (dist(e2) < 5e-5 && K[i] < 0.1 && rxA[i] < 0.025) {
        coords.push_back(ci::get_coordinate<SPACE>(v));
      }
      i++;
    }
    std::cout << "splatting: " << coords.size() << " " << vertices.size() << " "
              << double(coords.size()) / double(vertices.size()) << std::endl;
#if 1
    i = 0;
    int N = coords.size();
    std::uniform_real_distribution<> dist3(0.5, 1.5);
    for (auto cj : coords) {
      i = 0;
      real r = dist3(e2) * _max;
      bool flop = dist(e2) > 0.5;
      int card = idist(e2);
      double F = flop ? 1.0 : 0.0;
      for (auto v : vertices) {
        coordinate_type ci = ci::get_coordinate<SPACE>(v);
        coordinate_type dc = cj - ci;
        real dist = va::norm(dc);
        // std::cout << ci.transpose() << " " << cj.transpose() << std::endl;
        real C = r;
        T d2 = dist * dist * dist;
        T l2 = C * C * C;

        T kappa = exp(-0.5 * d2 / l2);

        if (card == 0) {
          rxB[i] = (1.0 - kappa) * rxB[i];
          rxA[i] += kappa;
        } else if (card == 1) {
          rxD[i] = (1.0 - kappa) * rxD[i];
          rxC[i] += kappa;
        } else {
          rxF[i] = (1.0 - kappa) * rxF[i];
          rxE[i] += kappa;
        }

        i++;
      }
    }
#endif
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
    std::vector<typename SPACE::real> rxC(vertices.size());
    std::vector<typename SPACE::real> rxD(vertices.size());
    std::vector<typename SPACE::real> rxE(vertices.size());
    std::vector<typename SPACE::real> rxF(vertices.size());

    int i = 0;
    coordinate_array coords;
    for (auto v : vertices) {
      rxA[i] = 0.0;
      rxC[i] = 0.0;
      rxE[i] = 0.0;

      if (dist(e2) > 0.75)
        rxA[i] = 1.0;
      rxB[i] = 1.0 - rxA[i];

      if (dist(e2) > 0.75)
        rxC[i] = 1.0;
      rxD[i] = 1.0 - rxC[i];

      if (dist(e2) > 0.75)
        rxE[i] = 1.0;
      rxF[i] = 1.0 - rxE[i];

      i++;
    }

    setRx(surf, rxA, SPACE::vertex_index::RXA);
    setRx(surf, rxB, SPACE::vertex_index::RXB);
    setRx(surf, rxC, SPACE::vertex_index::RXC);
    setRx(surf, rxD, SPACE::vertex_index::RXD);
    setRx(surf, rxE, SPACE::vertex_index::RXE);
    setRx(surf, rxF, SPACE::vertex_index::RXF);
  }

  void reset() { this->initRxMesh(_meshGraph); }

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

  virtual void calcGreyScott(std::vector<rxd3::real> &rxA, //
                             std::vector<rxd3::real> &rxB, //
                             std::vector<rxd3::real> &rxC, //
                             std::vector<rxd3::real> &rxD, //
                             std::vector<rxd3::real> &rxE, //
                             std::vector<rxd3::real> &rxF, //
                             std::vector<rxd3::real> &fA,  //
                             std::vector<rxd3::real> &fB,  //
                             std::vector<rxd3::real> &fC,
                             std::vector<rxd3::real> &fD, //
                             std::vector<rxd3::real> &fE, //
                             std::vector<rxd3::real> &fF) {
    rxParams.Dc = 1e-4;
    rxParams.Da = 0.5 * rxParams.Dc;
    rxParams.Db = 0.5 * rxParams.Da;

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    auto v2 = [](double v, double f) {
      double v2 = v * v;
      return v2 / (1.0 + f * v2);
    };

    double //

        F = 0.075,
        k = 0.056; //

    for (int i = 0; i < rxA.size(); i++) {

      double a = rxA[i];
      double b = rxB[i];
      double c = rxC[i];
      double a2 = v2(a, 0.0);
      // gray-scott
      double ap = a2 * c - (F + k) * a;
      double bp = 0.0;
      double cp = -a2 * c + F * (1 - c);

      fA[i] = ap;
      fB[i] = bp;
      fC[i] = cp;
      fD[i] = 0.0;
      fE[i] = 0.0;
      fF[i] = 0.0;
    }
  }

  virtual void calc3ScaleGreyScott(std::vector<rxd3::real> &rxA, //
                                   std::vector<rxd3::real> &rxB, //
                                   std::vector<rxd3::real> &rxC, //
                                   std::vector<rxd3::real> &rxD, //
                                   std::vector<rxd3::real> &rxE, //
                                   std::vector<rxd3::real> &rxF, //
                                   std::vector<rxd3::real> &fA,  //
                                   std::vector<rxd3::real> &fB,  //
                                   std::vector<rxd3::real> &fC,
                                   std::vector<rxd3::real> &fD, //
                                   std::vector<rxd3::real> &fE, //
                                   std::vector<rxd3::real> &fF) {
    double base = 5e-4;
    rxParams.Df = base;
    rxParams.De = 0.56 * rxParams.Df;

    rxParams.Dd = 0.56 * rxParams.De;
    rxParams.Dc = 0.56 * rxParams.Dd;

    rxParams.Db = 0.56 * rxParams.Dc;
    rxParams.Da = 0.56 * rxParams.Db;

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist11(0, 0.025);
    std::uniform_real_distribution<> dist01(0, 0.5);

    auto rxAB = [](double &ap, double &bp,           //
                   const double &a, const double &b, //
                   const double &ra, const double &rb) {
      double a2b = a * a * b;
      ap = a2b - ra * a;
      bp = -a2b + rb * (1 - b);
    };

    double //

        F0 = 0.06,
        k0 = 0.0633, //

        F1 = 0.065,
        k1 = 0.0633, //

        F2 = 0.067,
        k2 = 0.0633; //

    for (int i = 0; i < rxA.size(); i++) {

      double a = rxA[i];
      double b = rxB[i];
      double c = rxC[i];
      double d = rxD[i];
      double e = rxE[i];
      double f = rxF[i];

      double                      //
          Cac = 1.0 * dist11(e2), //
          Cce = 1.0 * dist11(e2), //
          Cfa = 1.0 * dist11(e2);
      double                      //
          Cbd = 1.0 * dist01(e2), //
          Cdf = 1.0 * dist01(e2), //
          Cfb = 1.0 * dist01(e2);

      double ra = F0 + k0 + Cac * c;
      double rb = F0 / (1.0 - Cbd * d);
      double rc = F1 + k1 + Cce * e;
      double rd = F1 / (1.0 - Cdf * f);
      double re = F2 + k2 + Cfa * a;
      double rf = F2 / (1.0 - Cfb * b);
      double ap, bp, cp, dp, ep, fp;
      rxAB(ap, bp, a, b, ra, rb);
      rxAB(cp, dp, c, d, rc, rd);
      rxAB(ep, fp, e, f, re, rf);
      fA[i] = ap;
      fB[i] = bp;
      fC[i] = cp;
      fD[i] = dp;
      fE[i] = ep;
      fF[i] = fp;
      /*
      fA[i] = ap;
      fB[i] = bp;
      fC[i] = cp;
      fD[i] = dp;
      fE[i] = ep;
      fF[i] = fp;
      */
    }
  }

  virtual void calcActivatorDepletion(std::vector<rxd3::real> &rxA, //
                                      std::vector<rxd3::real> &rxB, //
                                      std::vector<rxd3::real> &rxC, //
                                      std::vector<rxd3::real> &rxD, //
                                      std::vector<rxd3::real> &rxE, //
                                      std::vector<rxd3::real> &rxF, //
                                      std::vector<rxd3::real> &fA,  //
                                      std::vector<rxd3::real> &fB,  //
                                      std::vector<rxd3::real> &fC,
                                      std::vector<rxd3::real> &fD, //
                                      std::vector<rxd3::real> &fE, //
                                      std::vector<rxd3::real> &fF) {
    double base = 1e-4;
    rxParams.Db = 0.5 * base;
    rxParams.Da = 0.5 * rxParams.Db;
    rxParams.Dc = 0.1 * base;
    auto v2 = [](double v, double sa, double ba) {
      double v2 = v * v;
      return (v2 + ba) / (1.0 + sa * v2);
    };

    double                   //
        ba = 0.00,           //
        ra = 0.0753 + 0.075, //
        bb = 0.075,          //
        rb = bb,             //
        rc = 0.01,           //
        sa = 0.00, sb = 0.1, sc = 0.05, s = sb + sc;

    for (int i = 0; i < rxA.size(); i++) {

      double a = rxA[i];
      double b = rxB[i];
      double c = rxC[i];
      double a2 = v2(a, sa, ba);

      double ap = a2 * s * b / (sb + sc * c) - ra * a;
      double bp = -a2 * s * b / (sb + sc * c) + bb - rb * b;
      double at = a * (sin(a * 4.0 * M_PI));
      double cp = rc * (at - c);

      fA[i] = ap;
      fB[i] = bp;
      fC[i] = cp;
      fD[i] = 0;
      fE[i] = 0;
      fF[i] = 0;
    }
  }

  virtual void calcLateralActivator(std::vector<rxd3::real> &rxA, //
                                    std::vector<rxd3::real> &rxB, //
                                    std::vector<rxd3::real> &rxC, //
                                    std::vector<rxd3::real> &rxD, //
                                    std::vector<rxd3::real> &rxE, //
                                    std::vector<rxd3::real> &rxF, //
                                    std::vector<rxd3::real> &fA,  //
                                    std::vector<rxd3::real> &fB,  //
                                    std::vector<rxd3::real> &fC,
                                    std::vector<rxd3::real> &fD, //
                                    std::vector<rxd3::real> &fE, //
                                    std::vector<rxd3::real> &fF) {
    double base = 1e-4;
    rxParams.De = 0.5 * base;
    rxParams.Dd = 0.5 * base;
    rxParams.Dc = 0.00 * base;
    rxParams.Db = 1.0 * base;
    rxParams.Da = 1.00 * base;

    double            //
        c = 0.1,      //
        alpha = 0.04, //
        beta = 0.06,  //
        gamma = 0.04, //
        p0 = 0.00,
        p1 = 0.0; //

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0.98, 1.0);

    for (int i = 0; i < rxA.size(); i++) {

      double g1 = rxA[i];
      double g2 = rxB[i];
      double r = rxC[i];
      double s1 = rxD[i];
      double s2 = rxE[i];

      double s2g1 = c * s2 * g1 * g1;
      double s1g2 = c * s1 * g2 * g2;

      double g1p = s2g1 / r - alpha * g1 + p0;
      double g2p = s1g2 / r - alpha * g2 + p0;
      double rp = s2g1 + s1g2 - beta * r;
      double s1p = gamma * (g1 - s1) + p1;
      double s2p = gamma * (g2 - s2) + p1;

      fA[i] = g1p;
      fB[i] = g2p;
      fC[i] = rp;
      fD[i] = s1p;
      fE[i] = s2p;
      fF[i] = 0.0;
    }
  }

  virtual void calcThreeSystem(std::vector<rxd3::real> &rxA, //
                               std::vector<rxd3::real> &rxB, //
                               std::vector<rxd3::real> &rxC, //
                               std::vector<rxd3::real> &rxD, //
                               std::vector<rxd3::real> &rxE, //
                               std::vector<rxd3::real> &rxF, //
                               std::vector<rxd3::real> &fA,  //
                               std::vector<rxd3::real> &fB,  //
                               std::vector<rxd3::real> &fC,
                               std::vector<rxd3::real> &fD, //
                               std::vector<rxd3::real> &fE, //
                               std::vector<rxd3::real> &fF) {

    rxd3::real base = 5e-5;
    rxParams.Df = 0.5 * base;
    rxParams.De = 0.5 * rxParams.Df;

    rxParams.Dd = 0.0 * base;
    rxParams.Dc = 0.0 * rxParams.Dd;

    rxParams.Db = 1.0 * base;
    rxParams.Da = 0.50 * rxParams.Db;

    auto v2 = [](double v, double s, double b) {
      double v2 = v * v;
      return (v2 + b) / (1.0 + s * v2);
    };

    double //

        cbc = 0.6, //
        cbd = 0.2,
        cae = 0.4, //
        cea = 0.2, //

        cf = 0.00, //

        sa = 0.0,  //
        sd = 0.5,  //
        se = 0.00, //

        ba = 0.00,  //
        bb = 0.075, //
        bc = 0.01,  //
        bd = 0.06,  //
        be = 0.13,  //
        bf = 0.02,  //-2.0 * cea;

        ra = 0.13, //
        rb = bb,   //

        rc = bc + 0.05,      //
        rd = bd + rc + 0.02, //

        re = bf + 0.01, //
        rf = bf,        // rf < re
        s = 1.0;        //

    for (int i = 0; i < rxA.size(); i++) {

      double a = rxA[i];
      double b = rxB[i];
      double c = rxC[i];
      double d = rxD[i];
      double e = rxE[i];
      double f = rxF[i];

      double sba2 = s * b * v2(a, sa, ba);
      double rfe2 = re * f * v2(e, se, be);

      double ap = sba2 - ra * a - cae * e * a;
      double bp = bb * (1 + cbc * c + cbd * d) - sba2 - rb * b;

      double cp = rc * a * (c * c + bc) / (sd + d) - rc * c;
      double dp = rc * a * (c * c + bc) - rd * d + bd;

      double ep = rfe2 - re * e - cea * e * a;
      double fp = bf * a + cf - rfe2 - rf * f;

      fA[i] = ap;
      fB[i] = bp;
      fC[i] = cp;
      fD[i] = dp;
      fE[i] = ep;
      fF[i] = fp;
    }
  }

  virtual void diffuse(std::vector<rxd3::real> &rxA, //
                       std::vector<rxd3::real> &rxB, //
                       std::vector<rxd3::real> &rxC, //
                       std::vector<rxd3::real> &rxD, //
                       std::vector<rxd3::real> &rxE, //
                       std::vector<rxd3::real> &rxF, //
                       std::vector<rxd3::real> &fA,  //
                       std::vector<rxd3::real> &fB,  //
                       std::vector<rxd3::real> &fC,
                       std::vector<rxd3::real> &fD, //
                       std::vector<rxd3::real> &fE, //
                       std::vector<rxd3::real> &fF) {

    m2::diffuse<rxd3> diff(_meshGraph);

    double Da = rxParams.Da;
    double Db = rxParams.Db;
    double Dc = rxParams.Dc;
    double Dd = rxParams.Dd;
    double De = rxParams.De;
    double Df = rxParams.Df;
    double dt = rxParams.dt;

    std::cout << "Da: " << Da << std::endl;
    std::cout << "Db: " << Db << std::endl;
    std::cout << "Dc: " << Dc << std::endl;
    std::cout << "Dd: " << Dd << std::endl;
    std::cout << "De: " << De << std::endl;
    std::cout << "Df: " << Df << std::endl;

    rxA = diff.second_order(rxA, fA, dt, Da);
    rxB = diff.second_order(rxB, fB, dt, Db);
    rxC = diff.second_order(rxC, fC, dt, Dc);
    rxD = diff.second_order(rxD, fD, dt, Dd);
    rxE = diff.second_order(rxE, fE, dt, De);
    rxF = diff.second_order(rxF, fF, dt, Df);
  }

  virtual void clampRx(std::vector<rxd3::real> &rxA, //
                       std::vector<rxd3::real> &rxB, //
                       std::vector<rxd3::real> &rxC, //
                       std::vector<rxd3::real> &rxD, //
                       std::vector<rxd3::real> &rxE, //
                       std::vector<rxd3::real> &rxF) {

    for (int i = 0; i < rxA.size(); i++) {
      rxA[i] = std::clamp(rxA[i], 0.0, 1.0);
      rxB[i] = std::clamp(rxB[i], 0.0, 1.0);
      rxC[i] = std::clamp(rxC[i], 0.0, 1.0);
      rxD[i] = std::clamp(rxD[i], 0.0, 1.0);
      rxE[i] = std::clamp(rxE[i], 0.0, 1.0);
      rxF[i] = std::clamp(rxF[i], 0.0, 1.0);
    }
  }

  virtual void onAnimate(int frame) {

    _meshGraph->update_all();
    _meshGraph->reset_flags();
    _meshGraph->pack();

    _integrator->integrate();

    std::cout << "frame: " << frame << std::endl;
    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<rxd3>(_meshGraph) << std::endl;
    // gg::fillBuffer(_meshGraph, _obj);

    std::cout << "get rx " << std::endl;
    std::vector<rxd3::real> rxA = getRx(_meshGraph, rxd3::vertex_index::RXA);
    std::vector<rxd3::real> rxB = getRx(_meshGraph, rxd3::vertex_index::RXB);
    std::vector<rxd3::real> rxC = getRx(_meshGraph, rxd3::vertex_index::RXC);
    std::vector<rxd3::real> rxD = getRx(_meshGraph, rxd3::vertex_index::RXD);
    std::vector<rxd3::real> rxE = getRx(_meshGraph, rxd3::vertex_index::RXE);
    std::vector<rxd3::real> rxF = getRx(_meshGraph, rxd3::vertex_index::RXF);

    std::vector<rxd3::real> fA(rxA);
    std::vector<rxd3::real> fB(rxB);
    std::vector<rxd3::real> fC(rxC);
    std::vector<rxd3::real> fD(rxD);
    std::vector<rxd3::real> fE(rxE);
    std::vector<rxd3::real> fF(rxF);

    // std::cout << "curve to rx " << std::endl;
    curvature_to_rx(_meshGraph, rxA, rxB, rxC, rxD, rxE, rxF);

    // double Dv = rxParams.Dv;

    double dt = rxParams.dt;
    double dtt = rxParams.dtt;

    std::cout << "====== " << std::endl;
#if 1

    double An = 0.0, Bn = 0.0, Cn = 0.0;
    double Dn = 0.0, En = 0.0, Fn = 0.0;
    for (int i = 0; i < rxA.size(); i++) {
      An += rxA[i];
      Bn += rxB[i];
      Cn += rxC[i];
      Dn += rxD[i];
      En += rxE[i];
      Fn += rxF[i];
    }

    // clampRx(rxA, rxB, rxC, rxD, rxE, rxF);

    std::cout << "====== " << std::endl;
    std::cout << "rx stage norm" << std::endl;
    std::cout << "   An: " << An / double(rxA.size()) << std::endl;
    std::cout << "   Bn: " << Bn / double(rxB.size()) << std::endl;
    std::cout << "   Cn: " << Cn / double(rxC.size()) << std::endl;
    std::cout << "   Dn: " << Dn / double(rxA.size()) << std::endl;
    std::cout << "   En: " << En / double(rxB.size()) << std::endl;
    std::cout << "   Fn: " << Fn / double(rxC.size()) << std::endl;

    // calcGreyScott(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD, fE, fF);
    calc3ScaleGreyScott(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD, fE, fF);
    // calcActivatorDepletion(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD,
    // fE,fF);
    //  calcLateralActivator(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD, fE,
    //  fF); //
    //  calcThreeSystem(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD, fE, fF);

    diffuse(rxA, rxB, rxC, rxD, rxE, rxF, fA, fB, fC, fD, fE, fF);

    setRx(_meshGraph, rxA, rxd3::vertex_index::RXA);
    setRx(_meshGraph, rxB, rxd3::vertex_index::RXB);
    setRx(_meshGraph, rxC, rxd3::vertex_index::RXC);
    setRx(_meshGraph, rxD, rxd3::vertex_index::RXD);
    setRx(_meshGraph, rxE, rxd3::vertex_index::RXE);
    setRx(_meshGraph, rxF, rxd3::vertex_index::RXF);

    std::cout << "====== " << std::endl;
    std::cout << "area: " << m2::ci::area<rxd3>(_meshGraph) << std::endl;
    std::cout << " mean edge length: "
              << m2::ci::geometric_mean_length<rxd3>(_meshGraph) << std::endl;

#endif
    setRxColor(_meshGraph);
    auto colors = getRxColors(_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);

#if 1

    std::vector<rxd3::vec3> normals = calcNormals<rxd3>();
    std::cout << "normal length: " << normals.size() << std::endl;
    std::cout << "rx length: " << rxA.size() << std::endl;
    std::vector<rxd3::vec3> positions =
        m2::ci::get_coordinates<rxd3>(_meshGraph);
#if 0
    // std::cout << "print vecs" << std::endl;

    // print_vecs<rxd3>(positions, normals);
    std::cout << "rendering debug" << std::endl;
    _debugLines->renderLines();
#endif

    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    for (int i = 0; i < normals.size(); i++) {
      Nn += m2::va::norm<rxd3::real>(normals[i]);
      NRxn += m2::va::norm<rxd3::real>(normals[i]);
      positions[i] += dtt * normals[i];
    }

    std::cout << "pos stage norm: " << Nn / double(rxA.size()) << " "
              << NRxn / double(rxB.size()) << std::endl;

    m2::ci::set_coordinates<rxd3>(positions, _meshGraph);
    this->dump_gaudi(frame);
    _meshGraph->print();

#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "rxdiff." << frame << ".obj";
    m2::write_obj<rxd3>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    m2::flattened_surf<rxd3> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    m2::flattened_surf<rxd3> fsurf(_meshGraph);
    fsurf.write("rxd." + std::to_string(frame) + ".gaudi");
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
    double dt = 1.0;
    // double dtt = 0.0003;
    double dtt = 0.0001;

    double Da = 0.0;
    double Db = 0.0;
    double Dc = 0.0;
    double Dd = 0.0;
    double De = 0.0;
    double Df = 0.0;

    // double k = 0.061, F = 0.070; //pretty good one
    // double k = 0.0613, F = 0.0725; //this is a good one.
    // double k = 0.0613, F = 0.047; // this is a good one.
    // double k = 0.0613, F = 0.082;
    //  double k = 0.06132758, F = 0.037;
    //  double k = 0.059, F = 0.03;
    //  double k = 0.056, F = 0.098; //interesting?
    //  double k = 0.0613, F = 0.06; // uskatish
    double k = 0.0613, F = 0.075;

    //  double k = 0.0628, F = 0.0567;
    //  double k = 0.061, F = 0.42;
    //  double k = 0.045, F = 0.01; // waves

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

    /*
    Slider *slider0 = new Slider(window);
    slider0->setValue(pow(this->scene->rxParams.Du, 1.0 / Kd));
    slider0->setFixedWidth(w);

    slider0->setCallback(
        [this, Kd](float value) { this->scene->rxParams.Du = pow(value, Kd); });
    */
    Slider *slider2 = new Slider(window);
    slider2->setValue(this->scene->rxParams.F * Kf);
    slider2->setFixedWidth(w);

    Slider *slider3 = new Slider(window);
    slider3->setValue(this->scene->rxParams.k * Kf);
    slider3->setFixedWidth(w);

    Slider *slider4 = new Slider(window);
    slider4->setValue(this->scene->rxParams.dt / 2.0);
    slider4->setFixedWidth(w);

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
