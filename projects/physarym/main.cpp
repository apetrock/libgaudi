//#include "nanoguiincludes.h"

#include "manifold/coordinate_interface.hpp"
#include "manifold/m2.hpp"
#include <algorithm>
#include <cmath>
#include <exception>
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

#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"

#include "GaudiGraphics/geometry_logger.h"
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

#include <Eigen/Geometry>

//#include "m2Operators.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

template <typename T> class physarym_space {
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

  enum class edge_index { MAXINDEX = 0 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    SMOOTH = 2,
    VELOCITY = 3,
    MAGNITUDE = 4,
    MAXINDEX = 5
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
    case vertex_index::SMOOTH:
      return storage_type::REAL;
    case vertex_index::VELOCITY:
      return storage_type::VEC3;
    case vertex_index::MAGNITUDE:
      return storage_type::REAL;
    default:
      return storage_type::SIZE;
    }
  }
};

typedef physarym_space<double> physarym;

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

    m2::obj_loader<physarym> load;
    m2::subdivide<physarym> sub;
    m2::make<physarym> mk;
    m2::convex_hull<physarym> ch;

    m2::construct<physarym> bevel;
    m2::affine<physarym> mod;
    std::string start_frame = "";
    // std::string start_frame = "physarym.46.gaudi";
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

      m2::remesh<physarym> rem;
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

    //_integrator = new m2::surf_integrator<physarym>(_meshGraph,
    // 0.2, 3.0, 3.0);
    _integrator = new m2::surf_integrator<physarym>(_meshGraph, 0.2, 3.0, 2.0);

    _integrator->add_default_vertex_policy<typename physarym::real>(
        physarym::vertex_index::SMOOTH);
    _integrator->add_default_vertex_policy<typename physarym::vec3>(
        physarym::vertex_index::VELOCITY);
    _integrator->add_default_vertex_policy<typename physarym::real>(
        physarym::vertex_index::MAGNITUDE);

    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  template <typename SPACE>
  vector<m2::colorRGB> getColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);
    auto hot = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::MAGNITUDE);

    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    std::vector<m2::colorRGB> vert_colors(smooth.size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = K[i];
      typename SPACE::real N = 0;
      typename SPACE::real s = smooth[i];
      typename SPACE::real h = hot[i];

      typename SPACE::coordinate_type colorH(0.5, 0.5, 0.0);
      typename SPACE::coordinate_type colorS(0.0, 0.75, 0.75);
      typename SPACE::coordinate_type colorC(0.56, 0.50, 0.60);

      typename SPACE::coordinate_type mx = m2::va::mix(s, colorS, colorC);
      mx = m2::va::mix(h, colorH, mx);
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
                  const std::vector<typename SPACE::vec3> &p1, double D = 0.05,
                  int col = 0) {
    double mx =
        std::accumulate(p1.begin(), p1.end(), 0.0, [](double a, auto &c) {
          return max(a, m2::va::norm(c));
        });

    for (int i = 0; i < p0.size(); i++) {

      const auto &p = p0[i];
      const auto &a = p1[i];

      auto pa = p + D * a / mx;
      // std::cout << a.transpose() << std::endl;
      auto c = grey(m2::va::norm(a));
      if (col == 1)
        c = red(m2::va::norm(a));
      if (col == 2)
        c = green(m2::va::norm(a));
      if (col == 3)
        c = blue(m2::va::norm(a));

      gg::geometry_logger::line(p, pa, Vec4d(c[0], c[1], c[2], 1.0));
    }
  }

  template <typename SPACE>
  void stats(const std::vector<typename SPACE::vec3> &u0) {
    M2_TYPEDEFS;

    real e = 0.0, stddev = 0, mn = 9999999, mx = 0, sum = 0.0;
    int i = 0;

    for (auto &v : u0) {
      real u = u0[i].norm();
      sum += u;
      mn = std::min(u, mn);
      mx = std::max(u, mx);
      e += u;
      i++;
    }

    i = 0;
    e /= real(u0.size());
    for (auto &v : u0) {
      real s = u0[i].norm() - e;
      stddev += (s * s);
      i++;
    }

    stddev /= real(u0.size());
    std::cout << " avg scaling, mean: " << e << ", min: " << mn
              << ", max: " << mx << ", std_dev: " << sqrt(stddev)
              << ", sum: " << sum << std::endl;
  }

  template <typename SPACE>
  std::vector<typename SPACE::vec3>
  make_smooth(m2::surf<physarym> *surf,
              const std::vector<typename SPACE::vec3> &u0,
              const std::vector<physarym::vec3> &positions, double reg) {
    M2_TYPEDEFS;

    m2::mesh_calculator<SPACE> calc;

    std::vector<typename SPACE::vec3> uf =
        m2::ci::verts_to_faces<SPACE, typename SPACE::vec3>(u0, surf, false);

    std::vector<typename SPACE::vec3> u1 =
        calc.harmonicAvg(surf, uf, positions, reg);
    stats<SPACE>(u0);
    stats<SPACE>(uf);
    stats<SPACE>(u1);

    /*
        real ec = e0 / e1;
        for (auto &v : u1) {
          v *= ec;
        }
    */
    return u1;
  }

#endif
  template <typename SPACE>
  std::vector<typename SPACE::vec3>
  randomize(const std::vector<typename SPACE::vec3> &u0,
            const std::vector<typename SPACE::vec3> &normals,
            const double &turning) {
    M2_TYPEDEFS;

    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    std::vector<typename SPACE::vec3> u1 = u0;
    for (int i = 0; i < u0.size(); i++) {
      auto vel0 = u0[i];
      auto vel1 = vel0;

      if (unif01(rng) > 1.0) {
        double t = 0.1;
        coordinate_type velt = m2::va::orthogonal_project(normals[i], vel0);
        // physarym::coordinate_type cp = velocities_smooth[i];
        vel1 = t * velt + (1.0 - t) * vel0;
      } else {
        coordinate_type u = vel0.normalized();
        coordinate_type N = normals[i];
        coordinate_type w = m2::va::cross(N, u);
        coordinate_type v = m2::va::cross(N, w);

        double theta = turning;

        if (unif01(rng) > 0.5)
          theta *= -1;

        vel1 = Eigen::AngleAxisd(theta, w) * vel0;

        if (unif01(rng) > 0.5)
          theta *= -1;

        vel1 = Eigen::AngleAxisd(theta, v) * vel1;
      }

      u1[i] = vel1.normalized();
    }

    return u1;
  }

  template <typename SPACE>
  std::vector<typename SPACE::vec3> calc_growth_0(double dt) {
    M2_TYPEDEFS;

    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    m2::mesh_calculator<SPACE> calc;

    std::vector<coordinate_type> positions =
        m2::ci::get_coordinates<SPACE>(_meshGraph);

    std::vector<coordinate_type> normals =
        m2::ci::get_vertex_normals<SPACE>(_meshGraph);

    std::vector<coordinate_type> u0 = m2::ci::get<SPACE, coordinate_type>(
        _meshGraph, SPACE::vertex_index::VELOCITY);

    std::vector<real> mag =
        m2::ci::get<SPACE, real>(_meshGraph, SPACE::vertex_index::MAGNITUDE);

    real turning = 0.15 * M_PI;

    std::vector<typename SPACE::vec3> u1 =
        randomize<SPACE>(u0, normals, turning);
    std::cout << "calc_covariance:" << std::endl;

    std::vector<typename SPACE::mat43> cov =
        calc.template covariance(_meshGraph, positions, 2.0 * _max);
    // std::vector<typename SPACE::mat43> cov =
    //     calc.template curvature(_meshGraph, positions, 4.0 * _max);

    int i = 0;
    int k = 0;
#if 0
    for (auto p : positions) {
      typename SPACE::mat43 US = cov[i];
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type s = US.row(3);

      coordinate_type u = U.col(k).transpose();
      U.block(0, 0, 3, 1) *= s[0];
      U.block(0, 1, 3, 1) *= s[1];
      U.block(0, 2, 3, 1) *= s[2];

      gg::geometry_logger::frame(U, p, 0.1);
      i++;
    }
#endif
    for (int i = 0; i < u0.size(); i++) {
      double dm, thet0, thet1;
      dm = thet0 = thet1 = 99999;
      auto u0i = u0[i];
      auto u1i = u1[i];
      auto p = positions[i];

      auto angle = [](const coordinate_type &p, const coordinate_type &N,
                      real thet, const std::vector<coordinate_type> &targets,
                      double &tmin, int &jmin) {
        auto Nn = N.normalized();
        real dmin = 99999;
        tmin = 2.0 * M_PI;
        jmin = -1;

        for (int j = 0; j < targets.size(); j++) {
          coordinate_type dp = targets[j] - p;
          double d = dp.norm();
          dp /= d;
          real ti = acos(dp.dot(Nn));
          if (ti < thet && d < dmin) {
            dmin = d;
            jmin = j;
            tmin = ti;
          }
        }
        return jmin > 0 ? dmin : -0.1;
      };
      typename SPACE::mat43 US = cov[i];
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type uc0 = U.col(0).transpose();
      coordinate_type uc1 = U.col(1).transpose();
      coordinate_type uc2 = U.col(2).transpose();
      coordinate_type s = US.row(3);

      coordinate_type N = normals[i];
      coordinate_type N0 = m2::va::sgn(m2::va::dot(uc0, u0i)) * uc0;
      coordinate_type N1 = m2::va::sgn(m2::va::dot(uc1, u0i)) * uc1;
      coordinate_type N2 = m2::va::sgn(m2::va::dot(uc2, N)) * uc0;

      double a0, a1, aN;
      int j0, j1, jN;
      double d0 = angle(p, u0i, 2.0 * turning, _targets, a0, j0);
      double d1 = angle(p, u1i, 2.0 * turning, _targets, a1, j1);
      double dN = angle(p, N0, 2.0 * turning, _targets, aN, jN);
#if 0
      if (jN > 0)
        gg::geometry_logger::line(p, _targets[jN], Vec4d(0.0, 1.0, 0.0, 1.0));
#endif

      real magi = 0.0;
#if 1
      // distilled
      u0[i] = N0;
      magi = 1.0 * (pow(s[0] / s[1], 4.0) - 4.0);
      mag[i] = magi * N.dot(u0[i]);
      mag[i] = std::min(1.0, mag[i]);
      mag[i] = std::max(-1.0, mag[i]);
#else
      double s0 = s[0] / (s[1] + s[2]);
      double C0 = pow(s0, 1.0) - 1.05;
      // double C0 = s[0] - s[1] - s[2] - 0.1;

      // double C0 = s[2] - 0.05;

      double C1 = pow(s[1] / s[0], 2.0);

      double Cn = N.dot(u0[i].normalized());
      if ((d0 > 0 && d1 < 0) || (d0 > 0 && d1 > 0 && a0 < a1)) {
        // do nothing
      } else if (unif01(rng) > 0.5) {
        u0[i] = u1i.normalized();
        mag[i] = 1.0 * Cn * C0;
      }

      mag[i] = std::min(1.0, mag[i]);
      mag[i] = std::max(-1.0, mag[i]);

#endif
    }

    for (int i = 0; i < u1.size(); i++) {
      u1[i] = mag[i] * u0[i];
    }

    std::cout << "make smooth:" << std::endl;
    std::vector<coordinate_type> us =
        make_smooth<SPACE>(_meshGraph, u1, positions, 0.75 * _max);

    for (int i = 0; i < u0.size(); i++) {
      // u0[i] += 0.025 * us[i].normalized();
      u0[i].normalize();
    }

    m2::ci::set<SPACE, real>(_meshGraph, mag,
                             physarym::vertex_index::MAGNITUDE);
    m2::ci::set<SPACE, coordinate_type>(_meshGraph, u0,
                                        physarym::vertex_index::VELOCITY);

#if 0
    // std::cout << "print vecs" << std::endl;

    print_vecs<SPACE>(positions, us, 0.1);

    for (int j = 0; j < _targets.size() - 1; j++) {
      auto p0 = _targets[j];
      auto p1 = _targets[j + 1];

      gg::geometry_logger::line(p0, p1, Vec4d(1.0, 0.0, 0.0, 1.0));
    }
    std::cout << "rendering debug" << std::endl;

#endif

    return us;
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
    using vertex_array = m2::surf<physarym>::vertex_array;
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
              << m2::ci::geometric_mean_length<physarym>(_meshGraph)
              << std::endl;

    auto colors = getColor(_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);
#if 1
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::uniform_real_distribution<double> unif01(0.0, 1.0);

    std::vector<physarym::vec3> normals =
        m2::ci::get_vertex_normals<physarym>(_meshGraph);

    int i = 0;
    if (!_vel_inited) {
      vertex_array &verts = _meshGraph->get_vertices();
      for (auto v : verts) {
        double c0 = unif(rng);
        double c1 = unif(rng);
        double c2 = unif(rng);
        physarym::coordinate_type c(c0, c1, c2);
        // c = normals[i] + 0.01 * c;
        // c = normals[i] + 0.01 * c;

        v->template set<physarym::coordinate_type>(
            physarym::vertex_index::VELOCITY, c);
        i++;
      }
      for (int i = 0; i < 32; i++) {
        double c0 = unif(rng);
        double c1 = unif(rng);
        double c2 = unif(rng);
        physarym::coordinate_type c(c0, c1, c2);

        _targets.push_back(2.0 * c);
      }
      _vel_inited = true;
    }

    std::vector<physarym::vec3> positions =
        m2::ci::get_coordinates<physarym>(_meshGraph);
    //    std::vector<physarym::vec3> normals =
    //        m2::ci::get_vertex_normals<physarym>(_meshGraph);

    double dt = 0.01;
    std::vector<physarym::vec3> u1 = calc_growth_0<physarym>(dt);

    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    for (int i = 0; i < u1.size(); i++) {
      Nn += m2::va::norm<physarym::real>(u1[i]);
      NRxn += m2::va::norm<physarym::real>(u1[i]);
      positions[i] += dt * u1[i];
    }

    m2::ci::set_coordinates<physarym>(positions, _meshGraph);
    this->dump_gaudi(frame);
    _meshGraph->print();
    gg::geometry_logger::render();

#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "physarym." << frame << ".obj";
    m2::write_obj<physarym>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    m2::flattened_surf<physarym> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    m2::flattened_surf<physarym> fsurf(_meshGraph);
    fsurf.write("physarym." + std::to_string(frame) + ".gaudi");
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
  double _max = 0.0;
  double _min = 0.0;
  bool _vel_inited = false;

  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  m2::surf<physarym> *_meshGraph;
  m2::surf_integrator<physarym> *_integrator;
  std::vector<physarym::coordinate_type> _targets;
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
