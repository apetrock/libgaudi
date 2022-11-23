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

  enum class edge_index { MAXINDEX = 0 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    SMOOTH = 2,
    MAXINDEX = 3
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

    asawa::obj_loader<stretch> load;
    asawa::subdivide<stretch> sub;
    asawa::make<stretch> mk;
    asawa::convex_hull<stretch> ch;

    asawa::construct<stretch> bevel;
    asawa::affine<stretch> mod;
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

      asawa::remesh<stretch> rem;
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
    //_integrator = new asawa::surf_integrator<stretch>(_meshGraph, 0.2, 3.0,
    // 0.35);
    _integrator =
        new asawa::surf_integrator<stretch>(_meshGraph, 0.1, 3.0, 0.5);
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
  vector<asawa::colorRGB> getColor(asawa::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth =
        asawa::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);

    asawa::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    std::vector<asawa::colorRGB> vert_colors(smooth.size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = K[i];
      typename SPACE::real N = 0;
      typename SPACE::real s = smooth[i];

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorC(0.56, 0.50, 0.60);

      typename SPACE::coordinate_type mx = va::mix(s, colorS, colorC);
      vert_colors[i] = asawa::colorRGB(mx[0], mx[1], mx[2], 1.0);
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
        std::accumulate(p1.begin(), p1.end(), 0.0,
                        [](double a, auto &c) { return max(a, va::norm(c)); });

    for (int i = 0; i < p0.size(); i++) {

      const auto &p = p0[i];
      const auto &a = p1[i];

      auto pa = p + D * a / mx;
      // std::cout << a.transpose() << std::endl;
      auto c = grey(va::norm(a));
      if (col == 1)
        c = red(va::norm(a));
      if (col == 2)
        c = green(va::norm(a));
      if (col == 3)
        c = blue(va::norm(a));

      _debugLines->pushLine(Vec4(p[0], p[1], p[2], 1.0),
                            Vec4(pa[0], pa[1], pa[2], 1.0),
                            Vec4(c[0], c[1], c[2], 1.0));
    }
  }

  template <typename SPACE>
  std::vector<typename SPACE::vec3>
  make_smooth(asawa::surf<stretch> *surf,
              const std::vector<typename SPACE::vec3> &field_vert,
              const std::vector<stretch::vec3> &positions, double reg) {
    M2_TYPEDEFS;

    asawa::mesh_calculator<SPACE> calc;

    std::vector<typename SPACE::vec3> field_face =
        asawa::ci::verts_to_faces<SPACE, typename SPACE::vec3>(field_vert,
                                                               surf);

    std::vector<typename SPACE::vec3> field_vert_smooth =
        calc.harmonicAvg(surf, field_face, positions, reg);
    return field_vert_smooth;
  }

  template <typename SPACE>
  std::vector<typename SPACE::vec3>
  make_div_free(asawa::surf<stretch> *surf,
                const std::vector<typename SPACE::vec3> &field_vert,
                const std::vector<stretch::vec3> &positions, double reg) {
    M2_TYPEDEFS;

    asawa::mesh_calculator<SPACE> calc;

    std::vector<typename SPACE::vec3> field_face =
        asawa::ci::verts_to_faces<SPACE, typename SPACE::vec3>(field_vert,
                                                               surf);

    std::vector<typename SPACE::vec3> field_vert_smooth =
        calc.harmonicAvg(surf, field_face, positions, reg);
    // return field_vert_smooth;
    std::vector<typename SPACE::real> div_vert =
        calc.template divergence<typename SPACE::real, typename SPACE::vec3>(
            surf, field_face, positions, reg);

    std::vector<typename SPACE::real> div_face =
        asawa::ci::verts_to_faces<SPACE, typename SPACE::real>(div_vert, surf);

    std::vector<typename SPACE::vec3> grads =
        calc.template gradient<typename SPACE::vec3, typename SPACE::real>(
            surf, div_face, positions, reg);

    for (int i = 0; i < field_vert_smooth.size(); i++) {
      field_vert_smooth[i] = field_vert_smooth[i] + grads[i];
    }

    return field_vert_smooth;
  }

  template <typename SPACE> std::vector<typename SPACE::vec3> calc_vel() {
    M2_TYPEDEFS;

    asawa::mesh_calculator<SPACE> calc;
    std::vector<stretch::vec3> positions =
        asawa::ci::get_coordinates<stretch>(_meshGraph);

    std::cout << "creating vecs" << std::endl;
    std::vector<asawa::face_triangle<SPACE>> triangles;
    auto &faces = _meshGraph->get_faces();

    std::cout << 1 << std::endl;
    for (int i = 0; i < faces.size(); i++) {
      if (!_meshGraph->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<asawa::face_triangle<SPACE>> tris =
          asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    double reg = 1.0 * _max;
    using ls_sphere_t = ls_sphere<real, coordinate_type>;
    double Rc = 64.0;
    std::vector<ls_sphere_t> spheres =
        calc.least_squares_sphere(_meshGraph, positions, Rc * reg);
    std::vector<typename SPACE::vec3> centers(positions.size());
    std::vector<typename SPACE::vec3> dcenters(positions.size());

    for (int i = 0; i < spheres.size(); i++) {
      coordinate_type c;
      double r;
      spheres[i].calc(c, r);
      // std::cout << r << " - " << positions[i].transpose() << std::endl;
      // centers[i] = positions[i] - c;
      dcenters[i] = c;
      centers[i] = positions[i] + c;
    }
    dcenters = make_smooth<stretch>(_meshGraph, dcenters, positions, Rc * reg);

    std::vector<typename SPACE::mat43> cov =
        calc.template covariance(_meshGraph, centers, reg);

    coordinate_type c = asawa::ci::center<SPACE>(_meshGraph);
    std::vector<typename SPACE::vec3> stretchV(positions.size());

    int i = 0;
    int k = 0;
    for (auto p : positions) {
      coordinate_type dp = p - c;
      typename SPACE::mat43 US = cov[i];
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type u = U.col(k).transpose();

      coordinate_type s = US.row(3);
      typename SPACE::real du = va::dot(u, dp);
      // std::cout << s.transpose() << std::endl;
      stretchV[i] = va::sgn(du) * s[k] * u;
      i++;
    }
    // std::vector<typename SPACE::vec3> smooth_stretch =
    //     make_div_free<stretch>(_meshGraph, stretchV, positions, reg);
    std::vector<typename SPACE::vec3> smooth_stretch =
        make_smooth<stretch>(_meshGraph, stretchV, positions, reg);

    print_vecs<stretch>(positions, smooth_stretch, 0.1, 3);

    return smooth_stretch;
  }

  template <typename SPACE> std::vector<typename SPACE::vec3> calc_vel_1() {
    M2_TYPEDEFS;

    asawa::mesh_calculator<SPACE> calc;
    std::vector<stretch::vec3> positions =
        asawa::ci::get_coordinates<stretch>(_meshGraph);

    std::cout << "creating vecs" << std::endl;
    std::vector<asawa::face_triangle<SPACE>> triangles;
    auto &faces = _meshGraph->get_faces();

    std::cout << 1 << std::endl;
    for (int i = 0; i < faces.size(); i++) {
      if (!_meshGraph->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<asawa::face_triangle<SPACE>> tris =
          asawa::ci::get_tris<SPACE>(faces[i]);
      triangles.insert(triangles.end(), tris.begin(), tris.end());
    }

    double reg = 1.0 * _max;
    using ls_sphere_t = ls_sphere<real, coordinate_type>;
    double Rc = 16.0 * reg;
    std::cout << "reg = " << Rc * reg << std::endl;
    std::vector<ls_sphere_t> spheres =
        calc.least_squares_sphere(_meshGraph, positions, Rc);
    std::vector<typename SPACE::vec3> centers(positions.size());
    std::vector<typename SPACE::vec3> dcenters(positions.size());

    for (int i = 0; i < spheres.size(); i++) {
      coordinate_type c;
      double r;
      spheres[i].calc(c, r);
      // std::cout << r << " - " << positions[i].transpose() << std::endl;
      // centers[i] = positions[i] - c;
      dcenters[i] = c;
      centers[i] = positions[i] + c;
    }
    // dcenters = make_smooth<stretch>(_meshGraph, dcenters, positions, Rc);

    std::vector<typename SPACE::mat43> cov =
        calc.template curvature(_meshGraph, centers, reg);

    coordinate_type c = asawa::ci::center<SPACE>(_meshGraph);
    std::vector<typename SPACE::vec3> stretchV(positions.size());

    int i = 0;
    int k = 1;
    for (auto p : positions) {
      coordinate_type dp = p - c;
      typename SPACE::mat43 US = cov[i];
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type u = U.col(k).transpose();

      coordinate_type s = US.row(3);
      typename SPACE::real du = va::dot(u, dp);
      // std::cout << s.transpose() << std::endl;
      stretchV[i] = va::sgn(du) * s[k] * u;
      i++;
    }

    std::vector<typename SPACE::vec3> stretchF =
        asawa::ci::verts_to_faces<SPACE, typename SPACE::vec3>(stretchV,
                                                               _meshGraph);

    std::vector<typename SPACE::vec3> field_vert_smooth =
        calc.curl(_meshGraph, stretchF, positions, reg);

    // std::vector<typename SPACE::vec3> smooth_stretch =
    //     make_smooth<stretch>(_meshGraph, stretchV, positions, reg);

    // print_vecs<stretch>(positions, field_vert_smooth, 0.05, 3);

    return field_vert_smooth;
  }

  /*
    template <typename SPACE>
    std::vector<typename SPACE::vec3> find_singularities(std::vector<typename
    SPACE::positions> positions) { M2_TYPEDEFS; using real = typename
    SPACE::real; using vec3 = typename SPACE::vec3;

      asawa::mesh_calculator<SPACE> calc;
      std::vector<typename SPACE::mat43> cov =
          calc.template covariance(positions, reg);

      std::vector<vec3> stretchV(positions.size());

      int i = 0;
      int k = 2;
      for (auto p : positions) {
        // coordinate_type dp = p - ccenters[i];
        coordinate_type dp = p - c;

        typename SPACE::mat43 US = cov[i];
        typename SPACE::mat3 U = US.block(0, 0, 3, 3);
        coordinate_type u = U.col(k).transpose();
        coordinate_type s = US.row(3);
        real du = va::dot(u, dp);
        real sk = s[k];
        stretchV[i] = 1.0 * va::sgn(du) * sk * u;
        i++;
      }
    }
  */

  template <typename SPACE> std::vector<typename SPACE::vec3> calc_vel_2() {
    M2_TYPEDEFS;
    using real = typename SPACE::real;
    using vec3 = typename SPACE::vec3;

    asawa::mesh_calculator<SPACE> calc;
    std::vector<vec3> positions = asawa::ci::get_coordinates<SPACE>(_meshGraph);

    std::cout << "creating vecs" << std::endl;
    auto &faces = _meshGraph->get_faces();

    std::vector<asawa::face_triangle<SPACE>> tris(faces.size());
    std::vector<vec3> centers(faces.size());
    std::vector<real> areaF(faces.size());

    std::cout << 1 << std::endl;
    for (int i = 0; i < faces.size(); i++) {
      if (!_meshGraph->has_face(i))
        continue;
      if (faces[i]->size() < 3)
        continue;
      std::vector<asawa::face_triangle<SPACE>> t =
          asawa::ci::get_tris<SPACE>(faces[i]);
      tris[i] = t[0];
      centers[i] = t[0].center();
      areaF[i] = 1.0;
    }

    double reg = 1.0 * _max;
    std::vector<real> areaV(positions.size(), 1.0);

    using ls_sphere_t = ls_sphere<real, coordinate_type>;
    double Rc = 64.0;
    std::vector<ls_sphere_t> spheres =
        calc.least_squares_sphere(_meshGraph, positions, Rc * reg);
    std::vector<vec3> ccenters(positions.size());
    std::vector<real> regs(positions.size());

    for (int i = 0; i < spheres.size(); i++) {
      coordinate_type c;
      double r;
      spheres[i].calc(c, r);
      // std::cout << r << " - " << positions[i].transpose() << std::endl;
      // centers[i] = positions[i] - c;
      regs[i] = r;
      ccenters[i] = positions[i] + c;
    }

    std::vector<typename SPACE::mat43> cov =
        calc.template curvature(_meshGraph, ccenters, regs);

    // std::vector<typename SPACE::mat43> cov =
    //     calc.template curvature(_meshGraph, positions, reg);

    // std::vector<typename SPACE::mat43> cov =
    //     calc.template curvature(_meshGraph, ccenters, regs);

    coordinate_type c = asawa::ci::center<SPACE>(_meshGraph);
    std::vector<vec3> stretchV(positions.size());

    int i = 0;
    int k = 2;
    for (auto p : positions) {
      // coordinate_type dp = p - ccenters[i];
      coordinate_type dp = p - c;

      typename SPACE::mat43 US = cov[i];
      typename SPACE::mat3 U = US.block(0, 0, 3, 3);
      coordinate_type u = U.col(k).transpose();
      coordinate_type s = US.row(3);
      real du = va::dot(u, dp);
      real sk = s[k];
      stretchV[i] = 1.0 * va::sgn(du) * sk * u;
      i++;
    }
    // std::vector<typename SPACE::vec3> smooth_stretch =
    //     make_div_free<stretch>(_meshGraph, stretchV, positions, reg);
    std::vector<vec3> smooth_stretch =
        make_smooth<stretch>(_meshGraph, stretchV, positions, reg);

    /*
    std::vector<vec3> forces =
        calc.template gravitation<vec3, real>(_meshGraph, areaF, positions,
    areaV, 8.0 * reg);

    i = 0;
    for (auto &v : smooth_stretch) {
      //std::cout << forces[i] << std::endl;
      v += 10.0 * forces[i++];;
    }
    */

    // print_vecs<SPACE>(positions, smooth_stretch, 0.1, 3);
    return smooth_stretch;
  }

#endif

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

    _integrator->integrate();

    std::cout << "frame: " << frame << std::endl;
    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << asawa::ci::geometric_mean_length<stretch>(_meshGraph)
              << std::endl;

    double dt = _params.dt;
    auto colors = getColor(_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);
#if 1

    // std::vector<stretch::vec3> normals = calc_vel_1<stretch>();
    std::vector<stretch::vec3> normals = calc_vel_2<stretch>();
    // std::vector<stretch::vec3> normals = calc_vel<stretch>();

    std::vector<stretch::vec3> positions =
        asawa::ci::get_coordinates<stretch>(_meshGraph);
#if 1
    // std::cout << "print vecs" << std::endl;

    // print_vecs<stretch>(positions, normals);
    std::cout << "rendering debug" << std::endl;
    _debugLines->renderLines();
#endif

    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    for (int i = 0; i < normals.size(); i++) {
      Nn += va::norm<stretch::real>(normals[i]);
      NRxn += va::norm<stretch::real>(normals[i]);
      positions[i] += dt * normals[i];
    }

    asawa::ci::set_coordinates<stretch>(positions, _meshGraph);
    this->dump_gaudi(frame);
    _meshGraph->print();

#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "stretch." << frame << ".obj";
    asawa::write_obj<stretch>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    asawa::flattened_surf<stretch> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    asawa::flattened_surf<stretch> fsurf(_meshGraph);
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
  asawa::surf<stretch> *_meshGraph;
  asawa::surf_integrator<stretch> *_integrator;
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
