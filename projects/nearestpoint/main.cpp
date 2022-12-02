//#include "nanoguiincludes.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
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

#include <complex>
#include <iostream>
#include <random>
#include <string>

#include "GaudiGraphics/buffers.hpp"
#include "GaudiGraphics/geometry_logger.h"
#include "GaudiGraphics/mesh_helper.hpp"
#include "GaudiGraphics/viewer.hpp"

#include "manifold/asawa/asawa.h"

#include "manifold/hepworth/constraints_init.hpp"
#include "manifold/hepworth/objective_function.hpp"
#include "manifold/hepworth/optimizer.hpp"

#include "manifold/bins.hpp"

#include "manifold/calder/harmonic_integrators.hpp"
#include "manifold/vec_addendum.h"

#define TRACKBALLSIZE (0.8f)
#define RENORMCOUNT 97

using std::cerr;
using std::cout;
using std::endl;

using namespace GaudiMath;

template <typename T> class nearest_point_space {
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

  enum class edge_index { PHI0 = 0, MAXINDEX = 1 };

  enum class vertex_index {
    COORDINATE = 0,
    COLOR = 1,
    SMOOTH = 2,
    ACT = 3,

    MAXINDEX = 4
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
    case vertex_index::ACT:
      return storage_type::VEC3;
    default:
      return storage_type::SIZE;
    }
  }
};

////////////////////////////////////////////////////////////////////////////
// ENERGY
////////////////////////////////////////////////////////////////////////////

template <typename SPACE>
vector<typename SPACE::real>
calcEnergy(asawa::surf<SPACE> *mesh,
           const std::vector<typename SPACE::coordinate_type> &evalPoints,
           typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  using Avg_Integrator =
      calder::Geometry_Integrator<SPACE, real, triangle_type, real>;

  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;
  std::vector<vec3> vertex_normals = ci::get_vertex_normals<SPACE>(mesh);
  std::vector<real> vertex_areas = ci::get_vertex_areas<SPACE>(mesh);
  std::vector<real> faceVals = std::vector<real>(mesh->get_faces().size(), 0.0);

  if (evalPoints.empty())
    return std::vector<real>();

  auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                        ATree &tree, real &netCharge, coordinate_type &avgPoint,
                        coordinate_type &avgNormal) -> void {
    avgPoint = coordinate_type(0, 0, 0);
    netCharge = 0.0;

    T netWeight = 0;

    for (int i = node.begin; i < node.begin + node.size; i++) {
      int ii = tree.permutation[i];
      triangle_type tri = tris[ii];
      T w = tri.area();
      avgPoint += w * tri.center();
      avgNormal += w * tri.normal();
      netCharge += w * faceVals[ii];
      netWeight += w;
    }

    avgPoint /= netWeight;
  };

  auto compute = [&vertex_normals, &vertex_areas, faceVals, regLength](
                     int i_c, const real &wq, const coordinate_type &pc,
                     const coordinate_type &pe, const coordinate_type &N,
                     const vector<triangle_type> &tris, ANode &node,
                     ATree &tree) -> real {
    real out = 0.0;

    coordinate_type Nv = vertex_normals[i_c];
    real wv = vertex_areas[i_c];

    auto computeTangentRadius = [](coordinate_type dp, coordinate_type N) {
      real dpn = dp.norm();
      real dpnProj = ((N * N.transpose()) * dp).norm();
      return dpn * dpn / dpnProj;
    };

    auto K = [](coordinate_type dp, coordinate_type N, double p) {
      mat3 P = N * N.transpose();
      coordinate_type Pdp = P * dp;
      real ndp = dp.norm();

      if (ndp == 0.0)
        return 0.0;

      real nPdp = (Pdp).norm();
      real k = pow(nPdp, p) / pow(ndp, 2.0 * p);

      return k;
    };

    if (node.isLeaf()) {
      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        real qi = faceVals[ii];
        auto tri = tris[ii];
        auto w = tri.area();
        auto c = tri.center();
        auto Nf = tri.normal();
        coordinate_type dp = c - pe;
        real k = K(dp, Nv, 3.0);
        out += w * k;
      }
    } else {
      coordinate_type dp = pc - pe;
      real w = N.norm();
      coordinate_type Nf = N / w;
      real k = K(dp, Nv, 6.0);
      out += w * k;
    }
    return out;
  };

  vector<face_ptr> &faces = mesh->get_faces();
  vector<coordinate_type> normals;
  vector<triangle_type> triangles;
  for (int i = 0; i < faces.size(); i++) {
    if (!mesh->has_face(i))
      continue;
    if (faces[i]->size() < 3)
      continue;
    std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
    triangles.insert(triangles.end(), tris.begin(), tris.end());
  }

  for (auto t : triangles) {
    normals.push_back(t.normal());
  }
  vector<real> u(evalPoints.size(), z::zero<real>());
  Avg_Integrator integrator;
  integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);

  return u;
}

template <typename SPACE>
vector<typename SPACE::vec3>
calcGradient(asawa::surf<SPACE> *mesh,
             const std::vector<typename SPACE::coordinate_type> &evalPoints,
             typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  using Avg_Integrator =
      calder::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;
  std::vector<vec3> vertex_normals = ci::get_vertex_normals<SPACE>(mesh);
  std::vector<real> vertex_areas = ci::get_vertex_areas<SPACE>(mesh);

  std::vector<vec3> faceVals = ci::verts_to_faces<SPACE>(vertex_normals, mesh);

  if (evalPoints.empty())
    return std::vector<vec3>();

  auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                        ATree &tree, vec3 &netCharge, coordinate_type &avgPoint,
                        coordinate_type &avgNormal) -> void {
    avgPoint = coordinate_type(0, 0, 0);
    netCharge = z::zero<vec3>();

    T netWeight = 0;

    for (int i = node.begin; i < node.begin + node.size; i++) {
      int ii = tree.permutation[i];
      triangle_type tri = tris[ii];
      T w = tri.area();
      avgPoint += w * tri.center();
      avgNormal += w * tri.normal();
      netCharge += w * faceVals[ii];
      netWeight += w;
    }

    avgPoint /= netWeight;
  };

  auto compute = [&vertex_normals, &vertex_areas, faceVals, regLength](
                     int i_c, const vec3 &wq, const coordinate_type &pc,
                     const coordinate_type &pe, const coordinate_type &N,
                     const vector<triangle_type> &tris, ANode &node,
                     ATree &tree) -> vec3 {
    vec3 out = z::zero<vec3>();

    coordinate_type Nv = vertex_normals[i_c];
    real av = vertex_areas[i_c];

    // real wv = vertex_areas[i_c];

    auto computeTangentRadius = [](coordinate_type dp, coordinate_type N) {
      real dpn = dp.norm();
      real dpnProj = ((N * N.transpose()) * dp).norm();
      return dpn * dpn / dpnProj;
    };

    auto dK = [](const coordinate_type &dp, const coordinate_type &N,
                 const double &p) -> coordinate_type {
      mat3 P = N * N.transpose();
      coordinate_type Pdp = P * dp;
      real ndp = dp.norm();

      if (ndp == 0.0)
        return coordinate_type::Zero();

      real nPdp = (Pdp).norm();
      real k = pow(nPdp, p) / (pow(ndp, 2.0 * p));

      coordinate_type dk0 = P * Pdp / nPdp / nPdp;
      coordinate_type dk1 = 2.0 * dp / ndp / ndp;

      coordinate_type dk = p * k * (dk0 - dk1);
      // std::cout << (dk0 - dk1).transpose() << std::endl;
      return dk;
    };

    if (node.isLeaf()) {
      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        vec3 qi = faceVals[ii];
        auto tri = tris[ii];
        auto w = tri.area();
        auto c = tri.center();
        auto Nf = tri.normal();
        coordinate_type dp = c - pe;
        coordinate_type dkv = dK(dp, Nv, 3.0);
        coordinate_type dkf = dK(-dp, Nf, 3.0);
        out += w * (dkv - dkf);
      }
    } else {
      coordinate_type dp = pc - pe;
      real w = N.norm();
      coordinate_type Nf = N / w;
      coordinate_type dkv = dK(dp, Nv, 6.0);
      coordinate_type dkf = dK(-dp, Nf, 6.0);
      out += w * (dkv - dkf);
    }
    return out;
  };

  vector<face_ptr> &faces = mesh->get_faces();
  vector<coordinate_type> normals;
  vector<triangle_type> triangles;
  for (int i = 0; i < faces.size(); i++) {
    if (!mesh->has_face(i))
      continue;
    if (faces[i]->size() < 3)
      continue;
    std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
    triangles.insert(triangles.end(), tris.begin(), tris.end());
  }

  for (auto t : triangles) {
    normals.push_back(t.normal());
  }
  vector<vec3> u(evalPoints.size(), z::zero<vec3>());
  Avg_Integrator integrator;
  integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);

  return u;
}

template <typename SPACE>
vector<typename SPACE::vec3>
calcFracLaplace(asawa::surf<SPACE> *mesh,
                const std::vector<typename SPACE::coordinate_type> &evalPoints,
                const std::vector<typename SPACE::coordinate_type> &vertvals,
                typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  using Avg_Integrator =
      calder::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;
  std::vector<vec3> vertex_normals = ci::get_vertex_normals<SPACE>(mesh);
  std::vector<real> vertex_areas = ci::get_vertex_areas<SPACE>(mesh);

  std::vector<vec3> faceVals = ci::verts_to_faces<SPACE>(vertvals, mesh);

  if (evalPoints.empty())
    return std::vector<vec3>();

  auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                        ATree &tree, vec3 &netCharge, coordinate_type &avgPoint,
                        coordinate_type &avgNormal) -> void {
    avgPoint = coordinate_type(0, 0, 0);
    netCharge = z::zero<vec3>();

    T netWeight = 0;

    for (int i = node.begin; i < node.begin + node.size; i++) {
      int ii = tree.permutation[i];
      triangle_type tri = tris[ii];
      T w = tri.area();
      avgPoint += w * tri.center();
      avgNormal += w * tri.normal();
      netCharge += w * faceVals[ii];
      netWeight += w;
    }

    avgPoint /= netWeight;
  };

  auto compute = [&vertex_normals, &vertex_areas, faceVals, regLength](
                     int i_c, const vec3 &wq, const coordinate_type &pc,
                     const coordinate_type &pe, const coordinate_type &N,
                     const vector<triangle_type> &tris, ANode &node,
                     ATree &tree) -> vec3 {
    vec3 out = z::zero<vec3>();

    coordinate_type Nv = vertex_normals[i_c];
    real av = vertex_areas[i_c];

    // real wv = vertex_areas[i_c];

    if (node.isLeaf()) {
      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        vec3 qi = faceVals[ii];
        auto tri = tris[ii];
        auto w = tri.area();
        auto c = tri.center();
        auto Nf = tri.normal();
        coordinate_type dp = c - pe;
        coordinate_type dkv = dK(dp, Nv, 3.0);
        coordinate_type dkf = dK(-dp, Nf, 3.0);
        out += w * (dkv - dkf);
      }
    } else {
      coordinate_type dp = pc - pe;
      real w = N.norm();
      coordinate_type Nf = N / w;
      coordinate_type dkv = dK(dp, Nv, 6.0);
      coordinate_type dkf = dK(-dp, Nf, 6.0);
      out += w * (dkv - dkf);
    }
    return out;
  };

  vector<face_ptr> &faces = mesh->get_faces();
  vector<coordinate_type> normals;
  vector<triangle_type> triangles;
  for (int i = 0; i < faces.size(); i++) {
    if (!mesh->has_face(i))
      continue;
    if (faces[i]->size() < 3)
      continue;
    std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
    triangles.insert(triangles.end(), tris.begin(), tris.end());
  }

  for (auto t : triangles) {
    normals.push_back(t.normal());
  }
  vector<vec3> u(evalPoints.size(), z::zero<vec3>());
  Avg_Integrator integrator;
  integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);

  return u;
}

using namespace asawa;
template <typename SPACE> class nearest_point_integrator {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<nearest_point_integrator<SPACE>> ptr;
  static ptr create(real s) {
    return std::make_shared<nearest_point_integrator<SPACE>>(s);
  }

  nearest_point_integrator(real scale) : _scale(scale) {}

  coordinate_array step(surf_ptr surf, const coordinate_array &p0) {
    return update_position(surf, p0);
  }

  coordinate_array update_position(surf_ptr surf, const coordinate_array &p0) {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);
    int i = 0;

    coordinate_array p1(p0);
#if 1

    std::vector<real> kN = calcEnergy(surf, p0, 0.1 * _scale);
    coordinate_array fN = calcGradient(surf, p0, 0.1 * _scale);
    calder::mesh_calculator<SPACE> calc;

    // fN = calc.harmonicAvg(surf, fN, p0, 4.0 * _scale);

    i = 0;
    real C = 0.0;
    real dpd = 0.0;
    real stepsize = 0.0;
    real h = 99999;
    for (auto f : fN) {
      // gg::geometry_logger::line(p0[i], p0[i] + 1e-8 * f,
      //                          vec4(1.0, 1.0, 0.0, 1.0));

      C += kN[i];
      real ftf = f.transpose() * f;
      h = std::min(h, sqrt(_scale * _scale / ftf));
      dpd += ftf;

      if (ftf > 0.0)
        stepsize += (kN[i] / ftf);

      i++;
    }
    stepsize /= real(fN.size());
    // stepsize = C / dpd;
    i = 0;
    std::cout << " stepsize: " << stepsize << " " << C / dpd << " " << h
              << std::endl;
    for (auto f : fN) {
      // gg::geometry_logger::line(p0[i], p0[i] + 1e-8 * f,
      //                           vec4(1.0, 1.0, 0.0, 1.0));
      p1[i] += _scale * stepsize * f;
      i++;
    }
#endif

    return p1;
  }

  real _scale = 0.1;
  std::vector<coordinate_type> _targets;
};

typedef nearest_point_space<double> nearest_point;

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

    asawa::obj_loader<nearest_point> load;
    asawa::subdivide<nearest_point> sub;
    asawa::make<nearest_point> mk;

    asawa::construct<nearest_point> bevel;
    asawa::affine<nearest_point> mod;
    std::string start_frame = "";
    // start_frame = "nearest.test.gaudi";
    start_frame = "stretch.340.gaudi";

    if (!start_frame.empty()) {
      FILE *file;
      file = fopen(start_frame.c_str(), "rb");
      if (file != NULL) {
        this->load_gaudi(start_frame);
      }
    } else {

      _meshGraph = &load("assets/bunny.obj");
      //_meshGraph = &load("assets/close.obj");

      //_meshGraph = &load("assets/icosahedron.obj");
      //_meshGraph = &load("assets/sphere.obj");

      //_meshGraph = &load("assets/messer.obj");
      //_meshGraph = &load("assets/tet.obj");
      // std::cout << "--make cube" << std::endl;
      //_meshGraph = mk.cube(1.0, 1.0, 1.0);
      //_meshGraph = mk.tet();

      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      asawa::remesh<nearest_point> rem;
      rem.triangulate(_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      _meshGraph->update_all();
      _meshGraph->pack();

      mod.centerGeometry(*_meshGraph);
    }

    int N = 0;
    init_phi(_meshGraph);
    _integrator =
        new asawa::surf_integrator<nearest_point>(_meshGraph, 0.5, 2.5, 0.75);
    //_integrator = new asawa::surf_integrator<nearest_point>(_meshGraph,
    // 0.1, 3.0,
    // 0.5);
    _integrator->add_default_vertex_policy<typename nearest_point::real>(
        nearest_point::vertex_index::SMOOTH);
    _max = _integrator->_max;
    _min = _integrator->_min;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  template <typename SPACE> void init_phi(asawa::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    std::vector<edge_ptr> edges = surf->get_edges();
    for (auto e : edges) {
      e->template set<real>(SPACE::edge_index::PHI0, -9999);
    }
  }

  template <typename SPACE>
  vector<asawa::colorRGB> getColor(asawa::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth =
        asawa::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);
    auto hot =
        asawa::ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);

    asawa::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    real mn = K[0];
    real mx = K[0];
    std::for_each(std::begin(K), std::end(K), [&](const double d) {
      mn = std::min(mn, d);
      mx = std::max(mx, d);
    });

    std::vector<asawa::colorRGB> vert_colors(smooth.size());
    int i = 0;
    for (auto v : surf->get_vertices()) {
      typename SPACE::real k = (K[i] - mn) / (mx - mn);
      typename SPACE::real N = 0;
      typename SPACE::real s = smooth[i];
      typename SPACE::real h = hot[i].norm();

      typename SPACE::coordinate_type colorS(1.0, 0.0, 0.0);
      typename SPACE::coordinate_type colorC(0.56, 0.60, 0.40);
      typename SPACE::coordinate_type colorH(0.5, 0.0, 0.5);
      typename SPACE::coordinate_type colorK(0.75, 0.75, 0.0);
      typename SPACE::coordinate_type mx = colorC;
      // mx = va::mix(k, colorK, mx);
      mx = va::mix(s, colorS, mx);
      mx = va::mix(h, colorH, mx);
      // mx = va::mix(k, colorC, mx);

      vert_colors[i] = asawa::colorRGB(mx[0], mx[1], mx[2], 1.0);
      i++;
    }
    return vert_colors;
  }

  virtual void onAnimate(int frame) {

    if (!_grower) {
      _grower = nearest_point_integrator<nearest_point>::create(_max);
    }

    _meshGraph->update_all();
    _meshGraph->reset_flags();

    _meshGraph->pack();
    // if (frame == 1)

    _integrator->integrate();
    //    if (frame > 1)
    //      return;

    std::cout << "frame: " << frame << std::endl;
    std::cout << "====== " << std::endl;
    std::cout << " verts: " << _meshGraph->get_vertices().size() << std::endl;
    std::cout << " edges: " << _meshGraph->get_edges().size() << std::endl;
    std::cout << " faces: " << _meshGraph->get_faces().size() << std::endl;
    std::cout << " mean edge length: "
              << asawa::ci::geometric_mean_length<nearest_point>(_meshGraph)
              << std::endl;

    double dt = _params.dt;
    auto colors = getColor(_meshGraph);
#if 1

    // build constraints to capture current config
    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    using edge_ptr = typename asawa::surf<nearest_point>::edge *;
    using vertex_ptr = typename asawa::surf<nearest_point>::vertex *;

#if 1
    {
      std::cout << "building constraints " << std::endl;
      hepworth::constraint_set<nearest_point>::ptr constraints =
          hepworth::constraint_set<nearest_point>::create(_meshGraph);

      std::cout << "adding constraints " << std::endl;
      init_stretch_constraints<nearest_point>(_meshGraph, constraints, 1.0e-1);
      init_willmore_constraints<nearest_point>(_meshGraph, constraints, 1.0e-7);
      init_mem_bend_constraints<nearest_point>(_meshGraph, constraints, 1.0e-6,
                                               5.0e-6);

      constraints->add_constraint(
          hepworth::internal_collisions<nearest_point>::create(_meshGraph,
                                                               _max));

      std::vector<nearest_point::vec3> p0 =
          asawa::ci::get_coordinates<nearest_point>(_meshGraph);
      std::vector<nearest_point::vec3> p1 = _grower->step(_meshGraph, p0);

      hepworth::velocity_optimizer<nearest_point> opt(p0, p1);
      opt.update(constraints);
      p1 = constraints->get_positions();

      asawa::ci::set_coordinates<nearest_point>(p1, _meshGraph);
    }
#endif

#if 1

    // std::cout << "print vecs" << std::endl;
    std::cout << "rendering debug" << std::endl;
    gg::geometry_logger::render();
#endif

    _meshGraph->print();

    gg::fillBuffer(_meshGraph, _obj, colors);

    if (frame % 10 == 0) {
      save(frame);
      dump_gaudi(frame);
    }

#endif
  }

  virtual void save(int frame) {
    _meshGraph->pack();
    std::stringstream ss;
    ss << "nearest_point." << frame << ".obj";
    asawa::write_obj<nearest_point>(*_meshGraph, ss.str());
  }

  virtual void load_gaudi(std::string file_name) {
    std::cout << " loading" << std::endl;
    asawa::flattened_surf<nearest_point> fsurf;
    fsurf.clear();
    fsurf.read(file_name);
    _meshGraph = fsurf.to_surf();
    //_integrator->set_mesh(_meshGraph);
  }

  virtual void dump_gaudi(int frame = 0) {
    std::cout << " dumping" << std::endl;
    asawa::flattened_surf<nearest_point> fsurf(_meshGraph);
    fsurf.write("nearest_point." + std::to_string(frame) + ".gaudi");
  }

  virtual void onDraw(gg::Viewer &viewer) {

    std::for_each(mSceneObjects.begin(), mSceneObjects.end(),
                  [&](gg::DrawablePtr obj) mutable {
                    if (obj->isVisible)
                      obj->draw(viewer.getProjection(), viewer.getModelView());
                  });

    gg::geometry_logger::clear();
  }

  struct {
    double dt = 0.01;
  } _params;

private:
  double _max = 0.0;
  double _min = 0.0;

  std::vector<gg::DrawablePtr> mSceneObjects;

  gg::BufferObjectPtr _obj = NULL;
  asawa::surf<nearest_point> *_meshGraph;
  asawa::surf_integrator<nearest_point> *_integrator;
  nearest_point_integrator<nearest_point>::ptr _grower;
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
