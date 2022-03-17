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
    double l0 = m2::ci::geometric_mean_length<rxd3>(_meshGraph);
    l0 *= 0.5;

    _min = 0.5 * l0;
    _max = 3.0 * _min;

    std::cout << "avg length: " << l0 << std::endl;
    std::cout << "--int rx" << std::endl;

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

  template <typename SPACE> void cacheBary(m2::surf<rxd3> *surf) {
    M2_TYPEDEFS;
    int i = 0;
    for (auto f : surf->get_faces()) {

      if (!surf->has_face(i++))
        continue;

      if (f->size() < 3) {
        std::cout << " bary_size: " << std::endl;
        std::cout << "   face: " << f->size() << std::endl;
        std::cout << "   vert: " << f->fbegin()->vertex()->size() << std::endl;
        std::cout << "   edge: " << f->fbegin()->edge() << std::endl;
        f->print();
        f->fbegin()->vertex()->print();
        continue;
      }

      typename SPACE::coordinate_type c = m2::ci::center<SPACE>(f);
      typename SPACE::coordinate_type l = m2::ci::point_to_bary<SPACE>(f, c);

      assert(!isnan(l[0]));
      assert(!isnan(l[1]));
      assert(!isnan(l[2]));

      int i = 0;
      m2::for_each_face<SPACE>(f, [l, i](face_vertex_ptr fv) mutable {
        typename SPACE::real li = 0.0;
        if (i < 3)
          li = l[i];
        fv->template set<typename SPACE::real>(SPACE::face_vertex_index::BARY,
                                               li);
        i++;
      });
    }
  }

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
  void smoothMesh(m2::surf<SPACE> *surf, typename SPACE::real C, int N) {
    using namespace m2;
    M2_TYPEDEFS;
    m2::area_laplacian<SPACE, coordinate_type> M(surf);
    int i = 0;

    for (int k = 0; k < N; k++) {
      coordinate_array coords = m2::ci::get_coordinates<SPACE>(surf);
      coordinate_array normals = m2::ci::get_vertex_normals<SPACE>(surf);
      coordinate_array ncoords = M.mult(coords);

      for (int i = 0; i < coords.size(); i++) {
        coordinate_type cp = va::orthogonal_project(normals[i], ncoords[i]);

        if (isnan(cp[0])) {
          std::cout << "cp: " << cp.transpose() << std::endl;
          std::cout << "no: " << normals[i].transpose() << std::endl;
          std::cout << "nc: " << ncoords[i].transpose() << std::endl;
          std::cout << " c: " << coords[i].transpose() << std::endl;
        }

        coords[i] = coords[i] + C * cp;
        assert(!isnan(coords[i][0]));
        assert(!isnan(coords[i][1]));
        assert(!isnan(coords[i][2]));
      }
      m2::ci::set_coordinates<SPACE>(coords, surf);
    }
  }

  template <typename SPACE>
  void splitEdges(m2::surf<SPACE> *surf,
                  const std::vector<typename SPACE::real> &rxA,
                  const std::vector<typename SPACE::real> &rxB) {

    using namespace m2;
    M2_TYPEDEFS;

    m2::edge_splitter<SPACE> splitter(surf, _max);

    splitter.reset_flags();

    edge_array edgesToSplit = splitter.get_edges_to_split();
    std::cout << "calculating new vals" << std::endl;

    std::vector<real> rxAv(edgesToSplit.size());
    std::vector<real> rxBv(edgesToSplit.size());
    std::vector<coordinate_type> coords(edgesToSplit.size());

    int i = 0;
    for (auto e : edgesToSplit) {

      vertex_ptr v0 = e->v1()->vertex();
      vertex_ptr v1 = e->v2()->vertex();

      real rA0 = v0->template get<real>(SPACE::vertex_index::RXA);
      real rB0 = v1->template get<real>(SPACE::vertex_index::RXB);
      real rA1 = v0->template get<real>(SPACE::vertex_index::RXA);
      real rB1 = v1->template get<real>(SPACE::vertex_index::RXB);

      rxAv[i] = 0.5 * rA0 + 0.5 * rA1;
      rxBv[i] = 0.5 * rB0 + 0.5 * rB1;
      coordinate_type c0 = m2::ci::get_coordinate<SPACE>(v0);
      coordinate_type c1 = m2::ci::get_coordinate<SPACE>(v1);

      coords[i] = 0.5 * c0 + 0.5 * c1;

      i++;
    }

    std::cout << "splitting: " << edgesToSplit.size() << " edges" << std::endl;
    splitter.split_edges(edgesToSplit);

    vertex_array vertices = surf->get_vertices();

    std::cout << "assigning calc'd vertex values" << std::endl;

    for (auto v : vertices) {
      if (!v)
        continue;
      int topologyChangeId = v->topologyChangeId;
      if (topologyChangeId < 0)
        continue;
      // std::cout << "split:" << topologyChangeId << " " <<
      // v->position_in_set()
      // << std::endl;
      v->template set<typename rxd3::real>(rxd3::vertex_index::RXA,
                                           rxAv[topologyChangeId]);
      v->template set<typename rxd3::real>(rxd3::vertex_index::RXB,
                                           rxBv[topologyChangeId]);
      m2::ci::set_coordinate<SPACE>(coords[topologyChangeId], v);
      v->topologyChangeId = -1;
    }
  }

  template <typename SPACE>
  void collapsEdges(m2::surf<SPACE> *surf,
                    const std::vector<typename SPACE::real> &rxA,
                    const std::vector<typename SPACE::real> &rxB) {

    using namespace m2;
    M2_TYPEDEFS;

    m2::edge_collapser<SPACE> collapser(surf, _min);

    collapser.reset_flags();

    edge_array edgeToCollapse = collapser.get_edges_to_collapse();
    std::cout << "calculating new vals" << std::endl;

    std::vector<real> rxAv(edgeToCollapse.size());
    std::vector<real> rxBv(edgeToCollapse.size());
    std::vector<coordinate_type> coords(edgeToCollapse.size());

    int i = 0;
    for (auto e : edgeToCollapse) {

      vertex_ptr v0 = e->v1()->vertex();
      vertex_ptr v1 = e->v2()->vertex();
      int i0 = v0->position_in_set();
      int i1 = v1->position_in_set();

      real rA0 = v0->template get<real>(SPACE::vertex_index::RXA);
      real rB0 = v1->template get<real>(SPACE::vertex_index::RXB);
      real rA1 = v0->template get<real>(SPACE::vertex_index::RXA);
      real rB1 = v1->template get<real>(SPACE::vertex_index::RXB);
      rxAv[i] = 0.5 * rA0 + 0.5 * rA1;
      rxBv[i] = 0.5 * rB0 + 0.5 * rB1;
      coordinate_type c0 = m2::ci::get_coordinate<SPACE>(v0);
      coordinate_type c1 = m2::ci::get_coordinate<SPACE>(v1);

      coords[i] = 0.5 * c0 + 0.5 * c1;

      i++;
    }

    std::cout << "collapsing: " << edgeToCollapse.size() << " edges"
              << std::endl;
    collapser.collapse_edges(edgeToCollapse);

    vertex_array vertices = surf->get_vertices();

    std::cout << "assigning calc'd vertex values" << std::endl;

    for (auto v : vertices) {
      if (!v)
        continue;
      int topologyChangeId = v->topologyChangeId;
      if (topologyChangeId < 0)
        continue;
      // std::cout << "split:" << topologyChangeId << " " <<
      // v->position_in_set()
      // << std::endl;
      v->template set<typename rxd3::real>(rxd3::vertex_index::RXA,
                                           rxAv[topologyChangeId]);
      v->template set<typename rxd3::real>(rxd3::vertex_index::RXB,
                                           rxBv[topologyChangeId]);
      m2::ci::set_coordinate<SPACE>(coords[topologyChangeId], v);
    }
  }

  template <typename SPACE>
  void mergeTris(m2::surf<SPACE> *surf,
                 const std::vector<typename SPACE::real> &rxA,
                 const std::vector<typename SPACE::real> &rxB) {
    using namespace m2;
    M2_TYPEDEFS;
    face_merger<SPACE> merger(surf, 2.0 * _min);
    bool merging = true;

    merging = merger.merge();
  }

  template <typename SPACE>
  void mergeEdges(m2::surf<SPACE> *surf,
                  const std::vector<typename SPACE::real> &rxA,
                  const std::vector<typename SPACE::real> &rxB) {
    using namespace m2;
    M2_TYPEDEFS;
    edge_merger<SPACE> merger(surf, 10.0*_min);
    bool merging = true;

    merging = merger.merge();
  }

  template <typename SPACE>
  void flipEdges(m2::surf<SPACE> *surf,
                 const std::vector<typename SPACE::real> &rxA,
                 const std::vector<typename SPACE::real> &rxB) {
    using namespace m2;
    M2_TYPEDEFS;
    edge_flipper<SPACE> flipper(surf);
    flipper.flip();
  }

  template <typename SPACE>
  void delete_degenerates(m2::surf<SPACE> *surf,
                          const std::vector<typename SPACE::real> &rxA,
                          const std::vector<typename SPACE::real> &rxB) {
    using namespace m2;
    M2_TYPEDEFS;
    m2::construct<SPACE> cons;

    bool testing = true;
    int i;

    while (testing) {
      i = 0;
      bool deleting = false;
      for (auto e : _meshGraph->get_edges()) {
        if (!_meshGraph->has_edge(i++))
          continue;
        // std::cout << i << std::endl;
        deleting |= cons.delete_degenerates(_meshGraph, e);
      }
      testing = deleting;

      i = 0;
      deleting = false;
      for (auto v : _meshGraph->get_vertices()) {
        if (!_meshGraph->has_vertex(i++))
          continue;
        deleting = cons.delete_degenerates(_meshGraph, v);
      }
      testing != deleting;
    }
    
    i = 0;
    for (auto f : _meshGraph->get_faces()) {
      if (!_meshGraph->has_face(i++))
        continue;
      face_vertex_ptr fv = f->get_front();
      if (f->is_null() && fv->vertex()->size() > 1) {
        surf->remove_face(f->position_in_set());
        delete (fv);
      }
    }

    m2::remesh<SPACE> rem;
    rem.triangulate(_meshGraph);
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  diffuse(m2::surf<SPACE> *surf, const std::vector<typename SPACE::real> &u,
          const std::vector<typename SPACE::real> &f, rxd3::real dt,
          typename SPACE::real C) {

    using namespace m2;
    M2_TYPEDEFS;
    m2::laplacian<SPACE, real> M(surf);

    auto diffMult = [&M, dt, C, surf](const std::vector<real> &X) {
      std::vector<real> MX = M.multM(X);
      std::vector<real> CX = M.multC(X);

      int i = 0;
      for (int i = 0; i < MX.size(); i++) {
        MX[i] = MX[i] - dt * C * CX[i];
      }
      return MX;
    };

    std::vector<typename SPACE::real> u_f(u);

    for (int i = 0; i < u_f.size(); i++) {
      u_f[i] = u[i] + dt * f[i];
    }

    std::vector<typename SPACE::real> au_f = M.multM(u_f);

    int its;
    m2::gradient_descent<SPACE> solver;
    std::vector<real> x = solver.solveConjugateGradient(au_f, diffMult, its);

    return x;
  }

  template <typename SPACE>
  std::vector<typename SPACE::real>
  diffuse_second_order(m2::surf<SPACE> *surf,
                       const std::vector<typename SPACE::real> &u,
                       const std::vector<typename SPACE::real> &f,
                       rxd3::real dt, typename SPACE::real C) {

    using namespace m2;
    M2_TYPEDEFS;
    m2::laplacian<SPACE, real> M(surf);

    auto diffMult = [&M, dt, C, surf](const std::vector<real> &X) {
      std::vector<real> MX = M.multM(X);
      std::vector<real> CX = M.multC(X);

      int i = 0;
      for (int i = 0; i < MX.size(); i++) {
        MX[i] = MX[i] - 0.5 * dt * C * CX[i];
      }
      return MX;
    };

    std::vector<typename SPACE::real> u_f(u);

    for (int i = 0; i < u_f.size(); i++) {
      u_f[i] = u[i] + dt * f[i];
    }

    std::vector<typename SPACE::real> au_f = M.multM(u_f);
    std::vector<typename SPACE::real> lu = M.multC(u);

    for (int i = 0; i < u_f.size(); i++) {
      au_f[i] += 0.5 * dt * C * lu[i];
    }

    int its;
    m2::gradient_descent<SPACE> solver;
    std::vector<real> x = solver.solveConjugateGradient(au_f, diffMult, its);

    return x;
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
      if (dist(e2) > 0.97 && K[i] < 0.1 && rxA[i] > 0.97) {
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
#if 0
    int i = 0;
    coordinate_array coords;
    for (auto v : vertices) {
      if (dist(e2) > 0.80) {
        coords.push_back(ci::get_coordinate<SPACE>(v));
      }
      rxA[i] = 0.0;
      rxB[i] = 0.0;
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
        i++;
      }
    }
    double mx = 0.0;
    for (int i = 0; i < rxB.size(); i++) {
      mx = std::max(mx, rxA[i]);
      mx = std::max(mx, rxB[i]);
    }
    std::cout << " max: " << mx << std::endl;
    for (int i = 0; i < rxB.size(); i++) {
      rxB[i] /= mx;
      rxA[i] += (1.0 - rxB[i]);
    }

    mx = 0.0;
    for (int i = 0; i < rxB.size(); i++) {
      mx = std::max(mx, rxA[i]);
      mx = std::max(mx, rxB[i]);
    }
    std::cout << " max: " << mx << std::endl;

#else
    int i = 0;
    coordinate_array coords;
    for (auto v : vertices) {
      rxA[i] = 1.0;
      rxB[i] = 0.0;
      /*
      if (dist(e2) > 0.85) {
        coords.push_back(ci::get_coordinate<SPACE>(v));
        typename SPACE::real t = dist(e2);
        rxA[i] = 0.0;
        rxB[i] = 1.0;
      } else {
        rxA[i] = 1.0;
        rxB[i] = 0.0;
      }
      */
      i++;
    }
/*
    std::vector<typename SPACE::real> rxAn(rxA);
    std::vector<typename SPACE::real> rxBn(rxB);
    for (int k = 0; k < 3; k++) {

      for (auto v : vertices) {
        typename SPACE::real mx = 0.0;
        m2::for_each_vertex<SPACE>(v, [&mx, &rxB](face_vertex_ptr fv) mutable {
          int j = fv->next()->vertex()->position_in_set();
          mx = max(rxB[j], mx);
        });

        int i = v->position_in_set();
        rxBn[i] = mx;
        rxAn[i] = 1.0 - mx;
      }
      std::swap(rxA, rxAn);
      std::swap(rxB, rxBn);
    }
*/
#endif

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
    cacheBary<rxd3>(_meshGraph);

    std::cout << "ugly smoothing " << std::endl;
    smoothMesh(_meshGraph, 0.1, 10);
    std::cout << "done " << std::endl;
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

    rxA = diffuse_second_order<rxd3>(_meshGraph, rxA, fA, dt, Du);
    rxB = diffuse_second_order<rxd3>(_meshGraph, rxB, fB, dt, Dv);

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
    setRxColor(_meshGraph);

    auto colors = getRxColors(_meshGraph);
    // mod.centerGeometry(*_meshGraph);
    gg::fillBuffer(_meshGraph, _obj, colors);
#endif

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

    // mergeTris(_meshGraph, rxA, rxB);
    mergeEdges(_meshGraph, rxA, rxB);

    delete_degenerates(_meshGraph, rxA, rxB);
    splitEdges(_meshGraph, rxA, rxB);
    delete_degenerates(_meshGraph, rxA, rxB);
    collapsEdges(_meshGraph, rxA, rxB);
    delete_degenerates(_meshGraph, rxA, rxB);
    flipEdges(_meshGraph, rxA, rxB);
    delete_degenerates(_meshGraph, rxA, rxB);
    // get rid of some of this degenerate handling :P
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
    double k = 0.061, F = 0.074;
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
