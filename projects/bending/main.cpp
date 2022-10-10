//#include "nanoguiincludes.h"

#include "manifold/coordinate_interface.hpp"

#include "manifold/hepworth/constraints_init.hpp"
#include "manifold/hepworth/objective_function.hpp"

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
#include "manifold/m2Includes.h"
#include "manifold/m2Operators.h"
#include "manifold/make.hpp"
#include "manifold/objloader.hpp"

#include "manifold/hepworth/optimizer.hpp"

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

template <typename SPACE>
class action_policy
    : public m2::vertex_policy_t<SPACE, typename SPACE::coordinate_type> {
public:
  M2_TYPEDEFS;
  action_policy(typename SPACE::vertex_index id)
      : vertex_policy_t<SPACE, coordinate_type>(id) {
    std::cout << " id: " << int(this->_id) << std::endl;
  }
  virtual void calc(int i, edge_ptr &e, m2::op_type op) {

    face_vertex_ptr v0 = e->v1();
    face_vertex_ptr v1 = v0->prev();
    face_vertex_ptr v2 = e->v2();
    face_vertex_ptr v3 = v2->prev();

    coordinate_type t0 = v0->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t1 = v1->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t2 = v2->vertex()->template get<coordinate_type>(this->_id);
    coordinate_type t3 = v3->vertex()->template get<coordinate_type>(this->_id);

    if (op == m2::op_type::split) {
      this->_vals[i] =
          coordinate_type::Zero(); // 0.0 * (t0 + t2) - 0.0 * (t1 + t3);
      //_vals[i] = 0.0;
    } else {
      coordinate_type mx = t0.norm() > t2.norm() ? t0 : t2;
      this->_vals[i] = coordinate_type::Zero();
    }

    //_vals[i] = -0.25 * (t0 + t2);
    return;
  }
};

template <typename SPACE>
class edge_init_policy : public m2::edge_policy_t<SPACE, typename SPACE::real> {
public:
  M2_TYPEDEFS;
  edge_init_policy(typename SPACE::edge_index id)
      : edge_policy_t<SPACE, real>(id) {}

  virtual void calc(int i, edge_ptr &e, m2::op_type op) {
    this->_vals[4 * i + 0] = -9999;
    this->_vals[4 * i + 1] = -9999;
    this->_vals[4 * i + 2] = -9999;
    this->_vals[4 * i + 3] = -9999;
    return;
  }
};

////////////////////////////////////////////////////////////////////////////
// AVG
////////////////////////////////////////////////////////////////////////////

template <typename SPACE>
vector<typename SPACE::vec3>
calcPotential(m2::surf<SPACE> *mesh, std::vector<typename SPACE::vec3> vertVals,
              std::vector<typename SPACE::coordinate_type> evalPoints,
              typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  using Avg_Integrator =
      m2::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;

  std::vector<vec3> faceVals = ci::verts_to_faces<SPACE>(vertVals, mesh);

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

  auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
    T dist3 = dist * dist * dist;
    T l3 = C * C * C;
    T kappa = 1.0 / (dist3 + l3);
    return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
    T dist3 = C / dist / dist / dist;
    T dist6 = dist * dist * dist * dist * dist * dist;

    T l2 = C * C;
    T kappa = 1.0 / (dist3 + dist6);
    return kappa / 4.0 / M_PI;
#elif 1
    T kappa = 1.0 / (dist + C);
    return kappa / 4.0 / M_PI;
#elif 0
    T dist2 = dist * dist;
    T dt = 0.5;
    T kappa = exp(-dist2 / 4.0 / dt);
    return kappa / pow(4.0 * M_PI * dt, 1.5);
#endif
  };
  coordinate_array vertex_normals = ci::get_vertex_normals<SPACE>(mesh);

  auto compute = [&vertex_normals, &vertVals, faceVals, regLength,
                  computeK](int i_c, const vec3 &wq, const coordinate_type &pc,
                            const coordinate_type &pe, const coordinate_type &N,
                            const vector<triangle_type> &tris, ANode &node,
                            ATree &tree) -> vec3 {
    vec3 out = z::zero<vec3>();

    coordinate_type vn = vertex_normals[i_c];
    coordinate_type qc = vertVals[i_c];
    real qcmag = qc.norm();

    auto computeKelvin = [](coordinate_type r, coordinate_type f, real eps) {
      real mu = 10.0;
      real v = 2e-1;
      real a = 1.0 / 4.0 / mu;
      real b = a / (4.0 * (1.0 - v));
      real eps2 = 1e-10 * eps * eps;
      real re = 1.0 / sqrt(r.dot(r) + eps2);
      real re3 = re * re * re;

      mat3 R = r * r.transpose();
      mat3 I;
      I.setIdentity();

      coordinate_type u =
          ((a - b) / re * I + b / re3 * R + 0.5 * a * eps2 / re3 * I) * f;

      return u;
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

        real m = qi.norm();
        // out += computeKelvin(dp, w * qi, regLength);

        float sg = va::sgn(qi.dot(Nf));
        T dist = m2::va::norm(dp);
        T k = computeK(dist, regLength);
        dp /= dist;
        coordinate_type u = va::cross(dp, Nf);
        // out += w * k * qcmag * u;
        out += w * k * (m * sg * Nf + 0.5 * qi);
      }
    } else {
      // out += computeK(dist, regLength) * va::cross(q, dp);
      coordinate_type dp = pc - pe;

      // out += computeKelvin(dp, wq, regLength);

      T dist = m2::va::norm(dp);
      T k = computeK(dist, regLength);
      dp /= dist;
      coordinate_type u = va::cross(dp, N.normalized());

      real m = wq.norm();
      float sg = va::sgn(wq.dot(N));
      // out += k * qcmag * u;
      out += k * (m * sg * N + 0.5 * wq);
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
    std::vector<triangle_type> tris = m2::ci::get_tris<SPACE>(faces[i]);
    triangles.insert(triangles.end(), tris.begin(), tris.end());
  }

  for (auto t : triangles) {
    normals.push_back(t.normal());
  }
  vector<vec3> u(evalPoints.size(), z::zero<vec3>());
  Avg_Integrator integrator;
  integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);
  int i = 0;

  return u;
}

using namespace m2;
template <typename SPACE> class dendritic_growth {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<dendritic_growth<SPACE>> ptr;
  static ptr create(real s) {
    return std::make_shared<dendritic_growth<SPACE>>(s);
  }

  dendritic_growth(real scale) : _scale(scale) { init_targets(); }

  void init_velocity(surf_ptr surf) {
    auto &verts = surf->get_vertices();
    coordinate_array act(verts.size(), coordinate_type::Zero());
    ci::set<SPACE, coordinate_type>(surf, act, SPACE::vertex_index::ACT);
  }

  void init_targets() {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    for (int i = 0; i < 32; i++) {
      double c0 = unif(rng);
      double c1 = unif(rng);
      double c2 = unif(rng);
      coordinate_type c(c0, c1, c2);
      _targets.push_back(2.0 * c);
    }
  }

  void calc_stdev(const std::vector<real> &K, real &mean, real &stdev) {

    real sum = std::accumulate(std::begin(K), std::end(K), 0.0);
    mean = sum / K.size();
    real accum = 0.0;
    std::for_each(std::begin(K), std::end(K),
                  [&](const double d) { accum += (d - mean) * (d - mean); });
    stdev = sqrt(accum / (K.size() - 1));

    std::cout << " K mean: " << mean << " stdev: " << stdev << std::endl;
  }

  void frame(vertex_ptr v, coordinate_type &s, mat3 &U) {

    U = mat3::Zero();
    coordinate_type p0 = ci::get_coordinate<SPACE>(v);
    for_each_vertex<SPACE>(v, [&U, p0](face_vertex_ptr fv) {
      coordinate_type p1 = ci::get_coordinate<SPACE>(fv->next()->vertex());
      coordinate_type N = ci::normal<SPACE>(fv->next()->face());
      real A = ci::area<SPACE>(fv->next()->face());

      coordinate_type dp = A * N;
      U += dp * dp.transpose();
    });

    Eigen::JacobiSVD<mat3> svd(U, Eigen::ComputeFullU);
    U = svd.matrixU();
    s = svd.singularValues();
  }

  real angle(const coordinate_type &a, const coordinate_type &b) {
    coordinate_type an = a.normalized();
    coordinate_type bn = b.normalized();
    coordinate_type x = an.cross(bn);
    return atan2(x.norm(), an.dot(bn));
  }

  bool inCone(const coordinate_type &p, const coordinate_type &heading,
              real thet, real r, const std::vector<coordinate_type> &targets) {
    auto H = heading.normalized();
    real dmin = 99999;
    bool inCone = false;
    for (int j = 0; j < targets.size(); j++) {
      coordinate_type dp = targets[j] - p;
      double d = dp.norm();
      dp /= d;
      real ti = angle(dp, H);
      if (ti < thet && d < r) {
        inCone = true;
        gg::geometry_logger::line(p, p + r * heading, vec4(0.0, 1.0, 1.0, 1.0));
      }
    }
    return inCone;
  }

  coordinate_array step(surf_ptr surf, const coordinate_array &p0) {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif01(0.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = m2::ci::get_vertex_normals<SPACE>(surf);
    coordinate_array act =
        ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);

    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    // std::for_each(K.begin(), K.end(), [](auto &k) { k = 1.0 / k; });
    real mean, stdev;
    calc_stdev(K, mean, stdev);

    real N = double(verts.size());
    real d0 = 0.0;
    real dp = 0.0;
    real dn = 0.0;

    int i = 0;
    for (auto &v : verts) {
      coordinate_type av = act[i];
      real a = av.norm();
      d0 += fabs(a) / N;
      i++;
    }

    std::cout << "current density: " << d0 << std::endl;
    double max_d = 0.005;
    double prob_d = 0.25 / N / max_d;
    double d1 = d0;
    std::cout << "max density: " << max_d << ", prob density: " << prob_d
              << std::endl;

    i = 0;
    for (auto &v : verts) {
      coordinate_type av = act[i];
      real a = av.norm();

      if (unif01(rng) < prob_d && d1 < max_d && a < 1e-4) {
        // if (K[i] > mean + 0.0 * stdev) {
        if (unif01(rng) < 0.65) {
          av = normals[i];
          a = 1.0;
        } else /*if (K[i] < mean - 1.0 * stdev)*/ {
          a = 1.0;
          av = -normals[i];
        }
      }

      mat3 F;
      coordinate_type s;
      frame(v, s, F);

      if (fabs(a) > 0) {
        // if (unif01(rng) > 0.99)
        //   a = 0.0;
        // double theta = M_PI / 6.0 * unif01(rng);
        double theta = M_PI / 24.0;

        real sg = va::sgn(av.dot(normals[i]));

        bool sees_target =
            inCone(p0[i], av, 2.0 * theta, 12.0 * _scale, _targets);

        a *= (1.0 - 0.01 * unif01(rng)); // cool a down

        if (!sees_target) {
          coordinate_type uu = av.normalized();
          coordinate_type vv = uu.cross(normals[i]);
          coordinate_type ww = uu.cross(vv);
          uu.normalize();
          vv.normalize();
          ww.normalize();

          vv = F.block(0, 1, 3, 1);
          ww = F.block(0, 2, 3, 1);
          // coordinate_type v = F.block(0, 2, 3, 1);

          if (unif01(rng) > 0.5)
            theta *= -1;

          if (unif01(rng) > 0.5)
            av = Eigen::AngleAxisd(theta, vv) * av;
          else {
            av = Eigen::AngleAxisd(theta, ww) * av;
          }
        }
      }

      d1 += a / N;
      act[i] = a * av.normalized();
      // std::cout << s.transpose() << std::endl;
      real cov = s[0] / (s[0] + s[1] + s[2]);
      cov = 1.0 - 1.0 / cov;
      // std::cout << cov << std::endl;
      // act[i] += 1e-2 * _scale * normals[i];
      // gg::geometry_logger::line(p0[i], p0[i] + 0.1 * _scale * cov *
      // normals[i],
      //                           vec4(0.0, 1.0, 1.0, 1.0));

      i++;
    }
#if 1
    ///////////////
    // diffusion
    ///////////////
    std::vector<real> f(act.size());
    i = 0;
    for (auto a : act) {
      f[i] = act[i].norm();
      i++;
    }

    m2::laplacian<SPACE, real> lap(surf);

    double dt = 10.0 * pow(_scale, 2);
    std::vector<real> u = lap.diffuse(f, dt);
#if 0
    i = 0;
    for (auto xi : u) {
      xi = log(xi);
      gg::geometry_logger::line(p0[i], p0[i] + 0.01 * xi * normals[i],
                                vec4(1.0, 0.0, 0.0, 1.0));
      i++;
    }
#endif
    ///////////////
    // gradient
    ///////////////
    face_array &faces = surf->get_faces();
    edge_array &edges = surf->get_edges();

    coordinate_array gradU(faces.size(), coordinate_type::Zero());
    for (auto e : edges) {
      coordinate_type c0 = ci::get_coordinate<SPACE>(e->v1());
      coordinate_type c1 = ci::get_coordinate<SPACE>(e->v2());
      real u0 = u[e->v1()->prev()->vertex()->position_in_set()];
      real u1 = u[e->v2()->prev()->vertex()->position_in_set()];

      real A0 = m2::ci::area<SPACE>(e->v1()->face());
      real A1 = m2::ci::area<SPACE>(e->v2()->face());

      coordinate_type N0 = m2::ci::normal<SPACE>(e->v1()->face());
      coordinate_type N1 = m2::ci::normal<SPACE>(e->v2()->face());
      coordinate_type dp = c1 - c0;

      coordinate_type M0 = dp.cross(N0);
      coordinate_type M1 = dp.cross(N1);

      gradU[e->v1()->face()->position_in_set()] += M0 * u0 / 2.0 / A0;
      gradU[e->v2()->face()->position_in_set()] -= M1 * u1 / 2.0 / A1;
    }

    i = 0;
    for (auto f : faces) {
      typename SPACE::coordinate_type c = m2::ci::center<SPACE>(f);
      gradU[i].normalize();
      gg::geometry_logger::line(c, c + 0.25 * _scale * gradU[i],
                                vec4(0.0, 1.0, 1.0, 1.0));
      i++;
    }

    ///////////////
    // divergence
    ///////////////
    std::vector<real> divu(act.size(), 0.0);

    for (auto e : edges) {
      coordinate_type c0 = ci::get_coordinate<SPACE>(e->v1()->next());
      coordinate_type c1 = ci::get_coordinate<SPACE>(e->v2()->next());

      coordinate_type g0 = gradU[e->v1()->face()->position_in_set()];
      coordinate_type g1 = gradU[e->v2()->face()->position_in_set()];
      coordinate_type dp = c1 - c0;

      real cot0 = ci::abs_cotan<SPACE>(e->v1()->prev());
      real cot1 = ci::abs_cotan<SPACE>(e->v2()->prev());
      divu[e->v1()->vertex()->position_in_set()] += 0.5 * cot0 * dp.dot(g0);
      divu[e->v1()->vertex()->position_in_set()] += 0.5 * cot1 * dp.dot(g1);
      divu[e->v2()->vertex()->position_in_set()] -= 0.5 * cot0 * dp.dot(g0);
      divu[e->v2()->vertex()->position_in_set()] -= 0.5 * cot1 * dp.dot(g1);
    }
#if 0
    i = 0;
    for (auto xi : divu) {
      gg::geometry_logger::line(p0[i], p0[i] + 1.0 * xi * normals[i],
                                vec4(1.0, 0.0, 0.0, 1.0));
      i++;
    }
#endif

    std::vector<real> x = lap.solve(divu);
    // std::vector<real> x = lap.diffuse(divu, dt);
#if 1
    i = 0;
    for (auto xi : x) {
      gg::geometry_logger::line(p0[i], p0[i] + 10.0 * xi * normals[i],
                                vec4(1.0, 0.0, 0.0, 1.0));
      i++;
    }
#endif
#endif
    ci::set<SPACE, coordinate_type>(surf, act, SPACE::vertex_index::ACT);
    render_targets();
    return update_position(surf, p0);
  }

  coordinate_array update_position(surf_ptr surf, const coordinate_array &p0) {
    std::mt19937_64 rng;
    std::uniform_real_distribution<real> unif(-1.0, 1.0);

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = m2::ci::get_vertex_normals<SPACE>(surf);
    coordinate_array act =
        ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);

    coordinate_array p1(p0);
    int i = 0;
    double C_N = 0.25 * _scale;

    for (auto v : verts) {
      // auto f = va::sgn(wN[i].dot(normals[i])) * wN[i].norm();
      coordinate_type vi = act[i];
      real sg = va::sgn(vi.dot(normals[i]));

      coordinate_type Ni = sg * normals[i];

      real vm = vi.norm();
      coordinate_type B = Ni.cross(vi.normalized());
      // double thet = atan2(vv.norm(), uu.dot(Ni));

      real thet = angle(vi, Ni);
      Eigen::AngleAxisd R = Eigen::AngleAxisd(thet, B);

      for_each_vertex<SPACE>(
          v, [i, &p0, &p1, &vm, sg, R, &unif, &rng, &normals,
              C_N](typename m2::surf<SPACE>::face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            coordinate_type Nj = sg * normals[j];
            Nj = R * Nj;
            p1[j] += 0.1 * C_N * unif(rng) * p1[j];
            p1[j] += 4.0 * C_N * vm * Nj;
            /*
                        gg::geometry_logger::line(p0[j], p0[j] + 0.1 * Nj,
                                                  vec4(1.0, 0.0, 1.0, 1.0));
            */
          });
      /*
            gg::geometry_logger::line(p1[i], p1[i] + 0.1 * vi,
                                      vec4(0.0, 1.0, 1.0, 1.0));
      */

      p1[i] += 0.1 * C_N * unif(rng) * p1[i];
      p1[i] -= 0.5 * C_N * act[i];

      i++;
    }
    return p1;
  }

  void render_targets() {

    // std::cout << "print vecs" << std::endl;
    for (int j = 0; j < _targets.size() - 1; j++) {
      auto p0 = _targets[j];
      auto p1 = _targets[j + 1];
      gg::geometry_logger::line(p0, p1, vec4(1.0, 0.0, 0.0, 1.0));
    }
    std::cout << "rendering debug" << std::endl;
  }
  real _scale = 0.1;
  std::vector<coordinate_type> _targets;
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
    // start_frame = "stretch.1170.gaudi";
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

      m2::remesh<stretch> rem;
      rem.triangulate(_meshGraph);
      //_meshGraph = &sub.subdivide_control(*_meshGraph);

      _meshGraph->update_all();
      _meshGraph->pack();

      mod.centerGeometry(*_meshGraph);
    }

    int N = 0;
    init_phi(_meshGraph);
    _integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.5, 2.5, 0.75);
    //_integrator = new m2::surf_integrator<stretch>(_meshGraph, 0.1, 3.0,
    // 0.5);
    _integrator->add_default_vertex_policy<typename stretch::real>(
        stretch::vertex_index::SMOOTH);

    _integrator->add_vertex_policy(
        new action_policy<stretch>(stretch::vertex_index::ACT));
    _integrator->add_edge_policy(
        new edge_init_policy<stretch>(stretch::edge_index::PHI0));
    _max = _integrator->_max;
    _min = _integrator->_min;

    std::cout << "--init rx" << std::endl;

    std::cout << "creating buffer" << std::endl;

    _obj = gg::BufferObject::create();
    _obj->init();
    mSceneObjects.push_back(_obj);

    mSceneObjects.push_back(gg::geometry_logger::get_instance().debugLines);
  }

  template <typename SPACE> void init_phi(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    std::vector<edge_ptr> edges = surf->get_edges();
    for (auto e : edges) {
      e->template set<real>(SPACE::edge_index::PHI0, -9999);
    }
  }

  template <typename SPACE>
  vector<m2::colorRGB> getColor(m2::surf<SPACE> *surf) {
    M2_TYPEDEFS;
    auto smooth = m2::ci::get<SPACE, real>(surf, SPACE::vertex_index::SMOOTH);
    auto hot =
        m2::ci::get<SPACE, coordinate_type>(surf, SPACE::vertex_index::ACT);

    m2::cotan_curvature<SPACE> curve(surf);
    std::vector<real> K = curve();
    real mn = K[0];
    real mx = K[0];
    std::for_each(std::begin(K), std::end(K), [&](const double d) {
      mn = std::min(mn, d);
      mx = std::max(mx, d);
    });

    std::vector<m2::colorRGB> vert_colors(smooth.size());
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
      // mx = m2::va::mix(k, colorK, mx);
      mx = m2::va::mix(s, colorS, mx);
      mx = m2::va::mix(h, colorH, mx);
      // mx = m2::va::mix(k, colorC, mx);

      vert_colors[i] = m2::colorRGB(mx[0], mx[1], mx[2], 1.0);
      i++;
    }
    return vert_colors;
  }

  virtual void onAnimate(int frame) {

    if (!_grower) {
      _grower = dendritic_growth<stretch>::create(_max);
      _grower->init_velocity(_meshGraph);
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
              << m2::ci::geometric_mean_length<stretch>(_meshGraph)
              << std::endl;

    double dt = _params.dt;
    auto colors = getColor(_meshGraph);
#if 1

    // build constraints to capture current config
    std::cout << "====== " << std::endl;
    std::cout << "integrating " << std::endl;
    double Nn = 0, NRxn = 0;
    using edge_ptr = typename m2::surf<stretch>::edge *;
    using vertex_ptr = typename m2::surf<stretch>::vertex *;

    std::cout << "building constraints " << std::endl;
    hepworth::constraint_set<stretch>::ptr constraints =
        hepworth::constraint_set<stretch>::create(_meshGraph);

    std::cout << "adding constraints " << std::endl;
    init_stretch_constraints<stretch>(_meshGraph, constraints, 3.0e-2);
    // init_cross_constraints<stretch>(_meshGraph, constraints, 1e-5, 1.2);
    //    init_bend_constraints<stretch>(_meshGraph, constraints, 1.0e-7);
    init_willmore_constraints<stretch>(_meshGraph, constraints, 1.0e-7);
    // init_mem_bend_constraints<stretch>(_meshGraph, constraints,
    // 1e-6, 1.0e-4);

    constraints->add_constraint(
        hepworth::internal_collisions<stretch>::create(_meshGraph, _max));

    std::vector<stretch::vec3> p0 =
        m2::ci::get_coordinates<stretch>(_meshGraph);
    std::vector<stretch::vec3> p1 = _grower->step(_meshGraph, p0);
    // std::vector<stretch::vec3> p1 = p0;
    int i = 0;

#if 1
    hepworth::velocity_optimizer<stretch> opt(p0, p0);
    opt.update(constraints);
    p1 = constraints->get_positions();
#endif

    m2::ci::set_coordinates<stretch>(p1, _meshGraph);
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
  m2::surf<stretch> *_meshGraph;
  m2::surf_integrator<stretch> *_integrator;
  dendritic_growth<stretch>::ptr _grower;
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
