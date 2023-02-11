#include "Eigen/src/Geometry/Scaling.h"
#include "gaudi/common.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/dynamic_surface.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/manifold.hpp"
#include "gaudi/asawa/objloader_refactor.hpp"

#include "gaudi/asawa/datum_x.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/primitive_operations.hpp"

#include "gaudi/calder/tree_code.hpp"

#include <array>

#include <random>

#include <cmath>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

void center(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  vec3 dl = (max - min);
  real maxl = dl[0];
  maxl = maxl > dl[1] ? maxl : dl[1];
  maxl = maxl > dl[2] ? maxl : dl[2];
  real s = 2.0 / maxl;
  vec3 cen = (0.5 * (max + min) - min) * s;
  std::cout << " scale: " << s << std::endl;
  cen -= min;
  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= cen;
  }
}

manifold::ptr build_bunny() {
  std::string file("assets/bunny.obj");
  // std::string file("assets/skeleton.obj");
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  loadObjfile(file, vertices, faces);
  // make_cube(vertices, faces);
  // faces.pop_back();
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  manifold::ptr M = manifold::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}

class fast_summation_test {
public:
  typedef std::shared_ptr<fast_summation_test> ptr;

  static ptr create() { return std::make_shared<fast_summation_test>(); }

  fast_summation_test() {
    //__M = build_cube();
    __M = build_bunny();

    triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    std::vector<vec3> &x = x_datum->data();
    center(x);
  };

  std::vector<vec3> createPoints(int N) {
    auto randNormalVec = [](real mean, real std) {
      auto randomFunc =
          [distribution_ = std::normal_distribution<double>(mean, std),
           random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
            return vec3(distribution_(random_engine_),
                        distribution_(random_engine_),
                        distribution_(random_engine_));
            ;
          };
      return randomFunc;
    };

    std::vector<vec3> points;
    std::generate_n(std::back_inserter(points), N, randNormalVec(0, 0.5));
    return points;
  }
  std::vector<vec3> get_random_points(const int &N, const ext::extents_t &ext_t,
                                      manifold &M, const std::vector<vec3> &x) {
    ext::extents_t ext_m = asawa::ext(M, x);
    std::uniform_real_distribution<real> dist(0.0, 1.0);
    std::mt19937_64 re;
    double a_random_double = unif(re);
    int i = 0;
    auto scaled_rand_vec = [dist, re, ext_t, ext_m, &i]() mutable {
      vec3 p(dist(re), dist(re), dist(re));
      if (i++ == 0) {
        p = vec3(0.6, 0.6, 0.0);
      }
      vec3 dt = ext_t[1] - ext_t[0];
      vec3 dm = ext_m[1] - ext_m[0];

      p = Eigen::Scaling(dt) * p + ext_t[0];
      p = Eigen::Scaling(dm) * p;
      p += ext_m[0];
      return p;
    };

    std::vector<vec3> tpoints;
    std::generate_n(std::back_inserter(tpoints), N, scaled_rand_vec);

    return tpoints;
  }

  void fast_winding(manifold &M, const std::vector<vec3> &x, real l0) {

    real zp = 0.5 + 0.5 * sin(M_PI * real(_frame) / 100.0);
    // real zp = 1.0;

    real zs = 0.01;
    std::vector<vec3> pov = get_random_points(
        10000, {vec3(0.0, 0.0, zp - zs), vec3(1.0, 1.0, zp + zs)}, M, x);

    std::vector<vec3> N = asawa::face_normals(*__M, x);
    std::vector<real> w = asawa::face_areas(*__M, x);
    std::vector<vec3> wN(N);

    for (int i = 0; i < wN.size(); i++)
      wN[i] *= w[i];

    std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
    std::vector<index_t> face_map = M.get_face_map();
    std::vector<index_t> face_ids = M.get_face_range();

    arp::T3::ptr face_tree = arp::T3::create(face_vert_ids, x, 12);
    // face_tree->debug_half();

    calder::fast_summation<arp::T3> sum(*face_tree);
    sum.bind<vec3>(face_ids, wN);

    auto computeK = [](real dist, real C) {
      real dist3 = dist * dist * dist;
      real l3 = C * C * C;
      // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
      real kappa = 0.5 / M_PI / dist3;

      return kappa;
    };

    std::vector<real> u = sum.calc<real>(
        pov,
        [&computeK, l0](const index_t &i, const index_t &j, const vec3 &pi,
                        const std::vector<calder::datum::ptr> &data,
                        const arp::T3::node &node,
                        const arp::T3 &tree) -> real {
          const calder::vec3_datum::ptr N_datum =
              static_pointer_cast<calder::vec3_datum>(data[0]);

          const vec3 &N = N_datum->leaf_data()[j];
          vec3 p0 = tree.vert(3 * j + 0);
          vec3 p1 = tree.vert(3 * j + 1);
          vec3 p2 = tree.vert(3 * j + 2);
          vec3 pj = 1.0 / 3.0 * (p0 + p1 + p2);
          vec3 dp = pj - pi;
          real dist = va::norm(dp);
          // gg::geometry_logger::line(pi, pj, vec4(0.1, 0.7, 0.2, 0.5));

          return va::solidAngle(pi, p0, p1, p2);
          // real kappa = computeK(dist, 0.5 * l0);
          // return kappa * va::dot(N, dp);
          // return 0.0;
        },
        [&computeK, l0](const index_t &i, const index_t &j, const vec3 &pi,
                        const std::vector<calder::datum::ptr> &data,
                        const arp::T3::node &node,
                        const arp::T3 &tree) -> real {
          const calder::vec3_datum::ptr N_datum =
              static_pointer_cast<calder::vec3_datum>(data[0]);
          const vec3 &N = N_datum->node_data()[j];
          vec3 pj = node.center();
          vec3 dp = pj - pi;
          real dist = va::norm(dp);
          real kappa = computeK(dist, 0.5 * l0);

          // return kappa * dist;
          return kappa * va::dot(N, dp);
          // return 0.0;
        });

    for (int i = 0; i < pov.size(); i++) {
      gg::geometry_logger::point(pov[i], vec4(u[i], 0.4, 0.95, 0.5));
    }
  }

  void fast_frame(manifold &M, const std::vector<vec3> &x, real l0) {

    real zp = 0.5 + 0.5 * sin(M_PI * real(_frame) / 100.0);
    // real zp = 1.0;

    real zs = 0.01;
    // std::vector<vec3> pov = get_random_points(
    //     10000, {vec3(0.0, 0.0, zp - zs), vec3(1.0, 1.0, zp + zs)}, M, x);

    std::vector<vec3> E = asawa::edge_dirs(M, x);
    std::vector<real> w = asawa::edge_cotan_weights(M, x);
    std::vector<real> wa = asawa::edge_areas(M, x);
    std::vector<vec3> N = asawa::vertex_normals(*__M, x);

    std::vector<vec3> wE(E);

    for (int i = 0; i < wE.size(); i++) {
      // wE[i] *= w[i]; // * wE[i].normalized();
      wE[i] = w[i] * wa[i] * wE[i].normalized();
    }

    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> edge_map = __M->get_edge_map();
    arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

    calder::fast_summation<arp::T2> sum(*edge_tree);
    sum.bind(calder::edge_frame_datum::create(edge_verts, wE));

    auto computeK = [](real dist, real C) {
      real dist3 = dist * dist * dist;
      real l3 = C * C * C;
      // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
      real kappa = 0.5 / M_PI / (dist3 + l3);

      return kappa;
    };

    std::vector<mat3> u = sum.calc<mat3>(
        x,
        [&computeK, &N, l0](const index_t &i, const index_t &j, const vec3 &pi,
                            const std::vector<calder::datum::ptr> &data,
                            const arp::T2::node &node,
                            const arp::T2 &tree) -> mat3 {
          const calder::edge_frame_datum::ptr F_datum =
              static_pointer_cast<calder::edge_frame_datum>(data[0]);

          const vec3 &e = F_datum->leaf_data()[j];
          vec3 p0 = tree.vert(2 * j + 0);
          vec3 p1 = tree.vert(2 * j + 1);
          vec3 pj = 1.0 / 2.0 * (p0 + p1);

          vec3 dp = pj - pi;
          real dist = va::norm(dp);
          real kappa = computeK(dist, 2.0 * l0);
          // if (i == 1926) {
          //   gg::geometry_logger::frame(E * E.transpose(), pj, 1.0);
          //   gg::geometry_logger::line(pi, pj, vec4(0.8, 0.3, 0.9, 1.0));
          //   gg::geometry_logger::line(pj, pj + E, vec4(0.0, 1.0, 1.0, 1.0));
          //   gg::geometry_logger::line(p0, p1, vec4(1.0, 0.0, 1.0, 1.0));
          // }
          //
          return kappa * e * e.transpose();
        },
        [&computeK, &N, l0](const index_t &i, const index_t &j, const vec3 &pi,
                            const std::vector<calder::datum::ptr> &data,
                            const arp::T2::node &node,
                            const arp::T2 &tree) -> mat3 {
          const calder::edge_frame_datum::ptr F_datum =
              static_pointer_cast<calder::edge_frame_datum>(data[0]);
          const mat3 &E = F_datum->node_data()[j];

          vec3 pj = node.center();
          vec3 dp = pj - pi;
          real dist = va::norm(dp);
          real kappa = computeK(dist, 2.0 * l0);
          // if (i == 1926) {
          //   gg::geometry_logger::frame(E, pj, 1.0);
          //   gg::geometry_logger::line(pi, pj, vec4(0.1, 0.7, 0.2, 1.0));
          // }

          return kappa * E;
        });
#if 1
    for (int i = 0; i < x.size(); i++) {
      const vec3 &Ni = N[i];
      mat3 R = va::rejection_matrix(Ni);
      Eigen::JacobiSVD<mat3> svd(R * u[i], Eigen::ComputeFullU);
      mat3 U = svd.matrixU();
      vec3 s = svd.singularValues();
      // U.col(0) = vec3(s[0] * U.col(0));
      // U.col(1) = vec3(s[1] * U.col(1));
      // U.col(2) = vec3(s[2] * U.col(2));
      U.array().rowwise() *= s.transpose().array();
      // gg::geometry_logger::frame(U, x[i], 0.0005);
      gg::geometry_logger::frame(U, x[i], 2.0);

      // gg::geometry_logger::frame(u[i], x[i], 0.001);
    }
#endif
  }

  void test_edge_pyramids(manifold &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> edge_ids = __M->get_edge_range();

    std::vector<index_t> edge_map = __M->get_edge_map();
    arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 16);

    calder::test_extents(*edge_tree, edge_verts, x);
  }

  void test_pyramids(manifold &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    std::vector<vec3> N = asawa::face_normals(M, x);
    std::vector<real> w = asawa::face_areas(M, x);

    std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
    std::vector<index_t> face_map = M.get_face_map();
    std::vector<index_t> face_ids = M.get_face_range();

    arp::T3::ptr face_tree = arp::T3::create(face_vert_ids, x, 12);
    face_tree->debug();

    // calder::test_extents(*face_tree, face_vert_ids, x);
    calder::test_pyramid(*face_tree, face_ids, N, w);
  }

  void step(int frame) {

    _frame = frame;

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    real l0 = 1.0 * asawa::avg_length(*__M, x);

    // test_pyramids(*__M);
    fast_winding(*__M, x, l0);
    // test_edge_pyramids(*__M);
    // fast_frame(*__M, x, l0);
  }
  int _frame;
  manifold::ptr __M;
};

} // namespace duchamp
} // namespace gaudi
#endif