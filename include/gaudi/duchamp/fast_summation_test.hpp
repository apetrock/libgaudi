#include "Eigen/src/Geometry/Scaling.h"
#include "gaudi/common.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/calder/integrators.hpp"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"

#include "gaudi/calder/tree_code.hpp"

#include <array>

#include <math.h>
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

class fast_summation_test {
public:
  typedef std::shared_ptr<fast_summation_test> ptr;

  static ptr create() { return std::make_shared<fast_summation_test>(); }

  fast_summation_test() {
    //__M = load_cube();
    __M = shell::load_messer();

    shell::triangulate(*__M);
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
                                      shell::shell &M,
                                      const std::vector<vec3> &x) {
    ext::extents_t ext_m = asawa::shell::ext(M, x);
    std::uniform_real_distribution<real> dist(0.0, 1.0);
    std::mt19937_64 re;

    int i = 0;
    auto scaled_rand_vec = [dist, re, ext_t, ext_m, &i]() mutable {
      vec3 p(dist(re), dist(re), dist(re));
      if (i++ == 0) {
        p = vec3(0.57, 0.63, 0.0);
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

  void test_fast_winding(shell::shell &M, const std::vector<vec3> &x, real l0) {

    real zp = 0.5 + 0.5 * sin(M_PI * real(_frame) / 100.0);
    // real zp = 1.0;

    real zs = 0.01;
    std::vector<vec3> pov = get_random_points(
        10000, {vec3(0.0, 0.0, zp - zs), vec3(1.0, 1.0, zp + zs)}, M, x);
    std::vector<real> u = calder::fast_winding(M, x, pov, l0);

    for (int i = 0; i < pov.size(); i++) {
      // gg::geometry_logger::point(pov[i], vec4(u[i], 0.4, 0.95, 0.5));
      // gg::geometry_logger::point(pov[i], vec4(u[i], 0.4, 0.95, 0.5));
      gg::geometry_logger::line(pov[i], pov[i] + 1e-4 * vec3(1.0, 0.0, 0.0),
                                gg::sdf4(5.0 * u[i]));
    }
  }

  void test_parallel_transport(shell::shell &M, const std::vector<vec3> &x,
                               const std::vector<vec3> &Nx, real l0) {

    std::vector<mat3> Us = calder::fast_frame(M, x, x, Nx, l0);

    for (int i = 0; i < x.size(); i++) {

      // U.array().rowwise() *= s.transpose().array();
      mat3 U = Us[i];
      real thet = real(_frame) * M_PI / 120;
      real cx = cos(thet);
      real cy = sin(thet);
      vec3 pf = vec3(cx, cy, 0.0);
      pf = U * pf;

      gg::geometry_logger::line(x[i] - 0.5 * pf, //
                                x[i] + pf,       //
                                vec4(0.0, 0.7, 1.0, 1.0));
    }
  }

  void test_edge_pyramids(shell::shell &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    std::vector<index_t> edge_verts = __M->get_edge_vert_ids();
    std::vector<index_t> edge_ids = __M->get_edge_range();

    std::vector<index_t> edge_map = __M->get_edge_map();
    arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 16);

    calder::test_extents(*edge_tree, edge_verts, x);
  }

  void test_pyramids(shell::shell &M) {
    vec3_datum::ptr x_datum = static_pointer_cast<vec3_datum>(M.get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    std::vector<vec3> N = asawa::shell::face_normals(M, x);
    std::vector<real> w = asawa::shell::face_areas(M, x);

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
    std::vector<vec3> Nx = asawa::shell::vertex_normals(*__M, x);
    real l0 = 1.0 * asawa::shell::avg_length(*__M, x);

    // test_pyramids(*__M);
    // test_edge_pyramids(*__M);
    test_fast_winding(*__M, x, l0);
    // test_parallel_transport(*__M, x, Nx, 16.0 * l0);
  }
  int _frame;
  shell::shell::ptr __M;
};

class transport_mesh_test {
public:
  typedef std::shared_ptr<transport_mesh_test> ptr;

  static ptr create() { return std::make_shared<transport_mesh_test>(); }

  transport_mesh_test() {
    //__M = load_cube();
    __M = asawa::shell::load_bunny();
    //__M = asawa::load_heart();
    //__M = asawa::load_skeleton();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    center(x);
    // transform(x, mat3(Eigen::AngleAxisd(-0.5 * M_PI, vec3::UnitX())));
    /////////
    // dynamic surface
    /////////

    real l0 = 1.0 * asawa::shell::avg_length(*__M, x);
    __surf = shell::dynamic::create(__M, 1.0 * l0, 3.0 * l0, 1.0 * l0);
  };

  void test_parallel_transport(int frame) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, x);
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, x);

    std::vector<vec3> N = asawa::shell::vertex_normals(*__M, x);

    real l0 = 12.0 * asawa::shell::avg_length(*__M, x);
    // TODO: convert this to do fast frame on face centers

    std::vector<mat3> Us = calder::fast_frame(*__M, x, xf, Nf, l0);
    std::vector<vec3> df(xf.size());

    for (int i = 0; i < xf.size(); i++) {
      mat3 U = Us[i];
      real thet = real(frame) * M_PI / 120;
      real cx = cos(thet);
      real cy = sin(thet);
      vec3 pf = vec3(cx, cy, 0.0);
      // vec3 pf = vec3(1.0, 0.5, 0.0);

      pf = 4.0 * U * pf;

      df[i] = pf;
#if 0

      gg::geometry_logger::frame(U, xf[i], 0.1);
      gg::geometry_logger::line(xf[i] - 0.1 * df[i], //
                                xf[i] + 0.1 * df[i], //
                                vec4(0.0, 0.7, 1.0, 1.0));
#endif
    }

    std::vector<real> div = asawa::shell::divergence(*__M, df, x);
    std::cout << df.size() << " " << x.size() << " " << div.size() << std::endl;

    bontecou::laplacian L(__M, x);
    std::vector<real> p = L.solve(div);
    std::vector<vec3> dp = asawa::shell::gradient(*__M, p, x);

    for (int i = 0; i < df.size(); i++) {
      df[i] -= 1.0 * dp[i];
#if 1
      gg::geometry_logger::line(xf[i] - 0.075 * df[i], //
                                xf[i] + 0.075 * df[i], //
                                vec4(0.0, 0.7, 1.0, 1.0));
#endif
    }
  }

  void smoothMesh(real C, int N) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();
    bontecou::laplacian3 M(__M, x);

    for (int k = 0; k < N; k++) {
      std::cout << "." << std::flush;
      M.init();
      x = M.smooth(x, C, C + 1e-6);
      int i = 0;
      for (auto xi : x) {
        if (!std::isfinite(xi.dot(xi))) {
          std::cout << xi.transpose() << std::endl;
          __M->vprintv(i);
          i++;
        }
      }
    }
    // x_datum->data() = x;
    std::cout << "done!" << std::endl;
  }

  void step(int frame) {

    std::cout << "frame: " << frame << std::endl;
    // test_twist();
    test_parallel_transport(frame);

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    __surf->step(0.1, dx);

    smoothMesh(0.001, 10);
    //   arp::aabb_tree<1> tree(vids, x);

    // debug_shell(*__M, x_datum->data());
  }
  real omega0 = 18.0 * M_PI, phi0 = 0.0;
  index_t _iw;
  index_t _io;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif