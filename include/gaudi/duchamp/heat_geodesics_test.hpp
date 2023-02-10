#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"

#include "gaudi/bontecou/laplacian.hpp"

#include "gaudi/asawa/dynamic_surface.hpp"
#include "gaudi/asawa/faceloader.hpp"
#include "gaudi/asawa/manifold.hpp"
#include "gaudi/asawa/objloader_refactor.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/primitive_operations.hpp"
#include "gaudi/common.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

void debug_manifold(asawa::manifold &M, const std::vector<vec3> verts) {
  for (int i = 0; i < M.__corners_next.size(); i += 2) {
    if (M.__corners_next[i] < 0)
      continue;
    int i0 = i;
    int i1 = M.other(i0);
    int v0 = M.vert(i0);
    int v1 = M.vert(i1);
    gg::geometry_logger::line(verts[v0], verts[v1], vec4(0.5, 0.5, 0.5, 1.0));
  }
}

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

std::array<vec3, 2> extents(std::vector<vec3> &coords) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }
  return {min, max};
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

manifold::ptr build_cube() {
  std::vector<vec3> vertices;
  std::vector<std::vector<int>> faces;
  make_cube(vertices, faces);
  std::vector<index_t> corners_next, corners_vert, corners_face;
  assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
  manifold::ptr M = manifold::create(corners_next, corners_vert, corners_face);
  datum_t<vec3>::ptr vdata = datum_t<vec3>::create(prim_type::VERTEX, vertices);
  M->insert_datum(vdata);

  return M;
}

class heat_geodesics_test {
public:
  typedef std::shared_ptr<heat_geodesics_test> ptr;

  static ptr create() { return std::make_shared<heat_geodesics_test>(); }

  heat_geodesics_test() {
    //__M = build_cube();
    __M = build_bunny();

    triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    vec3_datum::ptr c_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));

    std::vector<vec3> &coords = c_datum->data();
    center(coords);

    /////////
    // dynamic surface
    /////////

    real l0 = 1.0 * asawa::avg_length(*__M, c_datum->data());
    __surf = dynamic_surface::create(__M, 1.0 * l0, 3.0 * l0, 1.0 * l0);

    /////////
    // weights
    /////////

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    std::vector<real> weights(__M->vert_count(), 0.0);

    real mx = -999.9, mn = 999.0;
    index_t mxi = 0, mni = 0.0;
    for (int i = 0; i < x.size(); i++) {
      if (x[i][0] > mx) {
        mx = x[i][0];
        mxi = i;
      }
      if (x[i][0] < mn) {
        mn = x[i][0];
        mni = i;
      }
    }
    std::cout << " imx: " << mxi << std::endl;
    weights[mni] = 1.0;
    weights[mxi] = -1.0;

    datum_t<real>::ptr wdata =
        datum_t<real>::create(prim_type::VERTEX, weights);
    _iw = __M->insert_datum(wdata);
  };

  real fib(real thet, real f) {
    real c0 = 2.0;
    real c1 = 3.0;
    real c2 = 5.0;
    real c3 = 8.0;
    real c4 = 13.0;
    real c5 = 21.0;
    real cs0 = 1.0 / c0 * cos(f * M_PI * (c0 * thet + 0.5));
    real cs1 = 1.0 / c1 * cos(f * M_PI * (c1 * thet + 0.33));
    real cs2 = 1.0 / c2 * cos(f * M_PI * (c2 * thet + 0.20));
    real cs3 = 1.0 / c3 * cos(f * M_PI * (c3 * thet + 0.125));
    real cs4 = 1.0 / c4 * cos(f * M_PI * (c4 * thet + 1.0 / 13.0));
    real cs5 = 1.0 / c5 * cos(f * M_PI * (c5 * thet + 1.0 / 21.0));
    return cs0 + cs1 + cs2 + cs3 + cs4 + cs5;
  }

  real kolmogorov(real thet, real f) {
    real k = 5 / 6;
    real c0 = pow(2.0, -k * 1.0);
    real c1 = pow(2.0, -k * 2.0);
    real c2 = pow(2.0, -k * 3.0);
    real c3 = pow(2.0, -k * 4.0);
    real c4 = pow(2.0, -k * 5.0);
    real cs0 = c0 * sin(f * M_PI * (1.0 * thet + 1.0));
    real cs1 = c1 * sin(f * M_PI * (2.0 * thet + 0.5));
    real cs2 = c2 * sin(f * M_PI * (4.0 * thet + 0.25));
    real cs3 = c3 * sin(f * M_PI * (8.0 * thet + 0.125));
    real cs4 = c4 * sin(f * M_PI * (16.0 * thet + 0.0625));
    return cs0 + cs1 + cs2 + cs3 + cs4;
  }

  real square(real thet, real f) {
    real x = 0.0;
    for (int i = 0; i < 100; i++) {
      real k = real(i);
      real K = 2.0 * k - 1.0;
      x += 4.0 / M_PI * 1.0 / K * sin(M_PI * K * f * thet);
    }
    return x;
  }

  real simple(real thet, real f) { return sin(M_PI * thet * f); }

  real freq_biased(real thet, real f) {
    real N = 4;
    real z = sin(M_PI * thet * f);
    real zp0 = pow(z, 1.0 * N);
    real zp1 = pow(1.0 - z, 64.0 * N);

    return zp0 / (zp0 + zp1) - 0.5;
  }

  void test_heat() {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    real_datum::ptr w_datum =
        static_pointer_cast<real_datum>(__M->get_datum(_iw));
    std::vector<real> &w = w_datum->data();
    real mx = -999.9, mn = 999.0;
    index_t mxi = 0, mni = 0.0;
    for (int i = 0; i < w.size(); i++) {
      if (w[i] > mx) {
        mx = w[i];
        mxi = i;
      }
      if (w[i] < mn) {
        mn = w[i];
        mni = i;
      }
    }

    std::vector<real> f(__M->vert_count(), 0.0);
    f[mni] = 1.0;
    f[mxi] = -1.0;

    bontecou::laplacian L(__M, x);
    std::vector<real> d = L.heatDist(f, 1e-3);

    mx = 0.0;
    for (int i = 0; i < d.size(); i++) {
      mx = max(d[i], mx);
    }
#if 1
    real As = asawa::surface_area(*__M, x);
    real l = sqrt(As);
    bool nn = false;
    for (int i = 0; i < __M->vert_count(); i++) {
      vec3 N = vert_normal(*__M, i, x);
      vec3 xi = x[i];

      real sn = pow(sin(M_PI * d[i] / mx / 2.0), 2.0);

      real cn = simple(d[i], 14.0);
      // real cn = kolmogorov(d[i], 8.0);
      // real cn = 1.0 * fib(d[i], 12.0);
      // real cn = 1.0 * freq_biased(d[i], 16.0);
      vec3 dN = cn * N;
      dx[i] = 0.5 * dx[i] + 0.01 * dN;
      // vec3 dN = pow(sn, 2.0) * kolmogorov(d[i], 2.0) * N;
      // vec3 dN = pow(sn, 2.0) * square(d[i], 12.0) * N;
      // vec3 dN = pow(sn, 2.0) * simple(d[i], 12.0) * N;
      // vec3 dN = pow(sn, 2.0) * freq_biased(d[i], 12.0) * N;

      // gg::geometry_logger::line(xi, xi + 0.02 * dN, vec4(0.3, 0.4, 0.95,
      // 0.5));
    }
#endif
    w = d;
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
    test_heat();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    __surf->step(0.1, dx);

    smoothMesh(0.015, 5);
    //  arp::aabb_tree<1> tree(vids, x);

    // debug_manifold(*__M, x_datum->data());
  }

  index_t _iw;
  manifold::ptr __M;
  dynamic_surface::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif