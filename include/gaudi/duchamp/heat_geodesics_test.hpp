#include "Eigen/src/Core/util/Meta.h"
#include "Eigen/src/Geometry/AngleAxis.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <ostream>
#include <unsupported/Eigen/FFT>

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/arp/arp.h"
#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/vec_addendum.h"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/shell.hpp"
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
  return cs0 + cs1 + cs2 + cs3;
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

void debug_shell(shell::shell &M, const std::vector<vec3> verts) {
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
  std::cout << " scale: " << s << std::endl;
  std::cout << " min/max: "           //
            << min.transpose() << " " //
            << max.transpose() << std::endl;

  for (auto &c : coords) {
    c -= min;
    c = s * c;
    c -= 0.5 * (max - min) * s;
  }
}

void transform(std::vector<vec3> &coords, const mat3 &M) {
  real accum = 0.0;
  vec3 min = coords[0];
  vec3 max = coords[0];

  for (auto &c : coords) {
    min = va::min(c, min);
    max = va::max(c, max);
  }

  for (auto &c : coords) {
    c = M * c;
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

class heat_geodesics_test {
public:
  typedef std::shared_ptr<heat_geodesics_test> ptr;

  static ptr create() { return std::make_shared<heat_geodesics_test>(); }

  heat_geodesics_test() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = asawa::load_heart();
    //__M = asawa::load_skeleton();

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
    // transform(x, mat3(Eigen::AngleAxisd(-0.5 * M_PI, vec3::UnitX())));
    /////////
    // dynamic surface
    /////////

    real l0 = 0.5 * shell::avg_length(*__M, x);
    __surf = shell::dynamic::create(__M, 1.0 * l0, 3.0 * l0, 1.0 * l0);

    /////////
    // weights
    /////////

    std::vector<real> weights(__M->vert_count(), 0.0);

    real mx = -999.9, mn = 999.0;
    index_t mxi = 0, mni = 0.0;
    index_t axis = 0;
    for (int i = 0; i < x.size(); i++) {
      if (x[i][axis] > mx) {
        mx = x[i][axis];
        mxi = i;
      }
      if (x[i][axis] < mn) {
        mn = x[i][axis];
        mni = i;
      }
    }
    std::cout << " imx: " << mxi << std::endl;
    weights[mni] = 1.0;
    weights[mxi] = -1.0;

    _iw = __M->insert_datum(real_datum::create(prim_type::VERTEX, weights));
    _io = __M->insert_datum(real_datum::create(prim_type::VERTEX, weights));
  };

  real calc_objective(const std::vector<real> &offset,
                      const std::vector<real> &dist, const real &omega,
                      const real &phi, const real &omega_init, const real &Ks) {
    real C = 0.0;
    for (int i = 0; i < dist.size(); i++) {
      real thet = dist[i];
      real S = offset[i];
      real st = sin(omega * thet + phi);
      C += pow(st - S, 2);
    }
    C += Ks * pow(omega - omega_init, 2);
    return C;
  }

  std::array<real, 3> calc_gradient(const std::vector<real> &offset,
                                    const std::vector<real> &dist,
                                    const real &omega, const real &phi,
                                    const real &omega_init, const real &Ks) {

    real CtC = 0.0;
    real g0 = 0.0;
    real g1 = 0.0;

    for (int i = 0; i < dist.size(); i++) {
      real thet = dist[i];
      real S = offset[i];
      real st = sin(omega * thet + phi);
      real ct = cos(omega * thet + phi);

      real g0i = 2.0 * thet * (st - S) * ct;
      real g1i = 2.0 * (st - S) * ct;
      g0 += g0i;
      g1 += g1i;
    }
    g0 += 2.0 * Ks * (omega - omega_init);

    CtC += pow(g0, 2.0);
    CtC += pow(g1, 2.0);
    return {g0, g1, CtC};
  }

  std::vector<real> calc_offset(const std::vector<real> &offset,
                                const std::vector<real> &dist, real &omega,
                                real &phi) {
    real omega_p = 18.0 * M_PI;
    real phi_p = 0.0;
    // real omega_p = omega;
    // real phi_p = phi;
    real omega_init = 18.0 * M_PI;
    real Ks = 0.5;
    real d = 1.0;
    std::cout << "C: " << std::flush;
    for (int k = 0; k < 40; k++) {
      real C = calc_objective(offset, dist, omega_p, phi_p, omega_init, Ks);
      if (k % 20 == 0)
        std::cout << C << " -> " << std::flush;
      auto G = calc_gradient(offset, dist, omega_p, phi_p, omega_init, Ks);

      real g0 = G[0];
      real g1 = G[1];
      real CtC = G[2];

      real alpha = 0.1 * C / CtC;
      // clamp the search around 2.0 * PI
      if (fabs(alpha * g0) > 0.25 * M_PI) {
        alpha = 2.0 * M_PI / fabs(g0);
      }

      bool backtracing = true;
      int i = 0;
      while (backtracing) {
        real omega_t = omega_p - alpha * g0;
        real phi_t = phi_p - alpha * g1;

        real Ci = calc_objective(offset, dist, omega_t, phi_t, omega_init, Ks);

        if (Ci < C)
          break;
        if (i++ > 20) {
          break;
          alpha = 0.0;
        }
        alpha /= 2.0;
      }
      omega_p -= alpha * g0;
      phi_p -= alpha * g1;
    }

    std::cout << " calc'd: " << omega_p << " " << phi_p << std::endl;

    std::vector<real> offset_out(offset);
    for (int i = 0; i < offset.size(); i++) {
      real thet = dist[i];
      offset_out[i] = sin(omega_p * thet + phi_p);
    }

    omega0 = omega_p;
    phi0 = phi_p;
    return offset_out;
  }

  std::vector<real> calc_mag(const std::vector<vec3> &v,
                             const std::vector<vec3> &x) {
    std::vector<real> m(v.size());
    for (int i = 0; i < v.size(); i++) {
      vec3 N = vert_normal(*__M, i, x);
      m[i] = va::sgn(v[i].dot(N)) * v[i].norm();
    }

    auto res = std::minmax_element(begin(m), end(m));
    real mn = *res.first;
    real mx = *res.second;

    if (mx - mn > 1e-16) {
      std::for_each(begin(m), end(m),
                    [mn, mx](real &v) { v = (v - mn) / (mx - mn); });
    }
    std::cout << "dx min/max: " << mn << " " << mx << std::endl;
    return m;
  }

  void test_heat(int frame) {

    vec3_datum::ptr x_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(0));
    std::vector<vec3> &x = x_datum->data();

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    real_datum::ptr w_datum =
        static_pointer_cast<real_datum>(__M->get_datum(_iw));
    std::vector<real> &w = w_datum->data();

    real_datum::ptr o_datum =
        static_pointer_cast<real_datum>(__M->get_datum(_io));
    std::vector<real> &o = o_datum->data();

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

    std::cout << " min/max ids: " << mni << " " << mxi << " w: " << w[mni]
              << " " << w[mxi] << std::endl;

    bontecou::laplacian L(__M, x);
    std::vector<real> d = L.heatDist(f, 0.2);
    for (int i = 0; i < __M->vert_count(); i++) {
      vec3 N = vert_normal(*__M, i, x);
    }
    mx = 0.0;
    for (int i = 0; i < d.size(); i++) {
      mx = max(d[i], mx);
    }
#if 1

    real As = shell::surface_area(*__M, x);
    real l = sqrt(As);
    bool nn = false;

    if (frame == 0)
      for (int i = 0; i < __M->vert_count(); i++) {
        o[i] = sin(omega0 * d[i] + phi0);
      }
    else {
      std::vector<real> o_dx = calc_mag(dx, x);
      o = calc_offset(o_dx, d, omega0, phi0);
    }

    for (int i = 0; i < __M->vert_count(); i++) {
      vec3 N = vert_normal(*__M, i, x);
      // gg::geometry_logger::line(x[i], x[i] + 0.05 * w[i] * N,
      //                           vec4(1.0, 0.0, 0.0, 1.0));

      dx[i] = 0.5 * dx[i] + 0.01 * o[i] * N;
    }
#endif
    w_datum->data() = d;
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
    test_heat(frame);

    vec3_datum::ptr v_datum =
        static_pointer_cast<vec3_datum>(__M->get_datum(1));
    std::vector<vec3> &dx = v_datum->data();

    __surf->step(0.1, dx);

    smoothMesh(0.015, 10);
    //  arp::aabb_tree<1> tree(vids, x);

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