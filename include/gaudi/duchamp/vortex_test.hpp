#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"
#include "gaudi/define_create_func.h"

#include "gaudi/bontecou/laplacian.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __TEST_TEMPLATE__
#define __TEST_TEMPLATE__

namespace gaudi {
namespace duchamp {
using namespace asawa;

struct circulation_datum : public datum_t<real> {
public:
  DEFINE_CREATE_FUNC(circulation_datum)

  circulation_datum(prim_type type, const std::vector<real> &data)
      : datum_t<real>(type, data){};
  virtual ~circulation_datum(){};

  virtual void calc(const shell::shell &M, const index_t &i,
                    const index_t &prim_id, const real &C,
                    const std::vector<index_t> &vals) {
    return; // do nothing
    // real val = this->__data[i];
    // this->__data[i] = val;
  }

  virtual void subdivide(const shell::shell &M, //
                         const index_t &i0, const index_t &i1,
                         const index_t &i2, const index_t &i3, //
                         const real &s, const index_t &source_corner) {
    index_t c0 = source_corner;
    index_t c1 = M.other(c0);
    if (__type != EDGE)
      return;
    // std::cout << this->__data.size() << std::endl;
    real val = this->__data[source_corner / 2];
    this->__data[i0] = val;
    this->__data[i1] = val;
    this->__data[i2] = 0.0;
    this->__data[i3] = 0.0;
  };

  virtual void flip(const shell::shell &M, const index_t &c0) {
    if (__type != EDGE)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(M, 0);
    index_t c1 = M.other(c0);
    index_t c0p = M.prev(c0);
    index_t c1p = M.prev(c1);
    vec3 e0 = x[M.vert(c1)] - x[M.vert(c0)];
    vec3 e1 = x[M.vert(c1p)] - x[M.vert(c0p)];

    real val = this->__data[c0 / 2];
    val = e1.dot(e0) / e1.dot(e1) * val;
  };
};

class vortex_test {
public:
  typedef std::shared_ptr<vortex_test> ptr;

  static ptr create() { return std::make_shared<vortex_test>(); }

  vortex_test() { load_shell(); };

  void load_shell() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    real l0 = asawa::shell::avg_length(*__M, x);

    /////////
    // dynamic surface
    /////////
    real C = 0.5;
    // real C = 1.5;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, 0.25 * C * l0);

    std::vector<real> g0(__M->corner_count() / 2, 0.0);
    circulation_datum::ptr gdata =
        circulation_datum::create(prim_type::EDGE, g0);
    _ig0 = __M->insert_datum(gdata);

    std::vector<real> g_0v(__M->vert_count(), 0.0);
    circulation_datum::ptr gdata_v =
        circulation_datum::create(prim_type::VERTEX, g_0v);
    _ig1 = __M->insert_datum(gdata_v);

    for (int i = 0; i < 5; i++) {
      __surf->step();
    }

    _eps = asawa::shell::avg_length(*__M, x);
  }

  TYPEDEF_MAT_NM(6, 5)
  TYPEDEF_MAT_NM(5, 6)
  TYPEDEF_VEC(5)
  TYPEDEF_MAT(5)

  std::vector<vec3> calc_vorticity_change(asawa::shell::shell &M, real alpha) {

    vec3 gravity(0.0, 9.8, 0.0);

    std::vector<index_t> face_range = M.get_face_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nf = asawa::shell::face_normals(M, x);

    std::vector<vec3> dw(M.face_count(), vec3::Zero());

    for (auto f : face_range) {
      vec3 N = Nf[f];
      real A = asawa::shell::face_area(M, f, x);
      vec3 w = alpha * A * N.cross(gravity);
#if 0
      vec3 xi = asawa::shell::face_center(M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e3 * w, vec4(0.0, 1.0, 0.0, 1.0));
#endif
      dw[f] = w;
    }

    return dw;
  }

  std::vector<real>
  calc_gamma_from_vorticity_edge(asawa::shell::shell &M,
                                 const std::vector<vec3> &dw) {

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> g(M.corner_count() / 2, 0.0);

    for (auto c0 : edge_range) {
      // std::cout << c0 << std::endl;
      index_t c1 = M.other(c0);
      index_t c0p = M.prev(c0);
      index_t c1p = M.prev(c1);

      index_t c0n = M.next(c0);
      index_t c1n = M.next(c1);

      index_t f0 = M.face(c0);
      index_t f1 = M.face(c1);
      index_t i_e = c0 / 2;
      mat6 E = mat6::Zero();
      vec6 wi = vec6::Zero();
      wi.segment(0, 3) = dw[f0];
      wi.segment(3, 3) = dw[f1];
      vec3 e0 = asawa::shell::g_edge_tangent(M, c0, x);
      vec3 e1 = asawa::shell::g_edge_tangent(M, c1, x);
      vec3 e0p = asawa::shell::g_edge_tangent(M, c0p, x);
      vec3 e1p = asawa::shell::g_edge_tangent(M, c1p, x);
      vec3 e0n = asawa::shell::g_edge_tangent(M, c0n, x);
      vec3 e1n = asawa::shell::g_edge_tangent(M, c1n, x);

      E.block(0, 0, 3, 1) = e0;
      E.block(0, 1, 3, 1) = e0n;
      E.block(0, 2, 3, 1) = e0p;
      E.block(3, 0, 3, 1) = e0;
      E.block(3, 3, 3, 1) = e1n;
      E.block(3, 4, 3, 1) = e1p;
      E.block(5, 0, 1, 6) = vec6(1.0, 1.0, 1.0, 1.0, 1.0, 1.0).transpose();
      // std::cout << K << std::endl;
      mat6 EtE = E.transpose() * E + 1.0e-12 * mat6::Identity();
      vec6 gi = EtE.ldlt().solve(E.transpose() * wi);
      // std::cout << wi.transpose() << std::endl;
#if 0
      vec3 xf0 = asawa::shell::face_center(M, f0, x);
      vec3 xf1 = asawa::shell::face_center(M, f1, x);
      vec3 w0 = gi[0] * e0 + gi[1] * e0n + gi[2] * e0p;
      vec3 w1 = -gi[0] * e1 + gi[3] * e1n + gi[4] * e1p;
      // vec3 w0 = dw[f0];
      // vec3 w1 = dw[f1];
      // gg::geometry_logger::line(xf0, xf0 + 1.0e3 * w0,
      //                          vec4(0.0, 1.0, 1.0, 1.0));
      // gg::geometry_logger::line(xf1, xf1 + 1.0e3 * w1,
      //                          vec4(0.0, 1.0, 1.0, 1.0));

#endif
#if 0
      g[c0 / 2] += gi[0];
#else
      g[c0 / 2] += gi[0];
      g[c0n / 2] += gi[1];
      g[c0p / 2] += gi[2];
      g[c1n / 2] += gi[3];
      g[c1p / 2] += gi[4];
#endif
    }
    return g;
  }

  std::vector<real>
  calc_gamma_from_vorticity_face(asawa::shell::shell &M,
                                 const std::vector<vec3> &dw) {
    std::vector<index_t> face_range = M.get_face_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> g(M.corner_count() / 2, 0.0);

    for (auto f : face_range) {
      // std::cout << c0 << std::endl;
      index_t c0 = M.fbegin(f);
      index_t c0n = M.next(c0);
      index_t c0p = M.prev(c0);
      vec3 w = dw[f];
      vec4 w4(w[0], w[1], w[2], 0.0);
      mat4 E = mat4::Zero();
      E.block(0, 0, 3, 1) = w4;
      vec3 e0 = asawa::shell::g_edge_tangent(M, c0, x);
      vec3 e0n = asawa::shell::g_edge_tangent(M, c0n, x);
      vec3 e0p = asawa::shell::g_edge_tangent(M, c0p, x);

      if (e0.norm() < 1.0e-12)
        continue;
      if (e0n.norm() < 1.0e-12)
        continue;
      if (e0p.norm() < 1.0e-12)
        continue;

      E.block(0, 0, 3, 1) = e0;
      E.block(0, 1, 3, 1) = e0n;
      E.block(0, 2, 3, 1) = e0p;
      E.block(3, 0, 1, 4) = vec4(1.0, 1.0, 1.0, 1.0).transpose();

      // std::cout << K << std::endl;
      mat4 EtE = E.transpose() * E;
      // Eigen::PartialPivLU<mat4> lu = EtE.partialPivLu();
      // vec4 gi = lu.solve(E.transpose() * w4);
      vec4 gi = EtE.ldlt().solve(E.transpose() * w4);

      g[c0 / 2] += gi[0];
      g[c0n / 2] += gi[1];
      g[c0p / 2] += gi[2];
    }
    return g;
  }

  std::vector<real>
  calc_gamma_from_vorticity_full(asawa::shell::shell &M,
                                 const std::vector<vec3> &dw) {

    std::vector<index_t> edge_range = M.get_edge_range();
    std::vector<index_t> edge_range_inv(edge_range.size(), -1);
    std::vector<index_t> edge_map(M.corner_count() / 2, -1);
    for (int i = 0; i < edge_range_inv.size(); i++) {
      edge_range_inv[i] = i;
    }

    for (int i = 0; i < edge_range.size(); i++) {
      edge_map[edge_range[i] / 2] = i;
    }

    std::vector<index_t> face_range = M.get_face_range();
    std::vector<index_t> face_range_inv(face_range.size(), 0);
    for (int i = 0; i < face_range_inv.size(); i++) {
      face_range_inv[i] = i;
    }

    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> g(M.corner_count() / 2, 0.0);

    matS E(3 * face_range_inv.size(), edge_range_inv.size());
    vecX w(3 * face_range_inv.size());

    std::vector<trip> triplets;
    for (auto fi : face_range_inv) {
      index_t f = face_range_inv[fi];
      w.segment(3 * fi, 3) = dw[f];
      M.for_each_face(f, [&](index_t c0, asawa::shell::shell &M) {
        asawa::vec3 e0 = asawa::shell::g_edge_tangent(M, c0, x);
        index_t ei = edge_map[c0 / 2];

        triplets.push_back(trip(3 * fi + 0, ei, e0[0]));
        triplets.push_back(trip(3 * fi + 1, ei, e0[1]));
        triplets.push_back(trip(3 * fi + 2, ei, e0[2]));

        // w[f] += g[c0 / 2] * asawa::shell::g_edge_tangent(M, c0, x);
      });
    }
    std::cout << " setting from triplets" << std::endl;
    E.setFromTriplets(triplets.begin(), triplets.end());
    matS EtE = E.transpose() * E;
    m_solver S(EtE);
    vecX Etw = E.transpose() * w;
    std::cout << " solving" << std::endl;
    vecX g1 = S.solve(Etw);

    for (int i = 0; i < g1.size(); i++) {
      g[edge_range[i] / 2] = g1[i];
    }
    return g;
  }

  std::vector<vec3> calc_vorticity_from_gamma(asawa::shell::shell &M,
                                              const std::vector<real> &g) {
    std::vector<index_t> face_range = M.get_face_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nf = asawa::shell::face_normals(M, x);

    std::vector<vec3> w(M.face_count(), vec3::Zero());

    for (auto f : face_range) {

      M.for_each_face(f, [&](index_t c0, asawa::shell::shell &M) {
        w[f] += g[c0 / 2] * asawa::shell::g_edge_tangent(M, c0, x);
      });

      if (w[f].norm() > 0.1 * _eps)
        w[f] = vec3::Zero();
// std::cout << w.transpose() << std::endl;
#if 0
      vec3 xi = asawa::shell::face_center(M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e2 * w[f],
                                vec4(1.0, 0.0, 0.0, 1.0));
#endif
    }
    return w;
  }

  std::vector<vec3> calc_edge_vorticity(asawa::shell::shell &M, real alpha) {
    std::vector<index_t> vert_range = M.get_vert_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> &g0 = asawa::get_real_data(M, _ig0);

    std::vector<vec3> dw = calc_vorticity_change(M, 5e-3 * _h);

    std::vector<real> dg = calc_gamma_from_vorticity_face(M, dw);
    // std::cout << dg.size() << " " << g0.size() << std::endl;

    for (int i = 0; i < g0.size(); i++) {
      vec3 Ne = asawa::shell::edge_normal(M, 2 * i, x);
      vec3 ce = asawa::shell::edge_center(M, 2 * i, x);
      g0[i] += dg[i];
      // gg::geometry_logger::line(ce, ce + 1.0e2 * g0[i] * Ne,
      //                           vec4(0.0, 0.0, 1.0, 1.0));
    }

    return calc_vorticity_from_gamma(M, g0);
  }

  void _diffuse(const std::vector<vec3> &x, std::vector<real> &f, real dt) {

    std::cout << x.size() << std::endl;
    bontecou::laplacian L(__M, x);
    std::vector<real> f_comp =
        asawa::shell::compress_to_vert_range<real>(*__M, f);

    std::vector<real> d = L.diffuse2(f_comp, dt);
    std::vector<real> d_exp =
        asawa::shell::expand_from_vert_range<real>(*__M, d);

    for (int k = 0; k < f.size(); k++) {
      f[k] = d_exp[k];
    }
  }

  std::vector<vec3> calc_vert_vorticity(asawa::shell::shell &M, real alpha) {

    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<index_t> vert_range = M.get_vert_range();
    std::vector<real> &g0 = asawa::get_real_data(M, _ig1);
    std::vector<real> dg(g0.size(), 0.0);
    vec3 B(0.0, 9.8, 0.0);
    for (auto v : vert_range) {
      real A = asawa::shell::vert_area(M, v, x);
      vec3 N = asawa::shell::vert_normal(M, v, x);
      dg[v] = alpha * A * N.dot(B);
    }

    _diffuse(x, dg, 1e-1);
    for (auto v : vert_range) {
      g0[v] += dg[v];
#if 0
      vec3 N = asawa::shell::vert_normal(M, v, x);
      gg::geometry_logger::line(x[v], x[v] + 1.0e2 * dg[v] * N,
                                vec4(1.0, 0.0, 0.0, 1.0));
#endif
    }
    _diffuse(x, dg, 1e-3);
    return asawa::shell::circulation(M, g0, x);
  }

  void step_dynamics(int frame) {

    asawa::shell::shell &M = *__M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);

    std::vector<vec3> w = calc_vert_vorticity(M, 1e-4 * _h);
    // std::vector<vec3> w = calc_edge_vorticity(M, 1e-4 * _h);

#if 0
    std::vector<index_t> face_range = M.get_face_range();
    for (auto f : face_range) {
      vec3 xi = asawa::shell::face_center(M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e2 * w[f],
                                vec4(1.0, 1.0, 0.0, 1.0));
    }
#endif
#if 1
    std::vector<vec3> f = calder::vortex_force(M, x, w, 32.0 * _eps, 3.0);
    // std::vector<vec3> fc = asawa::shell::vert_to_face<vec3>(M, f);
    //  std::vector<vec3> fs = calder::mls_avg<vec3>(M, fc, x, 4.0 * _eps, 3.0);

    for (int i = 0; i < f.size(); i++) {
      x[i] += _h * f[i];
    }
#endif
    // asawa::center(x, 2.0);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;

    for (int k = 0; k < 1; k++) {
      __surf->step(true);
    }
    step_dynamics(frame);
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  real _h = 0.05;
  real _eps;
  index_t _ig0, _ig1;

  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
};

} // namespace duchamp
} // namespace gaudi
#endif