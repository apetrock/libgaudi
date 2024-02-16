
#ifndef __DUCHAMP_VORTEX_MODULE__
#define __DUCHAMP_VORTEX_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/common.h"
#include "module_base_shell.hpp"
#include <vector>

namespace gaudi {
namespace duchamp {

struct circulation_datum : public asawa::datum_t<real> {
public:
  DEFINE_CREATE_FUNC(circulation_datum)

  circulation_datum(asawa::prim_type type, const std::vector<real> &data)
      : datum_t<real>(type, data){};
  virtual ~circulation_datum(){};

  virtual void calc(const asawa::shell::shell &M, const index_t &i,
                    const index_t &prim_id, const real &C,
                    const std::vector<index_t> &vals) {
    return; // do nothing
    // real val = this->__data[i];
    // this->__data[i] = val;
  }

  virtual void subdivide(const asawa::shell::shell &M, //
                         const index_t &i0, const index_t &i1,
                         const index_t &i2, const index_t &i3, //
                         const real &s, const index_t &source_corner) {
    index_t c0 = source_corner;
    index_t c1 = M.other(c0);
    if (__type != asawa::EDGE)
      return;
    // std::cout << this->__data.size() << std::endl;
    real val = this->__data[source_corner / 2];
    this->__data[i0] = val;
    this->__data[i1] = val;
    this->__data[i2] = 0.0;
    this->__data[i3] = 0.0;
  };

  virtual void flip(const asawa::shell::shell &M, const index_t &c0) {
    if (__type != asawa::EDGE)
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

TYPEDEF_MAT_NM(6, 5)
TYPEDEF_MAT_NM(5, 6)
TYPEDEF_VEC(5)
TYPEDEF_MAT(5)
////////////////////////////////////////////

class vortex_edge : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(vortex_edge)
  vortex_edge(asawa::shell::shell::ptr M, real b)
      : module_base_shell(M), _b(b) {

    std::vector<real> g0(_M->corner_count() / 2, 0.0);
    circulation_datum::ptr gdata =
        circulation_datum::create(asawa::prim_type::EDGE, g0);
    _i = _M->insert_datum(gdata);

    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    _eps = asawa::shell::avg_length(*_M, x);
  };

  virtual ~vortex_edge(){

  };

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
      index_t c0 = M.fbegin(f);
      index_t c0n = M.next(c0);
      index_t c0p = M.prev(c0);
      vec3 w = dw[f];
      vec4 w4(w[0], w[1], w[2], 0.0);
      mat4 E = mat4::Zero();

      E.block(0, 0, 4, 1) = w4;

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

      if (w[f].norm() > 0.1 * _eps || w[f].hasNaN())
        w[f] = vec3::Zero();

#if 0
      vec3 xi = asawa::shell::face_center(M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e2 * w[f],
                                vec4(1.0, 0.0, 0.0, 1.0));
#endif
    }
    std::cout << "done" << std::endl;
    return w;
  }

  std::vector<vec3> calc_vorticity(asawa::shell::shell &M, real alpha) {
    std::vector<index_t> vert_range = M.get_vert_range();
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<real> &g0 = asawa::get_real_data(M, _i);
    std::vector<vec3> dw = calc_vorticity_change(M, alpha);

    // std::vector<real> dg = calc_gamma_from_vorticity_face(M, dw);
    std::vector<real> dg = calc_gamma_from_vorticity_edge(M, dw);
    std::cout << dg.size() << " " << g0.size() << std::endl;

    for (int i = 0; i < g0.size(); i++) {
      vec3 Ne = asawa::shell::edge_normal(M, 2 * i, x);
      vec3 ce = asawa::shell::edge_center(M, 2 * i, x);
      g0[i] += dg[i];
      // gg::geometry_logger::line(ce, ce + 1.0e2 * g0[i] * Ne,
      //                           vec4(0.0, 0.0, 1.0, 1.0));
    }
    return calc_vorticity_from_gamma(M, g0);
  }

  virtual void step(real h) {
    // i don't love this model, tbf
    _w = calc_vorticity(*_M, _b * h);
  }

  std::vector<real> &get_circulation() { return asawa::get_real_data(*_M, _i); }
  
  std::vector<vec3> &get_vorticity() { return _w; }
  void set_voriticity(const std::vector<vec3> &w) { _w = w; }

  index_t _i;

  real _b;
  real _eps;
  std::vector<vec3> _w;
};

/////////////////////////////////////////////////////////////

class vortex_vert : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(vortex_vert)
  vortex_vert(asawa::shell::shell::ptr M, real b, real d0, real d1)
      : module_base_shell(M), _b(b), _d0(d0), _d1(d1) {

    _i = gaudi::asawa::init_vert_datum<real>(*_M, 0.0);
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    _eps = asawa::shell::avg_length(*_M, x);
  };

  virtual ~vortex_vert(){};

  std::vector<vec3> calc_vorticity(asawa::shell::shell &M, real alpha) {

    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<index_t> vert_range = M.get_vert_range();
    std::vector<real> &g0 = asawa::get_real_data(M, _i);
    std::vector<real> dg(g0.size(), 0.0);
    vec3 B(0.0, 9.8, 0.0);
    for (auto v : vert_range) {
      real A = asawa::shell::vert_area(M, v, x);
      vec3 N = asawa::shell::vert_normal(M, v, x);
      dg[v] = alpha * A * N.dot(B);
    }

    diffuse(_M, x, dg, _d0);
    // std::vector<real> dgf = asawa::shell::vert_to_face<real>(M, dg);
    // std::vector<real> dga = calder::mls_avg<real>(M, dgf, x, 1.0 *
    // _eps, 4.0);

    for (auto v : vert_range) {
      g0[v] += dg[v];
#if 0
      vec3 N = asawa::shell::vert_normal(M, v, x);
      gg::geometry_logger::line(x[v], x[v] + 1.0e2 * dg[v] * N,
                                vec4(1.0, 0.0, 0.0, 1.0));
#endif
    }
    diffuse(_M, x, dg, _d1);
    std::vector<vec3> w = asawa::shell::circulation(M, g0, x);
    std::vector<index_t> face_range = M.get_face_range();

    for (auto f : face_range) {
      if (w[f].norm() > 0.1 * _eps || w[f].hasNaN())
        w[f] = vec3::Zero();
#if 0
      vec3 xi = asawa::shell::face_center(M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e2 * w[f],
                                vec4(1.0, 0.0, 0.0, 1.0));
#endif
    }
    return w;
  }

  virtual void step(real h) {
    // i don't love this model, tbf
    _w = calc_vorticity(*_M, _b * h);
#if 0
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    std::vector<index_t> face_range = _M->get_face_range();
    for (auto f : face_range) {
      vec3 xi = asawa::shell::face_center(*_M, f, x);
      gg::geometry_logger::line(xi, xi + 1.0e2 * _w[f],
                                vec4(1.0, 1.0, 0.0, 1.0));
    }
#endif
  }

  std::vector<real> &get_circulation() { return asawa::get_real_data(*_M, _i); }
  std::vector<vec3> &get_vorticity() { return _w; }
  index_t _i;
  real _eps;
  real _b, _d0, _d1;
  std::vector<vec3> _w;
};

} // namespace duchamp
} // namespace gaudi

#endif