#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/asawa/shell/walk.hpp"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"

#include "gaudi/hepworth/block/rod_constraints.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"

#include "gaudi/hepworth/block/shell_constraints.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"

#include "gaudi/hepworth/block/coupling_collisions_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver.hpp"

#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/shell_integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <vector>
#include <zlib.h>

#ifndef __M2REFACTOR_TEST__
#define __M2REFACTOR_TEST__

namespace gaudi {
namespace duchamp {

using namespace asawa;

std::vector<vec3> createPoints(int N, real std = 0.5) {
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
  std::generate_n(std::back_inserter(points), N, randNormalVec(0, std));
  return points;
}

std::vector<vec3> createPoints(int N, const std::vector<vec3> &x, real std) {
  auto randNormalVec = [](vec3 xi, real mean, real std) {
    auto randomFunc =
        [xi, distribution_ = std::normal_distribution<double>(mean, std),
         random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
          return xi + vec3(distribution_(random_engine_),
                           distribution_(random_engine_),
                           distribution_(random_engine_));
          ;
        };
    return randomFunc;
  };
  std::vector<vec3> rpts;
  for (int i = 0; i < x.size(); i++) {
    std::generate_n(std::back_inserter(rpts), N, randNormalVec(x[i], 0, std));
  }
  return rpts;
}

class shell_normal_constraints {
public:
  typedef std::shared_ptr<shell_normal_constraints> ptr;

  static ptr create() { return std::make_shared<shell_normal_constraints>(); }

  shell_normal_constraints() {
    //__M = load_cube();
    __M = shell::load_bunny();
    //__M = shell::load_crab();

    shell::triangulate(*__M);
    for (int i = 0; i < __M->face_count(); i++) {
      if (__M->fbegin(i) > 0) {
        assert(__M->fsize(i) == 3);
      }
    }

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x, 2.0);

    /////////
    // dynamic surface
    /////////
    real l0 = asawa::shell::avg_length(*__M, x);
    real C = 0.85;
    // real C = 1.75;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);
    __surf->set_flip_pred([&](asawa::shell::shell &M, const index_t &c0) {
      index_t c1 = M.other(c0);
      return _adjacent.find(c0) == _adjacent.end() &&
             _adjacent.find(c1) == _adjacent.end();
    });
    /////////
    // constraints setup
    /////////
    std::vector<real> lr = asawa::shell::edge_lengths(*__M, x);

    datum_t<real>::ptr ldata = datum_t<real>::create(prim_type::EDGE, lr);
    _il0 = __M->insert_datum(ldata);

    /////////////////////
    // Rod
    /////////////////////
    std::vector<vec3> x_w = walk(2.0 * l0);

    __R = rod::rod::create(x_w, false);
    //__R->_update_frames(normals);

    real lavg = 1.25 * l0;
    __R->_r = 0.015;
    __Rd = rod::dynamic::create(__R, 0.35 * lavg, 2.0 * lavg, 0.25 * lavg);

    for (int i = 0; i < 5; i++) {
      __surf->step();
      __Rd->step();
    }
  };

  void load_loop_rod() {
    std::uniform_real_distribution<real> dist(0.5, 1.0);
    std::mt19937_64 re;
    vec3 p0(dist(re), dist(re), dist(re));
    vec3 p1(dist(re), dist(re), dist(re));
    real r0 = p0.norm();
    real r1 = p1.norm();

    p0.normalize();
    p1.normalize();
    vec3 f2 = p0.cross(p1).normalized();
    vec3 f1 = p0.cross(f2).normalized();
    vec3 f0 = f1.cross(f2).normalized();

    std::cout << p0 << ::std::endl;
    std::cout << ::std::endl;
    std::cout << p1 << ::std::endl;

    std::vector<vec3> points;
    int N = 256;
    for (int i = 0; i < N; i++) {
      real thet = 2.0 * M_PI * real(i) / real(N);
      vec3 pi = r0 * cos(thet) * f0 + r1 * sin(thet) * f1;
      points.push_back(pi);
    }
    __R = rod::rod::create(points);
    update_rod_normals();

    real lavg = 0.5 * __R->lavg();
    real C = 0.5;
    __Rd = rod::dynamic::create(__R, C * lavg, 3.0 * C * lavg, 1.5 * C * lavg);
  }

  vec3 align_walk(const vec3 x0, const vec3 &d0, const vec3 &N0, real li,
                  const std::vector<vec3> walk,
                  const std::vector<vec3> &normals, real eps = 1e-1) {
    // dumb little test to see if I can align the walk to neighboring lines...
    // we'll do this N^2 for fun
    if (walk.size() < 16)
      return d0;

    mat3 M = mat3::Zero();
    mat3 R = va::rejection_matrix(N0);
    real w = 0.0;
    for (int i = 1; i < walk.size() - 1; i++) {
      vec3 xi1 = walk[i];
      vec3 xi0 = walk[i - 1];
      vec3 Ni = normals[i];
      Eigen::Quaterniond q;
      q.setFromTwoVectors(Ni.normalized(), N0.normalized());

      // gg::geometry_logger::line(x0, x0 + 0.1 * N0, vec4(0.0, 0.0, 1.0, 1.0));
      // gg::geometry_logger::line(x0, x0 + 0.1 * Ni, vec4(0.0, 0.5, 1.0, 1.0));

      vec3 di = xi1 - xi0;
      di = q * di;
      real d = (x0 - 0.5 * (xi1 + xi0)).norm();

      real c2 = std::exp(-d * d / eps / eps);
      real s = va::sgn(d0.dot(di));
      real wi = c2;
      w += wi;
      // M += wi * (d0 * di.transpose() + di * d0.transpose());
      M += wi * di * di.transpose();
    }
    M /= w;
    Eigen::SelfAdjointEigenSolver<mat3> es(M);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(M);
    if (solver.info() != Eigen::Success) {
      std::cout << "SelfAdjointEigenSolver failed!" << std::endl;
    }
    if (es.eigenvalues().minCoeff() == 0.0) {
      std::cout << "SelfAdjointEigenSolver singular!" << std::endl;
      // matrix is singular
    }
    real ssum = es.eigenvalues().sum();
    real l0 = es.eigenvalues()[2] / ssum;
    vec3 f0 = es.eigenvectors().col(2);
    real l1 = es.eigenvalues()[1] / ssum;
    vec3 f1 = es.eigenvectors().col(1);
    vec3 dir = va::sgn(d0, vec3(l0 * f0)) * l0 * f0;
    // vec3 dir = d0 + 1.0e-5 * va::sgn(d0, vec3(f0)) * l0 * f0;
    dir += 2.78e-0 * va::sgn(d0, vec3(f1)) * l1 * f1;

    // dir = R * dir;
    //  dir += 1e-10 * va::sgn(d0, vec3(l0 * f1)) * l1 * f1;
    //      dir -= 1.0e-2 * va::sgn(d0, f1) * l1 * f1;

    // dir = va::sgn(dir.dot(d0)) * dir;
    dir.normalize();
    // gg::geometry_logger::line(x0, x0 + 0.1 * dir, vec4(1.0, 0.0, 0.0, 1.0));
    return dir;
  }

  vec3 rotate_walk(const index_t &ci, const vec3 &d0, const vec3 &N0, real li) {
    // dumb little test to see if I can align the walk to neighboring lines...
    // we'll do this N^2 for fun

    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    real a = asawa::shell::angle(*__M, ci, x);
    return Eigen::AngleAxis<real>(1.1e-1 * cos(-2.0 * a) * li * M_PI, N0) * d0;
  }

  std::vector<vec3> walk(real eps) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    index_t test = 50;

    vec3 N = asawa::shell::edge_normal(*__M, test, x);
    vec3 T = asawa::shell::edge_tangent(*__M, test, x).normalized();
    vec3 B = N.cross(T).normalized();

    real thet = 0.85 * M_PI;

    vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
    std::vector<index_t> corners;
    std::vector<real> S;
    std::vector<vec3> points;
    std::vector<vec3> normals;
    real l = 0.0;

    asawa::shell::walk(*__M, x, test, dir, 0.5, 5000, 1e-8,
                       [&](const asawa::shell::shell &M,
                           const std::vector<vec3> &x, const index_t &corner,
                           const real &s, const real &accumulated_length,
                           vec3 &dir) {
                         asawa::shell::index_t ci = corner;
                         S.push_back(s);
                         corners.push_back(corner);
                         vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
                         vec3 Ni = asawa::shell::edge_normal(*__M, corner, x);
                         real li = 0.0;
                         if (points.size() > 0)
                           li = (pt - points.back()).norm();
                         dir = align_walk(pt, dir, Ni, li, points, normals);
                         // dir = rotate_walk(corner, dir, Ni, li);

                         points.push_back(pt);
                         normals.push_back(Ni);
                         return true;
                       });
    std::cout << "walk resulted in: " << points.size() << " points"
              << std::endl;
    return points;
  }

  std::vector<vec3> get_rod_normals(asawa::rod::rod &R, asawa::shell::shell &M,
                                    real eps) {
    std::vector<vec3> &xr = __R->__x;
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> Nf = asawa::shell::face_normals(M, x);
    std::vector<vec3> Nr = calder::mls_avg<vec3>(M, Nf, xr, eps, 3.0);
    for (int i = 0; i < Nr.size(); i++) {
      vec3 T = R.dir(i);
      vec3 B = Nr[i].cross(T);
      Nr[i] = T.cross(B);
    }
    return Nr;
  }

  void update_rod_normals() {
    real eps = __surf->_Cc;
    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 1.0 * eps);
    __R->_update_frames(Nr);
  }

  std::vector<real> calc_conservative_collisions(asawa::rod::rod &R,
                                                 asawa::shell::shell &M,
                                                 real eps, int N_spread = 4) {
    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);

    // TODO:, these assume that everything is tightly packed, this is wrong
    // assumption if knots get more complicated
    std::vector<vec3> xr = R.x_infill(N_spread);
    std::vector<index_t> rverts(xr.size());
    for (int i = 0; i < rverts.size(); i++) {
      rverts[i] = i;
    }
    vector<std::array<index_t, 2>> nearest =
        __surf->get_pnt_tri_collisions(rverts, rverts, xr, M, eps);

    std::vector<real> dist0(xr.size(), 0.0);

    std::cout << "x.size() " << R.x().size() << std::endl;

    std::cout << "xr.size() " << xr.size() << std::endl;
    std::cout << "rverts.size() " << rverts.size() << std::endl;
    real lavg = R.lavg();
    for (auto &c : nearest) {
      auto consec = R.consec(c[0] / N_spread);
      if (consec[2] < 0)
        continue;

      if (c[1] < 0) {
        dist0[c[0]] = 4.0 * lavg;
        continue;
      }
#if 0
      index_t ivr = c[0];
      index_t ifs = c[1];

      vec3 xri = xr[ivr];
      vec3 xf = asawa::shell::face_pnt(xri, M, ifs, x_s);
      gg::geometry_logger::line(xri, xf, vec4(1.0, 0.0, 1.0, 1.0));
#endif
    }

    std::vector<real> dist(R.x().size(), 0.0);
    for (int i = 0; i < dist.size(); i++) {
      for (int k = 0; k < N_spread; k++) {
        dist[i] = std::max(dist0[i * N_spread + k], dist[i]);
      }
    }

    return dist;
  }

  std::vector<real> calc_dist_1(asawa::rod::rod &R, asawa::shell::shell &M,
                                real eps, int N_spread = 4) {

    const std::vector<vec3> &xr = R.x();
    std::vector<real> dist = calc_conservative_collisions(R, M, eps, N_spread);

    auto assign = [](index_t ip, index_t im, const std::vector<vec3> &x,
                     std::vector<real> &dist) {
      dist[ip] = std::min(dist[ip], dist[im] + (x[ip] - x[im]).norm());
    };
    // sweep forward

    std::vector<index_t> rverts = R.get_ordered_verts();

    for (int k = 0; k < 2; k++) {
      // this should only need one it, why not working?

      for (int i = 0; i < rverts.size(); i++) {
        index_t ip = R.next(rverts[i]);
        index_t im = i;
        if (ip < 0)
          continue;
        assign(ip, im, xr, dist);
      }

      for (int i = rverts.size() - 1; i > -1; i--) {
        index_t ip = R.next(rverts[i]);
        index_t im = i;
        if (ip < 0)
          continue;
        assign(im, ip, xr, dist);
      }
    }

    for (int k = 0; k < 8; k++)
      for (int i = 0; i < rverts.size(); i++) {
        auto cons = R.consec(rverts[i]);
        if (cons[2] < 0)
          continue;
        index_t ip = cons[2];
        index_t i0 = cons[1];
        index_t im = cons[0];

        dist[i0] = 0.25 * dist[im] + 0.5 * dist[i0] + 0.25 * dist[ip];
      }
#if 0
    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * eps);

    for (int i : rverts) {
      auto cons = R.consec(i);
      index_t i0 = cons[1];
      if (cons[2] < 0)
        continue;

      gg::geometry_logger::line(xr[i0],
                                xr[i0] + 1.0 * dist[i0] * Nr[i0].normalized(),
                                vec4(0.0, 1.0, 1.0, 1.0));
    }
#endif
    return dist;
  }

  std::vector<real> calc_rod_dist_grad_1(asawa::rod::rod &R,
                                         asawa::shell::shell &M, real eps,
                                         int N_spread = 4) {
    std::vector<real> dist = calc_dist_1(R, M, eps, N_spread);
    std::vector<real> g_d(dist.size(), 0.0);
    real lavg = R.lavg();
    std::vector<vec3> xr = R.x();
    for (int i = 0; i < dist.size(); i++) {
      auto idx = R.consec(i);
      index_t im = idx[1];
      index_t ip = idx[2];
      vec3 xr1 = xr[ip];
      vec3 xr0 = xr[im];
      real di = (xr1 - xr0).norm();
      di = std::max(di, lavg);
      real ddi = dist[ip] - dist[im];
      ddi = va::sgn(ddi) * std::min(abs(ddi), lavg);
      g_d[i] = 4.0 * ddi / di;
      g_d[i] = va::sgn(g_d[i]) * std::min(abs(g_d[i]), 2.0);
    }
    // g_d = __R->vert_avg(g_d);
    return g_d;
  }

  void init_rod_shell_weld(
      asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
      asawa::shell::dynamic &shell_d, const std::vector<vec3> &fm,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &wr, const real &ws,
      std::vector<hepworth::sim_block::ptr> blocks) {
    _adjacent.clear();
    const std::vector<vec3> &x0 = R.x();
    std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

    std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_R = R.get_edge_map();
    real eps = 2.0 * shell_d._Cm;

    std::vector<vec3> Nr = get_rod_normals(R, M, eps);
    std::vector<vec3> Nr0 = R.N1();

#if 1
    std::vector<index_t> edges = M.get_edge_range();
    std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
    std::vector<index_t> edge_map_M = M.get_edge_map();

    // auto g_d = calc_rod_dist_grad(R, M, 2.0 * eps);

    auto g_d = calc_rod_dist_grad_1(R, M, 0.5 * eps, 4);

    vector<std::array<index_t, 4>> sr_collisions =
        rod_d.get_collisions(edge_verts_M, x1, eps);

    for (auto &c : sr_collisions) {
      if (c[0] < 0)
        continue;
      if (c[1] < 0)
        continue;
      index_t vs0 = c[0];
      index_t vs1 = c[1];
      index_t vr0 = c[2];
      index_t vr1 = c[3];
      // real gr = 3e1 * eps;
      // if (abs(g_d[vr0]) > gr || abs(g_d[vr1]) > gr)
      //   continue;

      if (vr1 < 0 || vr1 < 0 || vs0 < 0 || vs1 < 0)
        continue;

      vec3 xr0 = x0[vr0];
      vec3 xr1 = x0[vr1];
      vec3 xs0 = x1[vs0];
      vec3 xs1 = x1[vs1];

      if ((xr1 - xr0).norm() < 1e-8)
        continue;

      if ((xs1 - xs0).norm() < 1e-8)
        continue;

      std::array<real, 3> d = va::distance_Segment_Segment(xr0, xr1, xs0, xs1);

      real g_di = va::mix(d[1], g_d[vr0], g_d[vr1]);
      vec3 xr = va::mix(d[1], xr0, xr1);
      vec3 xs = va::mix(d[2], xs0, xs1);
      vec3 dr = xr1 - xr0;
      vec3 dx = xr - xs;
      vec3 Ns0 = asawa::shell::vert_normal(M, vs0, x1);
      vec3 Ns1 = asawa::shell::vert_normal(M, vs1, x1);
      vec3 Ns = va::mix(d[2], Ns0, Ns1);

      real is_perp = pow(Ns.dot(dx.normalized()), 2.0);

      if ((is_perp > 0.75 && d[0] < 1.5 * eps)) {

        // gg::geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

        real lr = (xr1 - xr0).norm();
        vec3 Nri0 = Nr[vr0].normalized();
        vec3 Nri1 = Nr[vr1].normalized();
        vec3 Nri = va::mix(d[1], Nri0, Nri1);

        vec3 Nr0i0 = Nr0[vr0].normalized();
        vec3 Nr0i1 = Nr0[vr1].normalized();
        vec3 Nr0i = va::mix(d[1], Nr0i0, Nr0i1);
        /*
        real lbound = 1e-8;
        real ubound = 1.0 - lbound;
        if (d[2] < lbound)
          continue;
        if (d[2] > ubound)
          continue;
    */
        // gg::geometry_logger::line(xs, xs + 1.0 * g_di * dr,
        //                           vec4(0.5, 0.5, 1.0, 1.0));
        _adjacent.insert(M.find_edge_from_verts(vs0, vs1));
        // gg::geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));

        hepworth::block::edge_edge_weld::ptr constraint =
            hepworth::block::edge_edge_weld::create(
                std::vector<index_t>({vr0, vr1, vs0, vs1}), wr, ws, blocks);
        constraint->set_slide(d[1] + g_di);

        // constraint->set_rotate_to(0.001, Nri, Nr0i);
        constraints.push_back(constraint);
      }
    }
#endif
  }

  void init_rod_shell_creep(
      asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
      asawa::shell::dynamic &shell_d, const std::vector<vec3> &fm,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {

    const std::vector<vec3> &xr = R.__x;
    std::vector<vec3> xv = asawa::get_vec_data(M, 0);

    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, xv);

    real eps = 1.0 * __surf->_Cc;
    real r = 1.0 * eps;
    std::vector<vec3> Nr = get_rod_normals(R, M, eps);
    std::vector<vec3> Nr0 = R.N1();
    auto g_d = calc_rod_dist_grad_1(R, M, 2.0 * eps);
    // for (int i = 0; i < x1.size(); i++)
    //   x1[i] += _h * fm[i];

    std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_R = R.get_edge_map();

#if 1
    std::vector<index_t> verts_M = M.get_vert_range();
    std::vector<index_t> edges = M.get_edge_range();
    std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
    std::vector<index_t> edge_map_M = M.get_edge_map();

    vector<std::array<index_t, 3>> sr_collisions =
        rod_d.get_vert_collisions(verts_M, xv, 999999.9);

    for (auto &c : sr_collisions) {

      index_t vs = c[0];
      index_t vr0 = c[1];
      index_t vr1 = c[2];

      if (vr1 < 0 || vr1 < 0 || vs < 0)
        continue;

      // if (R.length(vr0) < 1e-8)
      //   continue;

      vec3 xr0 = xr[vr0];
      vec3 xr1 = xr[vr1];
      vec3 dr = xr1 - xr0;

      if (dr.dot(dr) < 1e-8)
        continue;
      real lr = (xr1 - xr0).norm();
      vec3 xs = xv[vs];

      real t0 = (xs - xr0).dot(dr) / dr.dot(dr);
      vec3 Nri0 = Nr[vr0].normalized();
      vec3 Nri1 = Nr[vr1].normalized();
      vec3 Nri = va::mix(t0, Nri0, Nri1);

      vec3 Nr0i0 = Nr0[vr0].normalized();
      vec3 Nr0i1 = Nr0[vr1].normalized();
      vec3 Nr0i = va::mix(t0, Nr0i0, Nr0i1);

      real g_di = va::mix(t0, g_d[vr0], g_d[vr1]);

      vec3 Nsi = Nv[vs].normalized();
      vec3 Tr = dr.normalized();
      vec3 Bri = Nri.cross(Tr).normalized();
      Nri = Tr.cross(Bri).normalized();

      vec3 xr = va::mix(t0, xr0, xr1);

      vec3 dX = xs - xr;

      real d = dX.norm();

      real w_dx = exp(-1.5 * d * d / r / r);
      if (w_dx < 1e-4)
        continue;
        // gg::geometry_logger::line(xs, xs + g_di * dr, vec4(1.0, 0.0,
        // 0.5, 1.0));
        //  std::cout << Nri << std::endl;
        //  gg::geometry_logger::line(xr, xr + 0.01 *
        //  Nri, vec4(1.0, 0.0, 0.0, 1.0));

#if 0
      gg::geometry_logger::line(xs, xr, vec4(1.0, 0.0, 0.0, 1.0));

      if (g_di < 0.0) {
        gg::geometry_logger::line(xs, xr0 + t1 * dr, vec4(1.0, 0.25, 0.0, 1.0));
        gg::geometry_logger::line(xs, xs + delta * Tr,
                                  vec4(1.0, 0.5, 0.0, 1.0));
      } else
        gg::geometry_logger::line(xs, xr0 + t1 * dr, vec4(0.0, 1.0, 0.25, 1.0));
      gg::geometry_logger::line(xs, xs + delta * Tr, vec4(0.0, 1.0, 0.5, 1.0));
#endif
      hepworth::block::point_edge_creep::ptr constraint =
          hepworth::block::point_edge_creep::create(
              std::vector<index_t>({vs, vr0, vr1}), r, Nri0, Nri1, Nsi,
              w_dx * w, blocks);
      // constraint->set_rotate_to(0.001, Nri, Nr0i);
      if (std::isnan(t0 + w_dx * g_di))
        continue;
      constraint->set_slide(t0 + w_dx * g_di);

      constraints.push_back(constraint);
    }
#endif
  }

  std::vector<vec3> calc_bend_grad() {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    const std::vector<vec3> &x_r = R.x();

    real eps = 0.5 * Md._Cc;
    std::vector<index_t> rverts = R.get_edge_vert_ids();
    std::vector<index_t> rverts_map = R.get_edge_map();

    std::vector<vec3> Nr0 = R.N0();
    std::vector<vec3> Nr2 = R.N2();

    vector<std::array<index_t, 2>> nearest =
        Md.get_edge_edge_collisions(rverts, rverts_map, x_r, eps);
    std::vector<real> phi(x_r.size(), 0.0);

    // std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * eps);
    for (auto &c : nearest) {
      if (c[0] < 0)
        continue;
      if (c[1] < 0)
        continue;

      index_t ivr0 = c[0];
      index_t ivr1 = R.next(c[0]);
      index_t icm0 = M.vert(c[1]);
      index_t icm1 = M.vert(M.other(c[1]));
      vec3 xri0 = x_r[ivr0];
      vec3 xri1 = x_r[ivr1];
      vec3 xr = va::mix(0.5, xri0, xri1);

      vec3 Nr = Nr0[ivr0];
      vec3 Tr = Nr2[ivr0];
      vec3 Ns = asawa::shell::edge_normal(M, c[1], x_s);

      vec3 Nrs = (Nr + Ns).normalized();

      vec3 rotationAxis = Nr.cross(Ns);
      real cosAngle = Nr.dot(Ns);
      if (std::isnan(cosAngle))
        continue;
      real angleRad = std::atan2(rotationAxis.dot(Tr), cosAngle);
      /*
      gg::geometry_logger::line(xr, xr + 0.025 * Nr, vec4(1.0, 0.0,
      0.0, 1.0)); gg::geometry_logger::line(xr, xr + 0.025 * angleRad * Nrs,
                                vec4(1.0, 0.0, 1.0, 1.0));
      gg::geometry_logger::line(xr, xr + 0.025 * Ns, vec4(0.0,
      0.0, 1.0, 1.0));
      */
      phi[ivr0] = angleRad;
    }
    std::vector<vec3> f = calder::vortex_force(R, x_s, phi, eps, 3.0);
    for (int i = 0; i < f.size(); i++) {
      f[i] *= -1e-2;
      /*
      gg::geometry_logger::line(x_s[i], x_s[i] + f[i],
                                vec4(1.0, 1.0, 0.0, 0.0));
*/
    }

    return f;
  }

  void test_quadric_sdf() {
    const std::vector<vec3> &x_s = asawa::get_vec_data(*__M, 0);
    const std::vector<vec3> &x_r = __R->x();

    std::vector<vec3> x = createPoints(200000, 0.5);
    real eps = 0.5 * __surf->_Cc;

    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, eps);

    std::vector<real> sdf = calder::quadric_sdf(*__R, Nr, x, eps);

    for (int i = 0; i < x.size(); i++) {
      gg::geometry_logger::line(x[i], x[i] + 1e-4 * vec3(1.0, 1.0, 1.0),
                                1000.0 * gg::sdf4(sdf[i]));
    }
  }

  std::vector<vec3> calc_quadric_grad() {
    asawa::shell::shell &M = *__M;
    asawa::rod::rod &R = *__R;
    asawa::shell::dynamic &Md = *__surf;
    asawa::rod::dynamic &Rd = *__Rd;

    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<vec3> x_s_f = asawa::shell::face_centers(M, x_s);
    const std::vector<vec3> &x_r = R.x();

    real eps = Md._Cc;

    std::vector<vec3> Nr = get_rod_normals(R, M, eps);

    std::vector<real> Q = calder::quadric_sdf(R, Nr, x_s_f, eps);
    // std::vector<vec3> N_s = asawa::shell::vertex_normals(*__M, x_s);
    std::vector<vec3> N_s_f = asawa::shell::face_normals(*__M, x_s);

#if 1
    int i = 0;
    for (vec3 &N : N_s_f) {
      N *= -1.0 * Q[i];
      // gg::geometry_logger::line(x_s_f[i], x_s_f[i] + 1.0 * N,
      //                           vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    std::vector<vec3> Nss =
        calder::mls_avg<vec3>(*__M, N_s_f, x_s, 2.0 * __surf->_Cc, 3.0);
#if 0
    i = 0;
    for (vec3 &N : Nss) {
      gg::geometry_logger::line(x_s[i], x_s[i] + 1.0 * N,
                                vec4(1.0, 1.0, 0.0, 1.0));
      i++;
    }
#endif
    return Nss;
  }

#if 1
  std::vector<vec3> compute_tangent_point_gradient() {
    real eps = __Rd->_Cc;
    std::vector<vec3> &x = __R->x();
    std::vector<real> l = __R->l0();
    std::vector<vec3> T = __R->N2c();
    std::vector<vec3> xc = __R->xc();

    std::vector<vec3> g0 =
        calder::tangent_point_force(*__R, x, l, T, 1.0 * eps, 6.0);
    for (int i = 0; i < g0.size(); i++) {
      // gg::geometry_logger::line(x[i], x[i] + 1.0e-7 * g0[i],
      //                           vec4(0.6, 0.0, 0.8, 1.0));
      g0[i] *= -2.0e-5;
    }
    return g0;
  }
#endif

  std::vector<real> get_dist_rod(const std::vector<index_t> &vert_ids,
                                 const std::vector<vec3> &x,
                                 const asawa::rod::rod &R,
                                 const asawa::rod::dynamic &Rd) {
    std::vector<vec3> &xr = __R->x();
    std::vector<real> df(x.size(), 0.0);
    vector<std::array<index_t, 3>> sr_collisions =
        __Rd->get_vert_collisions(vert_ids, x, 99999.9);
    for (auto &c : sr_collisions) {
      index_t vs = c[0];
      index_t vr0 = c[1];
      index_t vr1 = c[2];

      if (vr1 < 0 || vr1 < 0 || vs < 0)
        continue;
      vec3 xr0 = xr[vr0];
      vec3 xr1 = xr[vr1];
      vec3 xs = x[vs];
      df[vs] = va::distance_from_line(xr0, xr1, xs);
    }
    return df;
  }

  void init_weighted_area(
      const asawa::shell::shell &M,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {

    std::vector<vec3> &xv = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, xv);
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, xv);
    std::vector<vec3> &xr = __R->x();
    std::vector<index_t> verts_F = M.get_face_range();

    real eps = 2.0 * __surf->_Cc;
    std::vector<real> df = get_dist_rod(verts_F, xf, *__R, *__Rd);

#if 0
    for (int i = 0; i < df.size(); i++) {
      gg::geometry_logger::line(xf[verts_F[i]],
                                xf[verts_F[i]] + 0.1 * df[verts_F[i]] * Nf[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
#endif
    std::transform(df.begin(), df.end(), df.begin(),
                   [w, eps](double x) { return max(x - 2.5 * eps, 0.0); });
    auto [min_it, max_it] = std::minmax_element(df.begin(), df.end());
    real vmin = *min_it;
    real vmax = *max_it;
    std::transform(df.begin(), df.end(), df.begin(), [w, vmin, vmax](real x) {
      return w * (x - vmin) / (vmax - vmin);
    });

    hepworth::block::init_area(M, constraints, xv, df, blocks);
    // return g;
  }

  void assert_nan(index_t k) {
    std::cout << k << std::endl;
    for (int i = 0; i < __R->__u.size(); ++i) {
      if (__R->__u[i].coeffs().hasNaN()) {
        std::cout << "nan at " << i << std::endl;
        exit(0);
      }
    }
  }
  void step_quadric(int frame) {
    std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> fs = calc_quadric_grad();
    for (int i = 0; i < fs.size(); i++) {
      xs[i] += _h * fs[i];
    }
  }

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &xs = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas(*__M, xs);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, xs);

    std::vector<real> &lr = __R->l0();
    std::vector<vec3> &xr = __R->x();

    __R->update_lengths();

#if 0
    for (auto &l : lr)
      l *= 1.01;
#endif

    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, xs);
    }

    std::vector<vec3> fs(xs.size(), vec3::Zero());
    std::vector<vec3> fr(xr.size(), vec3::Zero());

    std::vector<hepworth::projection_constraint::ptr> constraints;
    //      std::vector<vec3> f = compute_ribbon_charge();

    fr = compute_tangent_point_gradient();
    // fr = compute_coulomb_gradient();
    // fr = compute_tangent_coulomb_gradient();

    // fs = calc_bend_grad();
    // test_quadric_sdf();
    fs = calc_quadric_grad();
    // fs = calc_ribbon_sdf();

    hepworth::vec3_block::ptr Xs = hepworth::vec3_block::create(M, xs, v, fs);

    hepworth::vec3_block::ptr Xr =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, fr);
    hepworth::quat_block::ptr Ur =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    std::vector<hepworth::sim_block::ptr> blocks;
    blocks.push_back(Xs);
    blocks.push_back(Xr);
    blocks.push_back(Ur);

    std::cout << "l[0]: " << l0[0] << std::endl;
    std::cout << "init weld" << std::endl;
#if 1
    init_rod_shell_weld(*__R, *__Rd, //
                        *__M, *__surf,
                        fs, //
                        constraints, 1.0, 1.0, {Xr, Xs});
#endif

#if 0
    std::cout << "init creep" << std::endl;
    init_rod_shell_creep(*__R, *__Rd, //
                         *__M, *__surf,
                         fs, //
                         constraints, 0.05, {Xs, Xr});
#endif
    std::cout << "main constraints" << std::endl;
    hepworth::block::init_stretch_shear(*__R, constraints, lr, 1.0e-1,
                                        {Xr, Ur});
    hepworth::block::init_bend_twist(*__R, constraints, 1.0e-1, {Ur});
// hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.0, 1.0),
//                             1e-1 * M_PI, 1e-3, {Ur});
#if 0
    hepworth::block::init_angle(*__R, constraints, vec3(1.0, 0.0, 0.0),
                                0.25 * M_PI, 0.05, {Ur});
    // hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.1, 0.0),
    //                             0.05 * M_PI, 0.1, {Ur});
    hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.0, 1.0),
                                0.22 * M_PI, 0.05, {Ur});
#endif
    hepworth::block::init_pinned(*__R, constraints, xr, 1.0e-1, {Xr});
    std::cout << "main collisions" << std::endl;

    real offset = 1.0;
    offset = min(1.0 + 0.025 * real(frame), 8.0);
    // std::cout << "offset: " << offset << std::endl;
    // offset = 12.0;
    hepworth::block::init_collisions(*__R, *__Rd, constraints, 1.0, {Xr, Xr},
                                     offset);

    hepworth::block::init_edge_strain(*__M, constraints, xs, l0, 5.0e-2, {Xs});
    hepworth::block::init_bending(*__M, constraints, xs, 5.0e-2, {Xs});
    // hepworth::block::init_laplacian(*__M, constraints, x, 1, 1.0e-1, {Xs});

    hepworth::block::init_edge_willmore(*__M, constraints, 5e-1, blocks);
    std::cout << "init area" << std::endl;
    init_weighted_area(*__M, constraints, 1e-3, blocks);

    std::cout << "init pnt trie collisions" << std::endl;
    real eps = 1.0 * __surf->_Cc;
    hepworth::block::init_pnt_tri_collisions(
        *__M, *__surf, constraints, xs, 0.5 * eps, 0.5 * eps, 1.0, {Xs, Xs});

    solver.set_constraints(constraints);

    solver.step(blocks, _h, 0.9, 30);
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    std::cout << "  -surface" << std::endl;
    std::cout << "    -corners: " << __M->corner_count() << std::endl;
    std::cout << "    -verts: " << __M->vert_count() << std::endl;
    std::cout << "    -faces: " << __M->face_count() << std::endl;
    std::cout << "  -curve" << std::endl;
    std::cout << "    -verts: " << __R->x().size() << std::endl;

    // walk(__surf->_Cc);
    if (frame < 600) {
      // step_quadric(frame);
      step_dynamics(frame);
    }

    if (frame > 1200)
      exit(0);

    for (int k = 0; k < 3; k++) {
      __surf->step(true);
      __Rd->step();
    }
    // step_sdf(frame);
  }
  // std::map<index_t, index_t> _rod_adjacent_edges;
  std::set<index_t> _adjacent;
  real _h = 0.05;
  index_t _il0;
  index_t _iwalk;
  std::vector<index_t> _walk;
  shell::shell::ptr __M;
  shell::dynamic::ptr __surf;
  rod::rod::ptr __R;
  rod::dynamic::ptr __Rd;
};

} // namespace duchamp
} // namespace gaudi
#endif