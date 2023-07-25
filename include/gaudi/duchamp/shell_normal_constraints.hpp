#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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

#include "gaudi/calder/integrators.hpp"

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
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
    asawa::center(x);

    /////////
    // dynamic surface
    /////////
    real l0 = asawa::shell::avg_length(*__M, x);
    // real C = 0.65;
    real C = 0.6;
    __surf = shell::dynamic::create(__M, C * l0, 2.5 * C * l0, C * l0);
    __surf->set_flip_pred([&](asawa::shell::shell &M, const index_t &c0) {
      index_t c1 = M.other(c0);
      bool flip = _adjacent.find(c0) == _adjacent.end();
      flip &= _adjacent.find(c1) == _adjacent.end();
      return flip;
    });
    /////////
    // constraints setup
    /////////
    std::vector<real> lr = asawa::shell::edge_lengths(*__M, x);

    datum_t<real>::ptr ldata = datum_t<real>::create(prim_type::EDGE, lr);
    _il0 = __M->insert_datum(ldata);

    for (int i = 0; i < 3; i++)
      __surf->step();

    /////////////////////
    // Rod
    /////////////////////
    std::vector<vec3> x_w = walk(2.0 * l0);

    __R = rod::rod::create(x_w, false);
    //__R->_update_frames(normals);

    real lavg = 1.5 * __R->lavg();
    __R->_r = 0.015;
    __Rd = rod::dynamic::create(__R, 0.35 * lavg, 2.0 * lavg, 0.25 * lavg);

    for (int i = 0; i < 5; i++) {
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

  vec3 align_walk(const vec3 x0, const vec3 &d0, const vec3 &N0,
                  const std::vector<vec3> walk,
                  const std::vector<vec3> &normals, real eps = 1e-2) {
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
    // vec3 dir = va::sgn(d0, vec3(l0 * f0)) * l0 * f0;
    vec3 dir = d0 + 6.3e-1 * va::sgn(d0, vec3(f0)) * l0 * f0;
    dir += 3.77e1 * va::sgn(d0, vec3(f1)) * l1 * f1;

    dir = R * dir;
    // dir += 1e-10 * va::sgn(d0, vec3(l0 * f1)) * l1 * f1;
    //     dir -= 1.0e-2 * va::sgn(d0, f1) * l1 * f1;

    // dir = va::sgn(dir.dot(d0)) * dir;
    dir.normalize();
    // gg::geometry_logger::line(x0, x0 + 0.1 * dir, vec4(1.0, 0.0, 0.0, 1.0));
    return dir;
  }

  std::vector<vec3> walk(real eps) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    index_t test = 1000;

    vec3 N = asawa::shell::edge_normal(*__M, test, x);
    vec3 T = asawa::shell::edge_tangent(*__M, test, x).normalized();
    vec3 B = N.cross(T).normalized();

    real thet = 1.05 * M_PI;

    vec3 dir = std::cos(thet) * T + std::sin(thet) * B;
    std::vector<index_t> corners;
    std::vector<real> S;
    std::vector<vec3> points;
    std::vector<vec3> normals;

    asawa::shell::walk(*__M, x, test, dir, 0.5, 4000, 1e-2,
                       [&](const asawa::shell::shell &M,
                           const std::vector<vec3> &x, const index_t &corner,
                           const real &s, const real &accumulated_length,
                           vec3 &dir) {
                         asawa::shell::index_t ci = corner;
                         S.push_back(s);
                         corners.push_back(corner);
                         vec3 pt = asawa::shell::edge_vert(M, ci, s, x);
                         vec3 Ni = asawa::shell::edge_normal(*__M, corner, x);

                         dir = align_walk(pt, dir, Ni, points, normals);

                         points.push_back(pt);
                         normals.push_back(Ni);
                         return true;
                       });

    return points;
  }
  std::vector<real> calc_dist(asawa::rod::rod &R, asawa::shell::shell &M,
                              real eps, const std::vector<vec3> &xr) {
    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    std::vector<index_t> rverts = R.get_vert_range();
    std::vector<index_t> rverts_map = R.get_vert_map();
#if 0
    vector<std::array<index_t, 2>> nearest =
        __surf->get_pnt_tri_collisions(rverts, rverts_map, xr, M, 2.0 * eps);
#else
    vector<std::array<index_t, 2>> nearest =
        __surf->get_pnt_tri_collisions(rverts, rverts_map, xr, M, 9999.9);
#endif
    std::vector<real> dist(xr.size(), 2.0 * eps);

    /// std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * eps);
    for (auto &c : nearest) {
      if (c[0] < 0)
        continue;
      if (c[1] < 0)
        continue;

      index_t ivr = c[0];
      index_t ifs = c[1];
      vec3 xri = xr[ivr];
      // vec3 xf = asawa::shell::face_center(M, ifs, x_s);
      vec3 xf = asawa::shell::face_pnt(xri, M, ifs, x_s);
      real d = (xri - xf).norm();
      // vec3 Nri = Nr[ivr].normalized();

      // gg::geometry_logger::line(xri, xf, vec4(1.0, 0.0, 1.0, 1.0));
      //  gg::geometry_logger::line(xr, xr + d * Nri, vec4(1.0, 1.0,
      //  0.0, 1.0));
      dist[ivr] = d;
    }
    return dist;
  }

  void update_rod_normals() {
    real eps = __surf->_Cc;
    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 1.0 * eps);
    __R->_update_frames(Nr);
  }

  std::vector<real> calc_rod_dist_grad(asawa::rod::rod &R,
                                       asawa::shell::shell &M, real eps) {
    const std::vector<vec3> &x_r = R.__x;
    std::vector<vec3> x_s = asawa::get_vec_data(M, 0);

    std::vector<index_t> rverts = R.get_vert_range();
    std::vector<index_t> rverts_map = R.get_vert_map();

    std::vector<vec3> x_rc(x_r);
    for (int i : rverts) {
      auto consec = R.consec(i);
      if (consec[2] < 0)
        continue;
      x_rc[i] = 0.5 * (x_r[consec[1]] + x_r[consec[2]]);
    }

    // std::vector<vec3> x_s = asawa::get_vec_data(M, 0);
    // std::vector<real> dist0 = calder::mls_dist(M, x_rc, 0.1 * eps, 3.0);

    std::vector<real> dist0 = calc_dist(R, M, eps, x_r);
    std::vector<real> dist1 = calc_dist(R, M, eps, x_rc);

    std::vector<real> g_d(dist0.size());

    // std::transform(dist.begin(), dist.end(), dist.begin(),
    //                [eps](double x) { return max(x - 1.0 * eps, 0.0); });
    for (int i = 0; i < dist0.size(); i++) {
      auto idx = R.consec(i);
      vec3 T = R.dir(i);
      real dT = T.norm();
      // vec3 Nri = Nr[idx[1]].normalized();
      //  gg::geometry_logger::line(x_r[idx[1]], x_r[idx[1]] + d0 * Nri,
      //                            vec4(1.0, 1.0, 0.0, 1.0));
      real d0 = dist0[idx[1]];
      real d1 = dist1[idx[1]];
      // d0 = std::min(d0, 16.0 * eps);
      // d1 = std::min(d1, 16.0 * eps);
      real ddi = d1 - d0;
      // dT = std::max(dT, eps / 16.0);
      g_d[i] = ddi / 2.0 / dT;
      // gg::geometry_logger::line(x_r[idx[1]], x_r[idx[1]] + g_d[i] * T,
      //                           vec4(1.0, 0.0, 1.0, 1.0));
    }

    return g_d;
  }

  void init_rod_shell_weld(
      asawa::rod::rod &R, asawa::rod::dynamic &rod_d, asawa::shell::shell &M,
      asawa::shell::dynamic &shell_d, const std::vector<vec3> &fm,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &wr, const real &ws,
      std::vector<hepworth::sim_block::ptr> blocks) {
    _adjacent.clear();
    const std::vector<vec3> &x0 = R.__x;
    std::vector<vec3> x1 = asawa::get_vec_data(M, 0);

    std::vector<index_t> edge_verts_R = R.get_edge_vert_ids();
    std::vector<index_t> edge_map_R = R.get_edge_map();
    real eps = 2.0 * shell_d._Cm;

#if 1
    std::vector<index_t> edges = M.get_edge_range();
    std::vector<index_t> edge_verts_M = M.get_edge_vert_ids();
    std::vector<index_t> edge_map_M = M.get_edge_map();

    auto g_d = calc_rod_dist_grad(R, M, 2.0 * eps);
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

      if (d[0] > 5.0e-3 * eps)
        continue;

      // gg::geometry_logger::line(xs0, xs1, vec4(0.0, 0.0, 1.0, 1.0));
      real g_di = va::mix(d[1], g_d[vr0], g_d[vr1]);
      vec3 xr = va::mix(d[1], xr0, xr1);
      vec3 dr = xr1 - xr0;
      real lr = (xr1 - xr0).norm();

      /*
      real lbound = 1e-8;
      real ubound = 1.0 - lbound;
      if (d[2] < lbound)
        continue;
      if (d[2] > ubound)
        continue;
  */
      gg::geometry_logger::line(xr, xr + 0.1 * g_di * dr,
                                vec4(0.5, 0.5, 1.0, 1.0));
      _adjacent.insert(vs0);
      hepworth::block::edge_edge_weld::ptr constraint =
          hepworth::block::edge_edge_weld::create(
              std::vector<index_t>({vr0, vr1, vs0, vs1}), wr, ws, blocks);
      constraint->set_slide(d[1] + 1.0 * g_di);
      constraints.push_back(constraint);
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
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, xv);

    real eps = 1.0 * __surf->_Cc;
    real r = 5.0 * eps;
    std::vector<vec3> Nr = calder::mls_avg<vec3>(*__M, Nf, xr, 2.0 * eps, 3.0);
    auto g_d = calc_rod_dist_grad(R, M, 2.0 * eps);
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
        rod_d.get_vert_collisions(verts_M, xv, 6.0 * eps);
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
      vec3 Nri = va::mix(t0, Nri0, Nri1).normalized();
      real g_di = va::mix(t0, g_d[vr0], g_d[vr1]);
      vec3 Nsi = Nv[vs].normalized();
      vec3 Tr = dr.normalized();
      vec3 Bri = Nri.cross(Tr).normalized();
      Nri = Tr.cross(Bri).normalized();
      vec3 xr = xr0 + t0 * dr;

      vec3 dX = xs - xr;
      real Nsdx = Nri.dot(dX);

      real d = dX.norm();

      real w_dx = exp(-1.5 * d * d / r / r);
      if (w_dx < 1e-4)
        continue;

        // std::cout << Nri << std::endl;
        // gg::geometry_logger::line(xr, xr + 0.01 * Nri, vec4(1.0, 0.0,
        // 0.0, 1.0));

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
              std::vector<index_t>({vs, vr0, vr1}), r, xr0, xr1, Nri0, Nri1,
              Nsi, w_dx * w, blocks);
      constraints.push_back(constraint);
      constraint->set_slide(w_dx * (t0 + 1.0 * g_di));
    }
#endif
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

  void test_ribbon_sdf() {
    std::vector<vec3> &xm = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> x = createPoints(200000, 0.5);
    real eps = __surf->_Cc;
    real r = 3.0 * eps;
    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * eps);
    std::vector<real> sd = calder::ribbon_sdf(*__R, Nr, x, r, eps, 3.0, true);
    for (int i = 0; i < x.size(); i++) {
      gg::geometry_logger::line(x[i], x[i] + 1e-4 * vec3(1.0, 1.0, 1.0),
                                1000.0 * gg::sdf4(sd[i]));
    }
  }

  std::vector<vec3> compute_ribbon_sdf() {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);

    real eps = 1.0 * __surf->_Cc;
    real r = 3.0 * eps;

    std::vector<vec3> Nr = get_rod_normals(*__R, *__M, 4.0 * __surf->_Cm);

    // std::vector<real> g =
    //     calder::compute_cos2(*__R, Nv, x, eps, 6.0, 4.0, true);
    std::vector<real> sd =
        calder::ribbon_sdf(*__R, Nr, x, r, 2.0 * eps, 2.0, true);
    std::vector<real> df = calder::mls_dist(*__R, x, 2.0 * eps, 2.0);
    std::transform(df.begin(), df.end(), df.begin(),
                   [r](double x) { return (exp(-1.0 * x * x / r / r)); });

    std::vector<vec3> f(x.size(), vec3::Zero());
    for (int i = 0; i < f.size(); i++) {
      vec3 fi = 1.0 * sd[i] * df[i] * Nv[i];
      gg::geometry_logger::line(x[i], x[i] + 0.1 * fi,
                                vec4(0.0, 0.5, 1.0, 1.0));

      f[i] += fi;
    }
    return f;
  }

  std::vector<vec3> compute_tangent_point_gradient(
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      std::vector<hepworth::sim_block::ptr> blocks) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, x);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, x);

    const std::vector<vec3> &xr = __R->__x;
    real eps = 1.0 * __surf->_Cc;
    real r = 1.0 * eps;
    std::cout << "mls normals" << std::endl;
    std::vector<vec3> Nr = calder::mls_avg<vec3>(*__M, Nf, xr, 4.0 * eps, 3.0);
    for (int i = 0; i < Nr.size(); i++) {
      vec3 T = __R->dir(i);
      vec3 B = Nr[i].cross(T);
      Nr[i] = T.cross(B);
    }

    std::cout << "tp grad" << std::endl;
    std::vector<vec3> gf =
        calder::tangent_point_grad(*__R, Nr, xf, Nv, 4.0 * eps, 6.0);
#if 0
    for (int i = 0; i < gf.size(); i++) {
      gg::geometry_logger::line(xf[i], xf[i] + 1e-9 * gf[i],
                                vec4(0.5, 0.5, 0.0, 1.0));
    }
#endif
    std::cout << "mls gf" << std::endl;
    std::vector<vec3> g =
        calder::mls_avg<vec3>(*__M, gf, x, 4.0 * __surf->_Cc, 3.0);
    std::vector<vec3> xc(x);
    for (int i = 0; i < g.size(); i++) {
      // std::cout << g[i] << std::endl;
      // g[i] = compute_grad(dp[i], Nr[i], 2.0);
      g[i] *= 8e-1;
      xc[i] += g[i];
      gg::geometry_logger::line(x[i], x[i] + 0.1 * g[i],
                                vec4(0.6, 0.0, 0.8, 1.0));
    }
    hepworth::block::init_pinned(*__M, constraints, xc, 5.0e-2, blocks);
    return g;
  }

  void init_weighted_willmore(
      const asawa::shell::shell &M,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, x);

    real l = 2.0 * __surf->_Cc;
    std::vector<real> df = calder::mls_dist(*__R, x, 2.0 * l, 3.0);
    std::transform(df.begin(), df.end(), df.begin(),
                   [l, w](double x) { return w * (exp(-x * x / l / l)); });
    hepworth::block::init_willmore(M, constraints, x, df, blocks);
    // return g;
  }

  void init_weighted_edge_willmore(
      const asawa::shell::shell &M,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xe = asawa::shell::edge_centers(*__M, x);

    real l = 2.0 * __surf->_Cc;
    std::vector<real> df = calder::mls_dist(*__R, xe, l, 3.0);
    std::transform(df.begin(), df.end(), df.begin(),
                   [l, w](double x) { return w * (exp(-x * x / l / l)); });
    hepworth::block::init_edge_willmore(M, constraints, df, blocks);
    // return g;
  }

  void init_weighted_area(
      const asawa::shell::shell &M,
      std::vector<hepworth::projection_constraint::ptr> &constraints,
      const real &w, std::vector<hepworth::sim_block::ptr> blocks) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> xf = asawa::shell::face_centers(*__M, x);
    std::vector<vec3> Nf = asawa::shell::face_normals(*__M, x);

    real eps = 2.0 * __surf->_Cc;
    std::vector<real> df = calder::mls_dist(*__R, xf, eps, 3.0);

    std::transform(df.begin(), df.end(), df.begin(),
                   [w, eps](double x) { return max(x - 2.5 * eps, 0.0); });
    auto [min_it, max_it] = std::minmax_element(df.begin(), df.end());
    real vmin = *min_it;
    real vmax = *max_it;
    std::transform(df.begin(), df.end(), df.begin(), [w, vmin, vmax](real x) {
      return w * (x - vmin) / (vmax - vmin);
    });
    /*
        real mean = std::accumulate(df.begin(), df.end(), 0.0f) / std::size(df);
        real variance = std::transform_reduce(
                            df.begin(), df.end(), 0.0f, std::plus<>(),
                            [mean](real x) { return (x - mean) * (x - mean); })
       / std::size(df); real stddev = std::sqrt(variance);



        // normalize each element of the vector

    */
    /*
    for (int i = 0; i < df.size(); i++) {
      gg::geometry_logger::line(xf[i], xf[i] + 0.1 * df[i] * Nf[i],
                                vec4(0.0, 0.5, 1.0, 1.0));
    }
    */
    hepworth::block::init_area(M, constraints, x, df, blocks);
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

  void step_dynamics(int frame) {

    hepworth::block::projection_solver solver;

    std::vector<real> &l0 = asawa::get_real_data(*__M, _il0);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);

    std::vector<vec3> M = asawa::shell::vertex_areas(*__M, x);
    std::vector<real> li = asawa::shell::edge_lengths(*__M, x);
    std::vector<real> &lr = __R->__l0;
    __R->update_lengths();
    for (auto &l : lr)
      l *= 1.02;

    int i = 0;
    auto range = __M->get_edge_range();
    for (auto c0 : range) {
      if (std::isnan(l0[c0 / 2]))
        continue;
      l0[c0 / 2] = asawa::shell::edge_length(*__M, c0, x);
    }

    std::vector<vec3> Ns = asawa::shell::vertex_normals(*__M, x);
    std::vector<vec3> f(x.size(), vec3::Zero());

    std::vector<hepworth::projection_constraint::ptr> constraints;
    //      std::vector<vec3> f = compute_ribbon_charge();
    // test_ribbon_sdf();
    // f = compute_tangent_point_gradient(constraints, {Xs});
    // f = compute_ribbon_sdf();
    // f = compute_covariance(constraints, {Xs});
    // f = compute_curve_grad();
    hepworth::vec3_block::ptr Xs = hepworth::vec3_block::create(M, x, v, f);

    hepworth::vec3_block::ptr Xr =
        hepworth::vec3_block::create(__R->__M, __R->__x, __R->__v, f);
    hepworth::quat_block::ptr Ur =
        hepworth::quat_block::create(__R->__J, __R->__u, __R->__o);

    std::vector<hepworth::sim_block::ptr> blocks;
    blocks.push_back(Xs);
    blocks.push_back(Xr);
    blocks.push_back(Ur);

    std::cout << "l[0]: " << l0[0] << std::endl;
#if 1
    init_rod_shell_weld(*__R, *__Rd, //
                        *__M, *__surf,
                        f, //
                        constraints, 1.0, 1.0, {Xr, Xs});
#endif
#if 1
    init_rod_shell_creep(*__R, *__Rd, //
                         *__M, *__surf,
                         f, //
                         constraints, 1.25, {Xs, Xr});
#endif
    hepworth::block::init_stretch_shear(*__R, constraints, lr, 1.0e-1,
                                        {Xr, Ur});
    hepworth::block::init_bend_twist(*__R, constraints, 1.0e-1, {Ur});
    hepworth::block::init_angle(*__R, constraints, vec3(0.0, 0.0, 1.0),
                                1e-1 * M_PI, 1e-3, {Ur});
    hepworth::block::init_pinned(*__R, constraints, __R->__x, 15.0e0, {Xr});
    real offset = min(1.0 + 0.025 * real(frame), 30.0);
    offset = 8.0;
    hepworth::block::init_collisions(*__R, *__Rd, constraints, 0.5, {Xr, Xr},
                                     offset);

    hepworth::block::init_edge_strain(*__M, constraints, x, l0, 1e-1, {Xs});
    hepworth::block::init_bending(*__M, constraints, x, 1.0e-1, {Xs});
    // hepworth::block::init_edge_growth(*__M, constraints, x, l0, 0.999, 1e-3,
    //                                   {Xs});
    // hepworth::block::init_laplacian(*__M, constraints, x, 1, 1.0e-1, {Xs});
    //        hepworth::block::init_willmore(*__M, constraints, x, 5.0e0, {Xs});
    //     hepworth::block::init_edge_willmore(*__M, constraints, 1.0e0, {Xs});
    //   hepworth::block::init_area(*__M, constraints, x, 1.0e-1, {Xs});
    hepworth::block::init_edge_willmore(*__M, constraints, 6e-1, blocks);

    //     init_willmore_from_normals(*__M, constraints, 1e-1, blocks);
    //     hepworth::block::init_edge_willmore(*__M, constraints, 5e-1, blocks);

    // init_weighted_edge_growth(*__M, constraints, 1e-1, blocks);
    init_weighted_area(*__M, constraints, 1e-0, blocks);
    // init_weighted_edge_willmore(*__M, constraints, 1e-1, blocks);

    real eps = 1.0 * __surf->_Cc;

    hepworth::block::init_pnt_tri_collisions(
        *__M, *__surf, constraints, x, 0.5 * eps, 0.5 * eps, 1.0, {Xs, Xs});

    solver.set_constraints(constraints);

    solver.step(blocks, _h, 0.9, 30);
  }

  void step_sdf(int frame) {
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);

    for (int k = 0; k < 5; k++) {
      std::vector<vec3> f = compute_ribbon_sdf();
      for (int i = 0; i < x.size(); i++) {
        x[i] += _h * f[i];
      }
      __surf->step(true);
    }
  }

  void step(int frame) {
    std::cout << "frame: " << frame << std::endl;
    // walk(__surf->_Cc);
    step_dynamics(frame);
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