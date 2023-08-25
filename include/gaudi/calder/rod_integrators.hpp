//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __ROD_INTEGRATOR__
#define __ROD_INTEGRATOR__

#include "gaudi/common.h"
#include "integrators.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

using Rod_Tree_Type = arp::T2;
using Rod_Sum_Type = calder::fast_summation<Rod_Tree_Type>;

using Rod_Bind_Fcn =
    std::function<void(const std::vector<index_t> &, Rod_Sum_Type &)>;
template <typename Q>
using Rod_Compute_Fcn =
    std::function<Q(const index_t &i, const index_t &j, const vec3 &,
                    const vec3 &, const std::vector<datum::ptr> &,
                    Rod_Sum_Type::Node_Type &, const Rod_Tree_Type::node &,
                    const Rod_Tree_Type &)>;

template <typename T>
T get_data(Rod_Sum_Type::Node_Type node_type, index_t j, index_t data_id,
           const std::vector<calder::datum::ptr> &data) {
  const typename calder::datum_t<T>::ptr F_datum =
      static_pointer_cast<typename calder::datum_t<T>>(data[data_id]);
  if (node_type == Rod_Sum_Type::Node_Type::LEAF) {
    return F_datum->leaf_data()[j];
  } else {
    return F_datum->node_data()[j];
  }
}

#if 1
template <typename T>
std::vector<T> integrate_over_rod(asawa::rod::rod &R,
                                  const std::vector<vec3> &p_pov,
                                  Rod_Bind_Fcn bind_fcn = nullptr,
                                  Rod_Compute_Fcn<T> compute_fcn = nullptr) {

  // std::vector<vec3> &x = R.__x;
  std::vector<vec3> x = R.xc();
  vec3 u0 = vec3(0.0, 1.0, 0.0);

  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "summing" << std::endl;
  std::cout << " -n_faces: " << edge_ids.size() << std::endl;
  std::cout << " -create: " << std::endl;
  Rod_Tree_Type::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  Rod_Sum_Type sum(*edge_tree);
  bind_fcn(edge_ids, sum);
  std::cout << " -compute: " << std::endl;
  std::vector<T> us = sum.calc<T>(
      p_pov,
      [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                     const std::vector<calder::datum::ptr> &data,
                     Rod_Sum_Type::Node_Type node_type,
                     const Rod_Tree_Type::node &node,
                     const Rod_Tree_Type &tree) -> T {
        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        real l = (x1 - x0).norm();
        vec3 pj = va::project_on_line(x0, x1, pi);
        return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
      },
      [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                     const std::vector<calder::datum::ptr> &data,
                     Rod_Sum_Type::Node_Type node_type,
                     const Rod_Tree_Type::node &node,
                     const Rod_Tree_Type &tree) -> T {
        vec3 pj = node.center();
        return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
      },
      0.5, false);
  return us;
}
#endif

template <typename T>
std::vector<T> mls_avg(asawa::rod::rod &R, const std::vector<T> &v,
                       const std::vector<vec3> &p_pov, real l0, real p = 3.0) {

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<T> us = integrate_over_rod<T>(
      R, p_pov,
      [&v](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, v));
      },
      [l0, &sums, p](const index_t i, const index_t j, //
                     const vec3 &pi, const vec3 &pj,
                     const std::vector<calder::datum::ptr> &data,
                     Rod_Sum_Type::Node_Type node_type, //
                     const Rod_Sum_Type::Node &node,    //
                     const Rod_Sum_Type::Tree &tree) -> T {
        T e = get_data<T>(node_type, j, 0, data);
        real dist = (pj - pi).norm();
        real kappa = computeK(dist, l0, p);
        sums[i] += kappa;
        return kappa * e;
      });
#if 1
  int max_count = 0;
  int max_count_i = 0;
  for (int i = 0; i < p_pov.size(); i++) {
    us[i] /= sums[i];
  }
#endif
  return us;
}

std::vector<vec3> coulomb_force(asawa::rod::rod &R,
                                const std::vector<vec3> &p_pov, real l0,
                                real p = 3.0) {
  std::vector<real> weights = R.l0();
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<vec3> us = integrate_over_rod<vec3>(
      R, p_pov,
      [&weights](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
      },
      [l0, &sums, p](const index_t i, const index_t j, //
                     const vec3 &pi, const vec3 &pj,
                     const std::vector<calder::datum::ptr> &data,
                     Rod_Sum_Type::Node_Type node_type, //
                     const Rod_Sum_Type::Node &node,    //
                     const Rod_Sum_Type::Tree &tree) -> vec3 {
        real w = get_data<real>(node_type, j, 0, data);
        vec3 dp = pj - pi;
        real kappa = computeK(dp.norm(), l0, p);
        return w * kappa * dp;
      });
  return us;
}

std::vector<vec3> null_coulomb_force(asawa::rod::rod &R,
                                     const std::vector<vec3> &p_pov, real l0,
                                     real p = 3.0) {
  std::vector<real> weights = R.l0();

  std::vector<vec3> Tc = R.N2c();

  std::vector<vec3> us = integrate_over_rod<vec3>(
      R, p_pov,
      [&weights, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      },
      [l0, p](const index_t i, const index_t j, //
              const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Rod_Sum_Type::Node_Type node_type, //
              const Rod_Sum_Type::Node &node,    //
              const Rod_Sum_Type::Tree &tree) -> vec3 {
        real w = get_data<real>(node_type, j, 0, data);
        vec3 T = get_data<vec3>(node_type, j, 1, data);
        vec3 dp = pj - pi;
        T.normalize();
        vec3 N = va::rejection_matrix(T) * dp;
        real kappa = computeK(dp.norm(), l0, p);
        return w * kappa * N;
      });
  return us;
}

std::vector<vec3> tangent_point_force(asawa::rod::rod &R,
                                      const std::vector<vec3> &p_pov,
                                      const std::vector<real> &w_pov,
                                      const std::vector<vec3> &T_pov, real l0,
                                      real p = 3.0) {
  std::vector<real> weights = R.l0();
  std::vector<vec3> Tc = R.N2c();

  std::vector<vec3> us = integrate_over_rod<vec3>(
      R, p_pov,
      [&weights, &Tc](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::scalar_datum::create(edge_ids, weights));
        sum.bind(calder::vec3_datum::create(edge_ids, Tc));
      },
      [l0, &w_pov, &T_pov, p](const index_t i, const index_t j, //
                              const vec3 &pi, const vec3 &pj,
                              const std::vector<calder::datum::ptr> &data,
                              Rod_Sum_Type::Node_Type node_type, //
                              const Rod_Sum_Type::Node &node,    //
                              const Rod_Sum_Type::Tree &tree) -> vec3 {
        real wi = w_pov[i];
        vec3 Ti = T_pov[i];
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Tj = get_data<vec3>(node_type, j, 1, data);
        Ti.normalize();
        Tj.normalize();

        vec3 dp = pj - pi;
        // vec3 Ni = va::rejection_matrix(Ti) * dp;
        // vec3 Nj = va::rejection_matrix(Tj) * dp;
        vec3 Bi = Ti.cross(dp).normalized();
        vec3 Bj = Tj.cross(dp).normalized();
        vec3 Ni = Bi.cross(Ti).normalized();
        vec3 Nj = Bj.cross(Tj).normalized();

        real kappa = computeK(dp.norm(), l0, 3.0);
        vec3 gj = compute_tangent_point_radius_grad(dp, Nj, l0, p);
        vec3 gi = compute_tangent_point_radius_grad(-dp, Ni, l0, p);

        return wi * gi - wj * gj;
      });
  return us;
}

std::vector<vec10> quadric(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                           const std::vector<vec3> &p_pov, real l0,
                           real p = 4.0, bool normalize = true) {
  std::vector<real> weights = R.l0();

  auto mkA = [](vec3 dx) {
    mat410 A;
    real x = dx[0];
    real y = dx[1];
    real z = dx[2];
    real xx = x * x;
    real yy = y * y;
    real zz = z * z;
    real xy = 2.0 * x * y;
    real xz = 2.0 * x * z;
    real yz = 2.0 * y * z;

    A.col(0) = vec4(xx, 2.0 * x, 0.0 * y, 0.0 * z);
    A.col(1) = vec4(yy, 0.0 * x, 2.0 * y, 0.0 * z);
    A.col(2) = vec4(zz, 0.0 * x, 0.0 * y, 2.0 * z);

    A.col(3) = vec4(xy, 2.0 * y, 2.0 * x, 0.0 * z);
    A.col(4) = vec4(xz, 2.0 * z, 0.0 * y, 2.0 * x);
    A.col(5) = vec4(yz, 0.0 * x, 2.0 * z, 2.0 * y);

    A.col(6) = vec4(2.0 * x, 2.0, 0.0, 0.0);
    A.col(7) = vec4(2.0 * y, 0.0, 2.0, 0.0);
    A.col(8) = vec4(2.0 * z, 0.0, 0.0, 2.0);

    A.col(9) = vec4(1.0, 0.0, 0.0, 0.0);
    return A;
  };

  auto mkN = [](vec3 N) { return vec4(0.0, N[0], N[1], N[2]); };
  std::vector<mat10> As(p_pov.size(), mat10::Zero());
  std::vector<vec10> ns(p_pov.size(), vec10::Zero());

  std::vector<real> us = integrate_over_rod<real>(
      R, p_pov,
      [&Nr](const std::vector<index_t> &edge_ids, Rod_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, Nr));
      },
      [&](const index_t i, const index_t j, //
          const vec3 &pi, const vec3 &pj,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> real {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        vec3 Nr = get_data<vec3>(node_type, j, 0, data);
        real w = Nr.norm();
        Nr /= w;
        real dist = (pj - pi).norm();
        real kappa = computeK(dist, l0, p);
        mat410 A = mkA(pj - pi);
        vec4 n = mkN(Nr);
        ns[i] += w * kappa * A.transpose() * n;
        As[i] += w * kappa * A.transpose() * A;
        return 0.0;
      });

  std::vector<vec10> out(p_pov.size(), vec10::Zero());
  for (int i = 0; i < As.size(); i++) {
    mat10 &A = As[i];
    vec10 &n = ns[i];
    out[i] = A.colPivHouseholderQr().solve(n);
  }

  return out;
}

std::vector<real> quadric_sdf(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                              const std::vector<vec3> &p_pov, real l0,
                              real p = 4.0, bool normalize = true) {
  std::vector<vec10> Q = quadric(R, Nr, p_pov, l0, p, normalize);
  std::vector<real> out(p_pov.size(), 0.0);
  for (int i = 0; i < Q.size(); i++) {
    out[i] = Q[i][9]; //* q.segment(6, 3);
    if (std::isnan(out[i])) {
      out[i] = 0.0;
    }
  }
  return out;
}

std::vector<vec3> vortex_force(asawa::rod::rod &R,
                               const std::vector<vec3> &p_pov,
                               const std::vector<real> &phi, real l0 = 1e-2,
                               real p = 4.0) {

  std::vector<vec3> &x = R.x();
  std::vector<vec3> T = R.dirs();

  index_t i = 0;
  for (auto &t : T) {
    t *= phi[i++];
  }

  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  Rod_Tree_Type::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<Rod_Tree_Type> sum(*edge_tree);

  sum.bind<vec3>(edge_ids, T);

  std::vector<vec3> us = sum.calc<vec3>(
      p_pov,
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 T = F_datum->leaf_data()[j];

        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        real l = (x1 - x0).norm();
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real kappa = computeK(dp.norm(), l0, p);
        /*
        gg::geometry_logger::line(pj, pj + 0.1 * T, vec4(0.0, 0.0, 1.0, 1.0));
        gg::geometry_logger::line(pi, pj, vec4(0.0, 1.0, 0.0, 1.0));
        gg::geometry_logger::line(pi, pi + 1.0 * dp.cross(T),
                                  vec4(1.0, 0.0, 1.0, 1.0));
*/
        return kappa * dp.cross(T);
      },
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 T = F_datum->node_data()[j];
        // Nr.normalize();
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real kappa = computeK(dp.norm(), l0, p);

        return kappa * dp.cross(T);
      });

  return us;
}

std::vector<vec3> tangent_point_grad(asawa::rod::rod &R, real l0 = 1e-2,
                                     real p = 4.0) {

  std::vector<vec3> &x = R.x();
  std::vector<vec3> N = R.N2();

  std::vector<vec3> xc = R.xc();
  std::vector<vec3> Tc = R.N2c();
  std::vector<vec3> Nc = R.N2c();
  std::vector<vec3> Bc = R.N1c();
  /*
    for (int i = 0; i < x.size(); i++) {
      gg::geometry_logger::line(x[i], x[i] + 0.05 * N[i],
                                vec4(0.0, 1.0, 1.0, 1.0));
      gg::geometry_logger::line(xc[i], xc[i] + 0.05 * Tc[i],
                                vec4(1.0, 0.0, 0.0, 1.0));
      gg::geometry_logger::line(xc[i], xc[i] + 0.05 * Nc[i],
                                vec4(0.0, 1.0, 0.0, 1.0));
      gg::geometry_logger::line(xc[i], xc[i] + 0.05 * Bc[i],
                                vec4(0.0, 0.0, 1.0, 1.0));
    }
  */
  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  Rod_Tree_Type::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<Rod_Tree_Type> sum(*edge_tree);
  sum.bind(calder::vec3_datum::create(edge_ids, N));

  std::vector<real> sums(xc.size(), 0.0);
  std::vector<vec3> us = sum.calc<vec3>(
      xc,
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 Tr = F_datum->leaf_data()[j];
        const vec3 &Tv = Nc[i];
        real w = Tr.norm();
        if (w < 1e-8)
          return vec3::Zero();

        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        real l = (x1 - x0).norm();
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        vec3 Nr = Tr.cross(dp).normalized();
        vec3 Nv = Tv.cross(dp).normalized();

        real kappa = computeK(dp.norm(), l0, 3.0);
        vec3 gV = compute_tangent_point_radius_grad(dp, Nv, l0, p);
        vec3 gR = compute_tangent_point_radius_grad(-dp, Nr, l0, p);
        // gg::geometry_logger::line(pi, pi + dp, vec4(0.0, 0.8, 1.0, 1.0));
        // gg::geometry_logger::line(pj, pj + 1.0 * gR, vec4(0.0, 1.0,
        // 0.8, 1.0));
        //  return w * (gV - gR);
        sums[i] += w;
        // return w * kappa * mob_grad(dp, 2.0);
        return w * (gV - gR);
        // return -w * kappa * gR;
      },
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 Tr = F_datum->node_data()[j];
        const vec3 &Tv = Nc[i];
        real w = Tr.norm();
        if (w < 1e-8)
          return vec3::Zero();
        // Nr.normalize();
        Tr /= w;
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        vec3 Nr = Tr.cross(dp).normalized();
        vec3 Nv = Tv.cross(dp).normalized();

        real kappa = computeK(dp.norm(), l0, 3.0);
        vec3 gV = compute_tangent_point_radius_grad(dp, Nv, l0, p);
        vec3 gR = compute_tangent_point_radius_grad(-dp, Nr, l0, p);
        /*
                if (i == 0) {
                  gg::geometry_logger::line(pi, pi + dp, vec4(0.0,
        0.5, 1.0, 1.0)); gg::geometry_logger::line(pj, pj + 0.1 * kappa * gR,
                                            vec4(1.0, 0.5, 0.0, 1.0));
                  gg::geometry_logger::line(pj, pj + 0.1 * kappa * Nr,
                                            vec4(0.5, 1.0, 0.0, 1.0));
                }
        */
        sums[i] += w;
        // return w * kappa * mob_grad(dp, 2.0);

        return w * (gV - gR);
        // return -w * kappa * gR;
      });

  return us;
}

std::vector<mat3> covariance(asawa::rod::rod &R, const std::vector<vec3> &Nr,
                             const std::vector<vec3> &p_pov, real r, real l0,
                             real p = 4.0, bool normalize = true) {

  std::vector<vec3> &x = R.__x;
  std::vector<quat> &q = R.__u;

  vec3 u0 = vec3(0.0, 0.0, 1.0);
  std::vector<vec3> ue;
  std::vector<real> ls;
  ue.reserve(q.size());
  for (int i = 0; i < q.size(); i++) {
    ue[i] = vec3::Zero();
    if (R.next(i) == -1)
      continue;
    asawa::rod::consec_t ids = R.consec(i);
    real l = (x[ids[2]] - x[ids[1]]).norm();
    ls.push_back(l);
    vec3 ui = q[i].normalized() * (l * u0);
    ue.push_back(ui);
  }

  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  Rod_Tree_Type::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<Rod_Tree_Type> sum(*edge_tree);
  sum.bind(calder::edge_frame_datum::create(edge_ids, ue));
  sum.bind(calder::scalar_datum::create(edge_ids, ls));

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<mat3> u = sum.calc<mat3>(
      p_pov,
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);

        const vec3 &e = F_datum->leaf_data()[j];

        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        vec3 pj = va::project_on_line(x0, x1, pi);
        pj -= r * Nr[i];
        vec3 dp = pj - pi;
        real dist = va::norm(dp);

        vec3 dpN = dp / dist;

        real kappa = computeK(dist, l0, p);
        real w = (x0 - x1).norm();
        sums[i] += w * kappa;
        return kappa * w * e * e.transpose();
        // return kappa * w * dp * dp.transpose();
      },
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          Rod_Sum_Type::Node_Type node_type, //
          const Rod_Sum_Type::Node &node,    //
          const Rod_Sum_Type::Tree &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);
        const mat3 &E = F_datum->node_data()[j];
        const calder::scalar_datum::ptr R_datum =
            static_pointer_cast<calder::scalar_datum>(data[1]);
        const real &w = R_datum->node_data()[j];

        vec3 pj = node.center();
        pj -= r * Nr[i];
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0, p);
        vec3 dpN = dp / dist;

        sums[i] += w * kappa;
        return kappa * E;
        // return kappa * w * dp * dp.transpose();
      });
#if 1

  std::vector<mat3> Us(p_pov.size());
  for (int i = 0; i < p_pov.size(); i++) {
    // u[i] /= sums[i];
    Eigen::JacobiSVD<mat3> svd(u[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();
    vec3 s = svd.singularValues();
    mat3 S = mat3::Zero();

    S.col(0) = s[0] * U.col(0);
    S.col(1) = s[1] * U.col(1);
    S.col(2) = s[2] * U.col(2);
    Us[i] = S;
  }

  return Us;
#endif
}

} // namespace calder
} // namespace gaudi
#endif