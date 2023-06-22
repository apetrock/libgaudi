//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __M2HARMONIC_INTEGRATOR__
#define __M2HARMONIC_INTEGRATOR__

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"

#include "tree_code.hpp"
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

std::vector<real> fast_winding(asawa::shell::shell &M,
                               const std::vector<vec3> &x,
                               const std::vector<vec3> &pov, real l0) {

  std::vector<vec3> N = asawa::shell::face_normals(M, x);
  std::vector<real> w = asawa::shell::face_areas(M, x);
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
                      const arp::T3::node &node, const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        const vec3 &N = N_datum->leaf_data()[j];
        vec3 p0 = tree.vert(3 * j + 0);
        vec3 p1 = tree.vert(3 * j + 1);
        vec3 p2 = tree.vert(3 * j + 2);
        vec3 pj = 1.0 / 3.0 * (p0 + p1 + p2);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        if (i == 0)
          gg::geometry_logger::line(pi, pj, vec4(0.1, 0.7, 0.2, 0.5));

        return va::solidAngle(pi, p0, p1, p2);
        // real kappa = computeK(dist, 0.5 * l0);
        // return kappa * va::dot(N, dp);
        // return 0.0;
      },
      [&computeK, l0](const index_t &i, const index_t &j, const vec3 &pi,
                      const std::vector<calder::datum::ptr> &data,
                      const arp::T3::node &node, const arp::T3 &tree) -> real {
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

  return u;
}

std::vector<mat3> fast_frame(asawa::shell::shell &M, const std::vector<vec3> &x,
                             const std::vector<vec3> &p_pov,
                             const std::vector<vec3> &p_normals, real l0) {
  std::vector<vec3> E = asawa::shell::edge_tangents(M, x);
  std::vector<real> w = asawa::shell::edge_cotan_weights(M, x);
  std::vector<real> wa = asawa::shell::edge_areas(M, x);

  std::vector<vec3> wE(E);

  for (int i = 0; i < wE.size(); i++) {
    // wE[i] *= w[i]; // * wE[i].normalized();
    wE[i] = w[i] * wE[i];
  }

  std::vector<index_t> edge_verts = M.get_edge_vert_ids();
  std::vector<index_t> edge_map = M.get_edge_map();
  std::vector<index_t> edge_ids = M.get_edge_range_2();

  std::cout << "edge_ids: ";
  for (auto id : edge_ids)
    std::cout << id << " ";
  std::cout << std::endl;
  std::cout << "edge_ids.size() = "
            << " " << edge_verts.size() << " " << edge_ids.size() << std::endl;
  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  std::cout << " ==== wE.size(): " << wE.size() << std::endl;
  sum.bind(calder::edge_frame_datum::create(edge_ids, wE));

  auto computeK = [](real dist, real C) {
    real dist3 = dist * dist * dist;
    real l3 = C * C * C;
    // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
    real kappa = 0.5 / M_PI / (dist3 + l3);

    return kappa;
  };

  std::vector<mat3> u = sum.calc<mat3>(
      p_pov,
      [&computeK, l0](const index_t &i, const index_t &j, const vec3 &pi,
                      const std::vector<calder::datum::ptr> &data,
                      const arp::T2::node &node, const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);

        const vec3 &e = F_datum->leaf_data()[j];

        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);

        vec3 dpN = dp / dist;

        real kappa = computeK(dist, l0);
        // if (i == 1926) {
        //   gg::geometry_logger::frame(E * E.transpose(), pj, 1.0);
        //   gg::geometry_logger::line(pi, pj, vec4(0.8, 0.3, 0.9, 1.0));
        //   gg::geometry_logger::line(pj, pj + E, vec4(0.0, 1.0, 1.0, 1.0));
        //   gg::geometry_logger::line(p0, p1, vec4(1.0, 0.0, 1.0, 1.0));
        // }
        //
        // return kappa * dpN.asDiagonal() * e * e.transpose();
        return kappa * e * e.transpose();
      },
      [&computeK, l0](const index_t &i, const index_t &j, const vec3 &pi,
                      const std::vector<calder::datum::ptr> &data,
                      const arp::T2::node &node, const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);
        const mat3 &E = F_datum->node_data()[j];

        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0);
        vec3 dpN = dp / dist;
        // if (i == 1926) {
        //   gg::geometry_logger::frame(E, pj, 1.0);
        //   gg::geometry_logger::line(pi, pj, vec4(0.1, 0.7, 0.2, 1.0));
        // }

        // return kappa * dpN.asDiagonal() * E;
        return kappa * E;
      });
#if 1

  std::vector<mat3> Us(p_pov.size());
  for (int i = 0; i < p_pov.size(); i++) {
    const vec3 &Ni = p_normals[i];
    mat3 R = va::rejection_matrix(Ni);
    // gg::geometry_logger::frame(R * u[i], x[i], 10.0);
    Eigen::JacobiSVD<mat3> svd(R * u[i],
                               Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();

    vec3 s = svd.singularValues();

    vec3 Nu = U.col(0).cross(U.col(1));
    if (Nu.dot(Ni) < 0) {
      U.col(1) *= -1;
    }
    if (U.col(2).dot(Ni) < 0.0) {
      U.col(2) *= -1;
    }
    Us[i] = sqrt(s[0]) * U;
  }

  return Us;
#endif
}

real compute_tan_rad(const vec3 &dp0, const vec3 &N, const double &p, real l0) {
  real ndp = dp0.norm() + l0;
  vec3 dp = ndp * dp0.normalized();
  mat3 P = N * N.transpose();
  vec3 Pdp = P * dp;
  real nPdp = (Pdp).norm();
  if (nPdp <= 1e-8 || ndp <= 1e-8)
    return 0.0;
  real k = pow(nPdp, p) / (pow(ndp, 2.0 * p));
  return k;
};

std::vector<vec3> mls_edge_interp(asawa::rod::rod &R,
                                  const std::vector<vec3> &p_pov, real l0) {

  std::vector<vec3> &x = R.__x;
  std::vector<quat> &q = R.__u;
  vec3 u0 = vec3(0.0, 1.0, 0.0);
  std::vector<vec3> u;
  u.reserve(q.size());

  // todo::these need to be tightly packed...
  for (int i = 0; i < q.size(); i++) {
    u[i] = vec3::Zero();

    if (R.next(i) == -1)
      continue;
    asawa::rod::consec_t ids = R.consec(i);
    real l = (x[ids[2]] - x[ids[1]]).norm();
    // u[i] = q[i].normalized() * (l * u0);
    vec3 ui = q[i].normalized() * (l * u0);
    u.push_back(ui);
    gg::geometry_logger::line(x[i], x[i] + ui, vec4(1.0, 0.0, 0.5, 1.0));
  }

  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  sum.bind(calder::vec3_datum::create(edge_ids, u));

  auto computeK = [](real dist, real C) {
    real dist3 = dist * dist * dist;
    real l3 = C * C * C;
    // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
    real kappa = 0.5 / M_PI / (dist3 + l3);

    return kappa;
  };
  real max_dist = 0.0;
  int max_dist_i = 0;
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<index_t> counts(p_pov.size(), 0);
  std::vector<vec3> us = sum.calc<vec3>(
      p_pov,
      [&computeK, l0, &sums, &counts, &max_dist,
       &max_dist_i](const index_t &i, const index_t &j, const vec3 &pi,
                    const std::vector<calder::datum::ptr> &data,
                    const arp::T2::node &node, const arp::T2 &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        const vec3 &e = F_datum->leaf_data()[j];
        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        real l = (x1 - x0).norm();
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0);
        // if (i == 6664) {
        // gg::geometry_logger::line(pi, pj, vec4(0.0, 0.5, 1.0, 1.0));
        // gg::geometry_logger::line(pj, pj + 0.5 * e, vec4(1.0, 0.5,
        // 0.0, 1.0));
        //}
        counts[i]++;
        if (dist > max_dist) {
          max_dist = dist;
          max_dist_i = i;
        }

        sums[i] += kappa * e.norm();
        return kappa * e;
      },
      [&computeK, l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                             const std::vector<calder::datum::ptr> &data,
                             const arp::T2::node &node,
                             const arp::T2 &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 e = F_datum->node_data()[j];
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0);

        // std::cout << kappa << " " << e.norm() << std::endl;
        // if (i == 6664) {
        //  gg::geometry_logger::line(pi, pj, vec4(0.5, 0.0, 1.0, 1.0));
        //  gg::geometry_logger::line(pj, pj + 0.1 * kappa * e,
        //                           vec4(0.0, 1.0, 0.1, 1.0));
        //}
        sums[i] += kappa * e.norm();
        return kappa * e;
      });
#if 1
  int max_count = 0;
  int max_count_i = 0;
  for (int i = 0; i < p_pov.size(); i++) {
    us[i] /= sums[i];

    // gg::geometry_logger::line(p_pov[i], p_pov[i] + 0.05 * us[i],
    //                           vec4(1.0, 0.0, 0.2, 1.0));
    //
    if (counts[i] > max_count) {
      max_count = counts[i];
      max_count_i = i;
    }
  }
  std::cout << "max count " << max_count << " " << max_count_i << std::endl;
  std::cout << "max dist " << max_dist << " " << max_dist_i << std::endl;
  return us;
#endif
}

std::vector<vec3> tangent_point_grad(asawa::rod::rod &R,
                                     const std::vector<vec3> &p_pov,
                                     const std::vector<vec3> &n_pov, real l0) {

  std::vector<vec3> &x = R.__x;
  std::vector<quat> &q = R.__u;
  vec3 u0 = vec3(0.0, 1.0, 0.0);
  std::vector<vec3> u;
  u.reserve(q.size());

  // todo::these need to be tightly packed...
  for (int i = 0; i < q.size(); i++) {
    u[i] = vec3::Zero();

    if (R.next(i) == -1)
      continue;
    asawa::rod::consec_t ids = R.consec(i);
    real l = (x[ids[2]] - x[ids[1]]).norm();
    // u[i] = q[i].normalized() * (l * u0);
    vec3 ui = q[i].normalized() * (l * u0);
    u.push_back(ui);
    gg::geometry_logger::line(x[i], x[i] + ui, vec4(1.0, 0.0, 0.5, 1.0));
  }

  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  std::vector<index_t> edge_ids = R.get_vert_range();

  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  sum.bind(calder::vec3_datum::create(edge_ids, u));
  real p = 4.0;

  auto computeK = [](real dist, real C, real p) {
    real distp = pow(dist, p);
    real lp = pow(C, p);
    // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
    real kappa = 1.0 / (distp + lp);

    return kappa;
  };

  auto computeMobiusK = [&](real dist, real C, real p) -> real {
    real k0 = 1.0 / (pow(dist, p) + pow(C, p));
    real k1 = 1.0 / (p * pow(dist, p) + pow(C, p));
    return k0 - k1;
    // return Pdp;
  };

  auto compute_grad = [&](const vec3 &dp0, const vec3 &N,
                          const double &p) -> vec3 {
    real Ndp = dp0.dot(N);

    real ndp = dp0.norm();
    vec3 dp = ndp * dp0.normalized();

    mat3 P = N * N.transpose();
    vec3 Pdp = P * dp;

    Pdp *= 1.0 - Ndp;

    real nPdp = (Pdp).norm();

    real lp = pow(l0, p);

    real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);

    vec3 dk0 = P * Pdp / (nPdp * nPdp + l0 * l0);

    // vec3 dk0 = N / (nPdp * nPdp + l0 * l0);

    vec3 dk1 = 2.0 * dp / (ndp * ndp + l0 * l0);

    vec3 dk = p * k * (dk0 - dk1);
    // vec3 dk = p * k * (dk0);
    if (!std::isfinite(dk.norm()) || std::isnan(dk.norm())) {
      std::cout << "dump: " << nPdp << " " << ndp << " " << lp << " " << k
                << " " << dk.norm() << std::endl;
      std::flush(std::cout);
      exit(0);
    }
    return dk;
    // return (1.0 - N.dot(dp.normalized())) * dk;
    //    return (1.0 - N.dot(dp.normalized())) * dk;
    //     return -dp.normalized();

    // return Pdp;
  };

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<vec3> us = sum.calc<vec3>(
      p_pov,
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          const arp::T2::node &node, const arp::T2 &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 Nr = F_datum->leaf_data()[j];
        const vec3 &Nv = n_pov[i];
        real w = Nr.norm();
        if (w < 1e-8)
          return vec3::Zero();
        Nr /= w;
        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        real l = (x1 - x0).norm();
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real kappa = computeK(dp.norm(), l0, 2.0);
        vec3 gV = compute_grad(dp, Nv, p);
        vec3 gR = compute_grad(-dp, Nr, p);
        // gg::geometry_logger::line(pi, pi + dp, vec4(0.0, 0.8, 1.0, 1.0));
        // gg::geometry_logger::line(pj, pj + 1.0 * gR, vec4(0.0, 1.0,
        // 0.8, 1.0));
        //  return w * (gV - gR);
        sums[i] += w;
        // return w * kappa * mob_grad(dp, 2.0);
        // return w * kappa * (gV - gR);
        return -w * kappa * gR;
      },
      [&](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          const arp::T2::node &node, const arp::T2 &tree) -> vec3 {
        const calder::vec3_datum::ptr F_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);

        vec3 Nr = F_datum->node_data()[j];
        const vec3 &Nv = n_pov[i];
        real w = Nr.norm();
        if (w < 1e-8)
          return vec3::Zero();
        // Nr.normalize();
        Nr /= w;
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real kappa = computeK(dp.norm(), l0, 2.0);
        vec3 gV = compute_grad(dp, Nv, p);
        vec3 gR = compute_grad(-dp, Nr, p);
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

        // return w * kappa * (gV - gR);
        return -w * kappa * gR;
      });
#if 1
  for (int i = 0; i < p_pov.size(); i++) {
    std::cout << us[i].norm() << " " << sums[i] << std::endl;
    us[i] /= sums[i];
  }
#endif
  return us;
}

std::vector<real> mls_dist(asawa::rod::rod &R, const std::vector<vec3> &p_pov,
                           real l0) {
  std::cout << __PRETTY_FUNCTION__ << " 0" << std::endl;

  std::vector<vec3> &x = R.__x;

  std::cout << __PRETTY_FUNCTION__ << " 1" << std::endl;
  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);
  std::cout << __PRETTY_FUNCTION__ << " 2" << std::endl;

  calder::fast_summation<arp::T2> sum(*edge_tree);

  auto computeK = [](real dist, real C) {
    real dist3 = dist * dist * dist;
    real l3 = C * C * C;
    // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
    real kappa = 0.5 / M_PI / (dist3 + l3);

    return kappa;
  };
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<real> u = sum.calc<real>(
      p_pov,
      [&computeK, l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                             const std::vector<calder::datum::ptr> &data,
                             const arp::T2::node &node,
                             const arp::T2 &tree) -> real {
        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);

        real kappa = computeK(dist, l0);
        sums[i] += kappa;
        return kappa * dist;
      },
      [&computeK, l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                             const std::vector<calder::datum::ptr> &data,
                             const arp::T2::node &node,
                             const arp::T2 &tree) -> real {
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0);
        sums[i] += kappa;
        return kappa * dist;
      });
#if 1
  std::cout << __PRETTY_FUNCTION__ << " 3" << std::endl;

  for (int i = 0; i < p_pov.size(); i++) {
    u[i] /= sums[i];
  }

  return u;
#endif
}

std::vector<vec3> mls_dp(asawa::rod::rod &R, const std::vector<vec3> &p_pov,
                         real l0) {
  std::cout << __PRETTY_FUNCTION__ << " 0" << std::endl;

  std::vector<vec3> &x = R.__x;

  std::cout << __PRETTY_FUNCTION__ << " 1" << std::endl;
  std::vector<index_t> edge_verts = R.get_edge_vert_ids();
  std::vector<index_t> edge_map = R.get_edge_map();
  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);
  std::cout << __PRETTY_FUNCTION__ << " 2" << std::endl;

  calder::fast_summation<arp::T2> sum(*edge_tree);

  auto computeK = [](real dist, real C) {
    real dist3 = dist * dist * dist;
    real l3 = C * C * C;
    // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
    real kappa = 0.5 / M_PI / (dist3 + l3);

    return kappa;
  };
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<vec3> u = sum.calc<vec3>(
      p_pov,
      [&computeK, l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                             const std::vector<calder::datum::ptr> &data,
                             const arp::T2::node &node,
                             const arp::T2 &tree) -> vec3 {
        vec3 x0 = tree.vert(2 * j + 0);
        vec3 x1 = tree.vert(2 * j + 1);
        vec3 pj = va::project_on_line(x0, x1, pi);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);

        real kappa = computeK(dist, l0);
        sums[i] += kappa;
        return kappa * dp;
      },
      [&computeK, l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                             const std::vector<calder::datum::ptr> &data,
                             const arp::T2::node &node,
                             const arp::T2 &tree) -> vec3 {
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0);
        sums[i] += kappa;
        return kappa * dp;
      });
#if 1
  std::cout << __PRETTY_FUNCTION__ << " 3" << std::endl;

  for (int i = 0; i < p_pov.size(); i++) {
    u[i] /= sums[i];
  }

  return u;
#endif
}
} // namespace calder
} // namespace gaudi
#endif