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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

real computeK(real dist, real C, real p) {
  // std laplace kernel
  real distp = pow(dist, p);
  real lp = pow(C, p);
  real kappa = 1.0 / (distp + lp);
  return kappa;
};

real computeKm(real dist, real C, real p) {
  // mollified kernel
  real distp = pow(dist, p);
  real lp = pow(C, p);

  real kappa = (1.0 - exp(-distp / lp)) / distp;
  return kappa;
};

real compute3K(real dist, real C) {
  real dist3 = dist * dist * dist;
  real l3 = C * C * C;
  // T kappa = (1.0 - exp(-dist3 / l3)) / dist3;
  real kappa = 0.5 / M_PI / dist3;

  return kappa;
};

real compute_tangent_point_radius(const vec3 &dp, const vec3 &N) {
  real ndp = dp.squaredNorm();
  real nPdp = (N * N.transpose() * dp).norm();
  return 0.5 * ndp / nPdp;
};

real compute_tangent_point_inverse_radius(const vec3 &dp, const vec3 &N,
                                          const double &p, real l0) {
  // computes the inverse radius of the tangent point to power p
  real ndp = dp.norm();
  real nPdp = (N * N.transpose() * dp).norm();
  real k = pow(nPdp, p) / (pow(ndp + l0, 2.0 * p));
  return k;
};

vec3 compute_tangent_point_radius_grad(const vec3 &dp0, const vec3 &N,
                                       const double &l0, const double &p) {
  real Ndp = dp0.dot(N);

  real ndp = dp0.norm();
  vec3 dp = ndp * dp0.normalized();

  mat3 P = N * N.transpose();
  vec3 Pdp = P * dp;

  real nPdp = (Pdp).norm();
  real lp = pow(l0, p);
  real l2 = pow(l0, 2.0);

  real k = pow(nPdp, p) / (pow(ndp, 2.0 * p) + lp);
  vec3 dk0 = P * Pdp / (nPdp * nPdp + l2);
  vec3 dk1 = 2.0 * dp / (ndp * ndp + l2);
  vec3 dk = p * k * (dk0 - dk1);
  // vec3 dk = p * k * (dk0);

  return dk;
};

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

  arp::T3::ptr face_tree = arp::T3::create(face_vert_ids, x, 24);
  // face_tree->debug_half();

  calder::fast_summation<arp::T3> sum(*face_tree);
  sum.bind<vec3>(face_ids, wN);

  std::vector<real> u = sum.calc<real>(
      pov,
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
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
        if (i == 500)
          gg::geometry_logger::line(pj, pj + N, vec4(0.1, 0.7, 0.2, 0.5));

        return va::solidAngle(pi, p0, p1, p2);
        // real kappa = computeK(dist, 0.5 * l0);
        // return kappa * va::dot(N, dp);
        // return 0.0;
      },
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
           const arp::T3::node &node, const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = compute3K(dist, 0.5 * l0);
        if (i == 500) {
          gg::geometry_logger::line(pj, pj + N, vec4(0.1, 0.7, 0.2, 0.5));
          std::cout << "kappa ndp: " << kappa * va::dot(N, dp) << std::endl;
        }

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

  std::vector<mat3> u = sum.calc<mat3>(
      p_pov,
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T2>::Node_Type node_type,
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

        real kappa = compute3K(dist, l0);
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
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T2>::Node_Type node_type,
           const arp::T2::node &node, const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);
        const mat3 &E = F_datum->node_data()[j];

        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = compute3K(dist, l0);
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

} // namespace calder
} // namespace gaudi
#endif