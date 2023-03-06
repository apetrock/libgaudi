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

#include "gaudi/asawa/datum_x.hpp"
#include "gaudi/asawa/shell.hpp"
#include "gaudi/geometry_types.hpp"

#include "tree_code.hpp"
#include <vector>

namespace gaudi {

namespace calder {

std::vector<real> fast_winding(asawa::shell &M, const std::vector<vec3> &x,
                               const std::vector<vec3> &pov, real l0) {

  std::vector<vec3> N = asawa::face_normals(M, x);
  std::vector<real> w = asawa::face_areas(M, x);
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

std::vector<mat3> fast_frame(asawa::shell &M, const std::vector<vec3> &x,
                             const std::vector<vec3> &p_pov,
                             const std::vector<vec3> &p_normals, real l0) {
  std::vector<vec3> E = asawa::edge_dirs(M, x);
  std::vector<real> w = asawa::edge_cotan_weights(M, x);
  std::vector<real> wa = asawa::edge_areas(M, x);

  std::vector<vec3> wE(E);

  for (int i = 0; i < wE.size(); i++) {
    // wE[i] *= w[i]; // * wE[i].normalized();
    wE[i] = w[i] * wa[i] * wE[i].normalized();
  }

  std::vector<index_t> edge_verts = M.get_edge_vert_ids();
  std::vector<index_t> edge_map = M.get_edge_map();
  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  sum.bind(calder::edge_frame_datum::create(edge_verts, wE));

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
        vec3 p0 = tree.vert(2 * j + 0);
        vec3 p1 = tree.vert(2 * j + 1);
        vec3 pj = 1.0 / 2.0 * (p0 + p1);

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
} // namespace calder
} // namespace gaudi
#endif