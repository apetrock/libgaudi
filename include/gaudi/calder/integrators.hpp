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

// gaussian kernel
real computeKg(real dist, real C, real p) {
  real distp = pow(dist, p);
  real lp = pow(C, p);
  real kappa = exp(-distp / lp);
  return kappa;
};

real compute3K(real dist) {
  real dist3 = dist * dist * dist;
  real kappa = 1.0 / dist3;
  return kappa;
};

std::vector<real> fast_winding(const arp::T3::ptr &face_tree,
                               const std::vector<vec3> &pov,
                               real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->indices();
  std::vector<index_t> face_ids(face_vert_ids.size() / 3);

  for (int i = 0; i < face_ids.size(); i++) {
    face_ids[i] = i;
  }

  std::vector<vec3> N(face_vert_ids.size() / 3);
  real l0 = 0.0;
  for (int i = 0; i < face_ids.size(); i++) {
    vec3 x0 = x[face_vert_ids[i * 3 + 0]];
    vec3 x1 = x[face_vert_ids[i * 3 + 1]];
    vec3 x2 = x[face_vert_ids[i * 3 + 2]];
    N[i] = (x1 - x0).cross(x2 - x0);
    real A = N[i].norm();
    l0 += A;
  }

  l0 /= face_ids.size();
  l0 = pow(l0, 0.5);
  l0 *= spread;

  calder::fast_summation<arp::T3> sum(*face_tree);
  sum.bind<vec3>(face_ids, N);

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
        return 0.25 / M_PI * va::solidAngle(pi, p0, p1, p2);
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
        real kappa = computeK(dist, 0.0, 3.0);
        return 0.25 / M_PI * kappa * va::dot(N, dp);
      });

  return u;
}

std::vector<real> fast_dist(const arp::T3::ptr &face_tree,
                            const std::vector<vec3> &pov, real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->indices();
  std::vector<index_t> face_ids(face_vert_ids.size() / 3);

  for (int i = 0; i < face_ids.size(); i++) {
    face_ids[i] = i;
  }

  std::vector<vec3> N(face_vert_ids.size() / 3);
  real l0 = 0.0;
  for (int i = 0; i < face_ids.size(); i++) {
    vec3 x0 = x[face_vert_ids[i * 3 + 0]];
    vec3 x1 = x[face_vert_ids[i * 3 + 1]];
    vec3 x2 = x[face_vert_ids[i * 3 + 2]];
    N[i] = (x1 - x0).cross(x2 - x0);
    real A = N[i].norm();
    l0 += A;
  }

  l0 /= face_ids.size();
  l0 = pow(l0, 0.5);
  l0 *= spread;

  calder::fast_summation<arp::T3> sum(*face_tree);
  sum.bind<vec3>(face_ids, N);
  std::vector<real> min_dists(pov.size(),
                              std::numeric_limits<real>::infinity());
  std::vector<real> u = sum.calc<real>(
      pov,
      [l0, &min_dists](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       calder::fast_summation<arp::T3>::Node_Type node_type,
                       const arp::T3::node &node, const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->leaf_data()[j];

        vec3 p0 = tree.vert(3 * j + 0);
        vec3 p1 = tree.vert(3 * j + 1);
        vec3 p2 = tree.vert(3 * j + 2);

        std::array<real, 4> cp = va::closest_point({p0, p1, p2}, pi);
        min_dists[i] = std::min(min_dists[i], cp[0]);
        return 0.0;
      },
      [l0, &min_dists](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       calder::fast_summation<arp::T3>::Node_Type node_type,
                       const arp::T3::node &node, const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        // real dist = va::project_to_nullspace(dp, N);
        real dist = va::norm(dp);
        min_dists[i] = std::min(min_dists[i], dist);

        return 0.0;
      });

  return min_dists;
}

std::vector<real> fast_view(const arp::T3::ptr &face_tree,
                            const std::vector<vec3> &pov,
                            const std::vector<vec3> N_pov, //
                            real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->indices();
  std::vector<index_t> face_ids(face_vert_ids.size() / 3);

  for (int i = 0; i < face_ids.size(); i++) {
    face_ids[i] = i;
  }

  std::vector<vec3> N(face_vert_ids.size() / 3);
  real l0 = 0.0;
  for (int i = 0; i < face_ids.size(); i++) {
    vec3 x0 = x[face_vert_ids[i * 3 + 0]];
    vec3 x1 = x[face_vert_ids[i * 3 + 1]];
    vec3 x2 = x[face_vert_ids[i * 3 + 2]];
    N[i] = (x1 - x0).cross(x2 - x0);
    real A = N[i].norm();
    l0 += A;
  }

  l0 /= face_ids.size();
  l0 = pow(l0, 0.5);
  l0 *= spread;

  calder::fast_summation<arp::T3> sum(*face_tree);
  sum.bind<vec3>(face_ids, N);
  std::vector<real> dists(pov.size(), 0.0);
  std::vector<real> u = sum.calc<real>(
      pov,
      [l0, &dists, &N_pov](const index_t &i, const index_t &j, const vec3 &pi,
                           const std::vector<calder::datum::ptr> &data,
                           calder::fast_summation<arp::T3>::Node_Type node_type,
                           const arp::T3::node &node,
                           const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->leaf_data()[j];
        const vec3 &Ni = N_pov[i];
        vec3 p0 = tree.vert(3 * j + 0);
        vec3 p1 = tree.vert(3 * j + 1);
        vec3 p2 = tree.vert(3 * j + 2);

        std::array<real, 4> cp = va::closest_point({p0, p1, p2}, pi);
        dists[i] = std::max(dists[i], cp[0]);
        return 0.0;
      },
      [l0, &dists, &N_pov](const index_t &i, const index_t &j, const vec3 &pi,
                           const std::vector<calder::datum::ptr> &data,
                           calder::fast_summation<arp::T3>::Node_Type node_type,
                           const arp::T3::node &node,
                           const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        const vec3 &Ni = N_pov[i];
        vec3 pj = node.center();
        vec3 dp = pj - pi;

        if (Ni.dot(N) < 0.0)
          return 0.0;
        // real dist = va::project_to_nullspace(dp, N);
        real dist = va::norm(dp);
        dists[i] = std::max(dists[i], dist);

        return 0.0;
      });

  return dists;
}

std::vector<real> fast_winding(asawa::shell::shell &M,
                               const std::vector<vec3> &x,
                               const std::vector<vec3> &pov, real l0) {

  std::vector<vec3> N = asawa::shell::face_normals(M, x);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  std::vector<vec3> wN(N);
  real total_area = 0.0;
  for (int i = 0; i < wN.size(); i++)
    total_area += weights[i];
  std::cout << "sum = " << total_area << std::endl;

  for (int i = 0; i < wN.size(); i++)
    wN[i] *= weights[i];

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
        return 0.25 / M_PI * va::solidAngle(pi, p0, p1, p2);
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
        real kappa = computeK(dist, 0.0, 3.0);
        return 0.25 / M_PI * kappa * va::dot(N, dp);
      });

  return u;
}

std::vector<mat3> fast_frame(asawa::shell::shell &M, const std::vector<vec3> &x,
                             const std::vector<vec3> &p_pov,
                             const std::vector<vec3> &p_normals, real l0) {
  // this is potentially broken...
  // the summation iterates over an edge range, but the weights
  // aren't over the range, they are over the edges themselves
  // so the weights are not aligned with the edge range

  std::vector<vec3> E = asawa::shell::edge_tangents(M, x);
  std::vector<real> w = asawa::shell::edge_cotan_weights(M, x);
  std::vector<real> wa = asawa::shell::edge_areas(M, x);

  std::vector<vec3> wE(E);

  for (int i = 0; i < wE.size(); i++) {
    // wE[i] *= w[i]; // * wE[i].normalized();
    wE[i] = wa[i] * w[i] * E[i];
  }

  std::vector<index_t> edge_verts = M.get_edge_vert_ids();
  std::vector<index_t> edge_map = M.get_edge_map();
  std::vector<index_t> edge_ids = M.get_edge_range_2();

  arp::T2::ptr edge_tree = arp::aabb_tree<2>::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  std::cout << " ==== wE.size(): " << wE.size() << std::endl;
  sum.bind(calder::edge_frame_datum::create(edge_ids, wE));
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<mat3> u = sum.calc<mat3>(
      p_pov,
      [l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
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
        real kappa = computeK(dist, l0, 3.0);
        sums[i] += kappa;
        return kappa * e * e.transpose();
      },
      [l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                  const std::vector<calder::datum::ptr> &data,
                  calder::fast_summation<arp::T2>::Node_Type node_type,
                  const arp::T2::node &node, const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);
        const mat3 &E = F_datum->node_data()[j];

        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0, 3.0);
        sums[i] += kappa;
        return kappa * E;
      });
#if 1

  std::vector<mat3> Us(p_pov.size());
  for (int i = 0; i < p_pov.size(); i++) {
    const vec3 &Ni = p_normals[i];
    mat3 R = va::rejection_matrix(Ni);
    // gg::geometry_logger::frame(R * u[i], x[i], 10.0);
    mat3 Ui = 1.0 / sums[i] * R * u[i];
    Eigen::JacobiSVD<mat3> svd(Ui, Eigen::ComputeFullU | Eigen::ComputeFullV);
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