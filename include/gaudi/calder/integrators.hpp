#ifndef __M2HARMONIC_INTEGRATOR__
#define __M2HARMONIC_INTEGRATOR__

#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/arp/arp.h"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/datums.hpp"

#include "weight_functions.hpp"

#include "tree_code.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>
#include "gaudi/geometry_logger.hpp"

namespace gaudi {

namespace calder {

std::vector<real> fast_winding(const arp::T3::ptr &face_tree,
                               const std::vector<vec3> &pov,
                               real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->adjacency();
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
  std::vector<vec3> Nc = asawa::shell::compress_to_range<vec3>(face_ids, N);
  sum.bind<vec3>(face_ids, Nc);

  std::vector<real> u = sum.calc<real>(
      pov,
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
           const arp::T3 &tree) -> real {
        auto simplex = tree.leaf_simplex(j);
        return 0.25 / M_PI * va::solidAngle(pi, simplex[0], simplex[1], simplex[2]);
      },
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
           const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        const ext::extents_t &ext = tree.bvh_.internal[j];
        vec3 pj = 0.5 * (ext[0] + ext[1]);
        vec3 dp = pj - pi;
        real kappa = calc_inv_dist(dp, 0.0, 3.0);
        return 0.25 / M_PI * kappa * va::dot(N, dp);
      },
      0.25);

  return u;
}

std::vector<real> fast_dist(const arp::T3::ptr &face_tree,
                            const std::vector<vec3> &pov, real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->adjacency();
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
  std::vector<vec3> Nc = asawa::shell::compress_to_range<vec3>(face_ids, N);
  sum.bind<vec3>(face_ids, Nc);
  std::vector<real> min_dists(pov.size(),
                              std::numeric_limits<real>::max());
  std::vector<real> u = sum.calc<real>(
      pov,
      [l0, &min_dists](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       calder::fast_summation<arp::T3>::Node_Type node_type,
                       const arp::T3 &tree) -> real {
        auto simplex = tree.leaf_simplex(j);
        std::array<real, 4> cp = va::closest_point({simplex[0], simplex[1], simplex[2]}, pi);
        min_dists[i] = std::min(min_dists[i], cp[0]);
        return 0.0;
      },
      [l0, &min_dists](const index_t &i, const index_t &j, const vec3 &pi,
                       const std::vector<calder::datum::ptr> &data,
                       calder::fast_summation<arp::T3>::Node_Type node_type,
                       const arp::T3 &tree) -> real {
        const ext::extents_t &ext = tree.bvh_.internal[j];
        vec3 pj = 0.5 * (ext[0] + ext[1]);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        min_dists[i] = std::min(min_dists[i], dist);
        return 0.0;
      },
      0.25);

  return min_dists;
}

std::vector<vec3> fast_dist_gradient(const arp::T3::ptr &face_tree,
                                     const std::vector<vec3> &pov,
                                     real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->adjacency();
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
  std::vector<vec3> Nc = asawa::shell::compress_to_range<vec3>(face_ids, N);
  sum.bind<vec3>(face_ids, Nc);
  std::vector<real> dists(pov.size(), 999.9);
  std::vector<vec3> normals(pov.size(), vec3(0.0, 0.0, 0.0));
  std::vector<real> W(pov.size(), 0.0);
  std::vector<real> u = sum.calc<real>(
      pov,
      [l0, &W, &normals,
       &dists](const index_t &i, const index_t &j, const vec3 &pi,
               const std::vector<calder::datum::ptr> &data,
               calder::fast_summation<arp::T3>::Node_Type node_type,
               const arp::T3 &tree) -> real {
        auto simplex = tree.leaf_simplex(j);
        std::array<real, 4> cp = va::closest_point({simplex[0], simplex[1], simplex[2]}, pi);
        vec3 pT = cp[1] * simplex[0] + cp[2] * simplex[1] + cp[3] * simplex[2];
        real dist = cp[0];
        if (dist < dists[i]) {
          dists[i] = dist;
          normals[i] = (pT - pi).normalized();
        }
        return 0.0;
      },
      [l0, &W, &normals,
       &dists](const index_t &i, const index_t &j, const vec3 &pi,
               const std::vector<calder::datum::ptr> &data,
               calder::fast_summation<arp::T3>::Node_Type node_type,
               const arp::T3 &tree) -> real {
        const ext::extents_t &ext = tree.bvh_.internal[j];
        vec3 pj = 0.5 * (ext[0] + ext[1]);
        vec3 dp = pj - pi;
        real dist = va::norm(dp);

        if (dist < dists[i]) {
          dists[i] = dist;
          normals[i] = dp.normalized();
        }
        return 0.0;
      },
      0.25);
  for (int i = 0; i < normals.size(); i++) {
    if (W[i] > 0.0)
      normals[i] /= W[i];
    normals[i].normalize();
  }
  return normals;
}

std::vector<real> fast_view(const arp::T3::ptr &face_tree,
                            const std::vector<vec3> &pov,
                            const std::vector<vec3> N_pov,
                            real spread = 1.0) {

  const std::vector<vec3> x = face_tree->verts();
  const std::vector<index_t> &face_vert_ids = face_tree->adjacency();
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
  std::vector<vec3> Nc = asawa::shell::compress_to_range<vec3>(face_ids, N);
  sum.bind<vec3>(face_ids, Nc);
  std::vector<real> dists(pov.size(), 0.0);
  std::vector<real> u = sum.calc<real>(
      pov,
      [l0, &dists, &N_pov](const index_t &i, const index_t &j, const vec3 &pi,
                           const std::vector<calder::datum::ptr> &data,
                           calder::fast_summation<arp::T3>::Node_Type node_type,
                           const arp::T3 &tree) -> real {
        auto simplex = tree.leaf_simplex(j);
        std::array<real, 4> cp = va::closest_point({simplex[0], simplex[1], simplex[2]}, pi);
        dists[i] = std::max(dists[i], cp[0]);
        return 0.0;
      },
      [l0, &dists, &N_pov](const index_t &i, const index_t &j, const vec3 &pi,
                           const std::vector<calder::datum::ptr> &data,
                           calder::fast_summation<arp::T3>::Node_Type node_type,
                           const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        const vec3 &Ni = N_pov[i];
        const ext::extents_t &ext = tree.bvh_.internal[j];
        vec3 pj = 0.5 * (ext[0] + ext[1]);
        vec3 dp = pj - pi;

        if (Ni.dot(N) < 0.0)
          return 0.0;
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
  for (int i = 0; i < wN.size(); i++)
    wN[i] *= weights[i];

  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  std::vector<index_t> face_map = M.get_face_map();
  auto face_ids = M.get_face_range();
  std::vector<index_t> face_ix(face_ids.begin(), face_ids.end());

  arp::T3::ptr face_tree = arp::T3::create(face_vert_ids, x, 24);

  calder::fast_summation<arp::T3> sum(*face_tree);
  std::vector<vec3> Nc = asawa::shell::compress_to_range<vec3>(face_ix, wN);

  sum.bind<vec3>(face_ix, Nc);

  std::vector<vec3> centroids = asawa::shell::face_centers(M, x);
  auto com = calder::com_datum::create(face_ix, weights, centroids);
  com->pyramid(*face_tree);

  std::vector<real> u = sum.calc<real>(
      pov,
      [l0](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
           const arp::T3 &tree) -> real {
        auto simplex = tree.leaf_simplex(j);
        return 0.25 / M_PI * va::solidAngle(pi, simplex[0], simplex[1], simplex[2]);
      },
      [l0, &com](const index_t &i, const index_t &j, const vec3 &pi,
           const std::vector<calder::datum::ptr> &data,
           calder::fast_summation<arp::T3>::Node_Type node_type,
           const arp::T3 &tree) -> real {
        const calder::vec3_datum::ptr N_datum =
            static_pointer_cast<calder::vec3_datum>(data[0]);
        const vec3 &N = N_datum->node_data()[j];
        vec3 pj = com->get_node_com(j);
        vec3 dp = pj - pi;
        real kappa = calc_inv_dist(dp, 0.0, 3.0);
        return 0.25 / M_PI * kappa * va::dot(N, dp);
      },
      0.2);

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
    wE[i] = wa[i] * w[i] * E[i];
  }

  std::vector<index_t> edge_verts = M.get_edge_vert_ids();
  std::vector<index_t> edge_map = M.get_edge_map();
  std::vector<index_t> edge_ids = M.get_edge_range_2();

  arp::T2::ptr edge_tree = arp::T2::create(edge_verts, x, 12);

  calder::fast_summation<arp::T2> sum(*edge_tree);
  std::cout << " ==== wE.size(): " << wE.size() << std::endl;
  std::vector<vec3> Ec = asawa::shell::compress_to_range<vec3>(edge_ids, wE);
  sum.bind(calder::edge_frame_datum::create(edge_ids, Ec));
  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<mat3> u = sum.calc<mat3>(
      p_pov,
      [l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                  const std::vector<calder::datum::ptr> &data,
                  calder::fast_summation<arp::T2>::Node_Type node_type,
                  const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);

        const vec3 &e = F_datum->sorted_leaf_data()[j];
        auto simplex = tree.leaf_simplex(j);
        vec3 pj = va::project_on_line(simplex[0], simplex[1], pi);
        vec3 dp = pj - pi;
        real kappa = calc_inv_dist(dp, l0, 3.0);
        sums[i] += kappa;
        return kappa * e * e.transpose();
      },
      [l0, &sums](const index_t &i, const index_t &j, const vec3 &pi,
                  const std::vector<calder::datum::ptr> &data,
                  calder::fast_summation<arp::T2>::Node_Type node_type,
                  const arp::T2 &tree) -> mat3 {
        const calder::edge_frame_datum::ptr F_datum =
            static_pointer_cast<calder::edge_frame_datum>(data[0]);
        const mat3 &E = F_datum->node_data()[j];

        const ext::extents_t &ext = tree.bvh_.internal[j];
        vec3 pj = 0.5 * (ext[0] + ext[1]);
        vec3 dp = pj - pi;
        real kappa = calc_inv_dist(dp, l0, 3.0);
        sums[i] += kappa;
        return kappa * E;
      });

  std::vector<mat3> Us(p_pov.size());
  for (int i = 0; i < p_pov.size(); i++) {
    const vec3 &Ni = p_normals[i];
    mat3 R = va::rejection_matrix(Ni);
    mat3 Ui = 1.0 / sums[i] * R * u[i];
    Eigen::JacobiSVD<mat3> svd(Ui, Eigen::ComputeFullU | Eigen::ComputeFullV);
    mat3 U = svd.matrixU();

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
}

} // namespace calder
} // namespace gaudi
#endif
