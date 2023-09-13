//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __SHELL_INTEGRATOR__
#define __SHELL_INTEGRATOR__

#include "integrators.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <vector>

namespace gaudi {

namespace calder {

using Shell_Tree_Type = arp::T3;
using Shell_Sum_Type = calder::fast_summation<Shell_Tree_Type>;

using Shell_Bind_Fcn =
    std::function<void(const std::vector<index_t> &, Shell_Sum_Type &)>;

template <typename Q>
using Shell_Compute_Fcn =
    std::function<Q(const index_t &i, const index_t &j, const vec3 &,
                    const vec3 &, const std::vector<datum::ptr> &,
                    Shell_Sum_Type::Node_Type &, const Shell_Tree_Type::node &,
                    const Shell_Tree_Type &)>;

template <typename T>
T get_data(Shell_Sum_Type::Node_Type node_type, index_t j, index_t data_id,
           const std::vector<calder::datum::ptr> &data) {
  const typename calder::datum_t<T>::ptr F_datum =
      static_pointer_cast<typename calder::datum_t<T>>(data[data_id]);
  if (node_type == Shell_Sum_Type::Node_Type::LEAF) {
    return F_datum->leaf_data()[j];
  } else {
    return F_datum->node_data()[j];
  }
}

template <typename T>
std::vector<T>
integrate_over_shell(asawa::shell::shell &M, const std::vector<vec3> &p_pov,
                     Shell_Bind_Fcn bind_fcn = nullptr,
                     Shell_Compute_Fcn<T> compute_fcn = nullptr) {

  std::vector<vec3> &x = asawa::get_vec_data(M, 0);
  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  std::vector<index_t> face_map = M.get_face_map();
  std::vector<index_t> face_ids = M.get_face_range();
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "summing" << std::endl;
  std::cout << " -n_faces: " << face_ids.size() << std::endl;
  std::cout << " -create: " << std::endl;

  Shell_Tree_Type::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);
  std::cout << " -sum: " << std::endl;
  Shell_Sum_Type sum(*face_tree);
  bind_fcn(face_ids, sum);

  std::cout << " -compute: " << std::endl;
  std::vector<T> us = sum.calc<T>(
      p_pov,
      [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                     const std::vector<calder::datum::ptr> &data,
                     Shell_Sum_Type::Node_Type node_type,
                     const Shell_Tree_Type::node &node,
                     const Shell_Tree_Type &tree) -> T {
        vec3 x0 = tree.vert(3 * j + 0);
        vec3 x1 = tree.vert(3 * j + 1);
        vec3 x2 = tree.vert(3 * j + 2);
        vec3 pj;
        real dist = va::distance_from_triangle({x0, x1, x2}, pi, pj);
        pj = 0.333 * (x0 + x1 + x2);
        return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
      },
      [&compute_fcn](const index_t &i, const index_t &j, const vec3 &pi,
                     const std::vector<calder::datum::ptr> &data,
                     Shell_Sum_Type::Node_Type node_type,
                     const Shell_Tree_Type::node &node,
                     const Shell_Tree_Type &tree) -> T {
        vec3 pj = node.center();
        return compute_fcn(i, j, pi, pj, data, node_type, node, tree);
      },
      0.5, false);
  return us;
}

template <typename T>
std::vector<T> mls_avg(asawa::shell::shell &M, const std::vector<T> &v,
                       const std::vector<vec3> &p_pov, real l0, real p = 3.0) {

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<T> us = integrate_over_shell<T>(
      M, p_pov,
      [&v](const std::vector<index_t> &face_ids, Shell_Sum_Type &sum) {
        sum.bind(calder::datum_t<T>::create(face_ids, v));
      },
      [l0, &sums, p](const index_t i, const index_t j, //
                     const vec3 &pi, const vec3 &pj,
                     const std::vector<calder::datum::ptr> &data,
                     Shell_Sum_Type::Node_Type node_type, //
                     const Shell_Sum_Type::Node &node,    //
                     const Shell_Sum_Type::Tree &tree) -> T {
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
    if (sums[i] < 1e-6)
      continue;
    us[i] /= sums[i];
  }
#endif
  return us;
}

std::vector<vec3> vortex_force(asawa::shell::shell &M,
                               const std::vector<vec3> &p_pov,
                               const std::vector<vec3> &omega, real l0,
                               real p = 3.0) {

  std::vector<vec3> us = integrate_over_shell<vec3>(
      M, p_pov,
      [&omega](const std::vector<index_t> &edge_ids, Shell_Sum_Type &sum) {
        sum.bind(calder::vec3_datum::create(edge_ids, omega));
      },
      [l0, p](const index_t i, const index_t j, //
              const vec3 &pi, const vec3 &pj,
              const std::vector<calder::datum::ptr> &data,
              Shell_Sum_Type::Node_Type node_type, //
              const Shell_Sum_Type::Node &node,    //
              const Shell_Sum_Type::Tree &tree) -> vec3 {
        vec3 w = get_data<vec3>(node_type, j, 0, data);
        vec3 dp = pj - pi;
        real dist = dp.norm();
        // dist = std::max(dist, l0);
        real kappa = computeKm(dist, l0, p);

        return -kappa * dp.cross(w);
      });
  return us;
}

} // namespace calder
} // namespace gaudi
#endif