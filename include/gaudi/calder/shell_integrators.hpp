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

template <typename T>
std::vector<T> mls_avg(asawa::shell::shell &M, const std::vector<T> &v,
                       const std::vector<vec3> &p_pov, real l0, real p = 3.0) {

  std::vector<vec3> &x = asawa::get_vec_data(M, 0);
  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  std::vector<index_t> face_map = M.get_face_map();
  std::vector<index_t> face_ids = M.get_face_range();
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "summing" << std::endl;
  std::cout << " -n_faces: " << face_ids.size() << " " << v.size() << std::endl;

  std::cout << " -create: " << std::endl;

  arp::T3::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);
  std::cout << " -sum: " << std::endl;
  calder::fast_summation<arp::T3> sum(*face_tree);
  std::cout << " -bind: " << std::endl;

  sum.bind(calder::datum_t<T>::create(face_ids, v));

  real max_dist = 0.0;
  int max_dist_i = 0;
  std::cout << " -compute: " << std::endl;

  std::vector<real> sums(p_pov.size(), 0.0);
  std::vector<T> us = sum.calc<T>(
      p_pov,
      [l0, &sums, &max_dist, &max_dist_i,
       p](const index_t &i, const index_t &j, const vec3 &pi,
          const std::vector<calder::datum::ptr> &data,
          calder::fast_summation<arp::T3>::Node_Type node_type,
          const arp::T3::node &node, const arp::T3 &tree) -> T {
        const typename calder::datum_t<T>::ptr F_datum =
            static_pointer_cast<typename calder::datum_t<T>>(data[0]);

        T e = F_datum->leaf_data()[j];
        vec3 x0 = tree.vert(3 * j + 0);
        vec3 x1 = tree.vert(3 * j + 1);
        vec3 x2 = tree.vert(3 * j + 2);
        vec3 pj;

        real dist = va::distance_from_triangle({x0, x1, x2}, pi, pj);
        vec3 dp = pj - pi;
        real kappa = computeK(dist, l0, p);

        sums[i] += kappa * e.norm();
        if (std::isnan(sums[i])) {
          std::cout << "leaf: " << j << std::endl;
          std::cout << "nan" << std::endl;
          std::cout << "e: " << e << std::endl;
          std::cout << "kappa: " << kappa << std::endl;
          std::cout << "dist: " << dist << std::endl;
          std::cout << "dp: " << dp << std::endl;
          std::cout << "pi: " << pi << std::endl;
          std::cout << "pj: " << pj << std::endl;
          exit(0);
        }
        return kappa * e;
      },
      [l0, &sums, p](const index_t &i, const index_t &j, const vec3 &pi,
                     const std::vector<calder::datum::ptr> &data,
                     calder::fast_summation<arp::T3>::Node_Type node_type,
                     const arp::T3::node &node, const arp::T3 &tree) -> T {
        const typename calder::datum_t<T>::ptr F_datum =
            static_pointer_cast<typename calder::datum_t<T>>(data[0]);

        T e = F_datum->node_data()[j];
        vec3 pj = node.center();
        vec3 dp = pj - pi;
        real dist = va::norm(dp);
        real kappa = computeK(dist, l0, p);

        // std::cout << kappa << " " << e.norm() << std::endl;
        // if (i == 6664) {
        //  gg::geometry_logger::line(pi, pj, vec4(0.5, 0.0, 1.0, 1.0));
        //  gg::geometry_logger::line(pj, pj + 0.1 * kappa * e,
        //                           vec4(0.0, 1.0, 0.1, 1.0));
        //}
        sums[i] += kappa * e.norm();
        if (std::isnan(sums[i])) {

          std::cout << "node: " << j << std::endl;
          std::cout << "nan" << std::endl;
          std::cout << "e: " << e << std::endl;
          std::cout << "kappa: " << kappa << std::endl;
          std::cout << "dist: " << dist << std::endl;
          std::cout << "dp: " << dp << std::endl;
          std::cout << "pi: " << pi << std::endl;
          std::cout << "pj: " << pj << std::endl;
          exit(0);
        }
        return kappa * e;
      },
      0.5, false);
#if 1
  int max_count = 0;
  int max_count_i = 0;
  for (int i = 0; i < p_pov.size(); i++) {
    if (sums[i] < 1e-8) {
      us[i] = vec3::Zero();
      continue;
    }
    us[i] /= sums[i];
  }
#endif
  return us;
}

} // namespace calder
} // namespace gaudi
#endif