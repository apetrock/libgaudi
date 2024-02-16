
#ifndef __DUCHAMP_TUNNEL_MODULE__
#define __DUCHAMP_TUNNEL_MODULE__

#include "gaudi/asawa/datums.hpp"
#include "gaudi/calder/shell_integrators.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace gaudi {

// probably doesn't need to be a module, could go straight in test class...
/////////////////////////////////////////////////////////////
namespace calder {

std::vector<std::vector<vec3>> tunnel(asawa::shell::shell &M,
                                      const std::vector<vec3> &start_points,
                                      const std::vector<vec3> &Np, real l0,
                                      real p) {

  std::vector<vec3> &x = asawa::get_vec_data(M, 0);
  std::vector<index_t> face_vert_ids = M.get_face_vert_ids();
  std::vector<index_t> face_map = M.get_face_map();
  std::vector<index_t> face_ids = M.get_face_range();

  std::vector<vec3> Nc = asawa::shell::face_normals(M, x);
  std::vector<real> weights = asawa::shell::face_areas(M, x);
  real total_weight = 0.0;
  for (int i = 0; i < weights.size(); i++) {
    Nc[i] *= weights[i];
    total_weight += weights[i];
  }

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "summing" << std::endl;
  std::cout << " -n_faces: " << face_ids.size() << std::endl;
  std::cout << " -create: " << std::endl;

  Shell_Tree_Type::ptr face_tree = arp::T3::create(face_vert_ids, x, 16);
  // calder::test_extents(*face_tree);
  //  calder::test_pyramid_scalar(*face_tree, face_ids, weights);

  std::cout << " -sum: " << std::endl;
  Shell_Sum_Type sum(*face_tree);

  sum.bind(calder::scalar_datum::create(face_ids, weights));
  sum.bind(calder::vec3_datum::create(face_ids, Nc));

  std::vector<std::vector<vec3>> paths;
  for (int i = 0; i < start_points.size(); i++) {
    std::cout << " -compute path: " << std::endl;
    vec3 Ni = Np[i];
    vec3 xi = start_points[i];
    std::vector<vec3> path = {xi};
    int k = 0;
    bool iter = true;
    while (k < 400 && iter) {
      vec3 Nh = vec3::Zero();
      real w = 0.0;
      real ww = 0.0;
      real wa = 0.0;
      real wn = 0.0;
      auto compute = [&](const index_t i, const index_t j, //
                         const vec3 &pi, const vec3 &pj,
                         const std::vector<calder::datum::ptr> &data,
                         Shell_Sum_Type::Node_Type node_type, //
                         const Shell_Sum_Type::Node &node,    //
                         const Shell_Sum_Type::Tree &tree) -> vec3 {
        real wj = get_data<real>(node_type, j, 0, data);
        vec3 Nj = get_data<vec3>(node_type, j, 1, data);
        vec3 dp = pj - pi;

        ww += compute3K(dp.norm()) * dp.dot(Nj);

        wn += Nj.norm();

        Ni = Ni.normalized();
        Nj = Nj.normalized();
        real dpNi = dp.dot(Ni);
        real dpNj = dp.dot(Nj);

        if (dpNi < 0.0 || dpNj < 0.0) {
          return vec3::Zero();
        }

        real ndp = dp.squaredNorm();
        vec3 Ndp = (Ni * Ni.transpose() * dp);
        real nNdp = Ndp.norm();

        real r = 0.5 * ndp / nNdp;
        vec3 cen = pj - r * Nj;

        real kappa = calc_inv_dist(dp, l0, p);
        w += wj * kappa;
        Nh += wj * kappa * dp;
        // Nh += wj * kappa * (cen - pi);
        //     logger::line(pj, pj + Ndp, vec4(0.0, 1.0, 0.5, 1.0));
        //    logger::line(pi, pi + dp, vec4(1.0, 1.0, 0.0, 1.0));
        return vec3::Zero();
      };

      std::vector<vec3> us = sum.calc<vec3>(
          {xi},
          [&wa, &compute](const index_t &i, const index_t &j, const vec3 &pi,
                          const std::vector<calder::datum::ptr> &data,
                          Shell_Sum_Type::Node_Type node_type,
                          const Shell_Tree_Type::node &node,
                          const Shell_Tree_Type &tree) -> vec3 {
            vec3 x0 = tree.vert(3 * j + 0);
            vec3 x1 = tree.vert(3 * j + 1);
            vec3 x2 = tree.vert(3 * j + 2);
            vec3 pj;
            real dist = va::distance_from_triangle({x0, x1, x2}, pi, pj);
            pj = 0.333 * (x0 + x1 + x2);

            const calder::scalar_datum::ptr N_datum =
                static_pointer_cast<calder::scalar_datum>(data[0]);
            const real wj = N_datum->leaf_data()[j];
            wa += wj;
            return compute(i, j, pi, pj, data, node_type, node, tree);
          },
          [&wa, &compute](const index_t &i, const index_t &j, const vec3 &pi,
                          const std::vector<calder::datum::ptr> &data,
                          Shell_Sum_Type::Node_Type node_type,
                          const Shell_Tree_Type::node &node,
                          const Shell_Tree_Type &tree) -> vec3 {
            vec3 pj = node.center();

            const calder::scalar_datum::ptr N_datum =
                static_pointer_cast<calder::scalar_datum>(data[0]);
            const real &wj = N_datum->node_data()[j];
            wa += wj;
            return compute(i, j, pi, pj, data, node_type, node, tree);
          },
          0.5, false);

      // std::cout << "w: " << w << " winding: " << ww << " area: " << wa
      //           << " normal: " << wn << " total_area: " << total_weight
      //           << std::endl;
      Nh /= w;

      if (ww < 0.0) {
        iter = false;
      }

      Nh.normalize();
      xi += 0.01 * Nh;
      Ni += va::mix(0.95, Ni, Nh);
      Ni.normalize();
      path.push_back(xi);
      k++;
    }
    paths.push_back(path);
    for (int i = 0; i < path.size() - 1; i++) {
      logger::line(path[i], path[i + 1], vec4(1.0, 0.0, 1.0, 1.0));
    }
  }
  return paths;
}

} // namespace calder
namespace duchamp {

class tunnel : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(tunnel)
  tunnel(asawa::shell::shell::ptr M) : module_base_shell(M) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
    _eps = asawa::shell::avg_length(*M, x);
  };

  virtual ~tunnel(){

  };

  virtual void step(real h) {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
  }

  std::vector<std::vector<vec3>> drill(const std::vector<vec3> &start_points,
                                       const std::vector<vec3> &headings) {
    asawa::shell::shell &M = *_M;
    std::cout << "start_points: " << start_points.size() << std::endl;
    std::vector<std::vector<vec3>> paths =
        calder::tunnel(M, start_points, headings, _eps, 1.5);
    return paths;
  }

  real t = 0.0;
  real _eps;
};

} // namespace duchamp
} // namespace gaudi

#endif