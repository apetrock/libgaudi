//
//  m2Includes.h
//  Manifold
//
//  Created by John Delaney on 5/22/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
// includes files that are required for the data structure only, not derived
// files that USE the data structure
#ifndef __M2TREE_CODE__
#define __M2TREE_CODE__

#include "datums.hpp"
#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include <vector>

namespace gaudi {
namespace calder {

template <typename TREE>
void test_extents(const TREE &tree,                    //
                  const std::vector<index_t> &indices, //
                  const std::vector<vec3> &x) {

  std::vector<ext::extents_t> extents = arp::build_extents(tree, indices, x);
  for (const auto &ext : extents) {
    vec4 c(0.5, 0.5, 0.1, 1.0);
    gg::geometry_logger::ext(ext[0], ext[1], c);
  }
}

template <typename TREE>
void test_pyramid(const TREE &tree,                      //
                  const std::vector<index_t> &q_indices, //
                  const std::vector<vec3> &q,            //
                  const std::vector<real> &q_weights) {

  std::vector<vec3> wq(q);
  for (int i = 0; i < q.size(); i++)
    wq[i] *= q_weights[i];

  datum_t<vec3>::ptr x_datum = datum_t<vec3>::create(q_indices, wq);
  x_datum->pyramid(tree);
  for (int i = 0; i < tree.nodes.size(); i++) {
    vec3 cen = tree.nodes[i].center();
    vec3 N = x_datum->__tree_data[i];
    vec4 c(0.0, 0.5, 0.8, 1.0);
    gg::geometry_logger::line(cen, cen + N, c);
  }
}

template <typename TREE> class fast_summation {
public:
  typedef typename TREE::node NODE;

  fast_summation(const TREE &tree) : __tree(tree) {}

  void bind(const datum::ptr &x) { __data.push_back(x); }

  template <typename T>
  void bind(const std::vector<index_t> &indices, const std::vector<T> &x) {
    __data.push_back(datum_t<T>::create(indices, x));
  }

  void set_threshold(real t) { __thresh = t; }

  template <typename Q>
  using ComputeFcn =
      std::function<Q(const index_t &, const index_t &, const vec3 &,
                      const std::vector<datum::ptr> &, const NODE &,
                      const TREE &)>;

  template <typename Q>
  std::vector<Q> calc(const std::vector<vec3> &pov,
                      ComputeFcn<Q> leafComputeFcn,
                      ComputeFcn<Q> nodeComputeFcn) {
    for (int i = 0; i < __data.size(); i++) {
      __data[i]->pyramid(__tree);
    }

    vector<Q> u(pov.size(), z::zero<Q>());

    std::vector<ext::extents_t> extents =
        arp::build_extents(__tree, __tree.indices(), __tree.verts());

    for (int i = 0; i < pov.size(); i++) {
      vec3 pi = pov[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int j = stack1.top();
        stack1.pop();
        const NODE &pNode = __tree.nodes[j];
        vec3 pj = pNode.center();

        real dc = va::dist(pi, pj);
        // T sc = va::norm(pNode.half);
        // real sc = pNode.mag();
        ext::extents_t ext = extents[j];
        vec3 de = ext[1] - ext[0];
        real sc = min(de[0], min(de[1], de[2]));
        //    std::cout << sc << " " << dc << std::endl;
        //     int ii = octree.permutation[pNode.begin];

        if (sc / dc < __thresh || pNode.isLeaf()) {

          vec4 c(0.0, 0.5, 0.8, 1.0);
          if (i == 0) {
            gg::geometry_logger::ext(ext[0], ext[1], c);
            gg::geometry_logger::line(pi, pj, vec4(0.1, 0.7, 0.2, 0.5));
          }

          if (pNode.isLeaf()) {

            for (int jn = pNode.begin; jn < pNode.begin + pNode.size; jn++) {
              int jj = __tree.permutation[jn];
              Q ui = leafComputeFcn(i, jj, pi, __data, pNode, __tree);
              u[i] += ui;
            }

          } else {
            Q ui = nodeComputeFcn(i, j, pi, __data, pNode, __tree);
            u[i] += ui;
          }
        }

        else {
          for (int j = 0; j < pNode.getNumChildren(); j++) {
            if (pNode.children[j] > -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }
      }
    }

    // for(int i = 0; i < 5; i++){
    //  std::cout << u[i].transpose() << std::endl;
    //}

    return u;
  }
  real __thresh = 0.5;
  std::vector<datum::ptr> __data;
  const TREE &__tree;
};
} // namespace calder
} // namespace gaudi
#endif