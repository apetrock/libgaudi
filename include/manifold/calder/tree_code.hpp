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
#include "manifold/geometry_types.hpp"
#include <vector>

namespace calder {

template <typename TREE>
void test_extents(const TREE &tree,                    //
                  const std::vector<index_t> &indices, //
                  const std::vector<vec3> &x) {

  std::vector<arp::extents_t> extents = arp::build_extents(tree, indices, x);
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
  typedef typename TREE::node node;

  fast_summation(const TREE &tree) : __tree(tree) {}

  void add_datum(datum &x) { __data.push_back(x); }
  template <typename T> void add_datum(std::vector<T> &x) {
    __data.push_back(datum_t<T>::create(x));
  }

  void set_threshold(real t) { __thresh = t; }
  template <typename Q> std::vector<Q> calc(const std::vector<vec3> &pov) {
    for (int i = 0; i < __data.size(); i++) {
      __data[i].pyramid(__tree);
    }

    vector<Q> u(pov.size(), z::zero<Q>());

    std::stack<int> stack;
    stack.push(0);

    for (int i = 0; i < pov.size(); i++) {
      int count = 0;

      vec3 pi = pov[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int j = stack1.top();
        stack1.pop();
        node &pNode = __tree.nodes[j];
        vec3 pj = pNode.center();

        real dc = va::dist(pi, pj);
        // T sc = va::norm(pNode.half);
        real sc = pNode.mag();
        // int ii = octree.permutation[pNode.begin];

        if (sc / dc < __thresh || pNode.isLeaf()) {
          if (pNode.isLeaf()) {
            for (int jn = pNode.begin; jn < pNode.begin + pNode.size; jn++) {
              int jj = __tree.permutation[jn];
              // Q ui = leafComputeFcn(i, jj, pi, __data, pNode, __tree);
              // u[i] += ui;
            }

          } else {
            // Q ui = nodeComputeFcn(i, j, pi, __data, pNode, __tree);
            // u[i] += ui;
          }
        }

        else {
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
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
  std::vector<datum> __data;
  const TREE &__tree;
  const std::vector<vec3> &__vert_x;
};
} // namespace calder
#endif