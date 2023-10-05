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

template <typename TREE> void test_extents(const TREE &tree) {

  std::vector<ext::extents_t> extents =
      arp::build_extents(tree, tree.indices(), tree.verts());
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

template <typename TREE>
void test_pyramid_scalar(const TREE &tree,                      //
                         const std::vector<index_t> &q_indices, //
                         const std::vector<real> &q) {

  datum_t<real>::ptr x_datum = datum_t<real>::create(q_indices, q);
  x_datum->pyramid(tree);
  for (int i = 0; i < tree.nodes.size(); i++) {
    typename TREE::node pNode = tree.nodes[i];
    real w = x_datum->__tree_data[i];
    if (pNode.isLeaf()) {
      real wc = 0.0;
      for (int jn = pNode.begin; jn < pNode.begin + pNode.size; jn++) {
        int jj = tree.permutation[jn];
        // wc += x_datum->leaf_data()[jj];
        wc += x_datum->leaf_data()[jj];
      }
      if (w - wc > 1e-6) {
        std::cout << "leaf check: " << pNode.level << " " << pNode.size << " "
                  << w << " " << wc << std::endl;
      }
    } else {
      real wc = 0.0;
      for (int j = 0; j < pNode.getNumChildren(); j++) {
        if (pNode.children[j] > -1) {
          wc += x_datum->node_data()[pNode.children[j]];
        }
      }
      if (w - wc > 1e-6) {
        std::cout << "node check: " << pNode.level << " "
                  << pNode.getNumChildren() << " " << w << " " << wc << " "
                  << w - wc << " " << x_datum->node_data()[pNode.children[0]]
                  << " " << x_datum->node_data()[pNode.children[1]]
                  << std::endl;
      }
    }
  }
}

template <typename TREE> class fast_summation {
public:
  typedef TREE Tree;
  typedef typename TREE::node Node;

  enum Node_Type {
    LEAF,
    BRANCH,
  };

  fast_summation(const TREE &tree) : __tree(tree) {}

  void bind(const datum::ptr &x) { __data.push_back(x); }

  template <typename T>
  void bind(const std::vector<index_t> &indices, const std::vector<T> &x) {
    __data.push_back(datum_t<T>::create(indices, x));
  }

  template <typename Q>
  using ComputeFcn =
      std::function<Q(const index_t &, const index_t &, const vec3 &,
                      const std::vector<datum::ptr> &, Node_Type, const Node &,
                      const Tree &)>;

  template <typename Q>
  std::vector<Q>
  calc(const std::vector<vec3> &pov, ComputeFcn<Q> leafComputeFcn,
       ComputeFcn<Q> nodeComputeFcn, real eps = 0.5, bool debug = false) {
    for (int i = 0; i < __data.size(); i++) {
      __data[i]->pyramid(__tree);
    }

    vector<Q> u(pov.size(), z::zero<Q>());

    std::vector<ext::extents_t> extents =
        arp::build_extents(__tree, __tree.indices(), __tree.verts());
#if 1
    if (debug) {
      for (auto ext : extents) {
        vec4 c(0.5, 0.5, 0.1, 1.0);
        // std::cout << ext[0].transpose() << " " << ext[1].transpose() <<
        // std::endl;
        gg::geometry_logger::ext(ext[0], ext[1], c);
      }
    }
#endif

    int total_count = 0;
    int leaf_count = 0;
    int node_count = 0;
    std::vector<int> counts(__tree.nodes.size(), 0);
    std::vector<int> pov_counts(pov.size(), 0);
#pragma omp parallel for

    for (int i = 0; i < pov.size(); i++) {
      vec3 pi = pov[i];

      std::stack<int> stack1;
      stack1.push(0);
      vec4 c(0.5, 0.5, 0.1, 1.0);

      while (stack1.size() > 0) {
        total_count++;
        int j = stack1.top();
        stack1.pop();
        const Node &pNode = __tree.nodes[j];
        vec3 pj = pNode.center();

        real dc = va::dist(pi, pj);

        ext::extents_t ext = extents[j];
        vec3 de = ext[1] - ext[0];
        real V = de[0] * de[1] * de[2];
        real sc = 1.0 * pow(0.75 * V / M_PI, 1.0 / 3.0);
        // if (pNode.isLeaf()) {
        //   real sc = pNode.mag();
        // }
        // T sc = va::norm(pNode.half);
        //  real sc = 0.75 * pNode.mag();

        if (sc < dc * eps || pNode.isLeaf()) {
          pov_counts[i]++;
          if (pNode.isLeaf()) {
            c = vec4(0.8, 0.5, 1.0, 1.0);

            for (int jn = pNode.begin; jn < pNode.begin + pNode.size; jn++) {
              leaf_count++;
              int jj = __tree.permutation[jn];
              counts[j]++;
              u[i] += leafComputeFcn(i, jj, pi, __data, LEAF, pNode, __tree);
            }

          } else {
            node_count++;
            c = vec4(0.0, 0.5, 0.8, 1.0);
            u[i] += nodeComputeFcn(i, j, pi, __data, BRANCH, pNode, __tree);
            ;
          }

        }

        else {
          for (int j = 0; j < pNode.getNumChildren(); j++) {
            if (pNode.children[j] > -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }
#if 0
        if (i == 500 && true) {
          gg::geometry_logger::line(pi, pj, c);
          // gg::geometry_logger::ext(pj - vec3(sc, sc, sc), pj + vec3(sc, sc,
          // sc),
          //                          vec4(0.8, 0.0, 0.6, 0.5));
          gg::geometry_logger::ext(ext[0], ext[1], c);
          std::cout << " u[" << i << "]: " << u[i] << std::endl;
        }
#endif
      }
    }

#if 0
    for (int i = 0; i < counts.size(); i++) {
      if (counts[i] > 0)
        std::cout << i << " " << counts[i] << std::endl;
    
    }
#endif
#if 0
    double mean = std::accumulate(pov_counts.begin(), pov_counts.end(), 0.0) /
                  pov_counts.size();

    // Compute the variance
    double variance = std::accumulate(pov_counts.begin(), pov_counts.end(), 0.0,
                                      [mean](double acc, double x) {
                                        return acc + std::pow(x - mean, 2);
                                      }) /
                      pov_counts.size();

    // Compute the standard deviation
    double stddev = std::sqrt(variance);
    std::cout << "==== counts ==== " << std::endl;
    std::cout << " -  pov count: " << pov.size() << std::endl;
    std::cout << " -total count: " << total_count << std::endl;
    std::cout << " - node count: " << node_count << std::endl;
    std::cout << " - leaf count: " << leaf_count << std::endl;
    std::cout << " - compute count: " << node_count + leaf_count << std::endl;
    std::cout << " - tree leaf nodes:      " << __tree.leafNodes.size()
              << std::endl;
    std::cout << " - nlogn:           " << pov.size() * log(pov.size())
              << std::endl;
    std ::cout << " - pov / compute:  "
               << float(pov.size()) / float(node_count + leaf_count)
               << std::endl;
    std ::cout << " - compute / pov:  "
               << float(node_count + leaf_count) / float(pov.size())
               << std::endl;
    std::cout << " - node/total: " << float(node_count) / float(total_count)
              << std::endl;
    std::cout << " - leaf/total: " << float(leaf_count) / float(total_count)
              << std::endl;
    std::cout << " - log(tree leaf nodes): " << log(__tree.leafNodes.size())
              << std::endl;
    std::cout << "mean/std: " << mean << " " << sqrt(variance) << std::endl;
    std::cout << "================ " << std::endl;
#endif

    // for(int i = 0; i < 5; i++){
    //  std::cout << u[i].transpose() << std::endl;
    //}

    return u;
  }
  std::vector<datum::ptr> __data;
  const TREE &__tree;
};
} // namespace calder
} // namespace gaudi
#endif