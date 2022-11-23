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
#include "manifold/asawa/geometry_types.hpp"
#include "manifold/bins.hpp"

//#include <execution>
#include <pstl/glue_execution_defs.h>

/*
stages
*/
using namespace bins;

namespace calder {

template <typename SPACE, typename CHARGE, typename OUTPUT>
class Simple_BarnesHutt {
  M2_TYPEDEFS;

public:
  using Tree = pole_tree<SPACE>;
  using Node = pole_node<SPACE>;
  using ComputeAccumulation = std::function<void()>;
  using ComputeCharge = std::function<void()>;

  using PreComputeFcn = std::function<void(
      const vector<CHARGE> &charges, const vector<T> &weights,
      const vector<coordinate_type> &points, Node &node, Tree &tree,
      const vector<int> &permutation, coordinate_type &avgPoint,
      CHARGE &netCharge)>;

  using ComputeFcn = std::function<OUTPUT(size_t i,                           //
                                          const CHARGE &ci, const CHARGE &cj, //
                                          const coordinate_type &pi,
                                          const coordinate_type &pj, //
                                          Node &node, Tree &tree)>;

  vector<OUTPUT> integrate(vector<CHARGE> &charges, vector<T> &weights,
                           vector<coordinate_type> &chargePoints,
                           vector<coordinate_type> &evalPoints,
                           PreComputeFcn preComputeFcn, ComputeFcn computeFcn) {

    vector<OUTPUT> u(evalPoints.size(), coordinate_type(0, 0, 0));

    Tree octree(chargePoints);

    vector<CHARGE> nodeCharges;
    vector<coordinate_type> nodePositions;
    nodeCharges.resize(octree.nodes.size());
    nodePositions.resize(octree.nodes.size());

    std::stack<int> stack;
    stack.push(0);
    T netWeight = 0;

    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = octree.nodes[pId];
      coordinate_type avgPoint(0, 0, 0);
      CHARGE netCharge = z::zero<CHARGE>();
      T netChargeMag = 0;
      int N = pNode.size;
      int beg = pNode.begin;

      preComputeFcn(charges, weights, chargePoints, pNode, octree, netCharge,
                    avgPoint);

      nodeCharges[pId] = netCharge;
      nodePositions[pId] = avgPoint;

      for (int j = 0; j < 8; j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }

    T thresh = 0.5;

    for (int i = 0; i < evalPoints.size(); i++) {
      int count = 0;

      coordinate_type pi = evalPoints[i];
      CHARGE ci = charges[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pj = nodePositions[pId];

        T dc = va::dist(pi, pj);
        // T sc = va::norm(pNode.half);
        T sc = pNode.half.maxCoeff();
        // int ii = octree.permutation[pNode.begin];

        if (sc / dc < thresh || pNode.isLeaf()) {
          if (pNode.isLeaf()) {
            for (int i = pNode.begin; i < pNode.begin + pNode.size; i++) {
              int ii = octree.permutation[i];
              CHARGE cj = charges[ii];
              OUTPUT ui = computeFcn(pId, ci, cj, pi, pj, pNode, octree);
              u[i] += ui;
            }

          } else {
            CHARGE cj = nodeCharges[pId];
            OUTPUT ui = computeFcn(pId, ci, cj, pi, pj, pNode, octree);
            u[i] += ui;
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
  } // integrator
};  // simple_barnes hutt

template <typename SPACE, typename CHARGE, typename PRIMITIVE, typename OUTPUT>
class Geometry_Integrator {
  M2_TYPEDEFS;

public:
  using Tree = aabb_tree<SPACE, PRIMITIVE>;
  using Node = aabb_node<SPACE>;

  using ComputeAccumulation = std::function<void()>;
  using ComputeCharge = std::function<void()>;

  using PreComputeFcn = std::function<void(
      const vector<PRIMITIVE> &points, Node &node, Tree &tree,
      CHARGE &netCharge, coordinate_type &avgPoint,
      coordinate_type &avgNormal)>;

  using ComputeFcn = std::function<OUTPUT(
      int i, const CHARGE &q, const coordinate_type &pc,
      const coordinate_type &pe, const coordinate_type &N,
      const vector<PRIMITIVE> &points, Node &node, Tree &tree)>;

  void integrate(vector<CHARGE> &charges, vector<PRIMITIVE> &chargePrimitives,
                 const vector<coordinate_type> &evalPoints, vector<OUTPUT> &u,
                 PreComputeFcn preComputeFcn, ComputeFcn computeFcn) {

    Tree tree(chargePrimitives);

    vector<CHARGE> nodeCharges;
    vector<coordinate_type> nodePositions;
    vector<coordinate_type> nodeNormals;
    nodeCharges.resize(tree.nodes.size());
    nodePositions.resize(tree.nodes.size());
    nodeNormals.resize(tree.nodes.size());

    std::stack<int> stack;
    stack.push(0);
    T netWeight = 0;

    while (stack.size() > 0) {
      int pId = stack.top();
      stack.pop();
      Node &pNode = tree.nodes[pId];
      coordinate_type avgPoint(0, 0, 0);
      coordinate_type avgNormal(0, 0, 0);

      CHARGE netCharge = z::zero<CHARGE>();
      T netChargeMag = 0;
      preComputeFcn(chargePrimitives, pNode, tree, netCharge, avgPoint,
                    avgNormal);

      nodeCharges[pId] = netCharge;
      nodePositions[pId] = avgPoint;
      nodeNormals[pId] = avgNormal;

      for (int j = 0; j < pNode.getNumChildren(); j++) {
        if (pNode.children[j] != -1)
          stack.push(pNode.children[j]);
      }
    }

    T thresh = 0.5;
    // std::vector<size_t> counter(evalPoints.size());
    // std::iota(std::begin(counter), std::end(counter), 0);
    //#pragma omp parallel for
    for (size_t i = 0; i < evalPoints.size(); i++) {
      // std::for_each(/*std::execution::seq, */counter.begin(), counter.end(),
      // [&](auto &&i) {
      int count = 0;
      coordinate_type pi = evalPoints[i];

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = tree.nodes[pId];
        coordinate_type pj = nodePositions[pId];
        coordinate_type Nj = nodeNormals[pId];

        coordinate_type dp = pi - pj;

        T dc = va::norm(dp);
        T sc = pNode.bbox.half.maxCoeff();
        // T sc = va::norm(pNode.bbox.half);

        if (sc / dc <= thresh || pNode.isLeaf()) {
          CHARGE cj = nodeCharges[pId];
          OUTPUT ui =
              computeFcn(i, cj, pj, pi, Nj, chargePrimitives, pNode, tree);
          u[i] += ui;
          // if(sc/dc > thresh){
        } else {
          for (int j = 0; j < pNode.getNumChildren(); j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }
      }
    } //);
  }   // integrator

}; // geometry_integrator

} // namespace calder
#endif