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

#include "bins.hpp"
#include "geometry_types.hpp"

/*
stages
*/

namespace m2 {

template <typename SPACE, typename CHARGE, typename OUTPUT>
class Simple_BarnesHutt {
  M2_TYPEDEFS;

public:
  using ComputeAccumulation = std::function<void()>;
  using ComputeCharge = std::function<void()>;

  using PreComputeFcn = std::function<void(
      const vector<CHARGE> &charges, const vector<T> &weights,
      const vector<coordinate_type> &points, int begin, int N,
      const vector<int> &permutation, coordinate_type &avgPoint,
      CHARGE &netCharge)>;

  using ComputeFcn = std::function<OUTPUT(const CHARGE &charge,
                                       const coordinate_type &chargePoint,
                                       const coordinate_type &evalPoint)>;

  vector<OUTPUT> 
  integrate(vector<CHARGE> &charges, vector<T> &weights,
                   vector<coordinate_type> &chargePoints,
                   vector<coordinate_type> &evalPoints,
                   PreComputeFcn preComputeFcn, ComputeFcn computeFcn) {
    vector<coordinate_type> u(evalPoints.size(), coordinate_type(0, 0, 0, 1.0));
    
    typedef pole_tree<SPACE> Tree;
    typedef pole_node<SPACE> Node;

    Tree octree(chargePoints);

    vector<coordinate_type> nodeCharges;
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
      coordinate_type avgPoint(0, 0, 0, 0.0);
      coordinate_type netCharge(0, 0, 0, 0.0);
      T netChargeMag = 0;
      int N = pNode.size;
      int beg = pNode.begin;

      preComputeFcn(charges, weights, chargePoints, beg, N,
                    octree.permutation, netCharge, avgPoint);


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

      std::stack<int> stack1;
      stack1.push(0);
      while (stack1.size() > 0) {

        int pId = stack1.top();
        stack1.pop();
        Node &pNode = octree.nodes[pId];
        coordinate_type pj = nodePositions[pId];

        coordinate_type dp = pi - pj;

        T dc = norm(dp);
        T sc = norm(pNode.half);

        // int ii = octree.permutation[pNode.begin];

        if (sc / dc > thresh) {
          // if(sc/dc > thresh){
          for (int j = 0; j < 8; j++) {
            if (pNode.children[j] != -1) {
              stack1.push(pNode.children[j]);
            }
          }
        }

        else {
          coordinate_type ci = nodeCharges[pId];
          OUTPUT ui = computeFcn(ci, pj, pi);
          u[i] += ui;
        }
      }
    }

    //for(int i = 0; i < 5; i++){
    //  std::cout << u[i].transpose() << std::endl;
    //}


    return u;
  } // integrator
#endif
}; // simple_barnes hutt
} // m2
