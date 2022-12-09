#include "m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "manifold/vec_addendum.h"

#include "GaudiGraphics/geometry_logger.h"

#include "dynamic_surface.hpp"
#include "m2_refactor.hpp"
#include "primitive_objects.hpp"
#include "primitive_operations.hpp"
#include <vector>
#include <zlib.h>

#ifndef __AAABBB__
#define __AAABBB__
namespace arp {

template <typename SPACE> struct aabb_node {
public:
  M2_TYPEDEFS;
  int dim;
  int id;
  int begin;
  int size;
  int level;
  int parent;
  int children[2];
  box_type bbox;
  // coordinate_type centerOfMass;
  // int neighbors[6]; to be implemented later

  aabb_node() {
    dim = 0;
    id = -1;
    begin = -1;
    size = -1;
    level = -1;
    parent = -1;
    children[0] = -1;
    children[1] = -1;
  }

  ~aabb_node() {}

  aabb_node(const aabb_node &rhs) {
    bbox = rhs.bbox;
    // centerOfMass     = rhs.centerOfMass;
    dim = rhs.dim;
    id = rhs.id;
    begin = rhs.begin;
    size = rhs.size;
    size = rhs.size;
    parent = rhs.parent;
    level = rhs.level;

    children[0] = rhs.children[0];
    children[1] = rhs.children[1];
  }

  aabb_node &operator=(const aabb_node &rhs) {
    // this = new aabb_node();
    if (this != &rhs) {

      bbox = rhs.bbox;
      // centerOfMass     = rhs.centerOfMass;
      dim = rhs.dim;
      id = rhs.id;
      begin = rhs.begin;
      size = rhs.size;
      parent = rhs.parent;
      level = rhs.level;

      children[0] = rhs.children[0];
      children[1] = rhs.children[1];
    }
    return *this;
  }

  int getNumChildren() { return 2; }

  bool isLeaf() const { return children[0] < 0 && children[1] < 0; }
};

template <typename SPACE, typename PRIMITIVE> struct aabb_tree {

  M2_TYPEDEFS;

public:
  typedef aabb_node<SPACE> node_type;

  aabb_tree() {}
  ~aabb_tree() {}

  aabb_tree(const aabb_tree &other) { *this = other; }

  aabb_tree(vector<PRIMITIVE> &points) { this->build(points, 24); }

  aabb_tree &operator=(const aabb_tree &rhs) {
    if (this != &rhs) {
      nodes = rhs.nodes;
      leafNodes = rhs.leafNodes;
      permutation = rhs.permutation;
    }
    return *this;
  }

  void calcHalfCenter(coordinate_type &half, coordinate_type &cen,
                      vector<PRIMITIVE> &primitives,
                      const vector<int> &permutation, int beg, int N) {

    if (permutation.empty())
      return;

    box_type bb = primitives[permutation[beg]].bbox();
    for (int i = beg; i < beg + N; i++) {
      PRIMITIVE p = primitives[permutation[i]];
      bb.expandBy(p.bbox());
    }
    half = bb.half;
    cen = bb.center;
  }

  void build(vector<PRIMITIVE> &primitives, int maxLevel) {
    // TIMER function//TIMER(__FUNCTION__);

    // inititalize permutation
    permutation.resize(primitives.size());
    leafIds.resize(primitives.size());
    // nodes.reserve(points.size()*points.size());
    // permutation.reserve(points.size());

    for (int i = 0; i < permutation.size(); i++)
      permutation[i] = i;

    node_type root;
    root.begin = 0;
    root.level = 0;
    root.size = primitives.size();
    root.id = nodes.size();
    root.parent = -1;
    if (primitives.empty())
      return;
    calcHalfCenter(root.bbox.half, root.bbox.center, primitives, permutation,
                   root.begin, root.size);

    stack<int> stack;
    nodes.reserve(log(primitives.size()) * primitives.size());
    stack.push(nodes.size());
    nodes.push_back(root);
    while (stack.size() > 0) {
      int pNodeId = stack.top();
      stack.pop();
      node_type pNode = nodes[pNodeId];

      int beg = pNode.begin;
      int N = pNode.size;
      int dim = pNode.dim;

      coordinate_type center = pNode.bbox.center;

      int cN[2] = {0, 0}, cCounter[2] = {0, 0}, cAccum[2] = {0, 0};
      std::vector<int> lPerm(N);

      for (int i = beg; i < beg + N; i++) {
        PRIMITIVE p = primitives[permutation[i]];
        int bin = (p.center()[dim] < center[dim]) ? 0 : 1;
        cN[bin]++;
      }

      for (int j = 1; j < 2; j++)
        cAccum[j] = cAccum[j - 1] + cN[j - 1];

      for (int i = beg; i < beg + N; i++) {
        PRIMITIVE p = primitives[permutation[i]];
        int bin = (p.center()[dim] < center[dim]) ? 0 : 1;
        lPerm[cAccum[bin] + cCounter[bin]] = permutation[i];
        cCounter[bin]++;
      }

      int ii = 0;

      for (int i = beg; i < beg + N; i++) {
        // update the global permutation with the local permutation
        permutation[i] = lPerm[ii];
        ii++;
      }

      for (int j = 0; j < 2; j++) {
        if (cN[j] == 0)
          continue;

        int cNodeId = nodes.size();
        node_type cNode;

        cNode.level = nodes[pNodeId].level + 1;
        cNode.begin = beg + cAccum[j];
        cNode.size = cN[j];
        cNode.id = cNodeId;
        cNode.parent = pNodeId;
        cNode.dim = (dim + 1) % 3;

        nodes[pNodeId].children[j] = cNodeId;

        calcHalfCenter(cNode.bbox.half, cNode.bbox.center, primitives,
                       permutation, cNode.begin, cNode.size);

        nodes.push_back(cNode);

        if (cNode.size < pNode.size && cNode.level < maxLevel)
          stack.push(cNodeId);
        else if (cNode.size == pNode.size || cNode.size == 1 ||
                 cNode.level == maxLevel)
          leafNodes.push_back(cNodeId);

        // leafIds.push_back(cNodeId);
      }
      // delete lPerm;
    }
  }

  vector<node_type> nodes;
  vector<int> leafIds;
  vector<int> leafNodes;
  vector<int> permutation;
};

template <typename SPACE, typename PRIMITIVE_A, typename PRIMITIVE_B>
PRIMITIVE_A
getNearest(PRIMITIVE_B &primB, const aabb_tree<SPACE, PRIMITIVE_A> &faceTree,
           const vector<PRIMITIVE_A> &primitives,
           std::function<typename SPACE::double_type(const PRIMITIVE_A &a,
                                                     const PRIMITIVE_B &b)>
               testAB,
           typename SPACE::double_type tol) {
  M2_TYPEDEFS;
  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SPACE, PRIMITIVE_A> tree_type;
  typedef typename tree_type::node_type Node;

  PRIMITIVE_A primMin;
  T dmin = std::numeric_limits<T>::infinity();

  const Node &root = faceTree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  box_type boxB = primB.bbox();
  boxB.inflate(coordinate_type(tol, tol, tol));

  while (cstack.size() > 0) {
    line_tests<SPACE> test;
    int cId = cstack.top();
    cstack.pop();
    const Node &cnode = faceTree.nodes[cId];

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {

      for (int k = cnode.begin; k < cnode.begin + cnode.size; k++) {

        const PRIMITIVE_A &primA = primitives[faceTree.permutation[k]];
        const box_type &boxA = primA.bbox();

        if (!boxA.overlap(boxB)) {
          continue;
        }

        T dist = testAB(primA, primB);

        if (dist < dmin && dist < std::numeric_limits<T>::infinity()) {
          dmin = dist;
          primMin = primA;
        }
      }
    }

    for (int i = 0; i < 2; i++) {
      if (cnode.children[i] > -1) {
        box_type boxA = faceTree.nodes[cnode.children[i]].bbox;
        if (boxA.overlap(boxB)) {
          cstack.push(cnode.children[i]);
        } else {
          // std::cout << cId << ": nover " << cnode.children[i] << std::endl;
          continue;
        }
      }
    }
  }
  return primMin;
};
} // namespace arp
#endif