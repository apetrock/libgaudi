#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "GaudiGraphics/geometry_logger.h"
#include "manifold/geometry_types.hpp"
#include "manifold/vec_addendum.h"
#include <array>
#include <cstdio>
#include <vector>
#include <zlib.h>

#ifndef __AAABBB__
#define __AAABBB__
namespace arp {
using real = double;
using vec3 = Eigen::Matrix<real, 3, 1>;
using vec4 = Eigen::Matrix<real, 4, 1>;
using mat3 = Eigen::Matrix<real, 3, 3>;

using index_t = int;
using extents_t = std::array<vec3, 2>;

template <int S>
extents_t calc_extents(index_t i, const std::vector<index_t> &indices,
                       const std::vector<vec3> &vertices) {
  vec3 min = vertices[indices[S * i + 0]];
  vec3 max = min;
  // std::cout << "E: " << 2 * i << " " << indices.size() << std::endl;
  // std::cout << "  " << indices[2 * i + 0] << " " << indices[2 * i + 1] << " "
  //           << vertices.size() << std::endl;

  for (int k = 0; k < S; k++) {
    vec3 p = vertices[indices[S * i + k]];
    min = va::min(p, min);
    max = va::max(p, max);
  }
  //  std::cout << "1:  " << min.transpose() << " - " << max.transpose()
  //            << std::endl;
  return {min, max};
}

bool overlap(const extents_t &A, const extents_t &B) {
  if (va::greater_than(A[0], B[1]))
    return false;
  if (va::less_than(A[1], B[0]))
    return false;

  return true;
}

template <int S>
vec3 calc_center(index_t i, const std::vector<index_t> &indices,
                 const std::vector<vec3> &vertices) {
  vec3 cen = vec3::Zero();

  for (int k = 0; k < S; k++) {
    vec3 p = vertices[indices[S * i + k]];
    cen += p;
  }
  return cen /= real(S);
}

extents_t inflate(extents_t e, real eps) {
  vec3 deps(eps, eps, eps);
  e[0] -= deps;
  e[1] += deps;
  return e;
}

template <typename T, typename CTYPE> struct half_space {
public:
  CTYPE N;
  CTYPE cen;

  real d;
  half_space() : d(0), N(CTYPE(0, 0, 0)){};

  half_space(const CTYPE &cen_, const CTYPE &N_) { set(cen_, N_); };
  void set(const CTYPE &cen_, const CTYPE &N_) {
    d = N_.dot(cen_);
    N = N_;
    cen = cen_;
  }

  index_t left_right(const CTYPE &p) const { return int(N.dot(p) - d >= 0); }

  index_t intersect(const extents_t ext) const {

    real pmin = N.dot(ext[0]) - d;
    real pmax = N.dot(ext[1]) - d;

    pmin = va::sgn(pmin);
    pmax = va::sgn(pmax);

    //-2 0 2 => left itx right
    return index_t(pmin + pmax);
  }
};

using half_s = half_space<real, vec3>;

template <int S> struct aabb_node {

public:
  int id;
  int begin;
  int size;
  int level;
  int parent;
  int children[2];

  half_s half;
  // coordinate_type centerOfMass;
  // int neighbors[6]; to be implemented later

  aabb_node() {
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

    half = rhs.half;
    // centerOfMass     = rhs.centerOfMass;
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

      half = rhs.half;
      // centerOfMass     = rhs.centerOfMass;
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

  inline void calcSVD(mat3 &M, vec3 &s) {
    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU);
    const mat3 U = svd.matrixU();
    s = svd.singularValues();
    M = U;
  }

  void calcHalfCenter(const std::vector<index_t> &indices,
                      const std::vector<vec3> &vertices,
                      const std::vector<index_t> &permutation) {

    if (permutation.empty())
      return;

    vec3 c = vec3::Zero();
    // avg is weighted, whereas min/max is absolute
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        c += p;
      }
    }

    c /= real(S * this->size);

    mat3 M = mat3::Zero();
    real al = 0.0;
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        vec3 dp = p - c;
        al += dp.squaredNorm();
        M += dp * dp.transpose();
      }
    }

    vec3 s;
    calcSVD(M, s);

    vec3 N = M.block(0, 0, 3, 1);
    half.set(c, N);
    // half.N = vec3(0.0, 0.0, 1.0);
    //  gg::geometry_logger::line(half.center, half.center + 0.5 * half.N,
    //                           vec4(1.0, 1.0 / real(level + 1), 0.0, 1.0));
  }

  void debug(const std::vector<index_t> &indices,
             const std::vector<vec3> &vertices,
             const std::vector<index_t> &permutation) {

    if (permutation.empty())
      return;
    vec3 N = half.N;
    N.normalize();
    vec4 c(N[0], N[1], N[2], 1.0);

    for (int i = this->begin; i < this->begin + this->size; i++) {
      vec3 prim_cen(0.0, 0.0, 0.0);
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        prim_cen += p;
      }
      prim_cen /= real(S);
      // vec3 hc = half.d * half.N;
      vec3 hc = half.cen;
      gg::geometry_logger::line(hc, prim_cen, c);
    }
  }
};

// S = stride
template <int S> struct aabb_tree {

public:
  typedef std::shared_ptr<aabb_tree<S>> ptr;

  static ptr create(std::vector<index_t> &indices, std::vector<vec3> &vertices,
                    int lvl = 8) {
    return std::make_shared<aabb_tree<S>>(indices, vertices);
  }

  using node = aabb_node<S>;
  aabb_tree() {}
  ~aabb_tree() {}

  aabb_tree(const aabb_tree &other) { *this = other; }

  aabb_tree(std::vector<index_t> &indices, std::vector<vec3> &vertices,
            int lvl = 8) {
    this->build(indices, vertices, lvl);
  }

  aabb_tree &operator=(const aabb_tree &rhs) {
    if (this != &rhs) {
      nodes = rhs.nodes;
      leafNodes = rhs.leafNodes;
      permutation = rhs.permutation;
    }
    return *this;
  }

  void build(const std::vector<index_t> &indices,
             const std::vector<vec3> &vertices, int maxLevel) {
    // TIMER function//TIMER(__FUNCTION__);

    // inititalize permutation
    permutation.resize(indices.size() / S);
    leafIds.resize(indices.size() / S);
    // nodes.reserve(points.size()*points.size());
    // permutation.reserve(points.size());

    for (int i = 0; i < permutation.size(); i++)
      permutation[i] = i;

    node root;
    root.begin = 0;
    root.level = 0;
    root.size = permutation.size();
    root.id = nodes.size();
    root.parent = -1;
    if (indices.empty())
      return;

    root.calcHalfCenter(indices, vertices, permutation);

    stack<int> stack;
    nodes.reserve(log(indices.size()) * indices.size());
    stack.push(nodes.size());
    nodes.push_back(root);

    while (stack.size() > 0) {
      int pNodeId = stack.top();
      stack.pop();

      const node &pNode = nodes[pNodeId];

      int beg = pNode.begin;
      int N = pNode.size;

      int cN[2] = {0, 0}, cCounter[2] = {0, 0}, cAccum[2] = {0, 0};
      std::vector<int> lPerm(N);

      for (int i = beg; i < beg + N; i++) {
        vec3 cen = calc_center<S>(permutation[i], indices, vertices);
        int bin = pNode.half.left_right(cen);
        cN[bin]++;
      }

      for (int j = 1; j < 2; j++)
        cAccum[j] = cAccum[j - 1] + cN[j - 1];

      for (int i = beg; i < beg + N; i++) {
        vec3 cen = calc_center<S>(permutation[i], indices, vertices);
        int bin = pNode.half.left_right(cen);
        lPerm[cAccum[bin] + cCounter[bin]] = permutation[i];
        cCounter[bin]++;
      }

      int ii = 0;

      for (int i = 0; i < N; i++) {
        // update the global permutation with the local permutation
        permutation[beg + i] = lPerm[i];
      }

      for (int j = 0; j < 2; j++) {
        if (cN[j] == 0)
          continue;

        int cNodeId = nodes.size();
        node cNode;
        cNode.level = pNode.level + 1;
        cNode.begin = beg + cAccum[j];
        cNode.size = cN[j];
        cNode.id = cNodeId;
        cNode.parent = pNodeId;
        nodes[pNodeId].children[j] = cNodeId;

        cNode.calcHalfCenter(indices, vertices, permutation);
        nodes.push_back(cNode);

        if (cNode.size < pNode.size && cNode.level < maxLevel)
          stack.push(cNodeId);
        else if (cNode.size == pNode.size || cNode.size == 1 ||
                 cNode.level == maxLevel)
          leafNodes.push_back(cNodeId);

        // leafIds.push_back(cNodeId);
      }
    }
  }

  void debug(const std::vector<index_t> &indices,
             const std::vector<vec3> &vertices) {

    for (auto n : nodes)
      if (n.isLeaf()) {
        node &pn = nodes[n.parent];
        n.debug(indices, vertices, permutation);
      }
  }

  vector<node> nodes;
  vector<index_t> leafIds;
  vector<index_t> leafNodes;
  vector<index_t> permutation;
};

#if 1
template <int ST, int SS> // T=test, S=set... DOH! T could equal tree...
index_t getNearest(index_t &idT, const std::vector<index_t> &t_inds,
                   const vector<vec3> &t_verts, //
                   const aabb_tree<SS> &s_tree,
                   const std::vector<index_t> &s_inds,
                   const vector<vec3> &s_verts, //
                   real tol,
                   std::function<real(const index_t &idT, //
                                      const std::vector<index_t> &t_inds,
                                      const vector<vec3> &t_verts, //
                                      const index_t &idS,          //
                                      const std::vector<index_t> &s_inds,
                                      const vector<vec3> &s_verts)>
                       testAB) {

  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SS> tree_type;
  typedef typename tree_type::node Node;

  index_t idMin = -1;
  real dmin = std::numeric_limits<real>::infinity();

  const Node &root = s_tree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  extents_t extT = calc_extents<ST>(idT, t_inds, t_verts);
  inflate(extT, tol);

  while (cstack.size() > 0) {
    int cId = cstack.top();
    cstack.pop();
    const Node &cnode = s_tree.nodes[cId];

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {

      for (int k = cnode.begin; k < cnode.begin + cnode.size; k++) {

        const index_t &idS = s_tree.permutation[k];

        extents_t extS = calc_extents<SS>(idS, s_inds, s_verts);
        if (!overlap(extT, extS)) {
          continue;
        }

        real dist = testAB(idT, t_inds, t_verts, //
                           idS, s_inds, s_verts);

        if (dist < dmin && dist < std::numeric_limits<real>::infinity()) {
          dmin = dist;
          idMin = idS;
        }
      }
    }

    index_t itx = cnode.half.intersect(extT);

    if (itx <= 0 && cnode.children[0] > 0) {
      cstack.push(cnode.children[0]);
    }

    if (itx >= 0 && cnode.children[1] > 0) {
      cstack.push(cnode.children[1]);
    }
  }
  return idMin;
};
#endif

} // namespace arp
#endif