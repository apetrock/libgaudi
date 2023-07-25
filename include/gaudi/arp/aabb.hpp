#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "GaudiGraphics/geometry_logger.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stack>
#include <vector>
#include <zlib.h>

#ifndef __AAABBB__
#define __AAABBB__
namespace gaudi {
namespace arp {

using real = double;
using vec3 = Eigen::Matrix<real, 3, 1>;
using vec4 = Eigen::Matrix<real, 4, 1>;
using mat3 = Eigen::Matrix<real, 3, 3>;

using index_t = int;

template <int S>
ext::extents_t calc_extents(index_t i, const std::vector<index_t> &indices,
                            const std::vector<vec3> &vertices) {
  vec3 min = vertices[indices[S * i + 0]];
  vec3 max = min;
  for (int k = 0; k < S; k++) {
    vec3 p = vertices[indices[S * i + k]];
    min = va::min(p, min);
    max = va::max(p, max);
  }

  return {min, max};
}

template <int S> ext::extents_t calc_extents(std::array<vec3, S> &verts) {
  vec3 min = verts[0];
  vec3 max = min;
  for (int k = 0; k < S; k++) {
    vec3 p = verts[k];
    min = va::min(p, min);
    max = va::max(p, max);
  }
  return {min, max};
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

template <typename T, typename CTYPE> struct half_space {
public:
  CTYPE N;
  CTYPE cen;

  real d;
  real mag;
  half_space() : d(0), N(CTYPE(0, 0, 0)){};

  half_space(const CTYPE &cen_, const CTYPE &N_) { set(cen_, N_); };
  void set(const CTYPE &cen_, const CTYPE &N_) {
    mag = N_.norm();
    N = N_ / mag;
    d = N.dot(cen_);
    cen = cen_;
  }

  index_t left_right(const CTYPE &p) const { return int(N.dot(p) - d >= 0); }

  index_t intersect(const ext::extents_t ext) const {
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

  int getNumChildren() const { return 2; }

  bool isLeaf() const { return children[0] < 0 || children[1] < 0; }

  inline void calcSVD(mat3 &M, vec3 &s) {
    Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU);
    const mat3 U = svd.matrixU();
    s = svd.singularValues();
    M = U;
  }

  const vec3 &center() const { return half.cen; }
  const real &mag() const { return half.mag; }

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

#if 1
    vec3 mx_p;
    index_t mx_d = 0;
    vec3 var = vec3::Zero();
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        vec3 dp = p - c;
        var += dp.cwiseProduct(dp);
      }
    }
    var /= real(this->size);
    index_t mx = var[0] > var[1] ? 0 : 1;
    mx = var[mx] > var[2] ? mx : 2;

    vec3 N = vec3::Zero();
    N[mx] = 1.0 * sqrt(var[mx]);

    vec3 Nn = N.normalized();
    vec3 dc = vec3::Zero();
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        vec3 dp = p - c;
        dc += Nn.dot(dp) * Nn;
      }
    }
    dc /= real(S * this->size);
    c += dc;
    // gg::geometry_logger::line(c, c + N, vec4(1.0, 0.0, 0.0, 0.0));
    half.set(c, N);
#else
    mat3 U = mat3::Zero();
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        vec3 dp = (p - c).normalized();
        U += dp * dp.transpose();
      }
    }
    Eigen::JacobiSVD<mat3> svd(U, Eigen::ComputeFullU);
    U = svd.matrixU();
    vec3 s = svd.singularValues();
    vec3 N = U.col(0).transpose();
    N.normalize();
    real mx = 0.0;
    for (int i = this->begin; i < this->begin + this->size; i++) {
      for (int k = 0; k < S; k++) {

        vec3 p = vertices[indices[S * permutation[i] + k]];
        vec3 dp = (p - c).normalized();
        real ndp = abs(N.dot(dp));
        mx = std::max(mx, dp.norm());
      }
    }
    gg::geometry_logger::line(c, c + mx * N, vec4(1.0, 0.0, 0.0, 0.0));
    half.set(c, mx * N);
#endif
  }

  void debug(const std::vector<index_t> &indices,
             const std::vector<vec3> &vertices,
             const std::vector<index_t> &permutation) {

    if (permutation.empty())
      return;
    vec3 N = half.N;
    vec3 h_cen = half.cen;
    vec4 c(N[0], N[1], N[2], 1.0);

    for (int i = this->begin; i < this->begin + this->size; i++) {
      vec3 prim_cen(0.0, 0.0, 0.0);
      for (int k = 0; k < S; k++) {
        vec3 p = vertices[indices[S * permutation[i] + k]];
        prim_cen += p;
      }
      prim_cen /= real(S);
      // vec3 hc = half.d * half.N;

      gg::geometry_logger::line(h_cen, prim_cen, c);
    }
    gg::geometry_logger::line(h_cen, h_cen + half.mag * half.N, c);
  }

  void debug_half() {

    vec3 N = half.N;
    real t = level / 10;
    vec4 c(cos(t), cos(t + M_PI / 3.0), cos(t + M_PI / 6.0), 1.0);
    gg::geometry_logger::line(half.cen, half.cen + half.mag * half.N, c);
  }
};

// S = stride
template <int S> struct aabb_tree {

public:
  typedef std::shared_ptr<aabb_tree<S>> ptr;
  typedef aabb_node<S> node;

  static ptr create(const std::vector<index_t> &indices,
                    const std::vector<vec3> &vertices, int lvl = 8) {
    return std::make_shared<aabb_tree<S>>(indices, vertices, lvl);
  }

  aabb_tree() {}
  ~aabb_tree() {}

  aabb_tree(const aabb_tree &other) { *this = other; }

  aabb_tree(const std::vector<index_t> &indices,
            const std::vector<vec3> &vertices, int lvl = 8)
      : __x(vertices), __indices(indices) {
    this->build(__indices, __x, lvl);
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
    leafNodes.reserve(indices.size() / S);
    nodes.reserve(4 * log(indices.size()) * indices.size());
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

        if (cNode.size < 2 || cNode.level == maxLevel)
          leafNodes.push_back(cNodeId);
        else
          stack.push(cNodeId);

        // leafIds.push_back(cNodeId);
      }
    }
    // debug(indices, vertices);
  }

  const std::vector<index_t> &indices() const { return this->__indices; }
  const std::vector<vec3> &verts() const { return this->__x; }
  const vec3 &vert(index_t i) const { return __x[__indices[i]]; }

  void debug() {
    for (int i = 0; i < leafNodes.size(); i++) {
      node &n = nodes[leafNodes[i]];
      n.debug(__indices, __x, permutation);
    }
  }

  void debug_half() {
    for (int i = 0; i < nodes.size(); i++) {
      node &n = nodes[i];
      n.debug_half();
      // std::cout << n.level << " " << n.children[0] << " " << n.children[1]
      //          << " " << n.size << std::endl;
      if (n.isLeaf())
        continue;
      vec4 c(0.5, 0.5, 0.5, 1.0);

      if (n.children[0] > 0) {
        node &n0 = nodes[n.children[0]];
        gg::geometry_logger::line(n.half.cen, n0.half.cen, c);
      }
      if (n.children[1] > 0) {
        node &n1 = nodes[n.children[1]];
        gg::geometry_logger::line(n.half.cen, n1.half.cen, c);
      }
    }
  }

  vector<node> nodes;
  vector<index_t> leafNodes;
  vector<index_t> permutation;
  const std::vector<index_t> &__indices;
  const std::vector<vec3> &__x;
};

#if 1
template <int ST, int SS> // T=test, S=set... DOH! T could equal tree...
std::vector<index_t>
getNearest(index_t &idT, const std::vector<index_t> &t_inds,
           const vector<vec3> &t_verts, //
           const aabb_tree<SS> &s_tree, real tol,
           std::function<real(const index_t &idT, //
                              const std::vector<index_t> &t_inds,
                              const vector<vec3> &t_verts, //
                              const index_t &idS,          //
                              const std::vector<index_t> &s_inds,
                              const std::vector<vec3> &s_verts)>
               testAB) {

  // TIMER function//TIMER(__FUNCTION__);
  typedef aabb_tree<SS> tree_type;
  typedef typename tree_type::node Node;

  bool expanding_rad = tol > 999.9;

  index_t idMin = -1;
  real dmin = std::numeric_limits<real>::infinity();

  const Node &root = s_tree.nodes[0];
  std::stack<int> cstack;
  cstack.push(0);
  bool hit = false;
  // T tol = 0.05;
  ext::extents_t extT = calc_extents<ST>(idT, t_inds, t_verts);
  extT = ext::inflate(extT, tol);
  std::vector<index_t> collisions;
  while (cstack.size() > 0) {
    int cId = cstack.top();
    cstack.pop();
    const Node &cnode = s_tree.nodes[cId];

    if (expanding_rad) {
      real d = ext::distance(extT, cnode.half.cen);
      tol = std::min(tol, d);
      extT = calc_extents<ST>(idT, t_inds, t_verts);
      extT = ext::inflate(extT, tol);
    }

    if (cnode.children[0] == -1 && cnode.children[1] == -1) {

      for (int k = cnode.begin; k < cnode.begin + cnode.size; k++) {

        const index_t &idS = s_tree.permutation[k];

        ext::extents_t extS =
            calc_extents<SS>(idS, s_tree.indices(), s_tree.verts());

        if (!ext::overlap(extT, extS)) {

          continue;
        }

        if (idT == 1000 && 0) {
          vec3 cT = 0.5 * (extT[0] + extT[1]);
          vec3 cS = 0.5 * (extS[0] + extS[1]);
          gg::geometry_logger::ext(extS[0], extS[1], vec4(0.0, 1.0, 0.0, 1.0));
          gg::geometry_logger::ext(extT[0], extT[1], vec4(1.0, 0.0, 0.0, 1.0));
          gg::geometry_logger::line(cT, cS, vec4(1.0, 1.0, 0.0, 1.0));
        }

        real dist = testAB(idT, t_inds, t_verts, //
                           idS, s_tree.indices(), s_tree.verts());

        if (dist < dmin) {
          dmin = dist;
          idMin = idS;
        }
        if (dist < tol) {
          collisions.push_back(idS);
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

  if (dmin < tol) {
    collisions.push_back(idMin);
  } else
    collisions.push_back(-1);

  return collisions;
};
#endif

template <int TREE_S, int NODE_S>
void for_each(
    index_t i, const aabb_tree<TREE_S> &tree,
    const std::vector<index_t> &indices,
    std::function<void(index_t id, const aabb_tree<TREE_S> &tree)> func) {
  typedef typename aabb_tree<TREE_S>::node node;
  const node &cnode = tree.nodes[tree.leafNodes[i]];
  index_t beg = cnode.begin;
  index_t end = beg + cnode.size;

  for (int j = beg; j < end; j++) {
    const index_t &id = tree.permutation[j];
    for (int k = 0; k < NODE_S; k++) {

      func(indices[NODE_S * id + k], tree);
    }
  }
}

bool l_isnan(vec3 v) { return v.hasNaN(); }
bool l_isnan(mat3 v) { return v.hasNaN(); }
bool l_isnan(mat4 v) { return v.hasNaN(); }

bool l_isnan(std::array<vec3, 2> v) { return v[0].hasNaN() || v[1].hasNaN(); }

bool l_isnan(real v) { return std::isnan(v); }

template <int TREE_S, int NODE_S, typename Q0, typename Q1>
std::vector<Q1>
__build_pyramid(const aabb_tree<TREE_S> &tree,       //
                const std::vector<index_t> &indices, //
                const std::vector<Q0> &x, const Q1 &q_init,
                std::function<Q1(const Q0 &q0, const Q1 &q1)> lfunc,
                std::function<Q1(const Q1 &qc, const Q1 &qp)> nfunc) {
  typedef aabb_tree<TREE_S> tree_type;
  typedef typename tree_type::node node;
  std::vector<Q1> charges(tree.nodes.size(), q_init);
  // init the bounds

  for (int i = 0; i < tree.leafNodes.size(); i++) {
    Q1 q1 = q_init;
    for_each<TREE_S, NODE_S>(
        i, tree, indices, [&](index_t j, const aabb_tree<TREE_S> &tree) {
          Q0 q0 = x[j];
          q1 = lfunc(q0, q1);
          if (l_isnan(q1)) {
            std::cout << "q0: " << q0 << std::endl;
            std::cout << j << " " << x.size() << std::endl;
            throw std::runtime_error(std::string(__PRETTY_FUNCTION__) +
                                     std::string(": nan"));
          }
        });
    charges[tree.leafNodes[i]] = q1;
  }
  // return charges;
  std::deque<int> stack(tree.leafNodes.begin(), tree.leafNodes.end());
  std::vector<bool> in_queue(charges.size(), false);
  while (stack.size() > 0) {
    index_t cNodeId = stack.back();

    stack.pop_back();
    in_queue[cNodeId] = false;
    const node &cnode = tree.nodes[cNodeId];
    index_t pNodeId = cnode.parent;

    if (pNodeId < 0)
      continue;

    const node &pnode = tree.nodes[cnode.parent];
    const Q1 &extC = charges[cNodeId];
    const Q1 &extP = charges[pNodeId];

    charges[pNodeId] = nfunc(extC, extP);
    if (!in_queue[pNodeId]) {
      stack.push_front(pNodeId);
      in_queue[pNodeId] = true;
    }
  }

  return charges;
}

template <int S>
std::vector<ext::extents_t> build_extents(const aabb_tree<S> &tree,
                                          const std::vector<index_t> &indices,
                                          const std::vector<vec3> &x) {
  double mx = std::numeric_limits<double>::max();
  std::vector<ext::extents_t> extents =
      __build_pyramid<S, S, vec3, ext::extents_t>(
          tree, indices, x,
          {
              vec3(mx, mx, mx),
              vec3(-mx, -mx, -mx),
          }, //
          [](const vec3 &q0, const ext::extents_t &q1) {
            return ext::expand(q1, q0);
          },
          [](const ext::extents_t &qc, const ext::extents_t &qp) {
            return ext::expand(qp, qc);
          });

  return extents;
}

template <int TREE_S, int NODE_S, typename Q>
std::vector<Q> build_pyramid(const aabb_tree<TREE_S> &tree,
                             const std::vector<index_t> &indices,
                             const std::vector<Q> &x) {
  std::vector<Q> pyramid = __build_pyramid<TREE_S, NODE_S, Q, Q>(
      tree, indices, x,
      z::zero<Q>(), //
      [](const Q &p, const Q &q) { return p + q; },
      [](const Q &qc, const Q &qp) { return qp + qc; });

  return pyramid;
}

using T1 = arp::aabb_tree<1>;
using T2 = arp::aabb_tree<2>;
using T3 = arp::aabb_tree<3>;

} // namespace arp
} // namespace gaudi
#endif