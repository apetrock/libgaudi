#ifndef __ARP_DATUMS__
#define __ARP_DATUMS__

#include <cassert>
#include <cstddef>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include <iostream>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>

namespace gaudi {
namespace calder {

struct datum {
public:
  typedef std::shared_ptr<datum> ptr;

  datum(){};
  virtual ~datum(){};

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) = 0;

  template <int N, typename MortonT>
  void pyramid(const arp::bvh_tree<N, MortonT> &tree) {
    do_pyramid(tree.indices_, tree.internal_nodes_, tree.leaf_nodes_);
  }
};

template <typename TYPE> struct datum_t : public datum {
public:
  typedef std::shared_ptr<datum_t<TYPE>> ptr;

  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<TYPE> &data) {
    return std::make_shared<datum_t<TYPE>>(ind, data);
  }

  datum_t(const std::vector<index_t> &ind, const std::vector<TYPE> &data)
      : __leaf_indices(ind), __leaf_data(data){};
  virtual ~datum_t(){};

  const std::vector<TYPE> &leaf_data() const { return __leaf_data; }
  const std::vector<TYPE> &sorted_leaf_data() const { return __sorted_leaf_data; }
  std::vector<TYPE> &node_data() { return __tree_data; }
  const std::vector<TYPE> &node_data() const { return __tree_data; }

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) override {
    size_t num_leaves = leaf_nodes.size();
    __sorted_leaf_data.resize(num_leaves);
    for (size_t i = 0; i < num_leaves; i++) {
      index_t orig = indices[i];
      __sorted_leaf_data[i] = __leaf_data[__leaf_indices[orig]];
    }
    __tree_data = arp::build_pyramid(
        __sorted_leaf_data, internal_nodes, leaf_nodes,
        [](const TYPE &a, const TYPE &b) -> TYPE { return a + b; },
        z::zero<TYPE>());
  }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<TYPE> &__leaf_data;
  std::vector<TYPE> __sorted_leaf_data;
  std::vector<TYPE> __tree_data;
};

using scalar_datum = datum_t<real>;
using vec3_datum = datum_t<vec3>;
using mat3_datum = datum_t<mat3>;

struct edge_frame_datum : public datum {
public:
  typedef std::shared_ptr<edge_frame_datum> ptr;

  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<vec3> &data) {
    return std::make_shared<edge_frame_datum>(ind, data);
  }
  edge_frame_datum(const std::vector<index_t> &ind,
                   const std::vector<vec3> &data)
      : __leaf_indices(ind), __leaf_data(data){};
  virtual ~edge_frame_datum(){};

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) override {
    size_t num_leaves = leaf_nodes.size();
    __sorted_leaf_data.resize(num_leaves);
    for (size_t i = 0; i < num_leaves; i++) {
      index_t orig = indices[i];
      __sorted_leaf_data[i] = __leaf_data[__leaf_indices[orig]];
    }
    __tree_data = arp::build_pyramid<vec3, mat3>(
        __sorted_leaf_data, internal_nodes, leaf_nodes,
        [](const vec3 &e, const mat3 &F) -> mat3 {
          return F + e * e.transpose();
        },
        [](const mat3 &a, const mat3 &b) -> mat3 { return a + b; },
        mat3::Zero());
  }

  const std::vector<vec3> &leaf_data() const { return __leaf_data; }
  const std::vector<vec3> &sorted_leaf_data() const { return __sorted_leaf_data; }
  std::vector<mat3> &node_data() { return __tree_data; }
  const std::vector<mat3> &node_data() const { return __tree_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<vec3> &__leaf_data;
  std::vector<vec3> __sorted_leaf_data;
  std::vector<mat3> __tree_data;
};

struct quat_datum : public datum {
public:
  typedef std::shared_ptr<quat_datum> ptr;

  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<quat> &data) {
    return std::make_shared<quat_datum>(ind, data);
  }
  quat_datum(const std::vector<index_t> &ind, const std::vector<quat> &data)
      : __leaf_indices(ind), __leaf_data(data){};
  virtual ~quat_datum(){};

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) override {
    size_t num_leaves = leaf_nodes.size();
    __sorted_leaf_data.resize(num_leaves);
    for (size_t i = 0; i < num_leaves; i++) {
      index_t orig = indices[i];
      __sorted_leaf_data[i] = __leaf_data[__leaf_indices[orig]];
    }
    __node_data = arp::build_pyramid<quat, mat4>(
        __sorted_leaf_data, internal_nodes, leaf_nodes,
        [](const quat &q, const mat4 &Q) -> mat4 {
          return Q + q.coeffs() * q.coeffs().transpose();
        },
        [](const mat4 &a, const mat4 &b) -> mat4 { return a + b; },
        mat4::Zero());
  }

  const std::vector<quat> &leaf_data() const { return __leaf_data; }
  const std::vector<quat> &sorted_leaf_data() const { return __sorted_leaf_data; }
  std::vector<mat4> &node_data() { return __node_data; }
  const std::vector<mat4> &node_data() const { return __node_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<quat> &__leaf_data;
  std::vector<quat> __sorted_leaf_data;
  std::vector<mat4> __node_data;
};

struct extents_datum : public datum {
public:
  typedef std::shared_ptr<extents_datum> ptr;

  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<vec3> &data) {
    return std::make_shared<extents_datum>(ind, data);
  }
  extents_datum(const std::vector<index_t> &ind,
                const std::vector<vec3> &data)
      : __leaf_indices(ind), __leaf_data(data){};
  virtual ~extents_datum(){};

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) override {
    size_t num_leaves = leaf_nodes.size();
    __sorted_leaf_data.resize(num_leaves);
    for (size_t i = 0; i < num_leaves; i++) {
      index_t orig = indices[i];
      __sorted_leaf_data[i] = __leaf_data[__leaf_indices[orig]];
    }
    __tree_data = arp::build_pyramid<vec3, ext::extents_t>(
        __sorted_leaf_data, internal_nodes, leaf_nodes,
        [](const vec3 &pt, const ext::extents_t &e) -> ext::extents_t {
          return ext::expand(e, pt);
        },
        [](const ext::extents_t &a, const ext::extents_t &b) -> ext::extents_t {
          return ext::expand(b, a);
        },
        ext::init());
  }

  const std::vector<vec3> &leaf_data() const { return __leaf_data; }
  const std::vector<vec3> &sorted_leaf_data() const { return __sorted_leaf_data; }
  std::vector<ext::extents_t> &node_data() { return __tree_data; }
  const std::vector<ext::extents_t> &node_data() const { return __tree_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<vec3> &__leaf_data;
  std::vector<vec3> __sorted_leaf_data;
  std::vector<ext::extents_t> __tree_data;
};

struct com_datum : public datum {
public:
  typedef std::shared_ptr<com_datum> ptr;

  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<real> &areas,
                    const std::vector<vec3> &centroids) {
    return std::make_shared<com_datum>(ind, areas, centroids);
  }

  com_datum(const std::vector<index_t> &ind,
            const std::vector<real> &areas,
            const std::vector<vec3> &centroids)
      : __leaf_indices(ind), __areas(areas), __centroids(centroids){};
  virtual ~com_datum(){};

  virtual void do_pyramid(const std::vector<index_t> &indices,
                          const std::vector<arp::radix_tree_node> &internal_nodes,
                          const std::vector<arp::radix_tree_node> &leaf_nodes) override {
    size_t num_leaves = leaf_nodes.size();
    __sorted_leaf_data.resize(num_leaves);
    for (size_t i = 0; i < num_leaves; i++) {
      index_t orig = indices[i];
      index_t face = __leaf_indices[orig];
      real a = __areas[face];
      const vec3 &c = __centroids[face];
      __sorted_leaf_data[i] = vec4(a * c[0], a * c[1], a * c[2], a);
    }
    __tree_data = arp::build_pyramid(
        __sorted_leaf_data, internal_nodes, leaf_nodes,
        [](const vec4 &a, const vec4 &b) -> vec4 { return a + b; },
        z::zero<vec4>());
  }

  vec3 get_node_com(size_t j) const {
    return __tree_data[j].head<3>() / __tree_data[j][3];
  }

  vec3 get_leaf_com(size_t j) const {
    return __sorted_leaf_data[j].head<3>() / __sorted_leaf_data[j][3];
  }

  const std::vector<vec4> &sorted_leaf_data() const { return __sorted_leaf_data; }
  const std::vector<vec4> &node_data() const { return __tree_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<real> &__areas;
  const std::vector<vec3> &__centroids;
  std::vector<vec4> __sorted_leaf_data;
  std::vector<vec4> __tree_data;
};

} // namespace calder
} // namespace gaudi

#endif // __ARP_DATUMS__
