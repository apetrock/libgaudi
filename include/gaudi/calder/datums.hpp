#include <cassert>
#include <cstddef>
#include <cxxabi.h>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>
#include <zlib.h>

#include "gaudi/arp/aabb.hpp"

#ifndef __CALDER_DATUM__
#define __CALDER_DATUM__

namespace gaudi {
namespace calder {

struct datum {
public:
  typedef std::shared_ptr<datum> ptr;

  // static ptr create() { return std::make_shared<datum>(); }

  datum(){};
  virtual ~datum(){};

  void pyramid(const arp::T1 &tree) { do_pyramid(tree); }
  void pyramid(const arp::T2 &tree) { do_pyramid(tree); }
  void pyramid(const arp::T3 &tree) { do_pyramid(tree); }
  virtual void do_pyramid(const arp::T1 &tree) = 0;
  virtual void do_pyramid(const arp::T2 &tree) = 0;
  virtual void do_pyramid(const arp::T3 &tree) = 0;
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
  std::vector<TYPE> &node_data() { return __tree_data; }
  const std::vector<TYPE> &node_data() const { return __tree_data; }

  // unsigned long operator[](int i) const { return __data[i]; }
  // unsigned long &operator[](int i) { return __data[i]; }

  template <typename TREE, int TREE_S> void __do_pyramid(const TREE &tree) {
    __tree_data =
        arp::build_pyramid<TREE_S, 1, TYPE>(tree, __leaf_indices, __leaf_data);
  }

  virtual void do_pyramid(const arp::T1 &tree) {
    __do_pyramid<arp::T1, 1>(tree);
  };
  virtual void do_pyramid(const arp::T2 &tree) {
    __do_pyramid<arp::T2, 2>(tree);
  };
  virtual void do_pyramid(const arp::T3 &tree) {
    __do_pyramid<arp::T3, 3>(tree);
  };

  const std::vector<index_t> &__leaf_indices;
  const std::vector<TYPE> &__leaf_data;
  std::vector<TYPE> __tree_data;
};

using scalar_datum = datum_t<real>;
using vec3_datum = datum_t<vec3>;
using mat3_datum = datum_t<mat3>;

std::vector<mat3> build_edge_frame_pyramid(const arp::T2 &tree,
                                           const std::vector<index_t> &indices,
                                           const std::vector<vec3> &x) {
  std::vector<mat3> pyramid = arp::__build_pyramid<2, 1, vec3, mat3>(
      tree, indices, x,
      mat3::Zero(), //
      [](const vec3 &e, const mat3 &F) { return F + e * e.transpose(); },
      [](const mat3 &Fc, const mat3 &Fp) { return Fp + Fc; });

  return pyramid;
}

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

  void __do_pyramid(const arp::T2 &tree) {
    __tree_data = build_edge_frame_pyramid(tree, __leaf_indices, __leaf_data);
  }

  virtual void do_pyramid(const arp::T2 &tree) { __do_pyramid(tree); };
  virtual void do_pyramid(const arp::T1 &tree){
      // do_nothing
  };
  virtual void do_pyramid(const arp::T3 &tree){
      // do_nothing
  };

  const std::vector<vec3> &leaf_data() const { return __leaf_data; }
  std::vector<mat3> &node_data() { return __tree_data; }
  const std::vector<mat3> &node_data() const { return __tree_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<vec3> &__leaf_data;
  std::vector<mat3> __tree_data;
};

std::vector<mat4> build_quat_pyramid(const arp::T2 &tree,
                                     const std::vector<index_t> &indices,
                                     const std::vector<quat> &x) {
  std::vector<mat4> pyramid = arp::__build_pyramid<2, 1, quat, mat4>(
      tree, indices, x,
      mat4::Zero(), //
      [](const quat &q, const mat4 &Q) {
        return Q + q.coeffs() * q.coeffs().transpose();
      },
      [](const mat4 &Fc, const mat4 &Fp) { return Fp + Fc; });

  return pyramid;
}

struct quat_datum : public datum {
public:
  typedef std::shared_ptr<quat_datum> ptr;
  // avg a quat:http://www.acsu.buffalo.edu/~johnc/ave_quat07.pdf
  static ptr create(const std::vector<index_t> &ind,
                    const std::vector<quat> &data) {
    return std::make_shared<quat_datum>(ind, data);
  }
  quat_datum(const std::vector<index_t> &ind, const std::vector<quat> &data)
      : __leaf_indices(ind), __leaf_data(data){};
  virtual ~quat_datum(){};

  void __do_pyramid(const arp::T2 &tree) {
    std::vector<mat4> tree_data =
        build_quat_pyramid(tree, __leaf_indices, __leaf_data);
    // do stuff to convert mat4 to quat
  }

  virtual void do_pyramid(const arp::T2 &tree) { __do_pyramid(tree); };
  virtual void do_pyramid(const arp::T1 &tree){
      // do_nothing
  };
  virtual void do_pyramid(const arp::T3 &tree){
      // do_nothing
  };

  const std::vector<quat> &leaf_data() const { return __leaf_data; }
  std::vector<quat> &node_data() { return __tree_data; }
  const std::vector<quat> &node_data() const { return __tree_data; }

  const std::vector<index_t> &__leaf_indices;
  const std::vector<quat> &__leaf_data;
  std::vector<quat> __tree_data;
};
} // namespace calder
} // namespace gaudi

#endif
