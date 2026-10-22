#include <cassert>
#include <cstddef>
#include <type_traits>

#include <iostream>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "shell/datum_x.hpp"
#include "shell/shell.hpp"

#ifndef __ASAWA_DATUM__
#define __ASAWA_DATUM__
namespace gaudi {

namespace asawa {

namespace rod {
class rod;
}

enum prim_type {
  VERTEX,
  FACE,
  EDGE,
  CORNER,
};

struct datum {
public:
  typedef std::shared_ptr<datum> ptr;

  // static ptr create() { return std::make_shared<datum>(); }

  datum(prim_type type) : __type(type){};
  virtual ~datum(){};

  void alloc(size_t sz) { do_alloc(sz); }
  void resize(size_t sz) { do_resize(sz); }
  void map(index_t i, index_t it) { do_map(i, it); }
  void permute(const std::vector<index_t> &permute) { do_permute(permute); };

  virtual void do_alloc(const size_t &sz) = 0;
  virtual void do_resize(const size_t &sz) = 0;
  virtual size_t size() = 0;

  /*migrate these calcs to be less generic and specialized for each op*/

  virtual void subdivide(const shell::shell &M, //
                         const index_t &i0, const index_t &i1,
                         const index_t &i2, const index_t &i3, //
                         const real &s, const index_t &source_corner){};

  virtual void collapse(const shell::shell &M, const index_t &i){};

  virtual void merge(const shell::shell &M,                //
                     const index_t &i0, const index_t &i1, //
                     const index_t &i2, const index_t &i3, //
                     const index_t vA0, const index_t vA1, //
                     const index_t vB0, const index_t vB1) {}

  virtual void flip(const shell::shell &M, const index_t &i){};

  /// Called before `dissolve_vertex` removes vertex `v` and triangulates its
  /// ring. The new diagonals reuse the freed spoke and fan slots, which keep
  /// their old values (as `flip` does).
  virtual void dissolve(const shell::shell &M, const index_t &v) {}

  /// Rod topology hooks (corner-indexed). Default: no-op on base.
  virtual void rod_split(rod::rod & /*R*/, index_t /*c_prev*/, index_t /*c0*/,
                         index_t /*c1*/, index_t /*c_next*/,
                         index_t /*c_new*/) {}
  virtual void rod_collapse(rod::rod & /*R*/, index_t /*c_prev*/,
                            index_t /*c0*/, index_t /*c1*/,
                            index_t /*c_next*/) {}

  virtual void calc(const rod::rod &R, const index_t &i, const index_t &prim_id,
                    const real &C, const std::vector<index_t> &vals) = 0;

  virtual void do_map(const index_t i, const index_t it) = 0;
  virtual void do_permute(const std::vector<index_t> &permute) = 0;

  virtual void write(FILE *file) const = 0;

  virtual void compute_vert_split(size_t v0, size_t v1) {}

  virtual void read(FILE *file) = 0;
  virtual void clear() = 0;
  const prim_type &type() { return __type; };

  prim_type __type;
  /// Index in the owning GeometryChannels / topology pool. Set at insert.
  index_t slot_id = -1;
};

using datum_ptr = datum::ptr;

template <typename TYPE> struct datum_t : public datum {
public:
  typedef std::shared_ptr<datum_t<TYPE>> ptr;

  static ptr create(prim_type type, const std::vector<TYPE> &data) {
    return std::make_shared<datum_t<TYPE>>(type, data);
  }

  datum_t(prim_type type, const std::vector<TYPE> &data)
      : datum(type), __data(data){};
  virtual ~datum_t(){};

  std::vector<TYPE> &data() { return __data; }
  const std::vector<TYPE> &data() const { return __data; }
  std::vector<TYPE> &data_mut() { return __data; }
  virtual size_t size() { return __data.size(); };

  virtual void calc(const shell::shell &M, const index_t &i,
                    const index_t &prim_id, const real &C,
                    const std::vector<index_t> &vals) {
    if (__type == VERTEX) {
      // real iN = 1.0 / real(vals.size());
      // TYPE vavg = z::zero<TYPE>();
      // for (const auto &iv : vals) {
      //   vavg += iN * __data[iv];
      // }
      // std::cout << C
      this->__data[i] =
          va::mix(C, this->__data[vals[0]], this->__data[vals[1]]);
      // this->__data[i] = va::mix(C, this->__data[0], this->__data[1]);
      //__data[i] = vavg;
    }

    else if (__type == EDGE || __type == FACE) {
      if constexpr (std::is_same_v<TYPE, quat>) {
        __data[i] = __data[vals[0]];
      } else {
        TYPE v0 = __data[vals[0]];
        __data[i] = C * v0;
      }
    }
  }

  virtual void subdivide(const shell::shell &M, //
                         const index_t &i0, const index_t &i1,
                         const index_t &i2, const index_t &i3, //
                         const real &s, const index_t &source_corner) {
    shell::CornerId c0 = shell::corner_id(source_corner);
    shell::CornerId c1 = M.other(c0);
    if (__type == VERTEX) {
      index_t v0 = M.vert(c0);
      index_t v1 = M.vert(c1);
      this->__data[i0] = va::mix(s, this->__data[v0], this->__data[v1]);
    } else if (__type == EDGE) {
      shell::FaceId f0 = M.face(c0);
      shell::FaceId f1 = M.face(c1);
      (void)f0;
      (void)f1;
      /*not implemented*/
    } else if (__type == FACE) {
      index_t e0 = c0 / 2;
      (void)e0;
      /*not implemented*/
    }
  };

  /// Called before `collapse_edge(M, i)`. That keeps `vert(other(i))`, so
  /// the edge midpoint goes there.
  virtual void collapse(const shell::shell &M, const index_t &i) {
    shell::CornerId c0 = shell::corner_id(i);
    shell::CornerId c1 = M.other(c0);
    if (__type == VERTEX) {
      index_t v0 = M.vert(c0);
      index_t v1 = M.vert(c1);
      this->__data[v1] = va::mix(0.5, this->__data[v0], this->__data[v1]);
    } else if (__type == EDGE) {
      shell::FaceId f0 = M.face(c0);
      shell::FaceId f1 = M.face(c1);
      (void)f0;
      (void)f1;
      /*not implemented*/
    } else if (__type == FACE) {
      index_t e0 = c0 / 2;
      (void)e0;
      /*not implemented*/
    }
  };

  virtual void merge(const shell::shell &M,                //
                     const index_t &i0, const index_t &i1, //
                     const index_t &i2, const index_t &i3, //
                     const index_t vA0, const index_t vA1, //
                     const index_t vB0, const index_t vB1) {
    TYPE tvA0 = this->__data[vA0];
    TYPE tvA1 = this->__data[vA1];
    TYPE tvB0 = this->__data[vB0];
    TYPE tvB1 = this->__data[vB1];

    if (__type == VERTEX) {
      this->__data[i0] = va::mix(0.5, this->__data[vA0], this->__data[vB0]);
      this->__data[i1] = va::mix(0.5, this->__data[vA1], this->__data[vB1]);
      this->__data[i2] = va::mix(1.0, this->__data[vA0], this->__data[vB0]);
      this->__data[i3] = va::mix(1.0, this->__data[vA1], this->__data[vB1]);
    } else if (__type == EDGE) {
      // index_t f0 = M.face(c0);
      // index_t f1 = M.face(c1);
      /*not implemented*/
    }
  };

  virtual void calc(const rod::rod &R, const index_t &i, const index_t &prim_id,
                    const real &C, const std::vector<index_t> &vals) {
    real iN = 1.0 / real(vals.size());
    TYPE vavg = z::zero<TYPE>();
    for (const auto &iv : vals) {
      vavg += iN * __data[iv];
    }

    __data[i] = vavg;
  }

  virtual void rod_split(rod::rod & /*R*/, index_t /*c_prev*/, index_t c0,
                         index_t c1, index_t /*c_next*/, index_t c_new) override {
    ensure_index(c_new);
    if constexpr (std::is_same_v<TYPE, quat>) {
      this->__data[c_new] = va::slerp(this->__data[c0], this->__data[c1], 0.5f);
    } else {
      this->__data[c_new] =
          va::mix(real(0.5), this->__data[c0], this->__data[c1]);
    }
  }

  virtual void rod_collapse(rod::rod & /*R*/, index_t c_prev, index_t /*c0*/,
                            index_t c1, index_t /*c_next*/) override {
    if (c_prev < 0)
      return;
    ensure_index(c1);
    if constexpr (std::is_same_v<TYPE, quat>) {
      this->__data[c1] =
          va::slerp(this->__data[c_prev], this->__data[c1], 0.5f);
    } else {
      this->__data[c1] =
          va::mix(real(0.5), this->__data[c_prev], this->__data[c1]);
    }
  }

  void ensure_index(index_t i) {
    if (i >= static_cast<index_t>(__data.size()))
      __data.resize(static_cast<size_t>(i) + 1, z::zero<TYPE>());
  }

  virtual void flip(const shell::shell &M, const index_t &i) {}

  virtual void do_alloc(const size_t &sz) {
    constexpr size_t kMaxDatumAlloc = 1ull << 26;
    if (sz > kMaxDatumAlloc) {
      std::cerr << "datum::do_alloc absurd size sz=" << sz
                << " data.size()=" << __data.size() << std::endl;
      std::abort();
    }
    __tmp = std::vector<TYPE>(sz, z::zero<TYPE>());
    __data.resize(__data.size() + sz, z::zero<TYPE>());
  }

  virtual void do_resize(const size_t &sz) { __data.resize(sz); }
  virtual void do_map(const index_t i, const index_t it) {
    //__data[i] = __tmp[it];
  }

  virtual void do_permute(const std::vector<index_t> &permute) {
    std::vector<TYPE> n_data(__data);
    // std::cout << "permute.size(): " << permute.size() << std::endl;
    // std::cout << "__data.size(): " << __data.size() << std::endl;
    size_t n = std::min(__data.size(), permute.size());
    for (int i = 0; i < n; i++) {
      // std::cout << i << " " << permute[i] << " " << __data.size() << " "
      //           << n_data.size() << std::endl;
      // std::cout << i << " " << __data[permute[i]] << std::endl;

      // std::cout << "permute.size(): " << permute.size() << std::endl;
      if (permute[i] >= __data.size()) {
        std::cout << "permute[i]: " << permute[i] << std::endl;
        std::cout << "__data.size(): " << __data.size() << std::endl;
        std::cout << "n_data.size(): " << n_data.size() << std::endl;
      }

      n_data[i] = __data[permute[i]];
    }
    // std::cout << "permute.size(): " << permute.size() << std::endl;
    __data = n_data;
  };

  unsigned long operator[](int i) const { return __data[i]; }
  unsigned long &operator[](int i) { return __data[i]; }

  virtual void write(FILE *file) const {
    size_t nData = __data.size();
    size_t e;
    e = fwrite((void *)&nData, sizeof(size_t), 1, file);
    e = fwrite((void *)__data.data(), sizeof(TYPE), nData, file);
  }

  virtual void read(FILE *file) {
    size_t nData;
    size_t e;
    e = fread((void *)&nData, sizeof(size_t), 1, file);
    __data.resize(nData);
    e = fread((void *)__data.data(), sizeof(TYPE), nData, file);
  }

  virtual void clear() { __data.clear(); }

  std::vector<TYPE> __data;
  std::vector<TYPE> __tmp;
};

using scalar_datum = datum_t<real>;
using vec3_datum = datum_t<vec3>;
using vec4_datum = datum_t<vec4>;
using quat_datum = datum_t<quat>;

template <>
inline void datum_t<quat>::calc(const rod::rod & /*R*/, const index_t &i,
                                const index_t & /*prim_id*/, const real & /*C*/,
                                const std::vector<index_t> &vals) {
  if (vals.empty())
    return;
  ensure_index(i);
  __data[i] = __data[vals[0]];
}

std::vector<real> &get_real_data(shell::shell &M, index_t h) {
  return static_pointer_cast<scalar_datum>(M.get_datum(h))->data();
}

std::vector<vec3> &get_vec_data(shell::shell &M, index_t h) {
  return static_pointer_cast<vec3_datum>(M.get_datum(h))->data();
}
const std::vector<vec3> &const_get_vec_data(const shell::shell &M, index_t h) {
  return static_pointer_cast<vec3_datum>(M.const_get_datum(h))->data();
}

template <typename T> index_t init_vert_datum(shell::shell &M, T v0) {
  std::vector<T> rx(M.vert_count(), v0);
  typename datum_t<T>::ptr adata = datum_t<T>::create(prim_type::VERTEX, rx);
  return M.insert_datum(adata);
}

template <typename T>
std::vector<T> &channel_data(GeometryChannels &ch, index_t slot) {
  return static_pointer_cast<datum_t<T>>(ch.get(slot))->data();
}

template <typename T>
const std::vector<T> &channel_data(const GeometryChannels &ch, index_t slot) {
  return static_pointer_cast<datum_t<T>>(ch.get(slot))->data();
}

template <typename T>
std::vector<T> &channel_data(GeometryChannels &ch, const std::string &name) {
  return static_pointer_cast<datum_t<T>>(ch.get(name))->data();
}

inline std::vector<vec3> &velocity(GeometryChannels &ch) {
  return channel_data<vec3>(ch, ch.velocity_slot());
}

inline const std::vector<vec3> &velocity(const GeometryChannels &ch) {
  return channel_data<vec3>(ch, ch.velocity_slot());
}

} // namespace asawa

namespace asawa {
namespace shell {

inline index_t shell::insert_datum(datum_ptr datum) {
  const index_t id = __channels.insert(datum);
  if (datum)
    datum->slot_id = id;
  return id;
}

} // namespace shell
} // namespace asawa
} // namespace gaudi
#endif