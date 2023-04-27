#include <cassert>
#include <cstddef>
#include <cxxabi.h>

#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <ostream>
#include <stdio.h>
#include <vector>
#include <zlib.h>

#include "gaudi/common.h"
#include "gaudi/geometry_types.hpp"
#include "shell/datum_x.hpp"
#include "shell/shell.hpp"

#include "rod/rod.hpp"

#ifndef __ASAWA_DATUM__
#define __ASAWA_DATUM__
namespace gaudi {

namespace asawa {

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
  virtual void calc(const shell::shell &M, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) = 0;

  virtual void flip(const shell::shell &M, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) = 0;

  virtual void calc(const rod::rod &R, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) = 0;

  virtual void do_map(const index_t i, const index_t it) = 0;
  virtual void do_permute(const std::vector<index_t> &permute) = 0;

  virtual void write(FILE *file) const = 0;

  virtual void compute_vert_split(size_t v0, size_t v1) {}

  virtual void read(FILE *file) = 0;
  virtual void clear() = 0;
  const prim_type &type() { return __type; };

  prim_type __type;
};

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
  virtual size_t size() { return __data.size(); };

  virtual void calc(const shell::shell &M, const index_t &i, const real &C,
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
      TYPE v0 = __data[vals[0]];
      __data[i] = C * v0;
    }
  }

  virtual void calc(const rod::rod &R, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) {
    real iN = 1.0 / real(vals.size());
    TYPE vavg = z::zero<TYPE>();
    for (const auto &iv : vals) {
      vavg += iN * __data[iv];
    }

    __data[i] = vavg;
  }

  virtual void flip(const shell::shell &M, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) {}

  virtual void do_alloc(const size_t &sz) {
    __tmp = std::vector<TYPE>(sz, z::zero<TYPE>());
    __data.resize(__data.size() + sz, z::zero<TYPE>());
  }

  virtual void do_resize(const size_t &sz) { __data.resize(sz); }
  virtual void do_map(const index_t i, const index_t it) {
    //__data[i] = __tmp[it];
  }

  virtual void do_permute(const std::vector<index_t> &permute) {
    std::vector<TYPE> n_data(__data);
    for (int i = 0; i < n_data.size(); i++) {
      n_data[i] = __data[permute[i]];
    }
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

struct vec3_datum : public datum_t<vec3> {
public:
  typedef std::shared_ptr<vec3_datum> ptr;

  static ptr create(prim_type type, const std::vector<vec3> &data) {
    return std::make_shared<vec3_datum>(type, data);
  }

  vec3_datum(prim_type type, const std::vector<vec3> &data)
      : datum_t<vec3>(type, data){};
  virtual ~vec3_datum(){};

  virtual void calc(const shell::shell &M, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) {
    real w = 0.0;
    /*
    vec3 vavg = z::zero<vec3>();
    for (const auto &iv : vals) {
      real wi = vert_area(M, iv, this->__data);
      // real wi = vert_cotan_weight(M, iv, this->__data);
      w += wi;
      vavg += wi * this->__data[iv];
    }
    */
    this->__data[i] = va::mix(C, this->__data[vals[0]], this->__data[vals[1]]);
  }

  virtual void flip(const shell::shell &M, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) {}

  virtual void calc(const rod::rod &R, const index_t &i, const real &C,
                    const std::vector<index_t> &vals) {}
}; // namespace asawa

std::vector<real> &get_real_data(shell::shell &M, index_t h) {
  return static_pointer_cast<scalar_datum>(M.get_datum(h))->data();
}

std::vector<vec3> &get_vec_data(shell::shell &M, index_t h) {
  return static_pointer_cast<vec3_datum>(M.get_datum(h))->data();
}

} // namespace asawa
} // namespace gaudi
#endif