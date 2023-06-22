
#ifndef __HEP_CONSTRAINTS__
#define __HEP_CONSTRAINTS__

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <memory.h>
#include <numeric>
#include <ostream>
#include <stdio.h>
#include <type_traits>
#include <vector>
#include <zlib.h>

#include "gaudi/common.h"
#include "gaudi/define_create_func.h"

namespace gaudi {
namespace hepworth {

inline void block33(int i, int j, real v, std::vector<trip> &triplets) {
  for (int ax = 0; ax < 3; ax++)
    triplets.push_back(trip(3 * i + ax, 3 * j + ax, v));
}
inline vec3 from3(int i, const vecX &q) { return q.block(3 * i, 0, 3, 1); }
inline void to3(int i, vecX &p, vec3 val) { p.block(i, 0, 3, 1) = val; }

class projection_constraint {
public:
  typedef std::shared_ptr<projection_constraint> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w) {
    return std::make_shared<projection_constraint>(ids, w);
  }

  projection_constraint(const std::vector<index_t> &ids, const real &w)
      : _ids(ids), _w(w) {}
  virtual void project(const vecX &q, vecX &p){};
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets){};

  std::vector<index_t> _ids;
  index_t _id0;
  real _w;
};

} // namespace hepworth
} // namespace gaudi
#endif