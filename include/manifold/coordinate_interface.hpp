
#ifndef __TWOMANIFOLD_COORDINATE_INTERFACE__
#define __TWOMANIFOLD_COORDINATE_INTERFACE__
#include "m2Includes.h"

namespace m2 {

template <typename SPACE>
typename SPACE::coordinate_type get_coordinate(typename surf<SPACE>::vertex_ptr v) {
  return v->template get<typename SPACE::coordinate_type>(
      SPACE::vertex_index::COORDINATE);
}

template <typename SPACE>
void set_coordinate(typename SPACE::coordinate_type p,
                             typename surf<SPACE>::vertex_ptr v) {
  v->template set<typename SPACE::coordinate_type>(
      p, SPACE::vertex_index::COORDINATE);
}

} // namespace m2

#endif