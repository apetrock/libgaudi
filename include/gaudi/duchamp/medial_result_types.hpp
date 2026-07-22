#ifndef __GAUDI_DUCHAMP_MEDIAL_RESULT_TYPES__
#define __GAUDI_DUCHAMP_MEDIAL_RESULT_TYPES__

#include <vector>

#include "gaudi/common.h"

namespace gaudi {
namespace duchamp {

// Lean search result for helpers. Graph nodes emit parallel arrays:
//   output: field_datum<vec3>
//   mask:   field_datum<real>  (1 = accepted, 0 = rejected)
struct medial_point {
  vec3 point = vec3::Zero();
  bool accepted = false;
};

inline int count_medial_accepted(const std::vector<real> &mask) {
  int n = 0;
  for (real m : mask) {
    if (m > real(0.5)) {
      ++n;
    }
  }
  return n;
}

inline int count_medial_accepted(const std::vector<medial_point> &points) {
  int n = 0;
  for (const medial_point &p : points) {
    if (p.accepted) {
      ++n;
    }
  }
  return n;
}

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_MEDIAL_RESULT_TYPES__
