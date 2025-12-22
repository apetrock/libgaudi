#ifndef __GAUDI_ARP_BRUTE_FORCE__
#define __GAUDI_ARP_BRUTE_FORCE__

#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/common.h"
#include <limits>
#include <vector>

namespace gaudi {
namespace arp {

// Simple O(n*m) brute-force nearest search for validation
// Returns indices of data elements within tolerance of query, or nearest if
// contracting_rad mode (tol > 999.9)
template <int N, Vec3View PTYPE, Vec3View TTYPE>
std::vector<index_t> brute_force_nearest(const PTYPE &query, const TTYPE &data,
                                         real tol, auto &&testAB) {
  if (data.empty()) {
    return {};
  }

  bool contracting_rad = tol > 999.9;
  std::vector<index_t> collisions;
  index_t idMin = -1;
  real mMin = std::numeric_limits<real>::max();

  // Iterate over all data elements (stride by N)
  for (size_t i = 0; i < data.size(); i += N) {
    slice<N, TTYPE> datum(data, i / N);
    real dist = testAB(query, datum);

    if (contracting_rad) {
      // Find single nearest
      if (dist < mMin) {
        mMin = dist;
        idMin = static_cast<index_t>(i / N);
      }
    } else {
      // Collect all within tolerance
      if (dist < tol) {
        collisions.push_back(static_cast<index_t>(i / N));
      }
    }
  }

  if (contracting_rad) {
    collisions.push_back(idMin);
  } else {
    collisions.push_back(-1); // Sentinel to match BVH behavior
  }

  return collisions;
}

// Convenience wrapper that selects the appropriate test function from
// pairwise_tests.hpp
template <int Nquery, int Ndata, Vec3View PTYPE, Vec3View TTYPE>
std::vector<index_t> brute_force_nearest_auto(const PTYPE &query,
                                              const TTYPE &data, real tol) {
  if constexpr (Nquery == 1 && Ndata == 1) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_point_point(q, d);
        });
  } else if constexpr (Nquery == 1 && Ndata == 2) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_point_line(q, d);
        });
  } else if constexpr (Nquery == 1 && Ndata == 3) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_point_tri(q, d);
        });
  } else if constexpr (Nquery == 2 && Ndata == 2) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_line_line(q, d);
        });
  } else if constexpr (Nquery == 2 && Ndata == 3) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_line_tri(q, d);
        });
  } else if constexpr (Nquery == 3 && Ndata == 3) {
    return brute_force_nearest<Ndata>(
        query, data, tol,
        [](const PTYPE &q, const slice<Ndata, TTYPE> &d) {
          return test_tri_tri(q, d);
        });
  } else {
    return {};
  }
}

} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_BRUTE_FORCE__
