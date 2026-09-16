#ifndef GAUDI_TEST_ARP_QUERY_HARNESS_HPP
#define GAUDI_TEST_ARP_QUERY_HARNESS_HPP

#include "gaudi/arp/brute_force.hpp"
#include "gaudi/arp/pairwise_tests.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/common.h"
#include "gaudi/test/test.hpp"

#include <cmath>
#include <vector>

namespace gaudi {
namespace test {

template <int N, typename ViewT>
inline real arp_query_dist(const std::array<vec3, 1> &q, ViewT &view,
                           index_t id) {
  if constexpr (N == 1)
    return arp::test_point_point(q, slice<1, ViewT>(view, id));
  else if constexpr (N == 2)
    return arp::test_point_line(q, slice<2, ViewT>(view, id));
  else
    return arp::test_point_tri(q, slice<3, ViewT>(view, id));
}

/// Compare backend nearest to brute-force reference (distance primary, id tie-tolerant).
template <int N, typename BackendTree, typename ViewT>
inline void run_nearest_equivalence(BackendTree &tree, ViewT &view,
                                    const std::vector<vec3> &queries,
                                    real tol = 1000.0, real eps = 1e-6) {
  for (const auto &qp : queries) {
    const std::array<vec3, 1> query = {qp};

    auto backend_res = tree.get_nearest(query, tol);
    auto brute_res = arp::brute_force_nearest_auto<1, N>(query, view, tol);

    const index_t bid =
        backend_res.empty() ? -1 : backend_res.back();
    const index_t ref_id = brute_res.empty() ? -1 : brute_res.back();

    GAUDI_ASSERT(bid >= 0);
    GAUDI_ASSERT(ref_id >= 0);

    const real bd = arp_query_dist<N, ViewT>(query, view, bid);
    const real ref_d = arp_query_dist<N, ViewT>(query, view, ref_id);

    GAUDI_EXPECT(std::abs(bd - ref_d) < eps);
    GAUDI_EXPECT(bid == ref_id || std::abs(bd - ref_d) < eps);
  }
}

template <int N, typename BackendTree>
inline void run_nearest_equivalence_mesh(BackendTree &tree,
                                         std::vector<vec3> &vertices,
                                         std::vector<index_t> &adj,
                                         const std::vector<vec3> &queries,
                                         real tol = 1000.0, real eps = 1e-6) {
  adjacency_view<std::vector<vec3>, std::vector<index_t>> view(vertices, adj);
  run_nearest_equivalence<N, BackendTree, decltype(view)>(tree, view, queries,
                                                          tol, eps);
}

} // namespace test
} // namespace gaudi

#endif // GAUDI_TEST_ARP_QUERY_HARNESS_HPP
