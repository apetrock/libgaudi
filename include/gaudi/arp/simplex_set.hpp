#ifndef __GAUDI_ARP_SIMPLEX_SET__
#define __GAUDI_ARP_SIMPLEX_SET__

#include "gaudi/arp/morton.hpp"
#include "gaudi/common.h"
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace arp {

template <int N> class simplex_set {
public:
  using vertex_container = std::vector<vec3>;
  using adjacency_container = std::vector<index_t>;
  using view_type = adjacency_view<vertex_container, adjacency_container>;
  // Permuted simplex view type (returns tuples)
  using permuted_view_type =
      permuted_simplex_view<N, vertex_container, adjacency_container>;
  // Alias for clarity
  using simplex_view_type = permuted_view_type;

  simplex_set(const vertex_container &vertices,
              const adjacency_container &adjacency)
      : vertices_(&vertices), adjacency_(&adjacency), view_(vertices, adjacency) {
    refresh();
  }

  void refresh() {
    if (adjacency_->size() % N != 0) {
      throw std::runtime_error("adjacency size must be a multiple of " +
                               std::to_string(N));
    }
    coms_ = calc_com<N>(view_);
    centroids_.clear();
    centroids_.reserve(coms_.size());
    for (const auto &[mass, com] : coms_) {
      centroids_.push_back(com);
    }
  }

  const vertex_container &vertices() const { return *vertices_; }
  const adjacency_container &adjacency() const { return *adjacency_; }
  const view_type &view() const { return view_; }
  size_t size() const { return adjacency_->size() / N; }

  std::array<index_t, N> tuple_ids(index_t simplex_id) const {
    return view_.template get_tuple_ids<N>(simplex_id);
  }

  slice<N, view_type> simplex(index_t simplex_id) const {
    return slice<N, view_type>(view_, simplex_id);
  }

  slice<N, view_type> take(index_t simplex_id) const {
    return slice<N, view_type>(view_, simplex_id);
  }

  // Returns permuted simplex view (tuples of vec3)
  permuted_view_type permuted_view(
      const std::vector<index_t> &permutation) const {
    return permuted_view_type(*vertices_, *adjacency_, permutation);
  }

  // Alias for consistency
  simplex_view_type simplex_view(
      const std::vector<index_t> &permutation) const {
    return permuted_view(permutation);
  }

  const std::vector<MassPoint> &coms() const { return coms_; }
  const std::vector<vec3> &centroids() const { return centroids_; }

private:
  const vertex_container *vertices_;
  const adjacency_container *adjacency_;
  view_type view_;
  std::vector<MassPoint> coms_;
  std::vector<vec3> centroids_;
};

} // namespace arp
} // namespace gaudi

#endif // __GAUDI_ARP_SIMPLEX_SET__
