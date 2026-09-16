#ifndef GAUDI_ARP_BRUTE_FORCE_TREE_HPP
#define GAUDI_ARP_BRUTE_FORCE_TREE_HPP

#include "gaudi/arp/brute_force.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/ray_traversal.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/arp/tree_view.hpp"

#include <limits>
#include <memory>
#include <numeric>
#include <vector>

namespace gaudi {
namespace arp {

/// Drop-in backend: linear scan over all leaves (no internal culling).
/// Same public surface as bvh_tree / aabb_tree for tests and validation.
template <int N> class brute_force_tree {
public:
  static constexpr int kSimplexN = N;

  std::vector<index_t> indices_;
  std::vector<radix_tree_node> internal_nodes_;
  std::vector<radix_tree_node> leaf_nodes_;
  TreeResult<ext::extents_t> bvh_;
  std::vector<vec3> data_;
  std::vector<index_t> adjacency_;
  std::vector<MassPoint> coms_;

  using ptr = std::shared_ptr<brute_force_tree<N>>;
  using permutation_index_type = std::vector<index_t>;
  using permuted_view =
      permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>;
  using value_type = typename permuted_view::value_type;

  std::optional<permuted_view> permuted_data_view_;

  static ptr create(const std::vector<index_t> &adjacency,
                    const std::vector<vec3> &vertices) {
    return std::make_shared<brute_force_tree<N>>(vertices, adjacency);
  }

  static ptr create(const simplex_set<N> &set) {
    return std::make_shared<brute_force_tree<N>>(set);
  }

  brute_force_tree(const std::vector<vec3> &data,
                   const std::vector<index_t> &adjacency) {
    update(data, adjacency);
  }

  brute_force_tree(const simplex_set<N> &set) { update(set); }

  index_t get_index(size_t i) const { return indices_[i]; }

  const std::vector<index_t> &adjacency() const { return adjacency_; }
  const std::vector<vec3> &verts() const { return data_; }

  void update(const std::vector<vec3> &data,
              const std::vector<index_t> &adjacency) {
    data_ = data;
    adjacency_ = adjacency;
    internal_nodes_.clear();

    const size_t n = adjacency_.size() / static_cast<size_t>(N);
    indices_.resize(n);
    std::iota(indices_.begin(), indices_.end(), 0);

    leaf_nodes_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      leaf_nodes_[i].start = static_cast<index_t>(i);
      leaf_nodes_[i].end = static_cast<index_t>(i + 1);
      leaf_nodes_[i].split = UNULL;
      leaf_nodes_[i].parent = UNULL;
    }

    permuted_data_view_.emplace(data_, adjacency_, indices_);

    bvh_.leaf.clear();
    bvh_.internal.clear();
    if (permuted_data_view_) {
      bvh_.leaf = calc_extents(*permuted_data_view_);
      coms_ = calc_com(*permuted_data_view_);
    } else {
      coms_.clear();
    }
  }

  void update(const simplex_set<N> &set) {
    update(set.vertices(), set.adjacency());
  }

  std::array<index_t, N> get_tuple_ids(const index_t &i) {
    return permuted_data_view_->get_tuple_ids(i);
  }

  value_type get_simplex(index_t i) const { return (*permuted_data_view_)[i]; }

  std::array<vec3, N> leaf_simplex(index_t sorted_id) const {
    return (*permuted_data_view_)[sorted_id];
  }

  vec3 get_com(index_t i) const { return std::get<1>(coms_[i]); }

  tree_view view() const {
    return make_tree_view(*this);
  }

  template <int Nq> auto dispatch_test() const {
    if constexpr (Nq == 1 && N == 1)
      return [](const auto &q, const auto &d) {
        return test_point_point_tuple(q, d);
      };
    else if constexpr (Nq == 1 && N == 2)
      return [](const auto &q, const auto &d) {
        return test_point_line_tuple(q, d);
      };
    else if constexpr (Nq == 1 && N == 3)
      return [](const auto &q, const auto &d) {
        return test_point_tri_tuple(q, d);
      };
    else if constexpr (Nq == 2 && N == 2)
      return [](const auto &q, const auto &d) {
        return test_line_line_tuple(q, d);
      };
    else if constexpr (Nq == 2 && N == 3)
      return [](const auto &q, const auto &d) {
        return test_line_tri_tuple(q, d);
      };
    else if constexpr (Nq == 3 && N == 3)
      return [](const auto &q, const auto &d) {
        return test_tri_tri_tuple(q, d);
      };
  }

  template <int Nq> static auto make_query(const auto &view) {
    if constexpr (Nq == 1)
      return std::array<vec3, 1>{view[0]};
    else if constexpr (Nq == 2)
      return std::array<vec3, 2>{view[0], view[1]};
    else if constexpr (Nq == 3)
      return std::array<vec3, 3>{view[0], view[1], view[2]};
  }

  template <Vec3View PTYPE> index_t find_nearest(const PTYPE &query) {
    if (!permuted_data_view_ || permuted_data_view_->empty())
      return -1;
    constexpr int Nq = view_stride_v<PTYPE>;
    const auto q = make_query<Nq>(query);
    const auto test = dispatch_test<Nq>();

    index_t id_min = -1;
    real d_min = std::numeric_limits<real>::max();
    for (size_t i = 0; i < permuted_data_view_->size(); ++i) {
      const real d = test(q, (*permuted_data_view_)[i]);
      if (d < d_min) {
        d_min = d;
        id_min = permuted_data_view_->get_index(static_cast<index_t>(i));
      }
    }
    return id_min;
  }

  template <Vec3View PTYPE>
  std::vector<index_t> find_neighbors(const PTYPE &query, real tol) {
    if (!permuted_data_view_ || permuted_data_view_->empty())
      return {};
    constexpr int Nq = view_stride_v<PTYPE>;
    const auto q = make_query<Nq>(query);
    const auto test = dispatch_test<Nq>();

    std::vector<index_t> out;
    for (size_t i = 0; i < permuted_data_view_->size(); ++i) {
      const real d = test(q, (*permuted_data_view_)[i]);
      if (d < tol)
        out.push_back(permuted_data_view_->get_index(static_cast<index_t>(i)));
    }
    return out;
  }

  template <Vec3View PTYPE>
  std::vector<index_t> get_nearest(const PTYPE &query, real tol) {
    std::vector<index_t> result;
    if (tol > 999.9) {
      result.push_back(find_nearest<PTYPE>(query));
    } else {
      result = find_neighbors<PTYPE>(query, tol);
      result.push_back(-1);
    }
    return result;
  }

  std::optional<ray_hit<N>> ray_pierce(const Ray &ray) const {
    if (!permuted_data_view_)
      return std::nullopt;
    return ray_slop_pierce_nearest<N>(view(), *permuted_data_view_, ray, 0.0);
  }

  std::vector<ray_hit<N>> ray_collect(const Ray &ray, real tol) const {
    if (!permuted_data_view_)
      return {};
    return ray_collect_within_tol<N>(view(), *permuted_data_view_, ray, tol);
  }

  std::vector<ray_hit<N>> ray_slop_pierce(const Ray &ray, real slop) const {
    if (!permuted_data_view_)
      return {};
    return ray_slop_pierce_heap<N>(view(), *permuted_data_view_, ray, slop);
  }
};

template <int N> using BF_T = brute_force_tree<N>;

using BF1 = brute_force_tree<1>;
using BF2 = brute_force_tree<2>;
using BF3 = brute_force_tree<3>;

} // namespace arp
} // namespace gaudi

#endif // GAUDI_ARP_BRUTE_FORCE_TREE_HPP
