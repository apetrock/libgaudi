#ifndef __AAABBB__
#define __AAABBB__

// ---------------------------------------------------------------------------
// Legacy half-space AABB tree, modernized onto the SimplexView model.
// ---------------------------------------------------------------------------
// This used to be a bucket tree (1 leaf = many primitives, addressed by an
// `S`/`NODE_S` stride). It now produces the SAME canonical representation as
// the Morton/radix BVH in hash_tree.hpp:
//
//   indices_         sorted-leaf id -> original simplex id (permutation)
//   internal_nodes_  radix_tree_node[] (N-1 nodes, Karras child/parent links)
//   leaf_nodes_      radix_tree_node[] (N leaves, 1 leaf = 1 simplex)
//   bvh_ / coms_     per-node extents + centers of mass
//
// The only difference from the Morton backend is how the binary tree is built:
// here the leaf order and node topology come from recursive half-space splits
// on simplex centroids (median on the widest axis) instead of a Morton sort.
// Because the emitted representation is identical, every downstream consumer
// (arp::tree_view, build_pyramid, calder Barnes-Hut traversal + integrators,
// get_nearest / get_neighbors) works against this tree unchanged.

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/geometry_types.hpp"
#include "gaudi/vec_addendum.h"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <tuple>
#include <vector>

namespace gaudi {
namespace arp {

// ---------------------------------------------------------------------------
// Standalone primitive-distance callbacks (used by asawa shell dynamics).
// Kept verbatim from the legacy tree; not tied to the tree structure.
// ---------------------------------------------------------------------------
real pnt_tri_min(const index_t &idT, //
                 const std::vector<index_t> &t_inds,
                 const vector<vec3> &t_x, //
                 const index_t &idS,      //
                 const std::vector<index_t> &s_inds, const vector<vec3> &s_x) {
  index_t vT0 = t_inds[idT];
  index_t vS0 = s_inds[3 * idS + 0];
  index_t vS1 = s_inds[3 * idS + 1];
  index_t vS2 = s_inds[3 * idS + 2];
  if (vT0 == vS0)
    return std::numeric_limits<real>::max();
  if (vT0 == vS1)
    return std::numeric_limits<real>::max();
  if (vT0 == vS2)
    return std::numeric_limits<real>::max();

  const vec3 &x0 = t_x[vT0];
  const vec3 &xt0 = s_x[vS0];
  const vec3 &xt1 = s_x[vS1];
  const vec3 &xt2 = s_x[vS2];

  std::array<real, 4> cp = va::closest_point({xt0, xt1, xt2}, x0);
  vec3 xT = cp[1] * xt0 + cp[2] * xt1 + cp[3] * xt2;
  return cp[0];
};

real line_line_min(const index_t &idT, //
                   const std::vector<index_t> &t_inds,
                   const vector<vec3> &t_x, //
                   const index_t &idS,      //
                   const std::vector<index_t> &s_inds,
                   const vector<vec3> &s_x) {
  index_t vT0 = t_inds[2 * idT + 0];
  index_t vT1 = t_inds[2 * idT + 1];
  index_t vS0 = s_inds[2 * idS + 0];
  index_t vS1 = s_inds[2 * idS + 1];
  if (idT >= idS)
    return std::numeric_limits<real>::max();

  if (vT0 == vS0)
    return std::numeric_limits<real>::max();
  if (vT1 == vS1)
    return std::numeric_limits<real>::max();
  if (vT0 == vS1)
    return std::numeric_limits<real>::max();
  if (vT1 == vS0)
    return std::numeric_limits<real>::max();

  const vec3 &xA0 = t_x[t_inds[2 * idT + 0]];
  const vec3 &xA1 = t_x[t_inds[2 * idT + 1]];
  const vec3 &xB0 = s_x[s_inds[2 * idS + 0]];
  const vec3 &xB1 = s_x[s_inds[2 * idS + 1]];
  std::array<real, 3> d = va::distance_Segment_Segment(xA0, xA1, xB0, xB1);
  return d[0];
};

// ---------------------------------------------------------------------------
// Half-space radix builder.
// ---------------------------------------------------------------------------
// Partition order[s..e] (inclusive) about the median centroid on the widest
// axis, returning m = last index of the left sub-range (left = [s, m], right =
// [m+1, e]).
inline index_t aabb_split(index_t s, index_t e, std::vector<index_t> &order,
                          const std::vector<vec3> &cen) {
  vec3 lo = cen[order[s]];
  vec3 hi = lo;
  for (index_t i = s + 1; i <= e; i++) {
    const vec3 &c = cen[order[i]];
    lo = va::min(lo, c);
    hi = va::max(hi, c);
  }
  const vec3 d = hi - lo;
  int axis = (d[1] > d[0]) ? 1 : 0;
  if (d[2] > d[axis])
    axis = 2;

  const index_t mid = s + (e - s) / 2;
  std::nth_element(order.begin() + s, order.begin() + mid,
                   order.begin() + e + 1, [&](index_t a, index_t b) {
                     return cen[a][axis] < cen[b][axis];
                   });
  return mid;
}

// Build the canonical (indices, internal_nodes, leaf_nodes) triple from a set
// of per-simplex centroids. Node ids follow the Karras convention used by the
// Morton tree: an internal node covering a left sub-range [s, m] is labeled m,
// a right sub-range [m+1, e] is labeled m+1, and the root is 0 -- a bijection
// onto {0, ..., n-2} so the arrays line up 1:1 with the datum/bvh arrays.
inline std::tuple<std::vector<index_t>, std::vector<radix_tree_node>,
                  std::vector<radix_tree_node>>
build_aabb_radix(const std::vector<vec3> &centroids) {
  const index_t n = static_cast<index_t>(centroids.size());
  std::vector<index_t> order(n);
  std::iota(order.begin(), order.end(), 0);

  if (n == 0)
    return {std::move(order), {}, {}};

  std::vector<radix_tree_node> internal(n - 1);
  std::vector<radix_tree_node> leaf(n);
  for (index_t i = 0; i < n; i++) {
    leaf[i].start = i;
    leaf[i].end = i + 1;
    leaf[i].split = UNULL;
    leaf[i].parent = UNULL;
  }

  if (n == 1)
    return {std::move(order), std::move(internal), std::move(leaf)};

  struct Frame {
    index_t s, e, id, parent;
  };
  std::vector<Frame> stack;
  stack.reserve(64);
  stack.push_back({0, n - 1, 0, UNULL});

  while (!stack.empty()) {
    const Frame f = stack.back();
    stack.pop_back();

    const index_t m = aabb_split(f.s, f.e, order, centroids);
    internal[f.id].start = f.s;
    internal[f.id].end = f.e;
    internal[f.id].split = m;
    internal[f.id].parent = f.parent;

    // Left child covers [s, m]; if singleton it is leaf[s], else internal[m].
    if (m == f.s)
      leaf[f.s].parent = f.id;
    else
      stack.push_back({f.s, m, m, f.id});

    // Right child covers [m+1, e]; if singleton it is leaf[e], else
    // internal[m+1].
    if (m + 1 == f.e)
      leaf[f.e].parent = f.id;
    else
      stack.push_back({m + 1, f.e, m + 1, f.id});
  }

  return {std::move(order), std::move(internal), std::move(leaf)};
}

// ---------------------------------------------------------------------------
// aabb_tree<N>: drop-in alternative backend to bvh_tree<N>.
// ---------------------------------------------------------------------------
// Mirrors bvh_tree<N>'s public surface (members + create/update/get_index/
// leaf_simplex/get_com/find_nearest/find_neighbors/get_nearest) so it satisfies
// arp::tree_view and the calder FMM path identically; only the build differs.
template <int N> class aabb_tree {
public:
  static constexpr int kSimplexN = N;

  std::vector<index_t> indices_;
  std::vector<radix_tree_node> internal_nodes_;
  std::vector<radix_tree_node> leaf_nodes_;
  TreeResult<ext::extents_t> bvh_;
  std::vector<vec3> data_;
  std::vector<index_t> adjacency_;
  std::vector<MassPoint> coms_;

  using ptr = std::shared_ptr<aabb_tree<N>>;
  using permutation_index_type = std::vector<index_t>;
  using permuted_view =
      permuted_simplex_view<N, std::vector<vec3>, std::vector<index_t>>;
  using value_type = typename permuted_view::value_type; // std::array<vec3, N>

  std::optional<permuted_view> permuted_data_view_;

  static ptr create(const std::vector<index_t> &adjacency,
                    const std::vector<vec3> &vertices, int lvl = 8) {
    return std::make_shared<aabb_tree<N>>(vertices, adjacency);
  }

  static ptr create(const simplex_set<N> &set, int lvl = 8) {
    return std::make_shared<aabb_tree<N>>(set);
  }

  aabb_tree(const std::vector<vec3> &data,
            const std::vector<index_t> &adjacency) {
    update(data, adjacency);
  }

  aabb_tree(const simplex_set<N> &set) { update(set); }

  index_t get_index(size_t i) const { return indices_[i]; }

  const std::vector<index_t> &adjacency() const { return adjacency_; }
  const std::vector<vec3> &verts() const { return data_; }
  const vec3 &vert(index_t i) const { return data_[adjacency_[i]]; }

  // Per-simplex geometric centroid in the unpermuted (original) order.
  std::vector<vec3> simplex_centroids() const {
    std::vector<index_t> identity(adjacency_.size() / N);
    std::iota(identity.begin(), identity.end(), 0);
    permuted_view view(const_cast<std::vector<vec3> &>(data_),
                       const_cast<std::vector<index_t> &>(adjacency_), identity);
    std::vector<vec3> cents(view.size());
    for (size_t i = 0; i < view.size(); i++) {
      value_type s = view[i];
      vec3 c = vec3::Zero();
      for (int k = 0; k < N; k++)
        c += s[k];
      cents[i] = c / real(N);
    }
    return cents;
  }

  void update(const std::vector<vec3> &data,
              const std::vector<index_t> &adjacency) {
    data_ = data;
    adjacency_ = adjacency;

    std::vector<vec3> centroids = simplex_centroids();
    auto [indices, internal_nodes, leaf_nodes] = build_aabb_radix(centroids);

    indices_ = std::move(indices);
    internal_nodes_ = std::move(internal_nodes);
    leaf_nodes_ = std::move(leaf_nodes);

    permuted_data_view_.emplace(data_, adjacency_, indices_);

    bvh_ = make_bvh(*permuted_data_view_, internal_nodes_, leaf_nodes_);
    coms_ = calc_com(*permuted_data_view_);
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

  template <int Nq> auto dispatch_test() const {
    if constexpr (Nq == 1 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_point_tuple(q, d); };
    else if constexpr (Nq == 1 && N == 2)
      return [](const auto &q, const auto &d) { return test_point_line_tuple(q, d); };
    else if constexpr (Nq == 1 && N == 3)
      return [](const auto &q, const auto &d) { return test_point_tri_tuple(q, d); };
    else if constexpr (Nq == 2 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_line_tuple(d, q); };
    else if constexpr (Nq == 2 && N == 2)
      return [](const auto &q, const auto &d) { return test_line_line_tuple(q, d); };
    else if constexpr (Nq == 2 && N == 3)
      return [](const auto &q, const auto &d) { return test_line_tri_tuple(q, d); };
    else if constexpr (Nq == 3 && N == 1)
      return [](const auto &q, const auto &d) { return test_point_tri_tuple(d, q); };
    else if constexpr (Nq == 3 && N == 2)
      return [](const auto &q, const auto &d) { return test_line_tri_tuple(d, q); };
    else if constexpr (Nq == 3 && N == 3)
      return [](const auto &q, const auto &d) { return test_tri_tri_tuple(q, d); };
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
    constexpr int Nq = view_stride_v<PTYPE>;
    return arp::get_nearest<Singulus<Nq>, permuted_view>(
        make_query<Nq>(query), *permuted_data_view_, internal_nodes_,
        leaf_nodes_, bvh_, dispatch_test<Nq>());
  }

  template <Vec3View PTYPE>
  std::vector<index_t> find_neighbors(const PTYPE &query, real tol) {
    constexpr int Nq = view_stride_v<PTYPE>;
    return arp::get_neighbors<Singulus<Nq>, permuted_view>(
        make_query<Nq>(query), *permuted_data_view_, internal_nodes_,
        leaf_nodes_, bvh_, tol, dispatch_test<Nq>());
  }

  // Legacy wrapper: returns vector with idMin as last element.
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
};

template <int N> using AABB_T = aabb_tree<N>;

using A1 = aabb_tree<1>;
using A2 = aabb_tree<2>;
using A3 = aabb_tree<3>;

} // namespace arp
} // namespace gaudi
#endif
