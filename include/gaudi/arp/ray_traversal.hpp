#ifndef GAUDI_ARP_RAY_TRAVERSAL_HPP
#define GAUDI_ARP_RAY_TRAVERSAL_HPP

#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/ray_simplex.hpp"
#include "gaudi/arp/tree_view.hpp"
#include "gaudi/geometry_types.hpp"

#include <functional>
#include <optional>
#include <vector>

namespace gaudi {
namespace arp {

/// Min-heap on ray_hit::t_ray; on ties prefer smaller axis distance when ray is set.
template <int Ndata> struct ray_hit_heap_cmp {
  const Ray *ray = nullptr;

  bool operator()(const ray_hit<Ndata> &a, const ray_hit<Ndata> &b) const {
    constexpr real tie = 1e-10;
    if (a.t_ray > b.t_ray + tie) {
      return true;
    }
    if (b.t_ray > a.t_ray + tie) {
      return false;
    }
    if (!ray) {
      return false;
    }
    return ray_distance_to_point(*ray, a.point) >
           ray_distance_to_point(*ray, b.point);
  }
};

template <int Ndata> using ray_hit_t_min = ray_hit_heap_cmp<Ndata>;

template <int Ndata>
inline void ray_hits_make_min_heap(std::vector<ray_hit<Ndata>> &hits,
                                   const Ray *ray = nullptr) {
  std::make_heap(hits.begin(), hits.end(), ray_hit_heap_cmp<Ndata>{ray});
}

template <int Ndata>
inline bool ray_hits_heap_empty(const std::vector<ray_hit<Ndata>> &hits) {
  return hits.empty();
}

template <int Ndata>
inline ray_hit<Ndata> ray_hits_pop_min(std::vector<ray_hit<Ndata>> &hits,
                                       const Ray *ray = nullptr) {
  std::pop_heap(hits.begin(), hits.end(), ray_hit_heap_cmp<Ndata>{ray});
  ray_hit<Ndata> out = hits.back();
  hits.pop_back();
  return out;
}

inline std::optional<std::pair<real, real>>
ray_aabb_interval(const Ray &ray, const ext::extents_t &box) {
  const vec3 inv_dir(
      std::abs(ray.dir[0]) > 1e-12 ? real(1.0) / ray.dir[0] : real(0.0),
      std::abs(ray.dir[1]) > 1e-12 ? real(1.0) / ray.dir[1] : real(0.0),
      std::abs(ray.dir[2]) > 1e-12 ? real(1.0) / ray.dir[2] : real(0.0));

  real t0 = ray.t_min;
  real t1 = ray.t_max;

  for (int axis = 0; axis < 3; ++axis) {
    if (std::abs(ray.dir[axis]) <= 1e-12) {
      if (ray.origin[axis] < box[0][axis] - 1e-10 ||
          ray.origin[axis] > box[1][axis] + 1e-10) {
        return std::nullopt;
      }
      continue;
    }
    real t_near = (box[0][axis] - ray.origin[axis]) * inv_dir[axis];
    real t_far = (box[1][axis] - ray.origin[axis]) * inv_dir[axis];
    if (inv_dir[axis] < 0.0)
      std::swap(t_near, t_far);
    t0 = std::max(t0, t_near);
    t1 = std::min(t1, t_far);
    if (t0 > t1)
      return std::nullopt;
  }
  return std::make_pair(t0, t1);
}

inline bool ray_aabb_overlaps(const Ray &ray, const ext::extents_t &box) {
  return ray_aabb_interval(ray, box).has_value();
}

inline bool ray_slop_intersects_bbox(const Ray &ray, real slop,
                                     const ext::extents_t &box) {
  if (slop <= 0.0) {
    return ray_aabb_overlaps(ray, box);
  }
  return ray_aabb_interval(ray, ext::inflate(box, slop)).has_value();
}

template <typename LeafFn>
inline void for_each_leaf(const tree_view &view, LeafFn &&leaf_fn) {
  for (size_t i = 0; i < view.leaf_count(); ++i)
    leaf_fn(static_cast<index_t>(i));
}

template <typename LeafFn, typename CullFn>
inline void ray_slop_traverse_leaves(const tree_view &view, LeafFn &&leaf_fn,
                                     CullFn &&should_continue) {
  if (view.empty())
    return;

  if (view.internal_count() == 0 || !view.has_bbox()) {
    for_each_leaf(view, std::forward<LeafFn>(leaf_fn));
    return;
  }

  traverse_bfs(
      view.internal_nodes(), view.leaf_nodes(),
      [&](index_t, const radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const radix_tree_node &) {
        leaf_fn(leaf_id);
      },
      [&](index_t node_id, const radix_tree_node &) -> bool {
        return should_continue(view.internal_bbox(node_id));
      });
}

template <int Ndata, SimplexView DataView>
inline std::vector<ray_hit<Ndata>>
ray_slop_pierce_heap(const tree_view &view, const DataView &data,
                     const Ray &ray, real slop) {
  std::vector<ray_hit<Ndata>> hits;

  auto try_leaf = [&](index_t leaf_id) {
    const index_t pid = view.index(leaf_id);
    const auto simplex = data[leaf_id];
    if (auto hit = ray_slop_intersects_simplex<Ndata>(ray, slop, simplex, pid)) {
      hits.push_back(*hit);
    }
  };

  ray_slop_traverse_leaves(
      view, try_leaf,
      [&](const ext::extents_t &box) { return ray_slop_intersects_bbox(ray, slop, box); });

  if (!hits.empty()) {
    ray_hits_make_min_heap(hits, &ray);
  }
  return hits;
}

template <int Ndata, SimplexView DataView>
inline std::optional<ray_hit<Ndata>>
ray_slop_pierce_nearest(const tree_view &view, const DataView &data,
                        const Ray &ray, real slop,
                        const std::function<bool(index_t)> &reject = {}) {
  std::vector<ray_hit<Ndata>> hits =
      ray_slop_pierce_heap<Ndata>(view, data, ray, slop);
  while (!hits.empty()) {
    ray_hit<Ndata> hit = ray_hits_pop_min(hits, &ray);
    if (!reject || !reject(hit.id)) {
      return hit;
    }
  }
  return std::nullopt;
}

/// Tree wrapper: works with bvh_tree, aabb_tree, brute_force_tree.
template <typename TreeT>
inline auto ray_pierce_tree(const TreeT &tree, const Ray &ray) {
  constexpr int N = TreeT::kSimplexN;
  if (!tree.permuted_data_view_)
    return std::optional<ray_hit<N>>{};
  return ray_slop_pierce_nearest<N>(make_tree_view(tree),
                                  *tree.permuted_data_view_, ray, 0.0);
}

template <typename TreeT>
inline auto ray_slop_pierce_tree(const TreeT &tree, const Ray &ray, real slop) {
  constexpr int N = TreeT::kSimplexN;
  if (!tree.permuted_data_view_) {
    return std::vector<ray_hit<N>>{};
  }
  return ray_slop_pierce_heap<N>(make_tree_view(tree),
                                 *tree.permuted_data_view_, ray, slop);
}

template <typename TreeT>
inline auto ray_slop_pierce_tree_nearest(const TreeT &tree, const Ray &ray,
                                         real slop,
                                         const std::function<bool(index_t)> &reject =
                                             {}) {
  constexpr int N = TreeT::kSimplexN;
  if (!tree.permuted_data_view_) {
    return std::optional<ray_hit<N>>{};
  }
  return ray_slop_pierce_nearest<N>(make_tree_view(tree),
                                    *tree.permuted_data_view_, ray, slop, reject);
}

template <int Ndata, SimplexView DataView>
inline std::vector<ray_hit<Ndata>>
ray_collect_within_tol(const tree_view &view, const DataView &data,
                       const Ray &ray, real tol) {
  std::vector<ray_hit<Ndata>> hits;

  auto try_leaf = [&](index_t leaf_id) {
    const index_t pid = view.index(leaf_id);
    const auto simplex = data[leaf_id];
    ray_collect_simplex_within_tol<Ndata>(ray, simplex, tol, pid, hits);
  };

  if (view.internal_count() == 0 || !view.has_bbox()) {
    for_each_leaf(view, try_leaf);
    return hits;
  }

  traverse_bfs(
      view.internal_nodes(), view.leaf_nodes(),
      [&](index_t, const radix_tree_node &) {},
      [&](index_t, index_t leaf_id, const radix_tree_node &) {
        try_leaf(leaf_id);
      },
      [&](index_t node_id, const radix_tree_node &) -> bool {
        const ext::extents_t box = view.internal_bbox(node_id);
        const vec3 c = ext::center(box);
        const real r = ext::mag(box) * 0.5 + tol;
        return ray_distance_to_point(ray, c) <= r;
      });

  return hits;
}

template <typename TreeT>
inline auto ray_collect_tree(const TreeT &tree, const Ray &ray, real tol) {
  constexpr int N = TreeT::kSimplexN;
  if (!tree.permuted_data_view_)
    return std::vector<ray_hit<N>>{};
  return ray_collect_within_tol<N>(make_tree_view(tree), *tree.permuted_data_view_,
                                   ray, tol);
}

} // namespace arp
} // namespace gaudi

#endif // GAUDI_ARP_RAY_TRAVERSAL_HPP
