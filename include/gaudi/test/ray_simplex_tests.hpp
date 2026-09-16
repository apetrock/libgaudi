#ifndef GAUDI_TEST_RAY_SIMPLEX_TESTS_HPP
#define GAUDI_TEST_RAY_SIMPLEX_TESTS_HPP

#include "gaudi/arp/aabb.hpp"
#include "gaudi/arp/brute_force_tree.hpp"
#include "gaudi/arp/disk_triangle_hit.hpp"
#include "gaudi/arp/hash_tree.hpp"
#include "gaudi/arp/ray_simplex.hpp"
#include "gaudi/arp/ray_traversal.hpp"
#include "gaudi/arp/simplex_set.hpp"
#include "gaudi/test/test.hpp"
#include "gaudi/vec_addendum.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace gaudi {
namespace test {

GAUDI_TEST(disk_triangle_hit_2d_origin_inside) {
  const vec2 p0(0, 0);
  const vec2 p1(2, 0);
  const vec2 p2(0, 2);
  const arp::disk_triangle_hit_result hit =
      arp::disk_triangle_hit_2d(0.5, p0, p1, p2, 5.0, 5.0, 5.0, 0.0, 100.0);
  GAUDI_ASSERT(hit.hit);
  GAUDI_EXPECT(std::abs(hit.s - 4.5) < 1e-5);
  GAUDI_EXPECT(std::abs(hit.u) < 1e-5);
  GAUDI_EXPECT(std::abs(hit.w) < 1e-5);
}

GAUDI_TEST(disk_triangle_hit_2d_vertex_grazing) {
  const vec2 p0(1, 0);
  const vec2 p1(2, 0);
  const vec2 p2(1, 1);
  const arp::disk_triangle_hit_result hit =
      arp::disk_triangle_hit_2d(1.0, p0, p1, p2, 10.0, 10.0, 10.0, 0.0, 100.0);
  GAUDI_ASSERT(hit.hit);
  GAUDI_EXPECT(std::abs(hit.s - 10.0) < 1e-5);
  GAUDI_EXPECT(std::abs(hit.u - 1.0) < 1e-5);
  GAUDI_EXPECT(std::abs(hit.w) < 1e-5);
}

GAUDI_TEST(disk_triangle_hit_2d_constant_s_edge_interior) {
  const vec2 p0(1, 0);
  const vec2 p1(0, 1);
  const vec2 p2(3, 3);
  const real sr = std::sqrt(2.0);
  const arp::disk_triangle_hit_result hit =
      arp::disk_triangle_hit_2d(sr, p0, p1, p2, 5.0, 5.0, 5.0, 0.0, 100.0);
  GAUDI_ASSERT(hit.hit);
  GAUDI_EXPECT(std::abs(hit.u - 0.5) < 1e-4);
  GAUDI_EXPECT(std::abs(hit.w - 0.5) < 1e-4);
  GAUDI_EXPECT(hit.s < 4.0);
}

GAUDI_TEST(disk_triangle_hit_2d_miss) {
  const vec2 p0(5, 0);
  const vec2 p1(6, 0);
  const vec2 p2(5, 1);
  const arp::disk_triangle_hit_result hit =
      arp::disk_triangle_hit_2d(0.5, p0, p1, p2, 10.0, 10.0, 10.0, 0.0, 100.0);
  GAUDI_EXPECT(!hit.hit);
}

GAUDI_TEST(ray_slop_bbox_thin_hit) {
  const arp::Ray ray(vec3(0, 0, -1), vec3(0, 0, 1), 0, 100);
  ext::extents_t box;
  box[0] = vec3(-1, -1, 0);
  box[1] = vec3(1, 1, 1);
  GAUDI_EXPECT(arp::ray_slop_intersects_bbox(ray, 0.0, box));
}

GAUDI_TEST(ray_slop_bbox_slop_only_hit) {
  const arp::Ray ray(vec3(2, 0, 0.5), vec3(0, 0, 1), 0, 100);
  ext::extents_t box;
  box[0] = vec3(-1, -1, 0);
  box[1] = vec3(1, 1, 1);
  GAUDI_EXPECT(!arp::ray_aabb_overlaps(ray, box));
  GAUDI_EXPECT(arp::ray_slop_intersects_bbox(ray, 1.5, box));
}

GAUDI_TEST(ray_slop_bbox_parallel_miss) {
  const arp::Ray ray(vec3(0, 5, 0.5), vec3(1, 0, 0), 0, 100);
  ext::extents_t box;
  box[0] = vec3(-1, -1, 0);
  box[1] = vec3(1, 1, 1);
  GAUDI_EXPECT(!arp::ray_aabb_overlaps(ray, box));
  GAUDI_EXPECT(!arp::ray_slop_intersects_bbox(ray, 0.5, box));
}

GAUDI_TEST(ray_pierce_triangle_analytic) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(1, 0, 0);
  const vec3 v2(0, 1, 0);
  const arp::Ray ray(vec3(0.2, 0.2, 1), vec3(0, 0, -1));

  const auto hit = arp::ray_pierce_triangle(ray, v0, v1, v2);
  GAUDI_ASSERT(hit.has_value());
  GAUDI_EXPECT(std::abs(hit->t_ray - 1.0) < 1e-5);
  GAUDI_EXPECT(std::abs(hit->local.u - 0.2) < 1e-4);
  GAUDI_EXPECT(std::abs(hit->local.v - 0.2) < 1e-4);
  GAUDI_EXPECT(std::abs(hit->point.z()) < 1e-5);
}

GAUDI_TEST(ray_pierce_triangle_miss_parallel) {
  const arp::Ray ray(vec3(0, 0, 1), vec3(1, 0, 0));
  const auto hit = arp::ray_pierce_triangle(ray, vec3(0, 0, 0), vec3(1, 0, 0),
                                          vec3(0, 1, 0));
  GAUDI_EXPECT(!hit.has_value());
}

GAUDI_TEST(ray_point_within_tol_collect) {
  const arp::Ray ray(vec3(0, 0, 0), vec3(0, 0, 1));
  const auto hit = arp::ray_point_within_tol(ray, vec3(0.1, 0, 0.5), 0.2, 7);
  GAUDI_ASSERT(hit.has_value());
  GAUDI_EXPECT(hit->id == 7);
  GAUDI_EXPECT(std::abs(hit->t_ray - 0.5) < 1e-6);
}

namespace {

template <typename TreeA, typename TreeB>
void expect_same_pierce(const TreeA &a, const TreeB &b, const arp::Ray &ray,
                        real eps) {
  const auto hit_a = arp::ray_pierce_tree(a, ray);
  const auto hit_b = arp::ray_pierce_tree(b, ray);
  GAUDI_ASSERT(hit_a.has_value());
  GAUDI_ASSERT(hit_b.has_value());
  GAUDI_EXPECT(hit_a->id == hit_b->id);
  GAUDI_EXPECT(std::abs(hit_a->t_ray - hit_b->t_ray) < eps);
  if constexpr (TreeA::kSimplexN == 3) {
    GAUDI_EXPECT(std::abs(hit_a->local.u - hit_b->local.u) < eps);
    GAUDI_EXPECT(std::abs(hit_a->local.v - hit_b->local.v) < eps);
  }
}

inline void sort_slop_hits(std::vector<arp::ray_hit<3>> &hits) {
  std::sort(hits.begin(), hits.end(),
            [](const arp::ray_hit<3> &a, const arp::ray_hit<3> &b) {
              if (a.id != b.id) {
                return a.id < b.id;
              }
              return a.t_ray < b.t_ray;
            });
}

inline void expect_same_slop_hits(const std::vector<arp::ray_hit<3>> &a_in,
                                  const std::vector<arp::ray_hit<3>> &b_in,
                                  real eps) {
  std::vector<arp::ray_hit<3>> a = a_in;
  std::vector<arp::ray_hit<3>> b = b_in;
  sort_slop_hits(a);
  sort_slop_hits(b);
  GAUDI_EXPECT(a.size() == b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    GAUDI_EXPECT(a[i].id == b[i].id);
    GAUDI_EXPECT(std::abs(a[i].t_ray - b[i].t_ray) < eps);
  }
}

inline void expect_slop_contact(const arp::Ray &ray, real slop, const vec3 &v0,
                                const vec3 &v1, const vec3 &v2,
                                real slop_tol = 1e-3) {
  const auto hit =
      arp::ray_slop_intersects_triangle(ray, slop, v0, v1, v2, 0);
  GAUDI_ASSERT(hit.has_value());
  const real ray_gap = arp::ray_distance_to_point(ray, hit->point);
  GAUDI_EXPECT(ray_gap <= slop + slop_tol);
  vec3 on_tri;
  const std::array<vec3, 3> tri = {v0, v1, v2};
  const real claim_gap = va::distance_from_triangle(tri, hit->point, on_tri);
  GAUDI_EXPECT(claim_gap < 1e-2);
}

} // namespace

GAUDI_TEST(ray_pierce_tri_backend_equivalence) {
  std::vector<vec3> verts = {vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0)};
  std::vector<index_t> adj = {0, 1, 2};
  arp::simplex_set<3> tris(verts, adj);

  auto bf = arp::brute_force_tree<3>::create(tris);
  auto morton = arp::bvh_tree<3>::create(tris);
  auto legacy = arp::aabb_tree<3>::create(tris);

  const arp::Ray ray(vec3(0.25, 0.25, 2), vec3(0, 0, -1));

  expect_same_pierce(*bf, *morton, ray, 1e-5);
  expect_same_pierce(*bf, *legacy, ray, 1e-5);
  expect_same_pierce(*morton, *legacy, ray, 1e-5);

  const auto m_hit = arp::ray_pierce_tree(*morton, ray);
  GAUDI_ASSERT(m_hit.has_value());
  GAUDI_EXPECT(m_hit->id == 0);
}

GAUDI_TEST(ray_collect_points_backend_equivalence) {
  std::vector<vec3> verts = {
      vec3(0.05, 0, 0), vec3(-0.05, 0, 0), vec3(0, 0.05, 1),
      vec3(0, -0.05, 1), vec3(5, 5, 5)};
  std::vector<index_t> adj = {0, 1, 2, 3, 4};
  arp::simplex_set<1> pts(verts, adj);

  auto bf = arp::brute_force_tree<1>::create(pts);
  auto morton = arp::bvh_tree<1>::create(pts);

  const arp::Ray ray(vec3(0, 0, 0), vec3(0, 0, 1));
  const real tol = 0.08;

  const auto bf_hits = bf->ray_collect(ray, tol);
  const auto bv_hits = arp::ray_collect_tree(*morton, ray, tol);

  GAUDI_EXPECT(bf_hits.size() == bv_hits.size());
  GAUDI_EXPECT(bf_hits.size() >= 2);
}

GAUDI_TEST(ray_slop_pierce_triangle_grazing) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(1, 0, 0);
  const vec3 v2(0, 1, 0);
  const arp::Ray ray(vec3(0.2, 0.2, 1), vec3(0, 0, -1));

  const auto slop =
      arp::ray_slop_intersects_triangle(ray, 0.05, v0, v1, v2, 0);
  GAUDI_ASSERT(slop.has_value());
  GAUDI_EXPECT(std::abs(slop->t_ray - 1.0) < 0.1);
  GAUDI_EXPECT(slop->point.z() < 0.05);
  expect_slop_contact(ray, 0.05, v0, v1, v2);
}

GAUDI_TEST(ray_slop_pierce_triangle_large_tmax) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(1, 0, 0);
  const vec3 v2(0, 1, 0);
  const arp::Ray ray(vec3(0.2, 0.2, 1), vec3(0, 0, -1), 0, 1e6);

  const auto slop =
      arp::ray_slop_intersects_triangle(ray, 0.05, v0, v1, v2, 0);
  GAUDI_ASSERT(slop.has_value());
  GAUDI_EXPECT(std::abs(slop->t_ray - 1.0) < 0.1);
}

GAUDI_TEST(ray_slop_pierce_triangle_cylinder_entry) {
  const vec3 v0(0, 0, 5);
  const vec3 v1(1, 0, 5);
  const vec3 v2(0, 0, 10);
  const arp::Ray ray(vec3(0.5, 0, 15), vec3(0, 0, -1), 0, 100);

  const auto slop =
      arp::ray_slop_intersects_triangle(ray, 0.5, v0, v1, v2, 0);
  GAUDI_ASSERT(slop.has_value());
  GAUDI_EXPECT(slop->t_ray < 6.0);
  GAUDI_EXPECT(slop->point.z() > 9.0);
  expect_slop_contact(ray, 0.5, v0, v1, v2);
}

GAUDI_TEST(ray_slop_parallel_plane_hit) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(100, 0, 0);
  const vec3 v2(0, 100, 0);
  const arp::Ray ray(vec3(-50, 0.02, 0), vec3(1, 0, 0), 0, 200);
  const real slop = 0.05;

  GAUDI_EXPECT(!arp::ray_pierce_triangle(ray, v0, v1, v2).has_value());
  const auto hit = arp::ray_slop_intersects_triangle(ray, slop, v0, v1, v2);
  GAUDI_ASSERT(hit.has_value());
  expect_slop_contact(ray, slop, v0, v1, v2);
  GAUDI_EXPECT(std::abs(hit->point.z()) < 0.05);
}

GAUDI_TEST(ray_slop_null_basis_stable) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(1, 0, 0);
  const vec3 v2(0, 1, 0);
  const real slop = 0.05;
  const std::vector<vec3> dirs = {vec3(0, 0, 1), vec3(0, 0, -1), vec3(1, 0, 0),
                                  vec3(0, 1, 0), vec3(1, 1, 1).normalized()};
  for (const vec3 &d : dirs) {
    const arp::Ray ray(vec3(0.2, 0.2, 1), d, 0, 100);
    if (!arp::ray_slop_intersects_triangle(ray, slop, v0, v1, v2).has_value()) {
      continue;
    }
    expect_slop_contact(ray, slop, v0, v1, v2);
  }
}

GAUDI_TEST(ray_slop_pierce_triangle_matches_thin_at_zero_slop) {
  const vec3 v0(0, 0, 0);
  const vec3 v1(1, 0, 0);
  const vec3 v2(0, 1, 0);
  const arp::Ray ray(vec3(0.2, 0.2, 1), vec3(0, 0, -1));

  const auto thin = arp::ray_pierce_triangle(ray, v0, v1, v2);
  const auto slop = arp::ray_slop_intersects_triangle(ray, 0.0, v0, v1, v2);
  GAUDI_ASSERT(thin.has_value());
  GAUDI_ASSERT(slop.has_value());
  GAUDI_EXPECT(std::abs(thin->t_ray - slop->t_ray) < 1e-5);
}

GAUDI_TEST(ray_slop_collect_tri_backend_equivalence) {
  const std::vector<vec3> verts = {
      vec3(0, 0, 9),  vec3(1, 0, 9),  vec3(0, 1, 9),  vec3(0, 0, 1),
      vec3(1, 0, 1),  vec3(0, 1, 1),  vec3(2, 0, 5),  vec3(3, 0, 5),
      vec3(2, 1, 5),  vec3(2, 0, 0), vec3(3, 0, 0), vec3(2, 1, 0),
  };
  const std::vector<index_t> faces = {0,  1,  2,  3,  4,  5,  6,  7,  8,
                                      9, 10, 11, 6,  7,  8,  3,  4,  5};
  arp::simplex_set<3> tris(verts, faces);
  auto bf = arp::brute_force_tree<3>::create(tris);
  auto bvh = arp::bvh_tree<3>::create(tris);

  const std::vector<arp::Ray> rays = {
      arp::Ray(vec3(0.2, 0.2, 15), vec3(0, 0, -1), 0, 100),
      arp::Ray(vec3(0.5, 0.1, 12), vec3(-0.05, 0, -1), 0, 100),
      arp::Ray(vec3(2.4, 0.5, 8), vec3(0, 0, -1), 0, 50),
      arp::Ray(vec3(1.0, 0.5, 6), vec3(0.2, 0.1, -1), 0, 80),
  };
  const real slop = 0.35;

  for (const arp::Ray &ray : rays) {
    const auto bf_hits = bf->ray_slop_pierce(ray, slop);
    const auto bvh_hits = arp::ray_slop_pierce_tree(*bvh, ray, slop);
    expect_same_slop_hits(bf_hits, bvh_hits, 1e-3);
  }
}

GAUDI_TEST(ray_slop_pierce_min_heap_order) {
  const vec3 v_near0(0, 0, 9);
  const vec3 v_near1(1, 0, 9);
  const vec3 v_near2(0, 1, 9);
  const vec3 v_far0(0, 0, 1);
  const vec3 v_far1(1, 0, 1);
  const vec3 v_far2(0, 1, 1);
  const std::vector<vec3> verts = {v_near0, v_near1, v_near2,
                                   v_far0,  v_far1,  v_far2};
  const std::vector<index_t> faces = {0, 1, 2, 3, 4, 5};
  arp::simplex_set<3> tris(verts, faces);
  auto tree = arp::bvh_tree<3>::create(tris);

  const arp::Ray ray(vec3(0.2, 0.2, 15), vec3(0, 0, -1), 0, 100);
  auto hits = arp::ray_slop_pierce_tree(*tree, ray, 0.5);
  GAUDI_EXPECT(hits.size() >= 2);

  arp::ray_hits_make_min_heap(hits, &ray);
  const arp::ray_hit<3> first = arp::ray_hits_pop_min(hits, &ray);
  const arp::ray_hit<3> second = arp::ray_hits_pop_min(hits, &ray);
  GAUDI_EXPECT(first.t_ray < second.t_ray);
  GAUDI_EXPECT(first.point.z() > second.point.z());
}

} // namespace test
} // namespace gaudi

#endif // GAUDI_TEST_RAY_SIMPLEX_TESTS_HPP
