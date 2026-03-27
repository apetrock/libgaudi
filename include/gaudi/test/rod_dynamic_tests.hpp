#ifndef __GAUDI_ROD_DYNAMIC_TESTS_HPP__
#define __GAUDI_ROD_DYNAMIC_TESTS_HPP__

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/test/test.hpp"
#include <array>
#include <vector>

namespace gaudi {
namespace test {

inline asawa::rod::dynamic::ptr make_simple_rod_dynamic() {
  std::vector<vec3> verts = {
      vec3(0.0, 0.0, 0.0),
      vec3(1.0, 0.0, 0.0),
      vec3(2.0, 0.0, 0.0),
  };
  auto rod = asawa::rod::rod::create(verts, false);
  return asawa::rod::dynamic::create(rod, 0.25, 2.0, 0.0);
}

inline asawa::rod::dynamic::ptr make_two_segment_collision_dynamic() {
  auto rod = asawa::rod::rod::create();
  rod->insert_strand({vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)}, false);
  rod->insert_strand({vec3(0.5, 0.1, -0.2), vec3(0.5, 0.1, 0.2)}, false);
  rod->set_radius(0.4);
  return asawa::rod::dynamic::create(rod, 0.25, 2.0, 0.0);
}

GAUDI_TEST(rod_dynamic_low_tol_edge_no_hit_is_sentinel) {
  auto dyn = make_simple_rod_dynamic();

  std::vector<vec3> query_verts = {
      vec3(100.0, 100.0, 100.0),
      vec3(101.0, 100.0, 100.0),
  };
  std::vector<index_t> query_edges = {0, 1};
  adjacency_view<std::vector<vec3>, std::vector<index_t>> query_view(
      query_verts, query_edges);

  auto collisions = dyn->get_collisions(query_view, 1e-3);
  GAUDI_ASSERT(collisions.size() == 1);
  GAUDI_EXPECT(collisions[0][0] == -1);
  GAUDI_EXPECT(collisions[0][1] == -1);
}

GAUDI_TEST(rod_dynamic_low_tol_point_no_hit_is_sentinel) {
  auto dyn = make_simple_rod_dynamic();

  std::vector<vec3> query_points = {
      vec3(-200.0, 50.0, -75.0),
  };
  std::vector<index_t> point_ids = {0};
  adjacency_view<std::vector<vec3>, std::vector<index_t>> point_view(
      query_points, point_ids);

  auto collisions = dyn->get_vert_collisions(point_view, 1e-3);
  GAUDI_ASSERT(collisions.size() == 1);
  GAUDI_EXPECT(collisions[0][0] == -1);
  GAUDI_EXPECT(collisions[0][1] == -1);
}

GAUDI_TEST(rod_dynamic_collision_ids_invalid_input_is_safe_sentinel) {
  auto dyn = make_simple_rod_dynamic();
  auto ids = dyn->get_collision_ids({-1, -1});
  GAUDI_EXPECT(ids[0] == -1);
  GAUDI_EXPECT(ids[1] == -1);
  GAUDI_EXPECT(ids[2] == -1);
  GAUDI_EXPECT(ids[3] == -1);
}

GAUDI_TEST(rod_dynamic_internal_collision_detects_close_two_segments) {
  auto dyn = make_two_segment_collision_dynamic();
  auto collisions = dyn->get_internal_collisions(1.0);

  bool found_valid = false;
  for (const auto &pair : collisions) {
    if (pair[0] < 0 || pair[1] < 0) {
      continue;
    }
    auto ids = dyn->get_collision_ids(pair);
    if (ids[0] >= 0 && ids[1] >= 0 && ids[2] >= 0 && ids[3] >= 0) {
      found_valid = true;
      break;
    }
  }
  GAUDI_EXPECT(found_valid);
}

} // namespace test
} // namespace gaudi

#endif // __GAUDI_ROD_DYNAMIC_TESTS_HPP__
