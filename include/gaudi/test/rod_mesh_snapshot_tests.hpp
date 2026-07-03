#ifndef __GAUDI_TEST_ROD_MESH_SNAPSHOT_TESTS_HPP__
#define __GAUDI_TEST_ROD_MESH_SNAPSHOT_TESTS_HPP__

#include <cmath>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/duchamp/demo_snapshots.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

namespace {

int count_rod_edges(const asawa::rod::rod &rod) {
  int edges = 0;
  for (int i = 0; i < static_cast<int>(rod.corner_count()); ++i) {
    const auto ci = asawa::rod::corner_id(i);
    if (rod.next(ci) >= asawa::rod::corner_id(0)) {
      ++edges;
    }
  }
  return edges;
}

} // namespace

GAUDI_TEST(rod_mesh_snapshot_loop_topology) {
  std::vector<vec3> points;
  const int corners = 8;
  const int section_count = 16;
  for (int i = 0; i < corners; ++i) {
    const real t = 2.0 * M_PI * real(i) / real(corners);
    points.emplace_back(cos(t), sin(t), 0.0);
  }

  auto rod = asawa::rod::rod::create(points, true);
  const auto snapshot =
      duchamp::make_rod_mesh_snapshot(*rod, vec3(1.0, 0.2, 0.5), section_count);

  GAUDI_ASSERT(snapshot.positions.size() ==
               static_cast<size_t>(corners * section_count));
  const int edges = count_rod_edges(*rod);
  GAUDI_ASSERT(snapshot.indices.size() == static_cast<size_t>(edges * section_count * 6));

  for (uint32_t index : snapshot.indices) {
    GAUDI_ASSERT(index < snapshot.positions.size());
  }
}

GAUDI_TEST(rod_mesh_snapshot_not_contiguous_polyline) {
  std::vector<vec3> points = {vec3(0, 0, 0), vec3(1, 0, 0), vec3(2, 0, 0),
                              vec3(3, 0, 0)};
  auto rod = asawa::rod::rod::create(points, false);
  const int section_count = 8;
  const auto snapshot =
      duchamp::make_rod_mesh_snapshot(*rod, vec3(0.2, 0.8, 1.0), section_count);

  GAUDI_ASSERT(snapshot.positions.size() ==
               static_cast<size_t>(rod->corner_count() * section_count));
  GAUDI_ASSERT(snapshot.indices.size() ==
               static_cast<size_t>(count_rod_edges(*rod) * section_count * 6));
  GAUDI_ASSERT(snapshot.indices.size() != static_cast<size_t>((points.size() - 1) * 2));
}

} // namespace test
} // namespace gaudi

#endif
