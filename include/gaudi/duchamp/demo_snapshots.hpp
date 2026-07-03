#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/demo_trait.hpp"
#include "gaudi/vec_addendum.h"

namespace gaudi {
namespace duchamp {

inline mesh_snapshot make_shell_mesh_snapshot(asawa::shell::shell &mesh,
                                              const vec3 &color) {
  mesh_snapshot snapshot;
  snapshot.positions = asawa::get_vec_data(mesh, 0);

  const std::vector<index_t> face_vert_ids = mesh.get_face_vert_ids(true);
  snapshot.indices.reserve(face_vert_ids.size());
  std::transform(face_vert_ids.begin(), face_vert_ids.end(),
                 std::back_inserter(snapshot.indices),
                 [](index_t idx) { return static_cast<uint32_t>(idx); });

  snapshot.colors.assign(snapshot.positions.size(), color);
  return snapshot;
}

inline mesh_snapshot make_rod_mesh_snapshot(const asawa::rod::rod &rod,
                                            const vec3 &color,
                                            int section_count = 32) {
  mesh_snapshot snapshot;
  if (rod.corner_count() == 0 || section_count < 3) {
    return snapshot;
  }

  const std::vector<vec3> &x = rod.x();
  const std::vector<quat> &u = rod.u();
  const int Nc = section_count;
  const real radius = 0.5 * rod._r;

  matX section(3, Nc);
  for (int i = 0; i < Nc; ++i) {
    const real theta = 2.0 * M_PI * real(i) / real(Nc);
    section.col(i) << radius * cos(theta), radius * sin(theta), 0.0;
  }

  snapshot.positions.reserve(static_cast<size_t>(rod.corner_count()) *
                             static_cast<size_t>(Nc));
  for (int i = 0; i < static_cast<int>(rod.corner_count()); ++i) {
    const vec3 center = x[i];
    const quat frame = i < static_cast<int>(u.size()) ? u[i] : quat::Identity();
    const matX section_world = frame.toRotationMatrix() * section;
    for (int j = 0; j < Nc; ++j) {
      snapshot.positions.push_back(center + section_world.col(j));
    }
  }

  snapshot.indices.reserve(static_cast<size_t>(rod.corner_count()) *
                           static_cast<size_t>(Nc) * 6);
  for (int i0 = 0; i0 < static_cast<int>(rod.corner_count()); ++i0) {
    const auto ci0 = asawa::rod::corner_id(i0);
    const auto ci1 = rod.next(ci0);
    if (ci1 < asawa::rod::corner_id(0)) {
      continue;
    }
    const int i1 = static_cast<int>(ci1);
    for (int j0 = 0; j0 < Nc; ++j0) {
      const int j1 = (j0 + 1) % Nc;
      const uint32_t v00 = static_cast<uint32_t>(Nc * i0 + j0);
      const uint32_t v01 = static_cast<uint32_t>(Nc * i0 + j1);
      const uint32_t v10 = static_cast<uint32_t>(Nc * i1 + j0);
      const uint32_t v11 = static_cast<uint32_t>(Nc * i1 + j1);
      snapshot.indices.push_back(v00);
      snapshot.indices.push_back(v01);
      snapshot.indices.push_back(v11);
      snapshot.indices.push_back(v00);
      snapshot.indices.push_back(v11);
      snapshot.indices.push_back(v10);
    }
  }

  snapshot.colors.assign(snapshot.positions.size(), color);
  return snapshot;
}

inline rod_snapshot make_rod_snapshot(const asawa::rod::rod &rod,
                                      const vec3 &color) {
  rod_snapshot snapshot;
  snapshot.positions = rod.x();
  snapshot.colors.assign(snapshot.positions.size(), color);
  return snapshot;
}

} // namespace duchamp
} // namespace gaudi
