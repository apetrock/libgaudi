#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/demo_trait.hpp"

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

inline rod_snapshot make_rod_snapshot(const asawa::rod::rod &rod,
                                      const vec3 &color) {
  rod_snapshot snapshot;
  snapshot.positions = rod.x();
  snapshot.colors.assign(snapshot.positions.size(), color);
  return snapshot;
}

} // namespace duchamp
} // namespace gaudi
