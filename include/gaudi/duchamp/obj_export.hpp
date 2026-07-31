#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gaudi/duchamp/demo_trait.hpp"

namespace gaudi {
namespace duchamp {

struct named_mesh_snapshot {
  std::string name;
  const mesh_snapshot *mesh = nullptr;
};

inline vec3 snapshot_kd(const mesh_snapshot &mesh,
                        const vec3 &fallback = vec3(0.72, 0.74, 0.78)) {
  if (!mesh.colors.empty())
    return mesh.colors.front();
  return fallback;
}

// Writes path_stem.obj + path_stem.mtl. Returns true if at least one mesh was written.
inline bool write_mesh_snapshots_obj(
    const std::string &path_stem,
    const std::vector<named_mesh_snapshot> &meshes) {
  namespace fs = std::filesystem;

  const fs::path stem(path_stem);
  if (stem.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(stem.parent_path(), ec);
  }

  const std::string obj_path = path_stem + ".obj";
  const std::string mtl_path = path_stem + ".mtl";
  const std::string mtl_name = stem.filename().string() + ".mtl";

  std::ofstream obj(obj_path);
  std::ofstream mtl(mtl_path);
  if (!obj || !mtl) {
    std::cerr << "obj_export: failed to open " << obj_path << " / " << mtl_path
              << "\n";
    return false;
  }

  obj << "# gaudi duchamp mesh dump\n";
  obj << "mtllib " << mtl_name << "\n";
  mtl << "# gaudi duchamp materials\n";

  std::size_t v_offset = 0;
  int written = 0;
  for (const auto &entry : meshes) {
    if (!entry.mesh || entry.mesh->positions.empty() ||
        entry.mesh->indices.empty())
      continue;

    const mesh_snapshot &mesh = *entry.mesh;
    const std::string mat = entry.name + "_mat";
    const vec3 kd = snapshot_kd(mesh);

    mtl << "newmtl " << mat << "\n";
    mtl << "Kd " << kd.x() << " " << kd.y() << " " << kd.z() << "\n";
    mtl << "Ka " << kd.x() << " " << kd.y() << " " << kd.z() << "\n";
    mtl << "d 1.0\n\n";

    obj << "o " << entry.name << "\n";
    obj << "usemtl " << mat << "\n";
    for (const vec3 &p : mesh.positions)
      obj << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";

    for (std::size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
      const std::size_t i0 = v_offset + mesh.indices[i + 0] + 1;
      const std::size_t i1 = v_offset + mesh.indices[i + 1] + 1;
      const std::size_t i2 = v_offset + mesh.indices[i + 2] + 1;
      obj << "f " << i0 << " " << i1 << " " << i2 << "\n";
    }

    v_offset += mesh.positions.size();
    ++written;
  }

  if (written == 0) {
    std::cerr << "obj_export: no mesh geometry to write\n";
    return false;
  }

  std::cout << "obj_export: wrote " << obj_path << " (+ " << mtl_name << ")\n";
  return true;
}

inline bool write_scene_meshes_obj(
    const std::string &path_stem,
    const std::optional<mesh_snapshot> &shell,
    const std::optional<mesh_snapshot> &rod) {
  std::vector<named_mesh_snapshot> meshes;
  if (shell)
    meshes.push_back({"shell", &(*shell)});
  if (rod)
    meshes.push_back({"rod", &(*rod)});
  return write_mesh_snapshots_obj(path_stem, meshes);
}

} // namespace duchamp
} // namespace gaudi
