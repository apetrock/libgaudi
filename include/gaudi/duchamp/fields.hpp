#ifndef __GAUDI_DUCHAMP_FIELDS__
#define __GAUDI_DUCHAMP_FIELDS__

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "gaudi/common.h"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/shell.hpp"

namespace gaudi {
namespace duchamp {

// Number of primitives of a given kind on a mesh. Used by surface-following
// fields to validate that their owned storage tracks the mesh.
template <typename MeshType>
inline index_t prim_count(const MeshType &M, asawa::prim_type pt) {
  switch (pt) {
  case asawa::prim_type::VERTEX:
    return static_cast<index_t>(M.vert_count());
  case asawa::prim_type::FACE:
    return static_cast<index_t>(M.face_count());
  case asawa::prim_type::EDGE:
    return static_cast<index_t>(M.edge_count());
  case asawa::prim_type::CORNER:
    return static_cast<index_t>(M.corner_count());
  }
  return 0;
}

// Polymorphic field interface. Concrete fields differ only in where the
// std::vector<DataType> lives (mesh-owned vs. field-owned).
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct field_base {
  using mesh_ptr = std::shared_ptr<MeshType>;
  using ptr = std::shared_ptr<field_base<MeshType, DataType, PT>>;

  mesh_ptr mesh;

  field_base(mesh_ptr m = nullptr) : mesh(std::move(m)) {}
  virtual ~field_base() = default;

  virtual const std::vector<DataType> &get() const = 0;
  virtual std::vector<DataType> &get() = 0;

  index_t size() const { return static_cast<index_t>(get().size()); }

  const DataType &operator[](index_t i) const { return get()[i]; }
  DataType &operator[](index_t i) { return get()[i]; }
};

// Mesh-owned data, accessed by the datum index returned from insert_datum.
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct bound_field : public field_base<MeshType, DataType, PT> {
  using base = field_base<MeshType, DataType, PT>;
  using mesh_ptr = typename base::mesh_ptr;

  index_t data_idx;

  bound_field(mesh_ptr m, index_t idx) : base(std::move(m)), data_idx(idx) {}

  std::vector<DataType> &get() override {
    return std::static_pointer_cast<asawa::datum_t<DataType>>(
               this->mesh->get_datum(data_idx))
        ->data();
  }

  const std::vector<DataType> &get() const override {
    return std::static_pointer_cast<const asawa::datum_t<DataType>>(
               this->mesh->const_get_datum(data_idx))
        ->data();
  }
};

// Field-owned data that should track the mesh primitive count.
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct surface_following_field : public field_base<MeshType, DataType, PT> {
  using base = field_base<MeshType, DataType, PT>;
  using mesh_ptr = typename base::mesh_ptr;

  std::vector<DataType> data;

  surface_following_field(mesh_ptr m) : base(std::move(m)) {}

  bool valid() const {
    if (!this->mesh) {
      return false;
    }
    return static_cast<index_t>(data.size()) == prim_count(*this->mesh, PT);
  }

  void validate() const {
    if (!valid()) {
      throw std::runtime_error(
          "surface_following_field: size does not match mesh primitive count");
    }
  }

  std::vector<DataType> &get() override { return data; }
  const std::vector<DataType> &get() const override { return data; }

  void resize(index_t n) { data.resize(n); }
  void reserve(index_t n) { data.reserve(n); }
};

// Field-owned data, independent of the mesh (mesh pointer optional).
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct sampled_field : public field_base<MeshType, DataType, PT> {
  using base = field_base<MeshType, DataType, PT>;
  using mesh_ptr = typename base::mesh_ptr;

  std::vector<DataType> data;

  sampled_field(mesh_ptr m = nullptr) : base(std::move(m)) {}

  std::vector<DataType> &get() override { return data; }
  const std::vector<DataType> &get() const override { return data; }

  void resize(index_t n) { data.resize(n); }
  void reserve(index_t n) { data.reserve(n); }
};

// ---- type aliases ---------------------------------------------------------

using shell_mesh = asawa::shell::shell;

using shell_vert_positions =
    bound_field<shell_mesh, vec3, asawa::prim_type::VERTEX>;
using shell_vert_velocities =
    bound_field<shell_mesh, vec3, asawa::prim_type::VERTEX>;
using shell_vert_normals =
    surface_following_field<shell_mesh, vec3, asawa::prim_type::VERTEX>;
using shell_ridge_points =
    sampled_field<shell_mesh, vec3, asawa::prim_type::VERTEX>;

// ---- polymorphic algorithms ----------------------------------------------

template <typename MeshType, typename DataType, asawa::prim_type PT>
DataType compute_average(
    typename field_base<MeshType, DataType, PT>::ptr field) {
  const auto &data = field->get();
  DataType sum = DataType::Zero();
  for (const auto &val : data) {
    sum += val;
  }
  return sum / static_cast<real>(std::max<size_t>(data.size(), size_t(1)));
}

} // namespace duchamp
} // namespace gaudi

#endif // __GAUDI_DUCHAMP_FIELDS__
