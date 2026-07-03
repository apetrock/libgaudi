# Field Wrappers TDD

## Overview

Type-safe, polymorphic field system for gaudi that allows:
- Different field types (bound, surface-following, sampled)
- Mesh type parameterization (shell, rod, etc.)
- Data type parameterization (vec3, vec14, quat, etc.)
- Primitive type parameterization (VERTEX, EDGE, FACE, CORNER)
- Interchangeable use via `field_base<...>::ptr`

> **API note (reconciled with the real `asawa` API).** Earlier drafts of this
> doc referenced an idealized mesh API (`create_cube()`, `create_data<T>()`,
> `prim_type::VERT`). The actual code is:
>
> | Earlier draft | Real `asawa` API |
> |---|---|
> | `asawa::shell::create_cube()` | `asawa::shell::load_cube()` (`shell/asset_loader.hpp`) |
> | `mesh->create_data<vec3>(prim_type::VERT)` | `asawa::init_vert_datum<vec3>(*mesh, v0)` returns `index_t`, or `datum_t<vec3>::create(prim_type::VERTEX, data)` + `mesh->insert_datum(d)` (`datums.hpp`) |
> | `prim_type::VERT / EDGE / FACE` | `enum prim_type { VERTEX, FACE, EDGE, CORNER };` (`datums.hpp`) |
> | `vec14` | `albers::vec14` (`albers/darboux_cyclide.hpp`) |
>
> **Dependency correction.** The field wrappers themselves wrap mesh-owned
> `asawa::datum_t<T>` storage and do **not** depend on `ext/liblombardi/`. The
> liblombardi datum pool is only needed to run fields through the *node graph*
> (see "Node Integration" below).

## Design

### 1. Base Interface

```cpp
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct field_base {
    using mesh_ptr = std::shared_ptr<MeshType>;
    using ptr = std::shared_ptr<field_base<MeshType, DataType, PT>>;

    mesh_ptr mesh;

    virtual ~field_base() = default;
    virtual const std::vector<DataType> &get() const = 0;
    virtual std::vector<DataType> &get() = 0;

    index_t size() const { return static_cast<index_t>(get().size()); }

    DataType operator[](index_t i) const { return get()[i]; }
    DataType &operator[](index_t i) { return get()[i]; }
};
```

### 2. Bound Fields

Mesh-owned data accessed by the datum index returned from `insert_datum`.
Changes with mesh topology (the mesh owns and resizes the storage).

```cpp
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct bound_field : public field_base<MeshType, DataType, PT> {
    index_t data_idx;

    bound_field(mesh_ptr m, index_t idx);

    const std::vector<DataType> &get() const override; // casts mesh->get_datum(idx)
    std::vector<DataType> &get() override;
};
```

Access pattern:
`std::static_pointer_cast<asawa::datum_t<DataType>>(mesh->get_datum(data_idx))->data()`.

**Use cases**: Positions (datum index 0 after `load_*`), velocities, forces.

### 3. Surface-Following Fields

Owns data, validates size against mesh primitive count. Tracks the mesh but is
independent of mesh data storage.

```cpp
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct surface_following_field : public field_base<MeshType, DataType, PT> {
    std::vector<DataType> data;

    surface_following_field(mesh_ptr m);

    bool valid() const;     // data.size() == prim_count(mesh, PT)
    void validate() const;  // throws std::runtime_error on mismatch

    const std::vector<DataType> &get() const override;
    std::vector<DataType> &get() override;

    void resize(index_t n);
    void reserve(index_t n);
};
```

`prim_count(mesh, PT)` maps `VERTEX/FACE/EDGE/CORNER` to
`vert_count()/face_count()/edge_count()/corner_count()`.

**Use cases**: Cyclide parameters (`albers::vec14`), curvature, per-vertex normals.

### 4. Sampled Fields

Completely independent data, optional mesh pointer.

```cpp
template <typename MeshType, typename DataType, asawa::prim_type PT>
struct sampled_field : public field_base<MeshType, DataType, PT> {
    std::vector<DataType> data;

    sampled_field(mesh_ptr m = nullptr);

    const std::vector<DataType> &get() const override;
    std::vector<DataType> &get() override;

    void resize(index_t n);
    void reserve(index_t n);
};
```

**Use cases**: Ridge points, intermediate results.

## Type Aliases

```cpp
using shell = asawa::shell::shell;

// Shell vertex fields
using shell_vert_positions  = bound_field<shell, vec3, asawa::prim_type::VERTEX>;
using shell_vert_velocities = bound_field<shell, vec3, asawa::prim_type::VERTEX>;
using shell_vert_normals    = surface_following_field<shell, vec3, asawa::prim_type::VERTEX>;
using shell_vert_cyclide_params =
    surface_following_field<shell, albers::vec14, asawa::prim_type::VERTEX>;
using shell_ridge_points    = sampled_field<shell, vec3, asawa::prim_type::VERTEX>;
```

## Polymorphic Algorithms

Algorithms accept `field_base<...>::ptr` and work with any concrete type:

```cpp
template <typename MeshType, typename DataType, asawa::prim_type PT>
DataType compute_average(
    typename field_base<MeshType, DataType, PT>::ptr field) {
    const auto &data = field->get();
    DataType sum = DataType::Zero();
    for (const auto &val : data) sum += val;
    return sum / static_cast<real>(std::max<size_t>(data.size(), 1));
}
```

## Node Integration (liblombardi bridge)

To run field data through the `liblombardi` graph, a thin `Datum` adapter bridges
the two systems, plus a concrete single-step node that validates the whole stack.

```cpp
// field_datum<T> : liblombardi::Datum  -- a pool datum holding std::vector<T>.
// march_node     : liblombardi::Node   -- ONE input port, ONE output port.
//   compute(): out[i] = in[i] + step * normal[i]
//   the normal field is a constructor parameter passed as field_base<...>::ptr,
//   so the node consumes field polymorphism while staying 1-in / 1-out.
```

This is the smallest node that exercises fields + the node system together and
"solidifies" both APIs: a position buffer in, a marched position buffer out.

## File Structure

```
include/gaudi/duchamp/
├── fields.hpp        # field_base, bound/surface/sampled, aliases, compute_average
└── field_nodes.hpp   # field_datum<T> + march_node (liblombardi bridge)

include/gaudi/test/
└── field_tests.hpp   # GAUDI_TEST cases (included by tests/gaudi_tests.cpp)
```

## Tests

All tests use the `GAUDI_TEST` framework (`gaudi/test/test.hpp`) and run inside
the `gaudi_tests` binary. `EXPECT_DEATH`-style checks are replaced with explicit
`valid()` assertions since the framework has no death tests.

### Test 1: Bound Field
```cpp
GAUDI_TEST(field_bound_access) {
    auto M = asawa::shell::load_cube();          // 8 verts, positions at datum 0
    auto xs = std::make_shared<shell_vert_positions>(M, 0);
    GAUDI_ASSERT(xs->size() == 8);
    (*xs)[0] = vec3(0, 0, 0);
    GAUDI_ASSERT(((*xs)[0] - vec3(0, 0, 0)).norm() < 1e-12);
}
```

### Test 2: Surface-Following Validation
```cpp
GAUDI_TEST(field_surface_following_validation) {
    auto M = asawa::shell::load_cube();           // 8 verts
    auto params = std::make_shared<shell_vert_normals>(M);
    params->resize(8);
    GAUDI_ASSERT(params->valid());
    params->resize(7);
    GAUDI_ASSERT(!params->valid());
}
```

### Test 3: Polymorphic Algorithm
```cpp
GAUDI_TEST(field_polymorphic_algorithm) {
    auto M = asawa::shell::load_cube();
    auto xs = std::make_shared<shell_vert_positions>(M, 0);
    auto ridges = std::make_shared<shell_ridge_points>(M);
    ridges->resize(3);
    (*ridges)[0] = vec3(1, 0, 0);
    (*ridges)[1] = vec3(0, 1, 0);
    (*ridges)[2] = vec3(0, 0, 1);
    vec3 avg_xs = compute_average<shell, vec3, asawa::prim_type::VERTEX>(xs);
    vec3 avg_r  = compute_average<shell, vec3, asawa::prim_type::VERTEX>(ridges);
    GAUDI_ASSERT((avg_r - vec3(1, 1, 1) / 3.0).norm() < 1e-12);
}
```

### Test 4: Sampled Field Independence
```cpp
GAUDI_TEST(field_sampled_independence) {
    auto ridges = std::make_shared<shell_ridge_points>(nullptr);
    ridges->resize(10);
    GAUDI_ASSERT(ridges->size() == 10);
    GAUDI_ASSERT(ridges->mesh == nullptr);
}
```

### Test 5: One-Level March Node
```cpp
GAUDI_TEST(field_march_node_one_step) {
    auto M = asawa::shell::load_sphere(1.0, 24, 16);
    const auto &x = asawa::const_get_vec_data(*M, 0);

    auto normals = std::make_shared<shell_vert_normals>(M);
    normals->get() = asawa::vertex_normals(*M, x);

    const real step = 0.1;
    liblombardi::GraphContext ctx;
    auto node = ctx.create_node<march_node>(normals, step);

    node->get_datum<march_node::InputPortDef>()->data() = x; // seed input
    ctx.run();
    const auto &out = node->get_datum<march_node::OutputPortDef>()->data();

    GAUDI_ASSERT(out.size() == x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        const vec3 expect = x[i] + step * normals->get()[i];
        GAUDI_ASSERT((out[i] - expect).norm() < 1e-10);
    }
}
```

## Implementation Checklist

### Prerequisites (done)
- [x] Datum pool system in `ext/liblombardi/`
- [x] Port system with enum-based ports + `PortRef`
- [x] Node base class with pool management

### Field System
- [ ] `include/gaudi/duchamp/fields.hpp`: `field_base`, `bound_field`,
      `surface_following_field`, `sampled_field`, `prim_count`, aliases,
      `compute_average`
- [ ] `include/gaudi/duchamp/field_nodes.hpp`: `field_datum<T>`, `march_node`
- [ ] `include/gaudi/test/field_tests.hpp`: Tests 1-5
- [ ] `#include` field_tests in `tests/gaudi_tests.cpp`
- [ ] Build + run `gaudi_tests`
