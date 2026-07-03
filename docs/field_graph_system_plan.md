# Field Graph System Implementation Plan

## Overview
Create a type-safe, pool-based field graph system with automatic data flow, similar to the render graph system in libvermeer but for simulation fields instead of textures.

---

## Phase 1: Base Infrastructure

### 1.1 Create `lombardi/field_port.hpp`

**Purpose:** Define generic port and handle types for field data

**Content:**
```cpp
#pragma once
#include <map>
#include <vector>
#include <memory>

namespace lombardi {

// Generic port handle with stable ID
template <typename Tag>
struct PortHandle {
  size_t id = 0;

  bool is_valid() const { return id != 0; }
};

// Field input port
template <typename Tag>
struct FieldInput {
  PortHandle<Tag> handle;
  bool connected = false;

  void connect(const PortHandle<Tag>& h) {
    handle = h;
    connected = true;
  }
};

// Field output port
template <typename Tag>
struct FieldOutput {
  PortHandle<Tag> handle;

  void connect(FieldInput<Tag>* input) {
    // Add input to dependency list
  }
};

// Wire function for connecting outputs to inputs
template <typename OutTag, typename InTag>
inline void wire(PortHandle<OutTag>& out, FieldInput<InTag>& in) {
  static_assert(std::is_same_v<OutTag, InTag>, "Only same type ports can be wired");
  in.connect(&out);
}

} // namespace lombardi
```

---

### 1.2 Create `liblombardi/include/liblombardi/slab.hpp`

**Purpose:** Slab allocator for efficient field memory management

**Content:**
```cpp
#pragma once
#include <vector>
#include <memory>
#include <cstdint>

namespace liblombardi {

// Slab allocator - allocates memory in fixed-size blocks
template <typename T>
class SlabAllocator {
public:
  SlabAllocator(size_t block_size = 64) : _block_size(block_size) {}

  // Allocate a single element (O(1))
  size_t allocate() {
    // Try to find a free slot in existing blocks
    for (auto& block : _blocks) {
      if (!_free_slots.empty()) {
        size_t idx = _free_slots.back();
        _free_slots.pop_back();
        return block.size() * _block_size + idx;
      }
    }

    // Create new block if needed
    if (_blocks.empty() || _blocks.back().size() >= _block_size) {
      _blocks.emplace_back(_block_size);
    }

    size_t idx = _blocks.back().size();
    _blocks.back().push_back(T{});
    return _blocks.size() * _block_size + idx;
  }

  // Deallocate a single element (O(1))
  void deallocate(size_t index) {
    size_t block_id = index / _block_size;
    size_t local_idx = index % _block_size;

    if (block_id < _blocks.size()) {
      _blocks[block_id][local_idx].~T();
      _free_slots.push_back(local_idx);
    }
  }

  // Get raw pointer to data (for direct access)
  T* data(size_t index = 0) {
    size_t block_id = index / _block_size;
    size_t local_idx = index % _block_size;
    return &_blocks[block_id][local_idx];
  }

  // Get block size
  size_t block_size() const { return _block_size; }

  // Get number of blocks
  size_t num_blocks() const { return _blocks.size(); }

  // Clear all memory
  void clear() {
    _blocks.clear();
    _free_slots.clear();
  }

private:
  std::vector<std::vector<T>> _blocks;   // Multiple blocks
  std::vector<size_t> _free_slots;       // Free slots per block
  size_t _block_size;                    // Elements per block
};

} // namespace liblombardi
```

---

### 1.3 Create `liblombardi/include/liblombardi/field_pool.hpp`

**Purpose:** Field pool using slab allocators for different data types

**Content:**
```cpp
#pragma once
#include "slab.hpp"
#include <map>
#include <typeindex>

namespace liblombardi {

// Field pool using slab allocators
class FieldPool {
public:
  // Allocate memory for a field of type T
  template <typename T>
  size_t allocate() {
    auto& slab = get_slab<T>();
    return slab.allocate();
  }

  // Deallocate memory for a field of type T
  template <typename T>
  void deallocate(size_t index) {
    auto& slab = get_slab<T>();
    slab.deallocate(index);
  }

  // Get a pointer to field data (for direct access)
  template <typename T>
  T* get(size_t index) {
    return get_slab<T>().data(index);
  }

  // Get raw pointer to data
  template <typename T>
  T* data() {
    return get_slab<T>().data(0);
  }

  // Get the slab for a type (for iteration)
  template <typename T>
  SlabAllocator<T>& get_slab() {
    static std::type_index type_key(typeid(T));
    auto it = _slabs.find(type_key);
    if (it == _slabs.end()) {
      it = _slabs.emplace(type_key, SlabAllocator<T>(64)).first;
    }
    return it->second;
  }

  // Get the slab for a type (const)
  template <typename T>
  const SlabAllocator<T>& get_slab() const {
    return const_cast<FieldPool*>(this)->get_slab<T>();
  }

  // Clear all allocations
  void clear() {
    for (auto& [type, slab] : _slabs) {
      slab.clear();
    }
  }

private:
  std::map<std::type_index, std::any> _slabs;  // Slabs by type
};

} // namespace liblombardi
```

---

### 1.4 Create `liblombardi/include/liblombardi/field_node_base.hpp`

**Purpose:** Base node class with pool management

**Content:**
```cpp
#pragma once
#include <memory>
#include <vector>
#include "field_port.hpp"
#include "field_pool.hpp"

namespace liblombardi {

// Base class for all field nodes
template <typename MeshType>
class FieldNode {
public:
  using mesh_ptr = std::shared_ptr<MeshType>;

  FieldNode(mesh_ptr m) : _mesh(m) {}

  virtual ~FieldNode() = default;

  // Pool management
  void set_pool(FieldPool& pool) { _pool = &pool; }

  // Resource usage reporting
  virtual std::vector<ResourceUse> resource_uses() const { return _resource_uses; }

protected:
  mesh_ptr _mesh;
  FieldPool* _pool = nullptr;
  std::vector<ResourceUse> _resource_uses;

  // Helper to register resource usage
  void add_resource_use(ResourceUseKind kind, const PortHandle<Tag>& handle) {
    ResourceUse use{kind, handle.id};
    _resource_uses.push_back(use);
  }

  // Template for resolving from pool
  template <typename Tag>
  auto resolve(const PortHandle<Tag>& handle) -> decltype(pool_get<Tag>()) {
    if (!_pool || !handle.is_valid()) return {};
    return pool_get<Tag>((*_pool)[handle.id]);
  }
};

// Resource use kind
enum class ResourceUseKind { read, write, read_write };

struct ResourceUse {
  ResourceUseKind kind;
  size_t id;
};

} // namespace liblombardi
```

---

## Phase 2: Field Input/Output System

### 2.1 Extend `lombardi/field_port.hpp`

**Additions:**
- `FieldInput<Tag>` wrapper with handle
- `FieldOutput<Tag>` wrapper with handle
- `wire()` template function
- `resolve()` helper function

### 2.2 Create `liblombardi/include/liblombardi/field_binder.hpp`

**Purpose:** Helper functions for binding fields to pools and creating outputs

**Content:**
```cpp
#pragma once
#include "field_port.hpp"
#include "field_pool.hpp"

namespace liblombardi {

// Bind a field output to a pool slot
template <typename Tag>
inline void bind_field_output(FieldOutput<Tag>& output, FieldPool& pool) {
  output.handle = pool.allocate<Tag>();
}

// Create a field from mesh data
template <typename MeshType, typename DataType, prim_type PT>
FieldOutput<DataType> create_field_output(MeshType& mesh, const std::string& name) {
  FieldOutput<DataType> output;
  output.handle = field_pool.allocate<DataType>();
  return output;
}

} // namespace liblombardi
```

---

## Phase 3: Field Nodes

### 3.1 Create `lombardi/field_nodes.hpp`

**Purpose:** Define field node base class with enum ports

**Content:**
```cpp
#pragma once
#include <map>
#include "field_node_base.hpp"
#include "field_port.hpp"

namespace lombardi {

// Base field node with enum ports
template <typename MeshType, typename DataType, prim_type PT>
class FieldNodeBase : public FieldNode<MeshType> {
public:
  using field_ptr = std::shared_ptr<gaudi::duchamp::field_base<MeshType, DataType, PT>>;

  enum class Input { position, velocity, force, other };
  enum class Output { field };

  FieldNodeBase(mesh_ptr m, field_ptr field)
    : FieldNode<MeshType>(m), _field(field) {}

  // Set input field
  void set_input(Input port, field_ptr field) {
    _inputs[port] = field;
  }

  // Get output handle
  PortHandle<DataType> output_handle() const {
    return _output.handle;
  }

  // Get the field data
  const std::vector<DataType>& get_field() const {
    return _field->get();
  }

protected:
  field_ptr _field;
  std::map<Input, field_ptr> _inputs;
  FieldOutput<DataType> _output;
};

} // namespace lombardi
```

### 3.2 Implement field computation nodes

**Create:** `lombardi/computed_field_nodes.hpp`

**Nodes:**
- `PassthroughFieldNode` - just passes through input data
- `RidgeFieldNode` - computes ridge lines from positions
- `CurvatureFieldNode` - computes curvature values
- `DivergenceFieldNode` - computes divergence
- Other computation nodes as needed

```cpp
#pragma once
#include "field_nodes.hpp"

namespace lombardi {

// Passthrough node - doesn't modify data
template <typename MeshType, typename DataType, prim_type PT>
class PassthroughFieldNode : public FieldNodeBase<MeshType, DataType, PT> {
public:
  using base = FieldNodeBase<MeshType, DataType, PT>;
  using base::base;

protected:
  void compute_field() override {
    // Just return the input field
  }
};

// Ridge field node - computes ridge lines
template <typename MeshType>
class RidgeFieldNode : public FieldNodeBase<MeshType, vec3, prim_type::VERT> {
public:
  using base = FieldNodeBase<MeshType, vec3, prim_type::VERT>;
  using base::base;

protected:
  void compute_field() override {
    auto& input = base::_inputs[Input::position];
    auto& output = base::_field;

    for (size_t i = 0; i < input->size(); i++) {
      output->get()[i] = compute_ridge(input->get()[i]);
    }
  }

  vec3 compute_ridge(const vec3& pos) {
    // Ridge computation logic
    return pos;
  }
};

} // namespace lombardi
```

---

## Phase 4: Integration with Solver

### 4.1 Modify `solver_node.hpp`

**Changes:**
- Make solver nodes inherit from `FieldNode`
- Connect solver blocks to field ports
- Automatic data flow from fields to solver

### 4.2 Create field bundles

**Create:** `lombardi/field_bundles.hpp`

**Functions:**
```cpp
// Create a bundle of field-based constraints
inline FieldBundle make_field_bundle(FieldNodeBase<...>& node, real dt, real damping) {
  return {
    [&node, dt, damping](Solver& solver) {
      auto constraint = make_field_constraint(node, dt, damping);
      solver.add_constraint(constraint);
    }
  };
}
```

---

## Phase 5: Tests

### 5.1 Create `tests/field_graph_test.cpp`

**Test cases:**
1. Pool-based resolution
2. Port wiring
3. Field computation nodes
4. Solver integration

```cpp
TEST(FieldGraphTest, PoolResolution) {
  FieldPool pool;
  auto pos_handle = pool.allocate<vec3>();
  auto node = create_passthrough_node();

  // Store handle (stable)
  node->set_output_handle(pos_handle);

  // Resolve each frame (may change)
  auto pos = node->resolve<vec3>(pos_handle);
  // Each frame, pos gets the current value from pool
}

TEST(FieldGraphTest, PortWiring) {
  auto out_node = create_output_node();
  auto in_node = create_input_node();

  wire(out_node->output_handle(), in_node->input_handle());
  // Connections made at compile time
}

TEST(FieldGraphTest, FieldComputation) {
  auto mesh = create_test_mesh();
  auto pos_field = create_field(mesh);
  auto ridge_node = std::make_shared<RidgeFieldNode>(mesh, pos_field);

  ridge_node->set_input(Input::position, pos_field);
  ridge_node->compute();

  // Field automatically updated
}
```

---

## Implementation Checklist

### Phase 1
- [ ] Create `lombardi/field_port.hpp`
- [ ] Create `lombardi/field_node_base.hpp`
- [ ] Test port handles and wire function

### Phase 2
- [ ] Extend `field_port.hpp` with FieldInput/FieldOutput
- [ ] Create `lombardi/field_binder.hpp`
- [ ] Test pool binding and allocation

### Phase 3
- [ ] Create `lombardi/field_nodes.hpp`
- [ ] Implement `PassthroughFieldNode`
- [ ] Implement `RidgeFieldNode`
- [ ] Implement other computed field nodes

### Phase 4
- [ ] Modify `solver_node.hpp` to use field nodes
- [ ] Create `lombardi/field_bundles.hpp`
- [ ] Test solver integration

### Phase 5
- [ ] Create `tests/field_graph_test.cpp`
- [ ] Write test cases
- [ ] Verify all tests pass

---

## Key Design Decisions

1. **Handles are stable, pool resolves change each frame**
   - Store handles in nodes
   - Resolve from pool in `compute()` or `render()`

2. **Pool is shared across all nodes**
   - Managed by base class
   - Each node registers resource usage

3. **Enum-based ports**
   - Type-safe port identification
   - Maps to specific input/output types

4. **Template-based for genericity**
   - Works with any field type (vec3, real, vec14, etc.)
   - Works with any mesh type (shell, rod, etc.)

---

## Open Questions

1. **FieldPool implementation**
   - Use existing `lewitt::render_targets::target_pool`?
   - Create new generic pool?

2. **Field data types**
   - vec3 (positions, velocities, forces)
   - real (mass, curvature values)
   - vec14 (cyclide parameters)
   - quat (orientations)

3. **Performance optimization**
   - Cache resolved values?
   - Use immutable data when possible?
   - Parallel computation support?

4. **Error handling**
   - Invalid handles?
   - Missing inputs?
   - Pool not set?
