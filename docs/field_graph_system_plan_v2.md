# Field Graph System Implementation Plan (v2)

## Overview
Create a type-safe, pool-based field graph system with automatic data flow, similar to the render graph system in libvermeer but for simulation fields instead of textures.

**Updated Design Decisions**:
- Use slab allocator for efficient datum memory management (O(1) allocate/deallocate)
- Follow asawa datum pattern: base class + template specialization
- Heap allocation happens once per datum (not per index)
- Port system uses enums (not strings) for type safety
- Pool is shared across all nodes (base class concern)
- Handles are stable IDs that resolve from pool each frame

**Key Insight**: The existing `lombardi/render_graph.hpp` provides a working pattern for node graphs that we can adapt.

---

## Phase 1: Base Infrastructure

### 1.1 Create `ext/liblombardi/include/liblombardi/slab.hpp`

**Purpose:** Slab allocator for efficient datum memory management

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

### 1.2 Create `ext/liblombardi/include/liblombardi/datum_pool.hpp`

**Purpose:** Datum storage following asawa datum pattern

**Content:**
```cpp
#pragma once
#include "slab.hpp"
#include <map>
#include <memory>

namespace liblombardi {

// Base datum interface
class Datum {
public:
  virtual ~Datum() = default;
  virtual void resize(index_t n) = 0;
  virtual index_t size() const = 0;
  virtual void clear() = 0;
  virtual void* get_data() = 0;
  virtual const void* get_data() const = 0;
};

// Template specialization for each data type
template <typename TYPE>
struct DatumImpl : public Datum {
  std::vector<TYPE> data;

  void resize(index_t n) override {
    data.resize(n);
  }

  index_t size() const override {
    return data.size();
  }

  void clear() override {
    data.clear();
  }

  void* get_data() override {
    return &data[0];
  }

  const void* get_data() const override {
    return &data[0];
  }
};

// Datum pool using slab allocators
class DatumPool {
public:
  // Insert a new datum
  index_t insert_datum(std::shared_ptr<Datum> datum) {
    _datums.push_back(datum);
    return _datums.size() - 1;  // Stable index!
  }

  // Get datum by index
  std::shared_ptr<Datum>& get_datum(index_t i) {
    return _datums[i];
  }

  const std::shared_ptr<Datum>& get_datum(index_t i) const {
    return _datums[i];
  }

  // Get datum of specific type (type-safe)
  template <typename T>
  std::vector<T>* get_data(index_t datum_index) {
    auto* impl = static_cast<DatumImpl<T>*>(_datums[datum_index].get());
    return &impl->data;
  }

  template <typename T>
  const std::vector<T>* get_data(index_t datum_index) const {
    auto* impl = const_cast<DatumImpl<T>*>(static_cast<const DatumImpl<T>*>(_datums[datum_index].get()));
    return &impl->data;
  }

  // Allocate and get datum of specific type (convenience)
  template <typename T>
  index_t allocate_datum() {
    auto datum = std::make_shared<DatumImpl<T>>();
    return insert_datum(datum);
  }

  // Clear all datums
  void clear() {
    _datums.clear();
  }

private:
  std::vector<std::shared_ptr<Datum>> _datums;  // Managed by pool
};

using datum_index_t = index_t;
using datum_ptr = std::shared_ptr<Datum>;

} // namespace liblombardi
```

### 1.3 Create `ext/liblombardi/include/liblombardi/ports.hpp`

**Purpose:** Define enum-based ports following ssao_node pattern

**Content:**
```cpp
#pragma once
#include "datum_pool.hpp"

namespace liblombardi {

// Port handle with stable ID
template <typename Tag>
struct PortHandle {
  datum_index_t index = 0;

  bool is_valid() const { return index != 0; }
};

// Port direction
enum class PortDirection { Input, Output };

// Port wrapper
template <typename Tag, PortDirection Dir>
struct Port {
  PortHandle<Tag> handle;
  PortDirection direction = Dir;

  Port() = default;
  Port(datum_index_t idx) : handle({idx}) {}
};

// Input port type
template <typename Tag>
using InputPort = Port<Tag, PortDirection::Input>;

// Output port type
template <typename Tag>
using OutputPort = Port<Tag, PortDirection::Output>;

} // namespace liblombardi
```

### 1.4 Create `ext/liblombardi/include/liblombardi/node_base.hpp`

**Purpose:** Base node class with pool management

**Content:**
```cpp
#pragma once
#include "datum_pool.hpp"
#include "ports.hpp"

namespace liblombardi {

// Base class for all simulation nodes
class Node {
public:
  Node() = default;
  virtual ~Node() = default;

  // Pool management - all nodes share the same pool
  void set_pool(DatumPool& pool) { _pool = &pool; }
  DatumPool* pool() const { return _pool; }

  // Allocate datum of specific type
  template <typename T>
  datum_index_t allocate_datum() {
    return _pool->allocate_datum<T>();
  }

  // Get datum data
  template <typename T>
  std::vector<T>* get_datum_data(datum_index_t idx) {
    return _pool->get_data<T>(idx);
  }

  template <typename T>
  const std::vector<T>* get_datum_data(datum_index_t idx) const {
    return _pool->get_data<T>(idx);
  }

  // Get datum
  datum_ptr get_datum(datum_index_t idx) {
    return _pool->get_datum(idx);
  }

  const datum_ptr get_datum(datum_index_t idx) const {
    return _pool->get_datum(idx);
  }

  // Compute (called each frame)
  virtual void compute() = 0;

protected:
  DatumPool* _pool = nullptr;
};

} // namespace liblombardi
```

### 1.5 Create `ext/liblombardi/tests/liblombardi_test.cpp`

**Purpose:** Tests for datum pool and basic node system

**Content:**
```cpp
#pragma once
#include <gtest/gtest.h>
#include "liblombardi/datum_pool.hpp"
#include "liblombardi/test_nodes.hpp"

using namespace liblombardi;

TEST(DatumPoolTest, BasicAllocation) {
  DatumPool pool;
  auto idx = pool.allocate_datum<int>();
  auto* data = pool.get_data<int>(idx);

  EXPECT_NE(data, nullptr);
  EXPECT_EQ(data->size(), 0);

  data->push_back(10);
  EXPECT_EQ(data->size(), 1);
  EXPECT_EQ((*data)[0], 10);
}

TEST(DatumPoolTest, MultipleTypes) {
  DatumPool pool;
  auto int_idx = pool.allocate_datum<int>();
  auto vec_idx = pool.allocate_datum<vec3>();

  auto* int_data = pool.get_data<int>(int_idx);
  auto* vec_data = pool.get_data<vec3>(vec_idx);

  int_data->push_back(42);
  vec_data->push_back(vec3(1, 2, 3));

  EXPECT_EQ((*int_data)[0], 42);
  EXPECT_EQ((*vec_data)[0], vec3(1, 2, 3));
}

TEST(DatumPoolTest, PortHandles) {
  PortHandle<int> handle;
  EXPECT_FALSE(handle.is_valid());

  handle.index = 42;
  EXPECT_TRUE(handle.is_valid());
  EXPECT_EQ(handle.index, 42);
}
```

---

## Phase 2: Test Nodes for Generic Network

### 2.1 Create `ext/liblombardi/include/liblombardi/test_nodes.hpp`

**Purpose:** Demonstrate network with integer datum and compute nodes

**Content:**
```cpp
#pragma once
#include "node_base.hpp"
#include "ports.hpp"

namespace liblombardi {

// Test node: N-binned prefix sum
template <typename T>
class BinnedPrefixSumNode : public Node {
public:
  using Node::Node;

  enum class InputPort { source };
  enum class OutputPort { prefix_sum };

  BinnedPrefixSumNode(datum_index_t source_idx, datum_index_t output_idx, int bin_count)
    : _source_idx(source_idx), _output_idx(output_idx), _bin_count(bin_count) {}

  void compute() override {
    if (!_pool) return;

    auto* source_data = get_datum_data<T>(_source_idx);
    auto* output_data = get_datum_data<T>(_output_idx);

    if (!source_data || !output_data) return;

    output_data->resize(_bin_count);

    // Compute prefix sum (simplified - just sum all values)
    T sum = T{};
    for (size_t i = 0; i < source_data->size(); ++i) {
      sum += (*source_data)[i];
    }

    T bin_size = sum / _bin_count;
    for (int i = 0; i < _bin_count; ++i) {
      (*output_data)[i] = bin_size * (i + 1);
    }
  }

private:
  datum_index_t _source_idx;
  datum_index_t _output_idx;
  int _bin_count;
};

// Test node: Radix sort based on prefix sum
template <typename T>
class RadixSortNode : public Node {
public:
  using Node::Node;

  enum class InputPort { source, prefix_sum };
  enum class OutputPort { sorted };

  RadixSortNode(datum_index_t source_idx, datum_index_t prefix_sum_idx, datum_index_t output_idx)
    : _source_idx(source_idx), _prefix_sum_idx(prefix_sum_idx), _output_idx(output_idx) {}

  void compute() override {
    if (!_pool) return;

    auto* source_data = get_datum_data<T>(_source_idx);
    auto* prefix_data = get_datum_data<T>(_prefix_sum_idx);
    auto* output_data = get_datum_data<T>(_output_idx);

    if (!source_data || !prefix_data || !output_data) return;

    // Simple radix sort (for demonstration)
    *output_data = *source_data;  // Copy for now

    // In real implementation, use prefix sum to build histogram and sort
  }

private:
  datum_index_t _source_idx;
  datum_index_t _prefix_sum_idx;
  datum_index_t _output_idx;
};

} // namespace liblombardi
```

### 2.2 Update tests

Add to `ext/liblombardi/tests/liblombardi_test.cpp`:
```cpp
TEST(TestNodesTest, BinnedPrefixSum) {
  DatumPool pool;
  auto source_idx = pool.allocate_datum<int>();
  auto output_idx = pool.allocate_datum<int>();
  auto node = BinnedPrefixSumNode<int>(source_idx, output_idx, 4);
  node.set_pool(pool);

  auto* source_data = pool.get_data<int>(source_idx);
  source_data->resize(10);
  for (int i = 0; i < 10; ++i) {
    (*source_data)[i] = i;
  }

  node.compute();

  auto* output_data = pool.get_data<int>(output_idx);
  EXPECT_EQ(output_data->size(), 4);
}

TEST(TestNodesTest, RadixSort) {
  DatumPool pool;
  auto source_idx = pool.allocate_datum<int>();
  auto prefix_idx = pool.allocate_datum<int>();
  auto output_idx = pool.allocate_datum<int>();

  auto node = RadixSortNode<int>(source_idx, prefix_idx, output_idx);
  node.set_pool(pool);

  auto* source_data = pool.get_data<int>(source_idx);
  source_data->push_back(5);
  source_data->push_back(2);
  source_data->push_back(8);
  source_data->push_back(1);

  node.compute();

  auto* output_data = pool.get_data<int>(output_idx);
  EXPECT_EQ(output_data->size(), 4);
  // Verify sorted order
  EXPECT_EQ((*output_data)[0], 1);
  EXPECT_EQ((*output_data)[3], 8);
}
```

---

## Phase 3: Integration with solver

### 3.1 Update `solver_node.hpp`

**Changes:**
- Make solver nodes inherit from `Node`
- Connect solver blocks to datum ports
- Automatic data flow from datums to solver

### 3.2 Create datum bundles

**Create:** `ext/liblombardi/include/liblombardi/datum_bundles.hpp`

**Functions:**
```cpp
#pragma once
#include "node_base.hpp"

namespace liblombardi {

// Create a bundle of datum-based constraints
template <typename NodeT>
inline auto make_datum_bundle(NodeT& node, real dt, real damping) {
  return [
    &node, dt, damping
  ](asawa::hepworth::block::solver::ptr solver) {
    auto constraint = make_datum_constraint(node, dt, damping);
    solver->add_constraint(constraint);
  };
}

} // namespace liblombardi
```

---

## Phase 4: Migration to datum pool

### 4.1 Migrate lombardi nodes

**Goal**: Update existing render nodes (ssao_node, etc.) to use datum pool instead of texture pools

**Strategy:**
1. Create conversion layer from texture inputs to datum ports
2. Use the same node graph pattern
3. Test incrementally

**Example:**
```cpp
// From lombardi/ssao_node.hpp, change:
// OLD:
lewitt::resources::handle<lewitt::resources::ssao_attachment> _ao_target;
void set_input(input, const lewitt::resources::texture_input<lewitt::resources::position_sampled> &position);

// NEW:
datum_index_t _ao_target_idx;
void set_input(input, PortHandle<position_sampled> position);

// In compute():
// Resolve each frame:
auto* ao_data = get_datum_data<position_sampled>(_ao_target_idx);
auto* pos_data = get_datum_data<position_sampled>(position);
```

### 4.2 Update render graph

Update `lombardi/render_graph.hpp` to use datum pools instead of texture pools.

---

## Phase 5: Field System

### 5.1 Implement field wrappers

Follow `docs/field_wrappers_tdd.md` using datum pool.

### 5.2 Implement solver node composition

Follow `docs/solver_node_composition_tdd.md` using datum pool.

---

## Implementation Checklist

### Phase 1: Base Infrastructure
- [ ] Create `ext/liblombardi/include/liblombardi/slab.hpp`
- [ ] Create `ext/liblombardi/include/liblombardi/datum_pool.hpp`
- [ ] Create `ext/liblombardi/include/liblombardi/ports.hpp`
- [ ] Create `ext/liblombardi/include/liblombardi/node_base.hpp`
- [ ] Create `ext/liblombardi/tests/liblombardi_test.cpp`
- [ ] Write tests for datum pool
- [ ] Write tests for test nodes (prefix sum, radix sort)

### Phase 2: Test Nodes
- [ ] Create `ext/liblombardi/include/liblombardi/test_nodes.hpp`
- [ ] Implement `BinnedPrefixSumNode`
- [ ] Implement `RadixSortNode`
- [ ] Test multi-node computation chain

### Phase 3: Integration
- [ ] Update `solver_node.hpp` to use datum pool
- [ ] Create `ext/liblombardi/include/liblombardi/datum_bundles.hpp`
- [ ] Test solver integration with datum pool

### Phase 4: Migration
- [ ] Migrate lombardi nodes to datum pool
- [ ] Migrate lewitt nodes to datum pool
- [ ] Update render graph to use datum pools
- [ ] Verify all render functionality still works

### Phase 5: Field System
- [ ] Implement field wrappers (from field_wrappers_tdd.md)
- [ ] Implement solver node composition (from solver_node_composition_tdd.md)
- [ ] Integrate field system with datum pool
- [ ] Update all tests

### Phase 6: Final Integration
- [ ] Test full simulation pipeline
- [ ] Performance profiling
- [ ] Documentation

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
   - Works with any data type (vec3, int, quat, etc.)
   - Works with any mesh type (shell, rod, etc.)

5. **Slab allocator**
   - O(1) allocate/deallocate
   - No fragmentation
   - Multiple blocks for scalability

6. **Heap allocation once per datum**
   - Not per index
   - Follows asawa datum pattern

## Files to Create

```
ext/liblombardi/
├── include/liblombardi/
│   ├── slab.hpp
│   ├── datum_pool.hpp
│   ├── ports.hpp
│   ├── node_base.hpp
│   ├── test_nodes.hpp
│   └── datum_bundles.hpp
└── tests/
    └── liblombardi_test.cpp

include/gaudi/
├── duchamp/
│   └── fields.hpp
└── hepworth/
    └── nodes/
        └── solver_node.hpp

tests/
├── field_test.cpp
└── solver_node_test.cpp
```

## Design Evolution

### Original Plan (Discarded)
- Used `std::any` for type storage
- Created separate field system
- Complex resource management

### New Approach (Adopted)
- Use slab allocator + datum pool (like asawa)
- Enum-based ports (like lombardi)
- Generic node base class
- Test nodes to validate concept
- Migration plan for vermeer

**Why the change?**
- More consistent with existing codebase
- Better performance (no std::any overhead)
- Easier to migrate existing code
- Clearer separation of concerns

## Open Questions

1. **Datum pool performance**
   - Cache resolved values?
   - Use immutable data when possible?
   - Parallel computation support?

2. **Error handling**
   - Invalid handles?
   - Missing inputs?
   - Pool not set?

3. **Migration strategy**
   - Gradual migration vs. all-at-once?
   - Backward compatibility?
   - Testing during migration?

4. **Field system integration**
   - How to handle mesh topology changes?
   - How to manage bound fields vs. independent fields?
   - Resource lifetime management?
