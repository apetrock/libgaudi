# liblombardi Migration Plan

## Overview
Migrate render graph patterns from `ext/libvermeer/include/lombardi` to a new generic `ext/liblombardi` library with a data flow graph API.

## Current State

### ✅ Completed
- PortDef CRTP pattern with type-safe port definitions
- Node base class with `GraphContextBase*` pointer
- Basic test nodes: Generator, MapPlusOne, Add
- Test suite written in `tests/liblombardi_tests.cpp`

### ❌ Broken
1. **graph_context.hpp**: Duplicate `link()` methods (lines 53-69 and 133-149)
2. **Missing `run()` method**: Tests call `ctx.run()` but it doesn't exist
3. **test_nodes.hpp**: Duplicate port definitions and implementations
4. **Test references**: Use inconsistent port def names (`AddNode::OutputPortDef` vs `AddNodeBase::OutputPortDef`)

## Phase 1: Fix liblombardi Core

### 1.1 Fix graph_context.hpp
- Remove duplicate `link()` method (lines 133-149)
- Keep only the templated version with SFINAE (lines 53-69)

### 1.2 Implement run() method
```cpp
void run() {
    // Topological sort to determine execution order
    // Execute nodes in order
    // Track which nodes have been executed
}
```

### 1.3 Fix test_nodes.hpp
- Consolidate to single implementation per node
- Remove duplicate port definitions
- Ensure consistent PortDef usage

### 1.4 Fix test references
- Update `liblombardi_tests.cpp` to use consistent PortDef names
- Change `AddNode::OutputPortDef` → `AddNodeBase::OutputPortDef`
- Change `AddNode::InputPortDef0` → `AddNodeBase::InputPortDef0`

## Phase 2: Test Suite Execution

### 2.1 Build tests
```bash
cmake --build build --target gaudi_tests -j8
```

### 2.2 Run tests
```bash
./build/tests/gaudi_tests --gtest_filter=liblombardi_*
```

### 2.3 Expected test outputs
- **BasicGraphExecution**: `[4, 5, 6, 7, 8, 9, 10, 11, 12, 13]`
- **NodeExecutionOrder**: `[3, 4, 5, 6, 7]`
- **PoolAllocation**: Buffer management working
- **ErrorHandling**: Error handling works

## Phase 3: Architecture Documentation

### 3.1 Testing documentation
Create `docs/testing.md` documenting:
- Test locations (`tests/`)
- Testing framework (GAUDI_TEST)
- Test paradigm (data flow graph testing)
- How to run tests
- How to add new tests

### 3.2 Build documentation
Create `docs/building.md` documenting:
- CMake build system
- Build targets
- How to compile
- Dependencies

### 3.3 Remove old lombardi documentation
- Update or remove `ext/libvermeer/` lombardi references
- Point to new `ext/liblombardi/` location

## Phase 4: Migrate vermeer nodes (optional, future)

### 4.1 Port graph nodes
Migrate render nodes from `ext/libvermeer/include/lombardi/`:
- `g_buffer_node` → `ext/liblombardi/nodes/`
- `ssao_node` → `ext/liblombardi/nodes/`
- `deferred_lighting_node` → `ext/liblombardi/nodes/`

### 4.2 Migrate render_graph
Migrate render graph from `ext/libvermeer/include/lombardi/render_graph.hpp`:
- `render_graph` → `ext/liblombardi/graph.hpp`
- Add datum pool integration
- Use PortDef-based API

### 4.3 Add GPU-specific nodes
- Use WebGPU bindings
- Implement render passes
- Integrate with lewitt rendering system

## Phase 5: Apply to Original TDDs

### 5.1 Update field_wrappers_tdd.md
- Remove duplicate design sections
- Update to reference datum pool instead of implementing from scratch
- Add dependency on liblombardi

### 5.2 Update solver_node_composition_tdd.md
- Remove duplicate design sections
- Update to use graph context
- Show examples with node-based API

### 5.3 Implement field wrappers
- Create `include/gaudi/duchamp/fields.hpp`
- Use datum pool for storage
- Implement bound/surface_following/sampled field types

## Phase 6: Extend with Advanced Nodes

### 6.1 Sorting nodes (from conversation)
Implement radix sort and histogram nodes:
- `BinnedPrefixSumNode`: Prefix sum of histogram bins
- `RadixSortNode`: Full radix sort implementation

### 6.2 Test patterns
- Generator → Histogram → Prefix Sum
- Generator → Prefix Sum → Radix Sort

## Decision Points

### Datum Pool vs Generic Buffer
- **Current**: Using `Buffer<T>` class
- **Future**: Could migrate to datum pool pattern
- **Decision**: Keep generic buffer for now, migrate after Phase 1

### Port Definition API
- **Current**: CRTP with `PortDef<T, NodeT, ID>`
- **Alternative**: Enum-based with helper typedefs
- **Decision**: Keep CRTP for type safety, add enum helpers if needed

### Graph Execution Model
- **Current**: Topological sort in `run()`
- **Alternative**: Event-driven or pipeline-based
- **Decision**: Keep topological sort for simplicity

## Implementation Checklist

### Phase 1: Fix Core
- [ ] Fix graph_context.hpp (remove duplicate link method)
- [ ] Implement GraphContext::run() with topological sort
- [ ] Fix test_nodes.hpp (consolidate implementations)
- [ ] Fix test references (use consistent PortDef names)
- [ ] Build and run tests

### Phase 2: Documentation
- [ ] Create `docs/testing.md`
- [ ] Create `docs/building.md`
- [ ] Update field_wrappers_tdd.md (remove duplicates)
- [ ] Update solver_node_composition_tdd.md (remove duplicates)

### Phase 3: Advanced Testing (optional)
- [ ] Implement BinnedPrefixSumNode
- [ ] Implement RadixSortNode
- [ ] Add sorting tests

### Phase 4: Migration (future)
- [ ] Migrate g_buffer_node
- [ ] Migrate ssao_node
- [ ] Migrate deferred_lighting_node
- [ ] Migrate render_graph

## Key Design Decisions

### Why CRTP?
- Compile-time type safety
- No runtime type checking
- Clean separation of interface and implementation

### Why GraphContextBase*?
- Enables flexible buffer management
- Supports lazy allocation
- Decouples node implementation from context

### Why generic Buffer<T>?
- Simple, type-safe
- Easy to debug
- Can migrate to datum pool later if needed

## Risks

1. **Memory leaks**: Buffer allocation without cleanup
   - **Mitigation**: Use smart pointers, ensure RAII

2. **Topological sort complexity**: O(N+E) but could be slow
   - **Mitigation**: Cache execution order, optimize for common patterns

3. **Circular dependencies**: Nodes that depend on themselves
   - **Mitigation**: Detect cycles in run() and throw error

4. **Type safety**: Dynamic casts could fail
   - **Mitigation**: Use compile-time assertions, ensure consistent usage

## Success Criteria

- All 4 liblombardi tests pass
- Documentation is complete and discoverable
- Code is clean and well-organized
- Architecture is extensible for future work
