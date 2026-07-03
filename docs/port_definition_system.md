# Port Definition System - Implementation Plan

## Overview

We're building a compiler/linker for a node-based graph system where:
- Nodes are like functions/modules with typed I/O
- Ports are the interface between nodes
- The GraphContext is the linker that allocates buffers and resolves connections

## Key Design Principles

1. **Declarative ports**: Nodes declare ports via templates, not runtime registration
2. **Two-level indirection**:
   - Local IDs: Static IDs defined in each node (0, 1, 2, ...)
   - Global IDs: Integer IDs from the linker/allocator's "global heap"
3. **Generic get/set**: Nodes use `get<PortDef>()` and `set<PortDef>(value)`
4. **Buffer abstraction**: `Buffer<T>` wraps `std::vector<T>`, agnostic to element type
5. **Context as linker**: GraphContext is the allocator and connection resolver

## Core Components

### 1. Port Definition Template

```cpp
template <typename T, typename NodeT, auto PortId>
struct PortDef {
    using buffer_type = T;  // Element type (int, vec3, etc.)
    static constexpr auto id = PortId;  // Local ID within node
};
```

**Example usage:**
```cpp
enum class GenPort { output };              // = 0
enum class MapPort { input, output };        // = 0, 1
enum class AddPort { input0, input1, output }; // = 0, 1, 2

typedef PortDef<int, GeneratorNode, GenPort::output> GenOutputPort;
typedef PortDef<int, MapPlusOneNode, MapPort::input> MapInputPort;
```

### 2. Buffer Abstraction

```cpp
template <typename T>
class Buffer {
    std::vector<T> _data;
public:
    size_t size() const;
    void resize(size_t n);
    T& get();
    void push_back(const T& value);
    // ... other methods
};
```

### 3. Node Base Class

```cpp
template <typename NodeT, typename... PortDefs>
class Node {
protected:
    GraphContext* _ctx = nullptr;
    std::array<int, sizeof...(PortDefs)> _buff_ids;  // LocalId -> GlobalBuffId

public:
    Node() = default;

    // Generic get/set operations
    template <typename PortDef>
    Buffer<typename PortDef::buffer_type>& get_buffer();

    template <typename PortDef>
    void set_buffer(const auto& value);

    // Allocation helper
    template <typename T>
    int allocate_buffer(size_t size = 0);

protected:
    // Methods to set up node
    void set_context(GraphContext* ctx) { _ctx = ctx; }
    void set_buffer_id(int local_id, int global_id) {
        _buff_ids[local_id] = global_id;
    }
};
```

### 4. GraphContext as Linker/Allocator

```cpp
class GraphContext {
private:
    std::vector<std::unique_ptr<BufferBase>> _global_buffers;  // Global heap
    std::map<std::tuple<Node*, int>, int> _connections;        // Node + LocalId -> GlobalBuffId

public:
    // Node creation
    template <typename TNode, typename... Args>
    auto create_node(Args&&... args) -> std::shared_ptr<TNode>;

    // Port linking (declarative)
    template <typename PortDef1, typename PortDef2>
    void link(PortDef1 port1, PortDef2 port2);

    // Buffer allocation
    template <typename T>
    int allocate_buffer(size_t size = 0);

    // Buffer access
    template <typename T>
    Buffer<T>& get_buffer(int global_id);
};
```

## Execution Flow

### 1. Node Creation

```cpp
auto gen = ctx.create_node<GeneratorNode>(10);
auto map = ctx.create_node<MapPlusOneNode>();
auto add = ctx.create_node<AddNode>();
```

During creation:
- Node's constructor is called (no context pointer needed)
- Node registers ports via template
- Local IDs are assigned (from PortDef::id)
- Buffers are allocated via `allocate_buffer<T>(size)`

### 2. Port Linking

```cpp
ctx.link<GeneratorNode::OutputPort, MapPlusOneNode::InputPort>(
    gen, map);
```

During linking:
- GraphContext creates a new `Buffer<T>` where `T` comes from PortDef
- Returns a `GlobalBuffId` (int)
- Sets `map._buff_ids[MapPlusOneNode::InputPort::id] = global_id`
- Sets `gen._buff_ids[GeneratorNode::OutputPort::id] = global_id`

### 3. Node Execution

```cpp
void GeneratorNode::compute() override {
    auto& buf = get_buffer<OutputPort>();  // Returns Buffer<int>&
    buf.push_back(i + 1);
}

void MapPlusOneNode::compute() override {
    auto& in = get_buffer<InputPort>();
    auto& out = get_buffer<OutputPort>();

    in.resize(out.size());
    for (size_t i = 0; i < in.size(); i++) {
        out.push_back(in.get() + 1);
    }
}
```

During execution:
- `get_buffer<PortDef>()` resolves: `gen._buff_ids[OutputPort::id]` → global_id
- Calls `_ctx->get_buffer<T>(global_id)`
- Returns `Buffer<T>&`
- Direct access to buffer data

## File Structure

```
ext/liblombardi/include/liblombardi/
├── port_def.hpp           # PortDef template
├── buffer.hpp             # Buffer<T> template
├── node_base.hpp          # Node<T, PortDefs...> base class
├── graph_context.hpp      # GraphContext linker/allocator
└── test_nodes/            # Example nodes
    ├── generator_node.hpp
    ├── map_plus_one_node.hpp
    └── add_node.hpp
```

## Migration from Old Approach

### What to Remove

1. **Node base class methods:**
   - `register_input()` and `register_output()` (no longer needed)
   - `set_context()` (context passed through templates instead)
   - `shared_from_this()` (no longer needed)

2. **GraphContext methods:**
   - `_input_ports` and `_output_ports` maps (replaced with buffer ID mapping)
   - `register_input()` and `register_output()` (replaced with `link()`)
   - `get_outputs()` and `get_inputs()` (not needed with declarative linking)

3. **Test nodes:**
   - `register_input()` and `register_output()` overrides
   - Manual port registration

### What to Add

1. **New files:**
   - `port_def.hpp` - PortDef template definition
   - `buffer.hpp` - Buffer<T> abstraction
   - Update `node_base.hpp` to use generic PortDefs
   - Update `graph_context.hpp` to be a linker/allocator
   - Rewrite test nodes to use new API

2. **Updated node base:**
   - Template parameter: `template <typename NodeT, typename... PortDefs>`
   - `_buff_ids` array: `std::array<int, sizeof...(PortDefs)>`
   - Generic `get_buffer<PortDef>()` method
   - Generic `set_buffer<PortDef>(value)` method

3. **Updated graph context:**
   - `link<PortDef1, PortDef2>(Node& from, Node& to)` method
   - Buffer allocation via template `allocate_buffer<T>(size)`
   - Global buffer heap

## Test Plan

### Test 1: Basic Node Execution

```cpp
GAUDI_TEST(BasicNodeExecution) {
    GraphContext ctx;

    // Create nodes
    auto gen = ctx.create_node<GeneratorNode>(10);
    auto map = ctx.create_node<MapPlusOneNode>();

    // Link ports
    ctx.link<GeneratorNode::OutputPort, MapPlusOneNode::InputPort>(gen, map);

    // Execute
    ctx.run();

    // Verify output
    auto& output = map.get_buffer<MapPlusOneNode::OutputPort>();
    GAUDI_ASSERT(output.size() == 10);
}
```

### Test 2: Multi-Port Node

```cpp
GAUDI_TEST(MultiPortNode) {
    GraphContext ctx;

    auto gen0 = ctx.create_node<GeneratorNode>(5);
    auto gen1 = ctx.create_node<GeneratorNode>(5);
    auto add = ctx.create_node<AddNode>();

    // Link two generators to add's inputs
    ctx.link<GeneratorNode::OutputPort, AddNode::Input0Port>(gen0, add);
    ctx.link<GeneratorNode::OutputPort, AddNode::Input1Port>(gen1, add);

    ctx.run();

    auto& result = add.get_buffer<AddNode::OutputPort>();
    // Verify result = gen0 + gen1
}
```

### Test 3: Complex Graph

```cpp
GAUDI_TEST(ComplexGraph) {
    GraphContext ctx;

    // Create multiple nodes
    auto gen0 = ctx.create_node<GeneratorNode>(10);
    auto gen1 = ctx.create_node<GeneratorNode>(10);
    auto map0 = ctx.create_node<MapPlusOneNode>();
    auto map1 = ctx.create_node<MapPlusOneNode>();
    auto add = ctx.create_node<AddNode>();

    // Wire graph
    ctx.link<GeneratorNode::OutputPort, MapPlusOneNode::InputPort>(gen0, map0);
    ctx.link<GeneratorNode::OutputPort, MapPlusOneNode::InputPort>(gen1, map1);
    ctx.link<MapPlusOneNode::OutputPort, AddNode::Input0Port>(map0, add);
    ctx.link<MapPlusOneNode::OutputPort, AddNode::Input1Port>(map1, add);

    ctx.run();

    auto& result = add.get_buffer<AddNode::OutputPort>();
    // Verify: map(map(gen0)) + map(map(gen1))
}
```

## Benefits of New Design

1. **No runtime port registration** - Everything is compile-time
2. **Type-safe** - Port linking is checked at compile time
3. **No shared_from_this()** - No ownership issues
4. **Cleaner separation** - Nodes define ports, context manages connections
5. **Explicit allocation** - Buffers are allocated declaratively
6. **Easier debugging** - Port connections are visible in the code
7. **No context pointers in user API** - Context is managed by templates

## Next Steps

1. ✅ Design port definition system
2. ✅ Design buffer abstraction
3. ⏳ Implement port_def.hpp
4. ⏳ Implement buffer.hpp
5. ⏳ Rewrite node_base.hpp with generic PortDefs
6. ⏳ Rewrite graph_context.hpp as linker/allocator
7. ⏳ Rewrite test nodes to use new API
8. ⏳ Update tests to verify new design
9. ⏳ Verify all existing functionality still works
