# liblombardi Node System Design

## Overview

A type-safe, data-flow graph system where nodes communicate through typed ports that are linked via a GraphContext. Each port has a known buffer type (datum) and a local ID.

## Key Design Principles

1. **Typed Ports**: Each port has a specific buffer type (datum) and local ID
2. **Nested Port Definitions**: Ports are organized under `input` and `output` namespaces/structs
3. **Port Identity**: Each port has `id` and `expected_type` defined at compile time
4. **Generic Buffer Access**: Nodes use `get_buffer<PortDef>()` to access typed buffers
5. **Topological Execution**: GraphContext runs nodes in valid execution order

## Core Components

### 1. Datum System

Datums are heap-allocated data containers. They can be any type:
- Vector-based: `std::vector<int>`
- Scalar: `int`
- Complex structures: custom data types

Datums provide a uniform interface for buffer access.

```cpp
/// Base datum interface
class Datum {
public:
    Datum() = default;
    virtual ~Datum() = default;

    virtual void resize(size_t size) = 0;
    virtual size_t size() const = 0;
    virtual void clear() = 0;

    /// Get raw pointer to data
    virtual void* get_data() = 0;
    virtual const void* get_data() const = 0;

private:
    datum_index_t _index;
};

/// Example 1: Vector-based datum
struct vec_int_datum : public Datum {
public:
    std::vector<int> data;

    vec_int_datum() = default;
    vec_int_datum(const std::vector<int>& data) : data(data) {}

    virtual void resize(size_t size) override {
        data.resize(size);
    }

    virtual size_t size() const override {
        return data.size();
    }

    virtual void clear() override {
        data.clear();
    }

    virtual void* get_data() override {
        return data.data();
    }

    virtual const void* get_data() const override {
        return data.data();
    }
};

/// Example 2: Scalar datum
struct int_datum : public Datum {
public:
    int value;

    int_datum() = default;
    int_datum(int value) : value(value) {}

    virtual void resize(size_t size) override {
        // Scalar datums don't resize
        (void)size;
    }

    virtual size_t size() const override {
        return 1;
    }

    virtual void clear() override {
        value = 0;
    }

    virtual void* get_data() override {
        return &value;
    }

    virtual const void* get_data() const override {
        return &value;
    }
};

/// Example 3: Complex datum (e.g., for physics data)
struct vec3_datum : public Datum {
public:
    std::vector<vec3> data;

    vec3_datum() = default;
    vec3_datum(const std::vector<vec3>& data) : data(data) {}

    virtual void resize(size_t size) override {
        data.resize(size);
    }

    virtual size_t size() const override {
        return data.size();
    }

    virtual void clear() override {
        data.clear();
    }

    virtual void* get_data() override {
        return data.data();
    }

    virtual const void* get_data() const override {
        return data.data();
    }
};
```

### 2. Port Definition Template

PortDef binds a buffer type to a port identifier.

```cpp
template <typename BufferType, auto PortId>
struct PortDef {
    using buffer_type = BufferType;  // The datum type (e.g., vec_int_datum)
    static constexpr auto id = PortId;  // Local ID (e.g., input::input0::id)
};
```

**Example:**
```cpp
// Node has ports defined under input and output namespaces
enum class {
  input0 = 0;
input1 =1;
output = 2;
}


// PortDef uses the datum type and port struct as ID
// buffer_type is the specific datum class
using InputPort0Def = PortDef<vec_int_datum, input0>;  // vec_int_datum is the datum type
using InputPort1Def = PortDef<vec_int_datum, input1>;     // int_datum is a different datum type
using OutputPortDef = PortDef<vec_int_datum, output>;   // vec3_datum is another datum type
```

### 3. Base Port Struct

Simple base for port IDs.

```cpp
struct Port {
  static inline constexpr int id = 0;
}
```

### 4. Node Base Class

```cpp
class Node {
public:
    using ptr = std::shared_ptr<Node>;

    Node() = default;
    virtual ~Node() = default;

    // Compute method to be overridden by nodes
    virtual void compute() = 0;

    // Get a buffer by its PortDef
    template <typename PortDef>
    auto& get_buffer();

    // Get buffer size
    template <typename PortDef>
    size_t get_buffer_size() {
        return get_buffer<PortDef>().size();
    }

    // Get buffer pointer
    template <typename PortDef>
    void* get_buffer_ptr() {
        return get_buffer<PortDef>().get_data();
    }

    // Set context (called by GraphContext during construction)
    void set_context(GraphContext* ctx) { _ctx = ctx; }

    // Number of ports defined by the node
    virtual uint port_count() const = 0;

protected:
    GraphContext* _ctx = nullptr;
};
```

### 5. GraphContext

The linker that allocates datums and resolves port connections.

```cpp
class GraphContext {
public:
    // Create a node with given type and constructor args
    template <typename TNode, typename... Args>
    auto create_node(Args&&... args) -> std::shared_ptr<TNode>;

    // Link two ports
    // PortDef1::id and PortDef2::id determine which buffers are connected
    template <typename PortDef1, typename PortDef2>
    void link(NodeBase& from, NodeBase& to);

    // Run all nodes in topological order
    void run();

private:
    // Allocate a datum of given type
    template <typename T>
    datum_index_t allocate_datum(size_t size = 0);

    // Get a datum by index
    Datum* get_datum(datum_index_t index);
    const Datum* get_datum(datum_index_t index) const;

    // Track which nodes use which datums
    // Node + LocalPortId -> DatumIndex
    std::unordered_map<std::tuple<NodeBase*, int>, datum_index_t> _node_port_to_datum;

    // Allocate datums
    std::vector<std::unique_ptr<Datum>> _datums;
    datum_index_t _next_datum_index = 0;
};
```

## Node Implementation Example

```cpp
namespace liblombardi::test_nodes {

// Custom datum for vector<int>
struct vec_int_datum : public Datum {
public:
    std::vector<int> data;

    vec_int_datum() = default;
    vec_int_datum(const std::vector<int>& data) : data(data) {}

    virtual void resize(size_t size) override {
        data.resize(size);
    }

    virtual size_t size() const override {
        return data.size();
    }

    virtual void clear() override {
        data.clear();
    }

    virtual void* get_data() override {
        return data.data();
    }

    virtual const void* get_data() const override {
        return data.data();
    }
};

// Base Port struct
struct Port {
    static inline constexpr int id = 0;
};

// Add node that adds two input arrays element-wise
class AddNode : public Node {
public:
    using ptr = std::shared_ptr<AddNode>;

    // Port definitions
    enum class Port {
        input0 = 0,
        input1 = 1,
        output = 2
    };

    struct input {
        struct input0 : Port {
            static inline constexpr int id = 0;
            using expected_type = vec_int_datum;
        };

        struct input1 : Port {
            static inline constexpr int id = 1;
            using expected_type = vec_int_datum;
        };
    };

    struct output {
        struct output : Port {
            static inline constexpr int id = 2;
            using expected_type = vec_int_datum;
        };
    };

    // PortDef typedefs
    using InputPort0Def = PortDef<vec_int_datum, input::input0>;
    using InputPort1Def = PortDef<vec_int_datum, input::input1>;
    using OutputPortDef = PortDef<vec_int_datum, output::output>;

    AddNode() = default;

    // Compute: add two arrays element-wise
    void compute() override {
        auto& input0 = get_buffer<InputPort0Def>();
        auto& input1 = get_buffer<InputPort1Def>();
        auto& output = get_buffer<OutputPortDef>();

        // Ensure inputs have the same size
        size_t size = input0.data().size();
        if (input1.data().size() != size) {
            throw std::runtime_error("AddNode: inputs must have the same size");
        }

        // Resize output to match input size
        output.data().resize(size);

        // Element-wise addition
        for (size_t i = 0; i < size; ++i) {
            output.data()[i] = input0.data()[i] + input1.data()[i];
        }
    }

    uint port_count() const override {
        return 3;  // Has 2 inputs and 1 output port
    }
};

} // namespace liblombardi::test_nodes
```

## Port Linking

When two ports are linked, they share the same datum:

```cpp
auto ctx = std::make_shared<GraphContext>();

// Create nodes
auto gen = ctx->create_node<GeneratorNode>(10);
auto add = ctx->create_node<AddNode>();

// Link: generator.output -> add.input0
ctx->link<GeneratorNode::OutputPortDef, AddNode::InputPort0Def>(
    *gen, *add);

// Link: generator.output -> add.input1 (shares same datum)
ctx->link<GeneratorNode::OutputPortDef, AddNode::InputPort1Def>(
    *gen, *add);

// Run the graph
ctx->run();

// Verify: add.output contains element-wise sum of generator.output
auto& result = add->get_buffer<AddNode::OutputPortDef>();
```

During linking:
- GraphContext checks that both ports have compatible types
- Creates a new datum of the appropriate type (if not already created)
- Sets up a mapping: `(node, local_port_id) -> datum_index`
- Both ports now share the same datum

## Execution Flow

### 1. Node Creation

```cpp
auto ctx = std::make_shared<GraphContext>();
auto gen = ctx->create_node<GeneratorNode>(10);
```

During creation:
- Node is constructed
- `set_context(ctx)` is called
- For each port, GraphContext allocates a datum and stores the mapping

### 2. Port Linking

```cpp
ctx->link<GeneratorNode::OutputPortDef, AddNode::InputPort0Def>(gen, add);
```

During linking:
- Retrieves datum index for generator.output
- Stores mapping for add.input0 → same datum index
- Both ports now share the same datum instance (the datum object itself, not a copy)

### 3. Node Execution

```cpp
ctx->run();
```

During execution:
- Topological sort determines valid execution order
- Nodes are executed in order
- Each node's `compute()` is called
- Nodes access buffers via `get_buffer<PortDef>()` which returns a reference to the datum
- Datums are shared, so changes propagate through the graph
- Nodes access buffers via `get_buffer<PortDef>()` which returns a reference to the datum

## File Structure

```
ext/liblombardi/include/liblombardi/
├── datum_pool.hpp          # Datum and DatumImpl classes
├── port_def.hpp            # PortDef template
├── node_base.hpp           # Node base class
├── graph_context.hpp       # GraphContext linker/allocator
└── test_nodes/             # Example nodes
    ├── generator_node.hpp
    ├── map_plus_one_node.hpp
    └── add_node.hpp
```

## Test Plan

### Test 1: Single Node Execution

```cpp
GAUDI_TEST(SingleNodeExecution) {
    auto ctx = std::make_shared<GraphContext>();
    auto gen = ctx->create_node<GeneratorNode>(5);

    ctx->run();

    auto& output = gen->get_buffer<GeneratorNode::OutputPortDef>();
    GAUDI_ASSERT(output.data().size() == 5);
    for (size_t i = 0; i < output.data().size(); ++i) {
        GAUDI_ASSERT(output.data()[i] == static_cast<int>(i + 1));
    }
}
```

### Test 2: Generator → MapPlusOne → Add Chain

```cpp
GAUDI_TEST(GraphChain) {
    auto ctx = std::make_shared<GraphContext>();

    auto gen = ctx->create_node<GeneratorNode>(10);  // [1, 2, ..., 10]
    auto map1 = ctx->create_node<MapPlusOneNode>();    // [+1]
    auto map2 = ctx->create_node<MapPlusOneNode>();    // [+1]
    auto add = ctx->create_node<AddNode>();            // element-wise sum

    // Connect generator -> map1
    ctx->link<GeneratorNode::OutputPortDef, MapPlusOneNode::InputPortDef>(
        *gen, *map1);

    // Connect generator -> map2 (same input)
    ctx->link<GeneratorNode::OutputPortDef, MapPlusOneNode::InputPortDef>(
        *gen, *map2);

    // Connect map1 -> add.input0
    ctx->link<MapPlusOneNode::OutputPortDef, AddNode::InputPort0Def>(
        *map1, *add);

    // Connect map2 -> add.input1
    ctx->link<MapPlusOneNode::OutputPortDef, AddNode::InputPort1Def>(
        *map2, *add);

    ctx->run();

    // Verify: map1 and map2 share generator's datum
    // add.output = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    auto& result = add->get_buffer<AddNode::OutputPortDef>();
    for (size_t i = 0; i < result.data().size(); ++i) {
        GAUDI_ASSERT(result.data()[i] == static_cast<int>((i + 1) + 1 + (i + 1) + 1));
    }
}
```

### Test 3: Topological Execution Order

```cpp
GAUDI_TEST(TopologicalOrder) {
    auto ctx = std::make_shared<GraphContext>();

    auto gen = ctx->create_node<GeneratorNode>(5);
    auto map = ctx->create_node<MapPlusOneNode>();
    auto add = ctx->create_node<AddNode>();

    // Connect: gen -> map -> add
    ctx->link<GeneratorNode::OutputPortDef, MapPlusOneNode::InputPortDef>(gen, map);
    ctx->link<MapPlusOneNode::OutputPortDef, AddNode::InputPort0Def>(map, add);

    ctx->run();

    // Verify execution order: compute(gen) -> compute(map) -> compute(add)
    GAUDI_ASSERT(gen->_computed);
    GAUDI_ASSERT(map->_computed);
    GAUDI_ASSERT(add->_computed);
}
```

## Benefits

1. **Type Safety**: Port linking is checked at compile time via PortDef
2. **Explicit Buffer Types**: Clear what data flows through each port
3. **Reusability**: Same datum type can be used across many ports
4. **In-place Updates**: Buffers are shared, avoiding copies
5. **Easy Debugging**: Port connections are explicit in code
6. **No Runtime Registration**: Everything is compile-time
