# Radix Sort Node-Based Workflow Design (v2)

## Overview
Radix sort as a proper node-based workflow where:
- Nodes are created explicitly
- Wiring is done via set_input/set_output
- Execution is controlled at the graph context level
- No orchestrator node needed - each pass reuses the same nodes

---

## Architecture

```
GraphContext (manages pool, nodes, execution order)
├── DigitExtractionNode (configured per-pass)
├── HistogramNode (reused across passes)
├── PrefixSumNode (reused across passes)
├── ScatterNode (reused across passes)
└── 4-8 passes of: extract → histogram → prefix_sum → scatter
```

---

## Graph Context

**Purpose**: Manage node pool, handle pool, and execution order

```cpp
#pragma once
#include "datum_pool.hpp"
#include "node_base.hpp"
#include <vector>
#include <memory>
#include <map>

namespace liblombardi {

// Base node interface
class Node {
public:
  virtual ~Node() = default;
  virtual void compute() = 0;

  // Getters for input/output datum indices
  std::vector<datum_index_t> get_input_ports() const { return _input_ports; }
  std::vector<datum_index_t> get_output_ports() const { return _output_ports; }

protected:
  // Helper to register input/output
  void register_input(datum_index_t idx) { _input_ports.push_back(idx); }
  void register_output(datum_index_t idx) { _output_ports.push_back(idx); }

  std::vector<datum_index_t> _input_ports;
  std::vector<datum_index_t> _output_ports;
};

// Graph context that manages node pool and execution
class GraphContext {
public:
  GraphContext() = default;

  // Register a node with input/output handles
  template<typename TNode, typename... Args>
  std::shared_ptr<TNode> create_node(Args&&... args) {
    auto node = std::make_shared<TNode>(std::forward<Args>(args)...);
    _nodes.push_back(node);
    return node;
  }

  // Set input for a node
  template<typename TNode>
  void set_input(TNode* node,
                 typename TNode::InputPort port,
                 datum_index_t input_idx) {
    node->set_input(port, input_idx);
  }

  // Get output from a node
  template<typename TNode>
  datum_index_t get_output(TNode* node,
                           typename TNode::OutputPort port) {
    return node->get_output(port);
  }

  // Allocate data in the pool
  template<typename T>
  datum_index_t allocate() {
    return _pool.allocate_datum<T>();
  }

  // Get data pointer
  template<typename T>
  std::vector<T>* get_data(datum_index_t idx) {
    return _pool.get_data<T>(idx);
  }

  // Run all nodes in dependency order
  void run() {
    // Simple topological sort
    // In production: more sophisticated ordering
    for (auto& node : _nodes) {
      node->compute();
    }
  }

  // Clear all data
  void clear() {
    _pool.clear();
  }

private:
  DatumPool _pool;
  std::vector<std::shared_ptr<Node>> _nodes;
};

} // namespace liblombardi
```

---

## Radix Sort Nodes

### Module 1: DigitExtractionNode

**Purpose**: Extract specific byte position from integers

```cpp
#pragma once
#include "node_base.hpp"

namespace liblombardi {

class ByteExtractionNode : public Node {
public:
  enum class InputPort { source };
  enum class OutputPort { bytes };

  ByteExtractionNode(datum_index_t source_idx, datum_index_t bytes_idx)
    : _source_idx(source_idx), _bytes_idx(bytes_idx) {
    register_input(InputPort::source);
    register_output(OutputPort::bytes);
  }

  void compute() override {
    if (!_pool) return;

    auto* source = get_data<int>(_source_idx);
    auto* bytes = get_data<uint8_t>(_bytes_idx);
    bytes->clear();

    for (int value : *source) {
      uint8_t byte = (value >> (_byte_offset * 8)) & 0xFF;
      bytes->push_back(byte);
    }
  }

  void set_byte_offset(int offset) { _byte_offset = offset; }
  int get_byte_offset() const { return _byte_offset; }

private:
  datum_index_t _source_idx;
  datum_index_t _bytes_idx;
  int _byte_offset = 0;  // 0 = LSB, 1 = next byte, etc.
};

// Type alias for 32-bit integers
using RadixSortNode32 = ByteExtractionNode;

} // namespace liblombardi
```

---

### Module 2: HistogramNode

**Purpose**: Count occurrences of each byte value (0-255)

```cpp
#pragma once
#include "node_base.hpp"

namespace liblombardi {

class HistogramNode : public Node {
public:
  enum class InputPort { bytes };
  enum class OutputPort { histogram };

  HistogramNode(datum_index_t bytes_idx, datum_index_t histogram_idx)
    : _bytes_idx(bytes_idx), _histogram_idx(histogram_idx) {
    register_input(InputPort::bytes);
    register_output(OutputPort::histogram);
  }

  void compute() override {
    if (!_pool) return;

    auto* bytes = get_data<uint8_t>(_bytes_idx);
    auto* histogram = get_data<uint64_t>(_histogram_idx);

    std::fill(histogram->begin(), histogram->end(), 0);
    for (uint8_t byte : *bytes) {
      (*histogram)[byte]++;
    }
  }

private:
  datum_index_t _bytes_idx;
  datum_index_t _histogram_idx;
};

} // namespace liblombardi
```

---

### Module 3: PrefixSumNode

**Purpose**: Compute cumulative prefix sums for scatter positions

```cpp
#pragma once
#include "node_base.hpp"

namespace liblombardi {

class PrefixSumNode : public Node {
public:
  enum class InputPort { histogram };
  enum class OutputPort { prefix_sums };

  PrefixSumNode(datum_index_t histogram_idx, datum_index_t prefix_idx)
    : _histogram_idx(histogram_idx), _prefix_idx(prefix_idx) {
    register_input(InputPort::histogram);
    register_output(OutputPort::prefix_sums);
  }

  void compute() override {
    if (!_pool) return;

    auto* histogram = get_data<uint64_t>(_histogram_idx);
    auto* prefix = get_data<uint64_t>(_prefix_idx);

    std::fill(prefix->begin(), prefix->end(), 0);
    uint64_t offset = 0;
    for (int i = 0; i < 256; ++i) {
      (*prefix)[i] = offset;
      offset += (*histogram)[i];
    }
  }

private:
  datum_index_t _histogram_idx;
  datum_index_t _prefix_idx;
};

} // namespace liblombardi
```

---

### Module 4: ScatterNode

**Purpose**: Place elements into their digit bucket positions

```cpp
#pragma once
#include "node_base.hpp"

namespace liblombardi {

class ScatterNode : public Node {
public:
  enum class InputPort { source, bytes, prefix_sums };
  enum class OutputPort { scattered };

  ScatterNode(datum_index_t source_idx, datum_index_t bytes_idx,
              datum_index_t prefix_idx, datum_index_t scattered_idx)
    : _source_idx(source_idx), _bytes_idx(bytes_idx),
      _prefix_idx(prefix_idx), _scattered_idx(scattered_idx) {
    register_input(InputPort::source);
    register_input(InputPort::bytes);
    register_input(InputPort::prefix_sums);
    register_output(OutputPort::scattered);
  }

  void compute() override {
    if (!_pool) return;

    auto* source = get_data<int>(_source_idx);
    auto* bytes = get_data<uint8_t>(_bytes_idx);
    auto* prefix = get_data<uint64_t>(_prefix_idx);
    auto* scattered = get_data<int>(_scattered_idx);

    int n = source->size();
    std::vector<int> element_indices(n);

    for (int i = 0; i < n; ++i) {
      uint8_t byte = (*bytes)[i];
      uint64_t position = (*prefix)[byte] + element_indices[byte];
      scattered->push_back((*source)[i]);
      element_indices[byte]++;
    }
  }

private:
  datum_index_t _source_idx;
  datum_index_t _bytes_idx;
  datum_index_t _prefix_idx;
  datum_index_t _scattered_idx;
};

} // namespace liblombardi
```

---

## Usage Example

### Basic Radix Sort 32-bit Integers

```cpp
void radix_sort_32bit(GraphContext& ctx, datum_index_t input_idx, datum_index_t output_idx) {
  // Allocate all outputs upfront
  datum_index_t bytes_idx = ctx.allocate<uint8_t>();
  datum_index_t histogram_idx = ctx.allocate<uint64_t>();
  datum_index_t prefix_idx = ctx.allocate<uint64_t>();
  datum_index_t scatter_idx = ctx.allocate<int>();

  // Create nodes
  auto byte_extraction = ctx.create_node<ByteExtractionNode>(
    input_idx, bytes_idx
  );
  auto histogram = ctx.create_node<HistogramNode>(
    bytes_idx, histogram_idx
  );
  auto prefix_sum = ctx.create_node<PrefixSumNode>(
    histogram_idx, prefix_idx
  );
  auto scatter = ctx.create_node<ScatterNode>(
    input_idx, bytes_idx, prefix_idx, scatter_idx
  );

  // Wire inputs (data handles)
  ctx.set_input(byte_extraction.get(), ByteExtractionNode::InputPort::source, input_idx);
  ctx.set_input(histogram.get(), HistogramNode::InputPort::bytes, bytes_idx);
  ctx.set_input(prefix_sum.get(), PrefixSumNode::InputPort::histogram, histogram_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::source, input_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::bytes, bytes_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::prefix_sums, prefix_idx);

  // Execute 4 passes (for 32-bit integers)
  for (int pass = 0; pass < 4; ++pass) {
    byte_extraction->set_byte_offset(pass);
    ctx.run();
  }

  // Get sorted result
  auto* output = ctx.get_data<int>(output_idx);
  *output = std::move(*ctx.get_data<int>(scatter_idx));
}
```

### Alternative: Template-based Factory

```cpp
namespace liblombardi {

// Factory function for 32-bit radix sort
inline void radix_sort_32bit(GraphContext& ctx,
                            datum_index_t input_idx,
                            datum_index_t output_idx,
                            int num_passes = 4) {
  // Allocate outputs
  std::array<datum_index_t, 3> intermediate_idx = {
    ctx.allocate<uint8_t>(),
    ctx.allocate<uint64_t>(),
    ctx.allocate<uint64_t>()
  };

  // Create nodes
  auto byte_extraction = ctx.create_node<ByteExtractionNode>(
    input_idx, intermediate_idx[0]
  );
  auto histogram = ctx.create_node<HistogramNode>(
    intermediate_idx[0], intermediate_idx[1]
  );
  auto prefix_sum = ctx.create_node<PrefixSumNode>(
    intermediate_idx[1], intermediate_idx[2]
  );
  auto scatter = ctx.create_node<ScatterNode>(
    input_idx, intermediate_idx[0], intermediate_idx[2], output_idx
  );

  // Wire
  ctx.set_input(byte_extraction.get(), ByteExtractionNode::InputPort::source, input_idx);
  ctx.set_input(histogram.get(), HistogramNode::InputPort::bytes, intermediate_idx[0]);
  ctx.set_input(prefix_sum.get(), PrefixSumNode::InputPort::histogram, intermediate_idx[1]);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::source, input_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::bytes, intermediate_idx[0]);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::prefix_sums, intermediate_idx[2]);

  // Execute passes
  for (int pass = 0; pass < num_passes; ++pass) {
    byte_extraction->set_byte_offset(pass);
    ctx.run();
  }
}

} // namespace liblombardi
```

---

## Multi-Stage Workflow Example

### Visualization: Separate Passes with Debugging

```cpp
void radix_sort_debug(GraphContext& ctx,
                     datum_index_t input_idx,
                     datum_index_t sorted_idx) {
  // Allocate for each stage
  auto digit_idx = ctx.allocate<uint8_t>();
  auto histogram_idx = ctx.allocate<uint64_t>();
  auto prefix_idx = ctx.allocate<uint64_t>();
  auto intermediate_idx = ctx.allocate<int>();

  // Create nodes
  auto extract = ctx.create_node<ByteExtractionNode>(input_idx, digit_idx);
  auto histo = ctx.create_node<HistogramNode>(digit_idx, histogram_idx);
  auto prefix = ctx.create_node<PrefixSumNode>(histogram_idx, prefix_idx);
  auto scatter = ctx.create_node<ScatterNode>(
    input_idx, digit_idx, prefix_idx, intermediate_idx
  );

  // Wire
  ctx.set_input(extract.get(), ByteExtractionNode::InputPort::source, input_idx);
  ctx.set_input(histo.get(), HistogramNode::InputPort::bytes, digit_idx);
  ctx.set_input(prefix.get(), PrefixSumNode::InputPort::histogram, histogram_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::source, input_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::bytes, digit_idx);
  ctx.set_input(scatter.get(), ScatterNode::InputPort::prefix_sums, prefix_idx);

  // Execute passes
  for (int pass = 0; pass < 4; ++pass) {
    extract->set_byte_offset(pass);
    ctx.run();

    // DEBUG: Inspect histogram
    auto* histo_data = ctx.get_data<uint64_t>(histogram_idx);
    std::cout << "Pass " << pass << " histogram:" << std::endl;
    for (int i = 0; i < 256; ++i) {
      if ((*histo_data)[i] > 0) {
        std::cout << "  Byte " << i << ": " << (*histo_data)[i] << std::endl;
      }
    }
  }

  // Copy final result
  auto* sorted = ctx.get_data<int>(sorted_idx);
  *sorted = std::move(*ctx.get_data<int>(intermediate_idx));
}
```

---

## Benefits of Node-Based Design

1. **No orchestrator node** - Each pass reuses the same nodes
2. **Explicit wiring** - Clear input/output connections
3. **Reusable nodes** - Same nodes used across all passes
4. **Debugging** - Can inspect each stage independently
5. **Flexibility** - Easy to swap or add nodes
6. **Modular** - Each node is independent and testable
7. **Configuration** - Each pass configures the same nodes differently

---

## Comparison with Orchestrator Approach

### Old (Orchestrator) ❌
```cpp
// All logic in one node
class RadixSortNode : public Node {
  void compute() override {
    for (int pass = 0; pass < num_passes; ++pass) {
      // Mix of extract, histogram, prefix, scatter
      // Hard to debug
      // No reuse of nodes
    }
  }
};
```

### New (Node-based) ✅
```cpp
// Create nodes once
auto extract = ctx.create_node<ByteExtractionNode>(...);
auto histo = ctx.create_node<HistogramNode>(...);
auto prefix = ctx.create_node<PrefixSumNode>(...);
auto scatter = ctx.create_node<ScatterNode>(...);

// Reuse same nodes for all passes
for (int pass = 0; pass < num_passes; ++pass) {
  extract->set_byte_offset(pass);
  ctx.run();
}

// Each stage is independent and testable
```

---

## Updated Plan Integration

### In `field_graph_system_plan_v2.md`

Replace Phase 2 with:

```markdown
## Phase 2: Radix Sort Node-Based Workflow

### 2.1 Create `ext/liblombardi/include/liblombardi/node_base.hpp`

**Purpose:** Base node interface with input/output registration

**Content:**
```cpp
#pragma once
#include "datum_pool.hpp"
#include <vector>
#include <memory>
#include <map>

namespace liblombardi {

// Base node interface
class Node {
public:
  virtual ~Node() = default;
  virtual void compute() = 0;

  std::vector<datum_index_t> get_input_ports() const { return _input_ports; }
  std::vector<datum_index_t> get_output_ports() const { return _output_ports; }

protected:
  void register_input(datum_index_t idx) { _input_ports.push_back(idx); }
  void register_output(datum_index_t idx) { _output_ports.push_back(idx); }

  std::vector<datum_index_t> _input_ports;
  std::vector<datum_index_t> _output_ports;
};

// Graph context that manages node pool and execution
class GraphContext {
public:
  template<typename TNode, typename... Args>
  std::shared_ptr<TNode> create_node(Args&&... args) {
    auto node = std::make_shared<TNode>(std::forward<Args>(args)...);
    _nodes.push_back(node);
    return node;
  }

  template<typename TNode>
  void set_input(TNode* node,
                 typename TNode::InputPort port,
                 datum_index_t input_idx) {
    node->set_input(port, input_idx);
  }

  template<typename TNode>
  datum_index_t get_output(TNode* node,
                           typename TNode::OutputPort port) {
    return node->get_output(port);
  }

  template<typename T>
  datum_index_t allocate() {
    return _pool.allocate_datum<T>();
  }

  template<typename T>
  std::vector<T>* get_data(datum_index_t idx) {
    return _pool.get_data<T>(idx);
  }

  void run() {
    for (auto& node : _nodes) {
      node->compute();
    }
  }

  void clear() {
    _pool.clear();
  }

private:
  DatumPool _pool;
  std::vector<std::shared_ptr<Node>> _nodes;
};

} // namespace liblombardi
```

### 2.2 Create `ext/liblombardi/include/liblombardi/test_nodes.hpp`

**Purpose:** Radix sort nodes (ByteExtractionNode, HistogramNode, PrefixSumNode, ScatterNode)

**Content:** (Same as above in this document)

### 2.3 Update tests

Add to `ext/liblombardi/tests/liblombardi_test.cpp`:
```cpp
TEST(TestNodesTest, GraphContextBasic) {
  GraphContext ctx;

  datum_index_t int_idx = ctx.allocate<int>();
  datum_index_t bytes_idx = ctx.allocate<uint8_t>();
  datum_index_t sorted_idx = ctx.allocate<int>();

  auto node = ctx.create_node<ByteExtractionNode>(int_idx, bytes_idx);

  ctx.set_input(node.get(), ByteExtractionNode::InputPort::source, int_idx);

  auto* data = ctx.get_data<int>(int_idx);
  data->push_back(255);
  data->push_back(0);

  ctx.run();

  auto* bytes = ctx.get_data<uint8_t>(bytes_idx);
  EXPECT_EQ(bytes->size(), 2);
  EXPECT_EQ((*bytes)[0], 255);
  EXPECT_EQ((*bytes)[1], 0);
}

TEST(TestNodesTest, RadixSort32bit) {
  GraphContext ctx;

  datum_index_t input_idx = ctx.allocate<int>();
  datum_index_t sorted_idx = ctx.allocate<int>();

  auto* input = ctx.get_data<int>(input_idx);
  input->push_back(42);
  input->push_back(17);
  input->push_back(99);
  input->push_back(3);
  input->push_back(23);

  radix_sort_32bit(ctx, input_idx, sorted_idx, 4);

  auto* sorted = ctx.get_data<int>(sorted_idx);
  EXPECT_EQ(sorted->size(), 5);
  EXPECT_EQ((*sorted)[0], 3);
  EXPECT_EQ((*sorted)[1], 17);
  EXPECT_EQ((*sorted)[2], 23);
  EXPECT_EQ((*sorted)[3], 42);
  EXPECT_EQ((*sorted)[4], 99);
}
```
```

---

## Implementation Checklist

- [ ] Create `ext/liblombardi/include/liblombardi/node_base.hpp`
- [ ] Create `ext/liblombardi/include/liblombardi/test_nodes.hpp` with radix sort nodes
- [ ] Write tests for GraphContext
- [ ] Write tests for radix sort 32-bit
- [ ] Write tests for individual nodes (ByteExtraction, Histogram, PrefixSum, Scatter)
- [ ] Test multi-pass execution
- [ ] Test with various input distributions
