# Radix Sort Node-Based Workflow Design

## Overview
Radix sort as a node-based workflow with modular stages, enabling:
- Clear data flow between stages
- Intermediate output debugging
- Reusable digit extraction/combiner nodes
- Flexible histogram and prefix sum computation

---

## Node-Based Workflow Architecture

```
Original Array (int) → Digit Extraction → Histogram → Prefix Sum → Scatter → Gather → Sorted Array
                        ↓                 ↓            ↓           ↓
                    Digit Values    Histogram    Prefix Sum   Scattered
```

---

## Proposed Modules

### Module 1: Digit Extraction Node
**Purpose**: Extract specific digit position from integers

**Inputs**:
- Source array: `std::vector<int>` or `std::vector<T>`

**Outputs**:
- Extracted digit values: `std::vector<uint8_t>` (0-255 for each position)

**Data Types**:
```cpp
// Input/output data types
struct DigitExtractionInput {
  const std::vector<int>* source;
  int digit_position;  // 0 = LSB, 1 = next byte, etc.
};

struct DigitExtractionOutput {
  std::vector<uint8_t> digits;
  int digit_count;  // Total number of digit values extracted
};
```

**Algorithm**:
```cpp
std::vector<uint8_t> extract_digit(const std::vector<int>& source, int digit_position) {
  std::vector<uint8_t> digits;
  digits.reserve(source.size());

  for (int value : source) {
    // Shift right by digit_position * 8 bits
    int digit = (value >> (digit_position * 8)) & 0xFF;
    digits.push_back(static_cast<uint8_t>(digit));
  }

  return digits;
}
```

---

### Module 2: Histogram Construction Node
**Purpose**: Count occurrences of each digit value

**Inputs**:
- Digit values: `std::vector<uint8_t>`

**Outputs**:
- Histogram: `std::vector<uint64_t>` (counts for digits 0-255)

**Data Types**:
```cpp
struct HistogramInput {
  const std::vector<uint8_t>* digit_values;
};

struct HistogramOutput {
  std::vector<uint64_t> histogram;  // histogram[0] = count of 0, histogram[255] = count of 255
  uint64_t total_count;  // Sum of all counts
};
```

**Algorithm**:
```cpp
std::vector<uint64_t> build_histogram(const std::vector<uint8_t>& digit_values) {
  std::vector<uint64_t> histogram(256, 0);

  for (uint8_t digit : digit_values) {
    histogram[digit]++;
  }

  // Calculate total count
  uint64_t total = 0;
  for (uint64_t count : histogram) {
    total += count;
  }

  return histogram;
}
```

---

### Module 3: Prefix Sum Node
**Purpose**: Compute cumulative prefix sums for scatter positioning

**Inputs**:
- Histogram: `std::vector<uint64_t>`

**Outputs**:
- Prefix sums: `std::vector<uint64_t>` (positions where each digit group starts)

**Data Types**:
```cpp
struct PrefixSumInput {
  const std::vector<uint64_t>* histogram;
};

struct PrefixSumOutput {
  std::vector<uint64_t> prefix_sums;
  uint64_t total_elements;  // Total number of elements being sorted
};
```

**Algorithm**:
```cpp
std::vector<uint64_t> compute_prefix_sums(const std::vector<uint64_t>& histogram) {
  std::vector<uint64_t> prefix_sums(256);
  uint64_t offset = 0;

  for (int i = 0; i < 256; ++i) {
    prefix_sums[i] = offset;
    offset += histogram[i];
  }

  return prefix_sums;
}
```

**Use Case**: These prefix sums tell us where each digit group starts in the final sorted array.

---

### Module 4: Scatter Node
**Purpose**: Place elements into their digit bucket positions

**Inputs**:
- Original array: `std::vector<int>`
- Digit values: `std::vector<uint8_t>`
- Prefix sums: `std::vector<uint64_t>`

**Outputs**:
- Scattered array: `std::vector<int>` (elements grouped by digit value)

**Data Types**:
```cpp
struct ScatterInput {
  const std::vector<int>* source;
  const std::vector<uint8_t>* digit_values;
  const std::vector<uint64_t>* prefix_sums;
};

struct ScatterOutput {
  std::vector<int> scattered;
  std::vector<int> element_indices;  // Track which source element went to which position
};
```

**Algorithm**:
```cpp
std::vector<int> scatter(const std::vector<int>& source,
                        const std::vector<uint8_t>& digit_values,
                        const std::vector<uint64_t>& prefix_sums) {
  int n = source.size();
  std::vector<int> scattered(n);
  std::vector<int> element_indices(n);  // Track source indices

  for (int i = 0; i < n; ++i) {
    uint8_t digit = digit_values[i];
    uint64_t position = prefix_sums[digit] + element_indices[digit];
    scattered[position] = source[i];
    element_indices[digit]++;
  }

  return scattered;
}
```

---

### Module 5: Gather Node
**Purpose**: Collect scattered elements into sorted order

**Inputs**:
- Scattered array: `std::vector<int>`
- Prefix sums: `std::vector<uint64_t>` (same as used in scatter)

**Outputs**:
- Sorted array: `std::vector<int>`

**Data Types**:
```cpp
struct GatherInput {
  const std::vector<int>* scattered;
  const std::vector<uint64_t>* prefix_sums;
};

struct GatherOutput {
  std::vector<int> sorted;
};
```

**Algorithm**:
```cpp
std::vector<int> gather(const std::vector<int>& scattered,
                       const std::vector<uint64_t>& prefix_sums) {
  int n = scattered.size();
  std::vector<int> sorted(n);

  for (int i = 0; i < n; ++i) {
    uint8_t digit = (scattered[i] >> 24) & 0xFF;  // Extract most significant digit
    uint64_t position = prefix_sums[digit] + i;
    sorted[position] = scattered[i];
  }

  return sorted;
}
```

---

### Module 6: Multi-Pass Radix Sort Node
**Purpose**: Apply digit extraction and processing repeatedly

**Inputs**:
- Source array: `std::vector<int>`
- Number of passes: `int` (e.g., 4 for 32-bit integers)

**Outputs**:
- Sorted array: `std::vector<int>`

**Data Types**:
```cpp
struct RadixSortNodeInput {
  const std::vector<int>* source;
  int num_passes;  // Typically 4 for 32-bit integers
};

struct RadixSortNodeOutput {
  std::vector<int> sorted;
};
```

**Algorithm**:
```cpp
std::vector<int> radix_sort(const std::vector<int>& source, int num_passes = 4) {
  std::vector<int> current = source;
  std::vector<int> digit_values;
  std::vector<int> scattered;

  for (int pass = 0; pass < num_passes; ++pass) {
    // Extract digits
    digit_values = extract_digit(current, pass);

    // Build histogram
    auto histogram = build_histogram(digit_values);

    // Compute prefix sums
    auto prefix_sums = compute_prefix_sums(histogram);

    // Scatter
    scattered = scatter(current, digit_values, prefix_sums);

    // Gather (extract new digits from scattered for next pass)
    current = scatter;  // In practice, we'd use the prefix_sums to gather
  }

  return current;
}
```

---

## Intermediate Output Analysis

### Why Multiple Outputs Matter

1. **Debugging**: Can inspect histogram distribution, prefix sums, etc.
2. **Verification**: Check that scatter positions are correct
3. **Optimization**: Can reuse histograms across passes
4. **Alternative algorithms**: Can use histogram to implement different sorts (e.g., counting sort)

### Proposed Intermediate Data Types

| Module | Output | Data Type | Purpose |
|--------|--------|-----------|---------|
| 1 | Extracted digits | `std::vector<uint8_t>` | Verify digit extraction |
| 2 | Histogram | `std::vector<uint64_t>` | Check digit distribution |
| 3 | Prefix sums | `std::vector<uint64_t>` | Verify position calculations |
| 4 | Scattered array | `std::vector<int>` | Verify scatter placement |
| 5 | Sorted array | `std::vector<int>` | Final result |

---

## Node-Based Implementation Plan

### Node Class Hierarchy

```cpp
namespace liblombardi {

// Base node interface
class Node {
public:
  virtual void compute() = 0;
  virtual std::string name() const = 0;
};

// Digit extraction node
class DigitExtractionNode : public Node {
public:
  enum class InputPort { source };
  enum class OutputPort { digits };

  void compute() override {
    auto* source = get_input_data<int>(InputPort::source);
    int digit_position = _digit_position;

    std::vector<uint8_t> digits;
    for (int value : *source) {
      uint8_t digit = (value >> (digit_position * 8)) & 0xFF;
      digits.push_back(digit);
    }

    auto* output = get_output_data<uint8_t>(OutputPort::digits);
    *output = std::move(digits);
  }

private:
  int _digit_position = 0;
};

// Histogram construction node
class HistogramNode : public Node {
public:
  enum class InputPort { digit_values };
  enum class OutputPort { histogram };

  void compute() override {
    auto* digit_values = get_input_data<uint8_t>(InputPort::digit_values);

    std::vector<uint64_t> histogram(256, 0);
    for (uint8_t digit : *digit_values) {
      histogram[digit]++;
    }

    auto* output = get_output_data<uint64_t>(OutputPort::histogram);
    *output = std::move(histogram);
  }
};

// Prefix sum node
class PrefixSumNode : public Node {
public:
  enum class InputPort { histogram };
  enum class OutputPort { prefix_sums };

  void compute() override {
    auto* histogram = get_input_data<uint64_t>(InputPort::histogram);

    std::vector<uint64_t> prefix_sums(256);
    uint64_t offset = 0;
    for (int i = 0; i < 256; ++i) {
      prefix_sums[i] = offset;
      offset += (*histogram)[i];
    }

    auto* output = get_output_data<uint64_t>(OutputPort::prefix_sums);
    *output = std::move(prefix_sums);
  }
};

// Scatter node
class ScatterNode : public Node {
public:
  enum class InputPort { source, digit_values, prefix_sums };
  enum class OutputPort { scattered };

  void compute() override {
    auto* source = get_input_data<int>(InputPort::source);
    auto* digit_values = get_input_data<uint8_t>(InputPort::digit_values);
    auto* prefix_sums = get_input_data<uint64_t>(InputPort::prefix_sums);

    int n = source->size();
    std::vector<int> scattered(n);
    std::vector<int> element_indices(n);

    for (int i = 0; i < n; ++i) {
      uint8_t digit = (*digit_values)[i];
      uint64_t position = (*prefix_sums)[digit] + element_indices[digit];
      scattered[position] = (*source)[i];
      element_indices[digit]++;
    }

    auto* output = get_output_data<int>(OutputPort::scattered);
    *output = std::move(scattered);
  }
};

// Multi-pass radix sort node (orchestrates all stages)
class RadixSortNode : public Node {
public:
  enum class InputPort { source };
  enum class OutputPort { sorted };

  RadixSortNode(datum_index_t source_idx, datum_index_t sorted_idx, int num_passes)
    : _source_idx(source_idx), _sorted_idx(sorted_idx), _num_passes(num_passes) {}

  void compute() override {
    if (!_pool) return;

    // Copy source data
    auto* source = get_datum_data<int>(_source_idx);
    std::vector<int> current = *source;

    // Apply radix sort passes
    for (int pass = 0; pass < _num_passes; ++pass) {
      // Extract digits
      auto digit_idx = allocate_datum<uint8_t>();
      auto digit_data = get_datum_data<uint8_t>(digit_idx);
      for (int value : current) {
        uint8_t digit = (value >> (pass * 8)) & 0xFF;
        digit_data->push_back(digit);
      }

      // Build histogram
      auto histogram_idx = allocate_datum<uint64_t>();
      auto histogram_data = get_datum_data<uint64_t>(histogram_idx);
      std::vector<uint64_t> histogram(256, 0);
      for (uint8_t digit : *digit_data) {
        histogram[digit]++;
      }
      *histogram_data = std::move(histogram);

      // Compute prefix sums
      auto prefix_idx = allocate_datum<uint64_t>();
      auto prefix_data = get_datum_data<uint64_t>(prefix_idx);
      std::vector<uint64_t> prefix_sums(256);
      uint64_t offset = 0;
      for (int i = 0; i < 256; ++i) {
        prefix_sums[i] = offset;
        offset += histogram[i];
      }
      *prefix_data = std::move(prefix_sums);

      // Scatter
      auto scatter_idx = allocate_datum<int>();
      auto scatter_data = get_datum_data<int>(scatter_idx);
      int n = current.size();
      std::vector<int> scattered(n);
      std::vector<int> element_indices(n);

      for (int i = 0; i < n; ++i) {
        uint8_t digit = (*digit_data)[i];
        uint64_t position = (*prefix_data)[digit] + element_indices[digit];
        scattered[position] = current[i];
        element_indices[digit]++;
      }
      *scatter_data = std::move(scattered);

      current = std::move(scattered);
    }

    // Final gather to sorted position
    auto* sorted_data = get_datum_data<int>(_sorted_idx);
    *sorted_data = std::move(current);
  }

private:
  datum_index_t _source_idx;
  datum_index_t _sorted_idx;
  int _num_passes;
};

} // namespace liblombardi
```

---

## Updated Plan Integration

### In `field_graph_system_plan_v2.md`

Replace Phase 2 with:

```markdown
## Phase 2: Radix Sort Node-Based Workflow

### 2.1 Create `ext/liblombardi/include/liblombardi/test_nodes.hpp`

**Purpose:** Demonstrate node-based radix sort with modular stages

**Content:**
```cpp
#pragma once
#include "node_base.hpp"
#include "ports.hpp"

namespace liblombardi {

// Module 1: Digit extraction node
class DigitExtractionNode : public Node {
public:
  enum class InputPort { source };
  enum class OutputPort { digits };

  DigitExtractionNode(datum_index_t source_idx, datum_index_t digits_idx, int digit_position)
    : _source_idx(source_idx), _digits_idx(digits_idx), _digit_position(digit_position) {}

  void compute() override {
    if (!_pool) return;

    auto* source = get_datum_data<int>(_source_idx);
    auto* digits_data = get_datum_data<uint8_t>(_digits_idx);
    digits_data->clear();

    for (int value : *source) {
      uint8_t digit = (value >> (_digit_position * 8)) & 0xFF;
      digits_data->push_back(digit);
    }
  }

private:
  datum_index_t _source_idx;
  datum_index_t _digits_idx;
  int _digit_position;
};

// Module 2: Histogram construction node
class HistogramNode : public Node {
public:
  enum class InputPort { digit_values };
  enum class OutputPort { histogram };

  HistogramNode(datum_index_t digit_idx, datum_index_t histogram_idx)
    : _digit_idx(digit_idx), _histogram_idx(histogram_idx) {}

  void compute() override {
    if (!_pool) return;

    auto* digit_data = get_datum_data<uint8_t>(_digit_idx);
    auto* histogram_data = get_datum_data<uint64_t>(_histogram_idx);

    std::fill(histogram_data->begin(), histogram_data->end(), 0);
    for (uint8_t digit : *digit_data) {
      (*histogram_data)[digit]++;
    }
  }

private:
  datum_index_t _digit_idx;
  datum_index_t _histogram_idx;
};

// Module 3: Prefix sum node
class PrefixSumNode : public Node {
public:
  enum class InputPort { histogram };
  enum class OutputPort { prefix_sums };

  PrefixSumNode(datum_index_t histogram_idx, datum_index_t prefix_idx)
    : _histogram_idx(histogram_idx), _prefix_idx(prefix_idx) {}

  void compute() override {
    if (!_pool) return;

    auto* histogram_data = get_datum_data<uint64_t>(_histogram_idx);
    auto* prefix_data = get_datum_data<uint64_t>(_prefix_idx);

    std::fill(prefix_data->begin(), prefix_data->end(), 0);
    uint64_t offset = 0;
    for (int i = 0; i < 256; ++i) {
      (*prefix_data)[i] = offset;
      offset += (*histogram_data)[i];
    }
  }

private:
  datum_index_t _histogram_idx;
  datum_index_t _prefix_idx;
};

// Module 4: Scatter node
class ScatterNode : public Node {
public:
  enum class InputPort { source, digit_values, prefix_sums };
  enum class OutputPort { scattered };

  ScatterNode(datum_index_t source_idx, datum_index_t digit_idx,
              datum_index_t prefix_idx, datum_index_t scattered_idx)
    : _source_idx(source_idx), _digit_idx(digit_idx),
      _prefix_idx(prefix_idx), _scattered_idx(scattered_idx) {}

  void compute() override {
    if (!_pool) return;

    auto* source = get_datum_data<int>(_source_idx);
    auto* digit_data = get_datum_data<uint8_t>(_digit_idx);
    auto* prefix_data = get_datum_data<uint64_t>(_prefix_idx);
    auto* scattered_data = get_datum_data<int>(_scattered_idx);

    int n = source->size();
    std::vector<int> scattered(n);
    std::vector<int> element_indices(n);

    for (int i = 0; i < n; ++i) {
      uint8_t digit = (*digit_data)[i];
      uint64_t position = (*prefix_data)[digit] + element_indices[digit];
      scattered[position] = (*source)[i];
      element_indices[digit]++;
    }

    *scattered_data = std::move(scattered);
  }

private:
  datum_index_t _source_idx;
  datum_index_t _digit_idx;
  datum_index_t _prefix_idx;
  datum_index_t _scattered_idx;
};

// Module 5: Multi-pass radix sort node (orchestrator)
class RadixSortNode : public Node {
public:
  enum class InputPort { source };
  enum class OutputPort { sorted };

  RadixSortNode(datum_index_t source_idx, datum_index_t sorted_idx, int num_passes = 4)
    : _source_idx(source_idx), _sorted_idx(sorted_idx), _num_passes(num_passes) {}

  void compute() override {
    if (!_pool) return;

    // Copy source data
    auto* source = get_datum_data<int>(_source_idx);
    std::vector<int> current = *source;

    // Apply radix sort passes
    for (int pass = 0; pass < _num_passes; ++pass) {
      // Extract digits
      auto digit_idx = allocate_datum<uint8_t>();
      auto digit_data = get_datum_data<uint8_t>(digit_idx);
      digit_data->clear();
      for (int value : current) {
        uint8_t digit = (value >> (pass * 8)) & 0xFF;
        digit_data->push_back(digit);
      }

      // Build histogram
      auto histogram_idx = allocate_datum<uint64_t>();
      auto histogram_data = get_datum_data<uint64_t>(histogram_idx);
      std::fill(histogram_data->begin(), histogram_data->end(), 0);
      for (uint8_t digit : *digit_data) {
        (*histogram_data)[digit]++;
      }

      // Compute prefix sums
      auto prefix_idx = allocate_datum<uint64_t>();
      auto prefix_data = get_datum_data<uint64_t>(prefix_idx);
      std::fill(prefix_data->begin(), prefix_data->end(), 0);
      uint64_t offset = 0;
      for (int i = 0; i < 256; ++i) {
        (*prefix_data)[i] = offset;
        offset += (*histogram_data)[i];
      }

      // Scatter
      auto scatter_idx = allocate_datum<int>();
      auto scatter_data = get_datum_data<int>(scatter_idx);
      int n = current.size();
      std::vector<int> scattered(n);
      std::vector<int> element_indices(n);

      for (int i = 0; i < n; ++i) {
        uint8_t digit = (*digit_data)[i];
        uint64_t position = (*prefix_data)[digit] + element_indices[digit];
        scattered[position] = current[i];
        element_indices[digit]++;
      }
      *scatter_data = std::move(scattered);

      current = std::move(scattered);
    }

    // Final result
    auto* sorted_data = get_datum_data<int>(_sorted_idx);
    *sorted_data = std::move(current);
  }

private:
  datum_index_t _source_idx;
  datum_index_t _sorted_idx;
  int _num_passes;
};

// Original simple nodes for comparison
template <typename T>
class BinnedPrefixSumNode : public Node {
  // ... (keep as-is)
};

template <typename T>
class RadixSortNodeLegacy : public Node {
  // ... (keep as-is, for comparison)
};

} // namespace liblombardi
```

### 2.2 Update tests

Add to `ext/liblombardi/tests/liblombardi_test.cpp`:
```cpp
TEST(TestNodesTest, DigitExtraction) {
  DatumPool pool;
  auto source_idx = pool.allocate_datum<int>();
  auto digit_idx = pool.allocate_datum<uint8_t>();
  auto node = DigitExtractionNode(source_idx, digit_idx, 0);
  node.set_pool(pool);

  auto* source = pool.get_data<int>(source_idx);
  source->push_back(0xFF);      // 0xFF = 255 (0xFF)
  source->push_back(0x00);      // 0x00 = 0
  source->push_back(0xAB);      // 0xAB = 171

  node.compute();

  auto* digits = pool.get_data<uint8_t>(digit_idx);
  EXPECT_EQ(digits->size(), 3);
  EXPECT_EQ((*digits)[0], 0xFF);
  EXPECT_EQ((*digits)[1], 0x00);
  EXPECT_EQ((*digits)[2], 0xAB);
}

TEST(TestNodesTest, Histogram) {
  DatumPool pool;
  auto digit_idx = pool.allocate_datum<uint8_t>();
  auto histogram_idx = pool.allocate_datum<uint64_t>();
  auto node = HistogramNode(digit_idx, histogram_idx);
  node.set_pool(pool);

  auto* digit = pool.get_data<uint8_t>(digit_idx);
  digit->push_back(0x01);
  digit->push_back(0x02);
  digit->push_back(0x01);

  node.compute();

  auto* histogram = pool.get_data<uint64_t>(histogram_idx);
  EXPECT_EQ((*histogram)[0x00], 0);
  EXPECT_EQ((*histogram)[0x01], 2);
  EXPECT_EQ((*histogram)[0x02], 1);
  EXPECT_EQ((*histogram)[0x03], 0);
  // ... verify all 256 values
}

TEST(TestNodesTest, PrefixSum) {
  DatumPool pool;
  auto histogram_idx = pool.allocate_datum<uint64_t>();
  auto prefix_idx = pool.allocate_datum<uint64_t>();
  auto node = PrefixSumNode(histogram_idx, prefix_idx);
  node.set_pool(pool);

  auto* histogram = pool.get_data<uint64_t>(histogram_idx);
  // Set specific values
  (*histogram)[0] = 10;
  (*histogram)[1] = 5;
  (*histogram)[2] = 3;
  (*histogram)[3] = 0;
  // ... set other values to 0

  node.compute();

  auto* prefix = pool.get_data<uint64_t>(prefix_idx);
  EXPECT_EQ((*prefix)[0], 0);
  EXPECT_EQ((*prefix)[1], 10);
  EXPECT_EQ((*prefix)[2], 15);
  EXPECT_EQ((*prefix)[3], 18);
  // ... verify all 256 prefix sums
}

TEST(TestNodesTest, Scatter) {
  DatumPool pool;
  auto source_idx = pool.allocate_datum<int>();
  auto digit_idx = pool.allocate_datum<uint8_t>();
  auto prefix_idx = pool.allocate_datum<uint64_t>();
  auto scattered_idx = pool.allocate_datum<int>();
  auto node = ScatterNode(source_idx, digit_idx, prefix_idx, scattered_idx);
  node.set_pool(pool);

  auto* source = pool.get_data<int>(source_idx);
  source->push_back(10);
  source->push_back(5);
  source->push_back(20);
  source->push_back(5);

  auto* digit = pool.get_data<uint8_t>(digit_idx);
  digit->push_back(0x01);  // digit 1
  digit->push_back(0x05);  // digit 5
  digit->push_back(0x01);  // digit 1
  digit->push_back(0x05);  // digit 5

  auto* prefix = pool.get_data<uint64_t>(prefix_idx);
  (*prefix)[0x01] = 0;     // start of digit 1
  (*prefix)[0x05] = 2;     // start of digit 5 (after digit 1)

  node.compute();

  auto* scattered = pool.get_data<int>(scattered_idx);
  EXPECT_EQ((*scattered)[0], 5);
  EXPECT_EQ((*scattered)[1], 5);
  EXPECT_EQ((*scattered)[2], 10);
  EXPECT_EQ((*scattered)[3], 20);
}

TEST(TestNodesTest, MultiPassRadixSort) {
  DatumPool pool;
  auto source_idx = pool.allocate_datum<int>();
  auto sorted_idx = pool.allocate_datum<int>();
  auto node = RadixSortNode(source_idx, sorted_idx, 4);
  node.set_pool(pool);

  auto* source = pool.get_data<int>(source_idx);
  source->push_back(42);
  source->push_back(17);
  source->push_back(99);
  source->push_back(3);
  source->push_back(23);

  node.compute();

  auto* sorted = pool.get_data<int>(sorted_idx);
  EXPECT_EQ(sorted->size(), 5);
  // Verify sorted order
  EXPECT_EQ((*sorted)[0], 3);
  EXPECT_EQ((*sorted)[1], 17);
  EXPECT_EQ((*sorted)[2], 23);
  EXPECT_EQ((*sorted)[3], 42);
  EXPECT_EQ((*sorted)[4], 99);
}
```

---

## Benefits of Node-Based Design

1. **Modularity**: Each stage is independent and testable
2. **Debugging**: Can inspect intermediate outputs (histograms, prefix sums, etc.)
3. **Reusability**: Digit extraction, histogram, and prefix sum nodes can be reused
4. **Optimization**: Can swap out individual stages (e.g., use different histogram algorithms)
5. **Visualization**: Can create debug renderers for intermediate steps
6. **Verification**: Can check that scatter/gather positions are correct
