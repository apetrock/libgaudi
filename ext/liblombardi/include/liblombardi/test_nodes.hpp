#pragma once

#include "node_base.hpp"

namespace liblombardi {
namespace test_nodes {

// ==================== vec_int_datum ====================
// Custom datum for vector<int>

struct vec_int_datum : public Datum {
public:
    using value_type = int;  // Element type for type checking
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

    std::vector<int>& values() { return data; }
    const std::vector<int>& values() const { return data; }
};

// ==================== GeneratorNode ====================
// Generates a random array of integers [1, 2, 3, ..., N]

class GeneratorNode : public Node {
public:
    using ptr = std::shared_ptr<GeneratorNode>;

    // Port definitions using enum
    enum class PortId {
        Output = 0
    };

    // Port namespaces (not required, but can organize ports)
    struct input {};
    struct output {};

    // PortDef typedefs
    using OutputPortDef = PortDef<vec_int_datum, PortId::Output>;

    GeneratorNode() = default;

    // Initialize with size N (generates [1, 2, 3, ..., N])
    explicit GeneratorNode(size_t size) : _size(size) {}

    // Compute: generate the array
    void compute() override {
        auto datum = get_datum<OutputPortDef>();
        datum->resize(_size);
        for (size_t i = 0; i < _size; ++i) {
            datum->data[i] = static_cast<int>(i + 1);
        }
    }

    uint port_count() const override {
        return 1;  // Has 1 output port
    }

    PortRef<GeneratorNode, OutputPortDef> output() { return {*this}; }

    size_t get_size() const { return _size; }

private:
    size_t _size = 10;  // Default size
};

// ==================== MapPlusOneNode ====================
// Adds +1 to each element in the input array

class MapPlusOneNode : public Node {
public:
    using ptr = std::shared_ptr<MapPlusOneNode>;

    // Port definitions using enum
    enum class PortId {
        Input = 0,
        Output = 1
    };

    // Port namespaces
    struct input {};
    struct output {};

    // PortDef typedefs
    using InputPortDef = PortDef<vec_int_datum, PortId::Input>;
    using OutputPortDef = PortDef<vec_int_datum, PortId::Output>;

    MapPlusOneNode() = default;

    // Compute: map +1 to each element
    void compute() override {
        auto input_datum = get_datum<InputPortDef>();
        auto output_datum = get_datum<OutputPortDef>();
        output_datum->resize(input_datum->size());
        auto& input_data = input_datum->data;
        auto& output_data = output_datum->data;
        for (size_t i = 0; i < input_data.size(); ++i) {
            output_data[i] = input_data[i] + 1;
        }
    }

    uint port_count() const override {
        return 2;  // Has 1 input and 1 output port
    }

    PortRef<MapPlusOneNode, InputPortDef> input() { return {*this}; }
    PortRef<MapPlusOneNode, OutputPortDef> output() { return {*this}; }
};

// ==================== AddNode ====================
// Element-wise addition of two input arrays

class AddNode : public Node {
public:
    using ptr = std::shared_ptr<AddNode>;

    // Port definitions using enum
    enum class PortId {
        Input0 = 0,
        Input1 = 1,
        Output = 2
    };

    // Port namespaces
    struct input {};
    struct output {};

    // PortDef typedefs
    using InputPort0Def = PortDef<vec_int_datum, PortId::Input0>;
    using InputPort1Def = PortDef<vec_int_datum, PortId::Input1>;
    using OutputPortDef = PortDef<vec_int_datum, PortId::Output>;

    AddNode() = default;

    // Compute: add two arrays element-wise
    void compute() override {
        auto input0_datum = get_datum<InputPort0Def>();
        auto input1_datum = get_datum<InputPort1Def>();
        auto output_datum = get_datum<OutputPortDef>();

        size_t size = input0_datum->size();
        if (input1_datum->size() != size) {
            throw std::runtime_error("AddNode: inputs must have the same size");
        }

        output_datum->resize(size);
        auto& input0_data = input0_datum->data;
        auto& input1_data = input1_datum->data;
        auto& output_data = output_datum->data;
        for (size_t i = 0; i < size; ++i) {
            output_data[i] = input0_data[i] + input1_data[i];
        }
    }

    uint port_count() const override {
        return 3;  // Has 2 inputs and 1 output port
    }

    PortRef<AddNode, InputPort0Def> input0() { return {*this}; }
    PortRef<AddNode, InputPort1Def> input1() { return {*this}; }
    PortRef<AddNode, OutputPortDef> output() { return {*this}; }
};

} // namespace test_nodes
} // namespace liblombardi
