#include <iostream>
#include <cassert>

#include "liblombardi/composite_nodes.hpp"
#include "liblombardi/graph_context.hpp"
#include "liblombardi/test_nodes.hpp"

using namespace liblombardi;
using namespace liblombardi::test_nodes;

/// Test: Basic graph execution with Generator -> MapPlusOne -> Add chain

void test_basic_graph_execution() {
    std::cout << "Running test: BasicGraphExecution..." << std::endl;

    // Create graph context
    auto ctx = std::make_shared<GraphContext>();

    // Create nodes using ctx->create_node
    auto generator = ctx->create_node<GeneratorNode>(10);  // [1, 2, 3, ..., 10]
    auto map_plus_one1 = ctx->create_node<MapPlusOneNode>();
    auto map_plus_one2 = ctx->create_node<MapPlusOneNode>();
    auto add = ctx->create_node<AddNode>();

    // Set up connections using ctx->link
    ctx->link(generator->output(), map_plus_one1->input());
    ctx->link(generator->output(), map_plus_one2->input());
    ctx->link(map_plus_one1->output(), add->input0());
    ctx->link(map_plus_one2->output(), add->input1());

    // Run the graph
    ctx->run();

    // Verify outputs
    // Expected: map_plus_one1 output = [2, 3, 4, ..., 11]
    // Expected: map_plus_one2 output = [2, 3, 4, ..., 11]
    // Expected: add output = [4, 6, 8, ..., 22] (element-wise sum)

    auto output = generator->get_datum<GeneratorNode::OutputPortDef>();
    assert(output->size() == 10);

    std::cout << "Generator output: [";
    for (size_t i = 0; i < output->size(); ++i) {
        std::cout << output->data[i];
        if (i < output->size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // First element: 1
    assert(output->data[0] == 1);
    // Last element: 10
    assert(output->data[9] == 10);

    std::cout << "PASSED: BasicGraphExecution" << std::endl << std::endl;
}

/// Test: Single node execution (generator only)

void test_single_node_execution() {
    std::cout << "Running test: SingleNodeExecution..." << std::endl;

    auto ctx = std::make_shared<GraphContext>();
    auto generator = ctx->create_node<GeneratorNode>(5);  // [1, 2, 3, 4, 5]

    // Run the graph
    ctx->run();

    auto output = generator->get_datum<GeneratorNode::OutputPortDef>();
    assert(output->size() == 5);

    for (size_t i = 0; i < output->size(); ++i) {
        assert(output->data[i] == static_cast<int>(i + 1));
    }

    std::cout << "Generator output: [";
    for (size_t i = 0; i < output->size(); ++i) {
        std::cout << output->data[i];
        if (i < output->size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "PASSED: SingleNodeExecution" << std::endl << std::endl;
}

/// Test: MapPlusOne node alone

void test_map_plus_one_single_node() {
    std::cout << "Running test: MapPlusOneSingleNode..." << std::endl;

    auto ctx = std::make_shared<GraphContext>();
    auto map_plus_one = ctx->create_node<MapPlusOneNode>();

    // Run the graph (MapPlusOneNode has no input, so compute does nothing)
    ctx->run();

    auto output = map_plus_one->get_datum<MapPlusOneNode::OutputPortDef>();
    assert(output->size() == 0);

    std::cout << "MapPlusOne output is empty: " << (output->data.empty() ? "PASS" : "FAIL") << std::endl;

    std::cout << "PASSED: MapPlusOneSingleNode" << std::endl << std::endl;
}

/// Test: Chain with internal intermediate buffer sharing

void test_intermediate_buffer_sharing() {
    std::cout << "Running test: IntermediateBufferSharing..." << std::endl;

    auto ctx = std::make_shared<GraphContext>();

    // Create nodes
    auto generator = ctx->create_node<GeneratorNode>(5);  // [1, 2, 3, 4, 5]
    auto map1 = ctx->create_node<MapPlusOneNode>();
    auto map2 = ctx->create_node<MapPlusOneNode>();
    auto add = ctx->create_node<AddNode>();

    // Connect generator -> map1
    ctx->link(generator->output(), map1->input());

    // Connect generator -> map2 (same output)
    ctx->link(generator->output(), map2->input());

    // Connect map1 -> add.input0
    ctx->link(map1->output(), add->input0());

    // Connect map2 -> add.input1
    ctx->link(map2->output(), add->input1());

    // Run the graph
    ctx->run();

    // Verify: map1 and map2 should share the same buffer (generator's output)
    // map1 output: [2, 3, 4, 5, 6]
    // map2 output: [2, 3, 4, 5, 6]
    // add output: [4, 6, 8, 10, 12]

    auto add_output = add->get_datum<AddNode::OutputPortDef>();
    assert(add_output->size() == 5);

    // Verify element-wise sums: 2+2=4, 3+3=6, 4+4=8, 5+5=10, 6+6=12
    for (size_t i = 0; i < add_output->size(); ++i) {
        assert(add_output->data[i] == static_cast<int>((i + 1) + 1 + (i + 1) + 1));
    }

    std::cout << "Add output: [";
    for (size_t i = 0; i < add_output->size(); ++i) {
        std::cout << add_output->data[i];
        if (i < add_output->size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "PASSED: IntermediateBufferSharing" << std::endl << std::endl;
}

/// Test: two MapPlusOne nodes wrapped in a shared-context composite group

void test_map_plus_two_group() {
    std::cout << "Running test: MapPlusTwoGroup..." << std::endl;

    auto ctx = std::make_shared<GraphContext>();
    auto generator = ctx->create_node<GeneratorNode>(5);  // [1, 2, 3, 4, 5]
    map_plus_two_group group(*ctx);

    // Composite expands into the parent graph; it is not a scheduled node.
    assert(ctx->nodes().size() == 3);

    ctx->link(generator->output(), group.input());
    ctx->run();

    // Each +1 applied twice: [1..5] -> [3, 4, 5, 6, 7]
    auto output = group.output_datum();
    assert(output->size() == 5);
    for (size_t i = 0; i < output->size(); ++i) {
        assert(output->data[i] == static_cast<int>(i + 3));
    }

    std::cout << "Group output: [";
    for (size_t i = 0; i < output->size(); ++i) {
        std::cout << output->data[i];
        if (i < output->size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "PASSED: MapPlusTwoGroup" << std::endl << std::endl;
}

int main() {
    std::cout << "=== liblombardi Test Suite ===" << std::endl << std::endl;

    test_basic_graph_execution();
    test_single_node_execution();
    test_map_plus_one_single_node();
    test_intermediate_buffer_sharing();
    test_map_plus_two_group();

    std::cout << "=== All tests passed! ===" << std::endl;

    return 0;
}
