#include <cassert>

#include "gaudi/common.h"
#include "liblombardi/graph_context.hpp"
#include "liblombardi/test_nodes/test_nodes.hpp"

using namespace liblombardi;
using namespace liblombardi::test_nodes;

GAUDI_TEST(BasicGraphExecution) {
    // Create graph context
    GraphContext ctx;

    // Create 2 generator nodes (each produces [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    auto gen0 = ctx.create_node<GeneratorNode>(10);
    auto gen1 = ctx.create_node<GeneratorNode>(10);

    // Create 2 map+1 nodes
    auto map0 = ctx.create_node<MapPlusOneNode>();
    auto map1 = ctx.create_node<MapPlusOneNode>();

    // Create 1 add node
    auto add = ctx.create_node<AddNode>();

    // Get outputs from generators
    auto outputs = ctx.get_outputs(gen0);
    GAUDI_ASSERT(outputs.size() == 1);
    auto gen0_output = outputs[0];

    outputs = ctx.get_outputs(gen1);
    GAUDI_ASSERT(outputs.size() == 1);
    auto gen1_output = outputs[0];

    // Get outputs from map nodes
    outputs = ctx.get_outputs(map0);
    GAUDI_ASSERT(outputs.size() == 1);
    auto map0_output = outputs[0];

    outputs = ctx.get_outputs(map1);
    GAUDI_ASSERT(outputs.size() == 1);
    auto map1_output = outputs[0];

    // Get output from add node
    outputs = ctx.get_outputs(add);
    GAUDI_ASSERT(outputs.size() == 1);
    auto add_output = outputs[0];

    // Wire graph: generator -> map+1 -> add
    ctx.set_input(map0, gen0_output);
    ctx.set_input(map1, gen1_output);
    ctx.set_input(add, map0_output);
    ctx.set_input(add, map1_output);

    // Run the graph
    ctx.run();

    // Get final output data
    // Expected: add(gen0) + add(gen1) = [1+1, 2+1, 3+1, 4+1, 5+1, 6+1, 7+1, 8+1, 9+1, 10+1]
    //          = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11] + [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    //          = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
    const auto& add_result = ctx.pool().get_data<int>(add_output);

    // Verify output size
    GAUDI_ASSERT(add_result.size() == 10);

    // Verify exact values
    int expected[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    for (size_t i = 0; i < add_result.size(); ++i) {
        GAUDI_ASSERT(add_result[i] == expected[i]);
    }
}

GAUDI_TEST(NodeExecutionOrder) {
    GraphContext ctx;

    // Create nodes
    auto gen0 = ctx.create_node<GeneratorNode>(5);
    auto map0 = ctx.create_node<MapPlusOneNode>();
    auto add = ctx.create_node<AddNode>();

    // Wire graph
    auto gen0_output = ctx.get_outputs(gen0)[0];
    auto map0_output = ctx.get_outputs(map0)[0];
    auto add_output = ctx.get_outputs(add)[0];

    ctx.set_input(map0, gen0_output);
    ctx.set_input(add, map0_output);

    // Get execution order
    const auto& execution_order = ctx.execution_order();

    // Verify order: gen0 -> map0 -> add
    GAUDI_ASSERT(execution_order.size() == 3);
    GAUDI_ASSERT(execution_order[0].get() == gen0.get());
    GAUDI_ASSERT(execution_order[1].get() == map0.get());
    GAUDI_ASSERT(execution_order[2].get() == add.get());

    // Run and verify add node has correct output
    ctx.run();

    const auto& add_result = ctx.pool().get_data<int>(add_output);
    GAUDI_ASSERT(add_result.size() == 5);

    // Expected: [1+1, 2+1, 3+1, 4+1, 5+1] + [1, 2, 3, 4, 5] = [3, 4, 5, 6, 7]
    int expected[] = {3, 4, 5, 6, 7};
    for (size_t i = 0; i < add_result.size(); ++i) {
        GAUDI_ASSERT(add_result[i] == expected[i]);
    }
}

GAUDI_TEST(PoolAllocation) {
    GraphContext ctx;

    // Allocate some datums
    auto idx1 = ctx.allocate_datum<int>(10);
    auto idx2 = ctx.allocate_datum<float>(5);

    // Verify datums exist
    GAUDI_ASSERT(ctx.pool().size() == 2);
    GAUDI_ASSERT(ctx.pool().get_data<int>(idx1) != nullptr);
    GAUDI_ASSERT(ctx.pool().get_data<float>(idx2) != nullptr);

    // Get data and verify size
    const auto& int_data = ctx.pool().get_data<int>(idx1);
    const auto& float_data = ctx.pool().get_data<float>(idx2);

    GAUDI_ASSERT(int_data->size() == 10);
    GAUDI_ASSERT(float_data->size() == 5);
}

GAUDI_TEST(ErrorHandling) {
    GraphContext ctx;

    // Try to access invalid datum
    GAUDI_ASSERT_THROW(ctx.pool().get_data<int>(0));
}
