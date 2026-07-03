#include "gaudi/test/test.hpp"

#include <vector>

#include "liblombardi/composite_nodes.hpp"
#include "liblombardi/graph_context.hpp"
#include "liblombardi/junction_node.hpp"
#include "liblombardi/test_nodes.hpp"

using gaudi::index_t;

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

    // Wire graph: generator -> map+1 -> add
    ctx.link(gen0->output(), map0->input());
    ctx.link(gen1->output(), map1->input());
    ctx.link(map0->output(), add->input0());
    ctx.link(map1->output(), add->input1());

    // Run the graph
    ctx.run();

    // Get final output datum from add node's output
    // Expected: map0 produces [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    //          map1 produces [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    //          add adds them: [2+2, 3+3, 4+4, 5+5, 6+6, 7+7, 8+8, 9+9, 10+10, 11+11]
    //          = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    auto add_result = add->get_datum<AddNode::OutputPortDef>();

    // Verify output size
    GAUDI_ASSERT(add_result->size() == 10);

    // Debug: print actual and expected values
    std::cout << "\n=== BasicGraphExecution Test ===\n";
    std::cout << "Output size: " << add_result->size() << "\n";
    
    // Check if output is all zeros
    bool all_zero = true;
    for (size_t i = 0; i < add_result->size(); ++i) {
        if (add_result->data[i] != 0) {
            all_zero = false;
            break;
        }
    }
    std::cout << "All zeros: " << (all_zero ? "yes" : "no") << "\n";
    
    int expected[] = {4, 6, 8, 10, 12, 14, 16, 18, 20, 22};
    for (size_t i = 0; i < add_result->size(); ++i) {
        std::cout << "  [" << i << "] expected=" << expected[i] << " actual=" << add_result->data[i] << "\n";
        if (add_result->data[i] != expected[i]) {
            std::cout << "    ERROR: mismatch at index " << i << "\n";
        }
        GAUDI_ASSERT(add_result->data[i] == expected[i]);
    }
}

GAUDI_TEST(NodeExecutionOrder) {
    std::cout << "\n=== NodeExecutionOrder Test ===\n";
    
    GraphContext ctx;

    // Create nodes
    auto gen0 = ctx.create_node<GeneratorNode>(5);
    auto map0 = ctx.create_node<MapPlusOneNode>();
    auto add = ctx.create_node<AddNode>();

    // Wire graph: generator -> map+1 -> add
    // This creates: [1, 2, 3, 4, 5] -> [2, 3, 4, 5, 6] -> [4, 6, 8, 10, 12]
    // Note: link(from, to) where from is producer, to is consumer
    ctx.link(gen0->output(), map0->input());
    ctx.link(map0->output(), add->input0());
    ctx.link(map0->output(), add->input1());  // Use same output for both inputs

    // Run and verify add node has correct output
    ctx.run();

    // Expected: [2, 3, 4, 5, 6] + [2, 3, 4, 5, 6] = [4, 6, 8, 10, 12]
    auto add_result = add->get_datum<AddNode::OutputPortDef>();
    GAUDI_ASSERT(add_result->size() == 5);

    std::cout << "Output size: " << add_result->size() << "\n";
    
    // Check if output is all zeros
    bool all_zero = true;
    for (size_t i = 0; i < add_result->size(); ++i) {
        if (add_result->data[i] != 0) {
            all_zero = false;
            break;
        }
    }
    std::cout << "All zeros: " << (all_zero ? "yes" : "no") << "\n";
    
    int expected[] = {4, 6, 8, 10, 12};
    for (size_t i = 0; i < add_result->size(); ++i) {
        std::cout << "  [" << i << "] expected=" << expected[i] << " actual=" << add_result->data[i] << "\n";
        if (add_result->data[i] != expected[i]) {
            std::cout << "    ERROR: mismatch at index " << i << "\n";
        }
        GAUDI_ASSERT(add_result->data[i] == expected[i]);
    }
}

GAUDI_TEST(PoolAllocation) {
    GraphContext ctx;

    // Allocate some datums directly
    auto idx1 = ctx.allocate_datum<vec_int_datum>(10);
    auto idx2 = ctx.allocate_datum<vec_int_datum>(5);

    // Verify datums exist
    GAUDI_ASSERT(ctx.datum_count() == 2);

    // Get datum and verify size
    auto datum1 = ctx.get_datum(idx1);
    auto datum2 = ctx.get_datum(idx2);

    GAUDI_ASSERT(datum1->size() == 10);
    GAUDI_ASSERT(datum2->size() == 5);
}

namespace {
std::vector<int> g_execution_order;
}

GAUDI_TEST(ExplicitEdgeScheduling) {
    struct order_node : Node {
        int id = 0;
        void compute() override { g_execution_order.push_back(id); }
        uint port_count() const override { return 0; }
    };

    g_execution_order.clear();
    GraphContext ctx;

    auto a = ctx.create_node<order_node>();
    a->id = 1;
    auto b = ctx.create_node<order_node>();
    b->id = 2;
    auto c = ctx.create_node<order_node>();
    c->id = 3;

    ctx.record_edge(*a, *b);
    ctx.record_edge(*b, *c);
    ctx.run();

    GAUDI_ASSERT(g_execution_order.size() == 3);
    GAUDI_ASSERT(g_execution_order[0] == 1);
    GAUDI_ASSERT(g_execution_order[1] == 2);
    GAUDI_ASSERT(g_execution_order[2] == 3);
}

GAUDI_TEST(MapPlusTwoGroup) {
    GraphContext ctx;
    auto generator = ctx.create_node<GeneratorNode>(5);  // [1, 2, 3, 4, 5]
    map_plus_two_group group(ctx);

    GAUDI_ASSERT(ctx.nodes().size() == 3);

    ctx.link(generator->output(), group.input());
    ctx.run();

    auto output = group.output_datum();
    GAUDI_ASSERT(output->size() == 5);
    for (size_t i = 0; i < output->size(); ++i) {
        GAUDI_ASSERT(output->data[i] == static_cast<int>(i + 3));
    }
}

GAUDI_TEST(ErrorHandling) {
    GraphContext ctx;

    const datum_index_t invalid_id = 999;

    try {
        ctx.get_datum(invalid_id);
        GAUDI_ASSERT(false);
    } catch (const std::out_of_range&) {
        // Expected
    }
}

GAUDI_TEST(JunctionNodeTwoInputParity) {
    GraphContext ctx;
    auto gen0 = ctx.create_node<GeneratorNode>(5);
    auto gen1 = ctx.create_node<GeneratorNode>(5);
    auto add = ctx.create_node<junction_node<2, vec_int_datum, add_op<int>>>();

    ctx.link(gen0->output(), add->input<0>());
    ctx.link(gen1->output(), add->input<1>());
    ctx.run();

    auto out = add->get_datum<junction_node<2, vec_int_datum, add_op<int>>::OutputPortDef>();
    GAUDI_ASSERT(out->size() == 5);
    for (size_t i = 0; i < out->size(); ++i) {
        GAUDI_ASSERT(out->data[i] == static_cast<int>(2 * (i + 1)));
    }
}

GAUDI_TEST(JunctionNodeThreeInputs) {
    GraphContext ctx;
    auto gen0 = ctx.create_node<GeneratorNode>(4);
    auto gen1 = ctx.create_node<GeneratorNode>(4);
    auto gen2 = ctx.create_node<GeneratorNode>(4);
    auto add = ctx.create_node<junction_node<3, vec_int_datum, add_op<int>>>();

    ctx.link(gen0->output(), add->input<0>());
    ctx.link(gen1->output(), add->input<1>());
    ctx.link(gen2->output(), add->input<2>());
    ctx.run();

    auto out = add->get_datum<junction_node<3, vec_int_datum, add_op<int>>::OutputPortDef>();
    GAUDI_ASSERT(out->size() == 4);
    for (size_t i = 0; i < out->size(); ++i) {
        GAUDI_ASSERT(out->data[i] == static_cast<int>(3 * (i + 1)));
    }
}
