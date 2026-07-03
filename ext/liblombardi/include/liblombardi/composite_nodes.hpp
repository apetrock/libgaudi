#pragma once

#include "graph_context.hpp"
#include "test_nodes.hpp"

namespace liblombardi {

// Shared-context composite: input -> MapPlusOne -> MapPlusOne -> output.
//
// Does NOT own the graph context. Builds its children into the parent context
// at construction, owns only the child node handles, and forwards boundary
// ports so the group links like a single node. The composite itself is not
// scheduled; the parent graph runs the children inline.
class map_plus_two_group {
public:
  explicit map_plus_two_group(GraphContext &ctx)
      : _first(ctx.create_node<test_nodes::MapPlusOneNode>()),
        _second(ctx.create_node<test_nodes::MapPlusOneNode>()) {
    ctx.link(_first->output(), _second->input());
  }

  PortRef<test_nodes::MapPlusOneNode, test_nodes::MapPlusOneNode::InputPortDef>
  input() {
    return _first->input();
  }

  PortRef<test_nodes::MapPlusOneNode, test_nodes::MapPlusOneNode::OutputPortDef>
  output() {
    return _second->output();
  }

  std::shared_ptr<test_nodes::vec_int_datum> output_datum() {
    return _second->get_datum<test_nodes::MapPlusOneNode::OutputPortDef>();
  }

private:
  test_nodes::MapPlusOneNode::ptr _first;
  test_nodes::MapPlusOneNode::ptr _second;
};

} // namespace liblombardi
