#pragma once

#include <memory>
#include <vector>
#include <stdexcept>

#include "datum_pool.hpp"
#include "ports.hpp"
#include "port_def.hpp"
#include "graph_context_base.hpp"

namespace liblombardi {

/// Base class for all nodes in the graph
/// Provides virtual interface for compute()
class NodeBase {
public:
    NodeBase() = default;
    virtual ~NodeBase() = default;

    /// Compute node outputs - must be implemented by derived classes
    virtual void compute() = 0;

    /// Get the number of ports for this node
    /// Returns 0 by default, overridden by derived classes
    virtual uint port_count() const { return 0; }

    /// Get the graph context
    GraphContextBase* context() const { return _context; }

public:
    /// Set the graph context
    void set_context(GraphContextBase* ctx) {
        _context = ctx;
    }

    /// Stable id assigned at registration time (index into the graph's node list).
    /// Lets the scheduler map a node back to its slot without searching.
    void set_id(size_t id) { _id = id; }
    size_t id() const { return _id; }

protected:

private:
    /// Graph context for buffer access
    GraphContextBase* _context = nullptr;

    /// Registration index, set by GraphContext::create_node.
    static constexpr size_t k_invalid_id = static_cast<size_t>(-1);
    size_t _id = k_invalid_id;
};

/// Node class that provides PortDef support and buffer access
/// Inherited by specific node types (GeneratorNode, MapPlusOneNode, etc.)
/// Note: NOT templated on buffer type - nodes are agnostic to buffer types
struct Node : NodeBase {
    Node() = default;
    virtual ~Node() = default;

    /// Compute node outputs - must be implemented by derived classes
    virtual void compute() = 0;

    /// Get the number of ports for this node
    virtual uint port_count() const = 0;

    /// Get datum by port definition.
    /// If the port is not linked to a datum yet (e.g. an output port that has no
    /// downstream consumer), lazily allocate one so it can be written/read.
    template <typename PortDef>
    std::shared_ptr<typename PortDef::buffer_type> get_datum() {
        const uint port = port_index_from_def_v<PortDef>;
        const datum_index_t global_id =
            context()->ensure_datum_for<typename PortDef::buffer_type>(id(), port);
        return std::dynamic_pointer_cast<typename PortDef::buffer_type>(
            context()->get_datum(global_id));
    }

    /// Get datum by port definition (const version)
    template <typename PortDef>
    std::shared_ptr<const typename PortDef::buffer_type> get_datum() const {
        const uint port = port_index_from_def_v<PortDef>;
        const auto global_id = context()->binding_for(id(), port);
        if (!global_id) {
            throw std::out_of_range("Port is not bound to a datum");
        }
        return std::dynamic_pointer_cast<const typename PortDef::buffer_type>(
            context()->get_datum(*global_id));
    }

    /// Set datum value by port definition (for simple types like int, float, etc.)
    template <typename PortDef>
    void set_datum(const std::vector<typename PortDef::buffer_type::value_type>& value) {
        static_assert(port_index_from_def_v<PortDef> < port_count(), "Port ID out of range");
        const uint port = port_index_from_def_v<PortDef>;
        const datum_index_t global_id =
            context()->ensure_datum_for<typename PortDef::buffer_type>(id(), port);
        auto datum = context()->get_datum(global_id);
        auto vec_datum = std::dynamic_pointer_cast<typename PortDef::buffer_type>(datum);
        if (vec_datum) {
            vec_datum->data() = value;
        }
    }

    /// Allocate a new datum of type T
    template <typename T>
    datum_index_t allocate_datum(size_t size = 0) {
        return context()->allocate_datum<T>(size);
    }
};

} // namespace liblombardi
