#pragma once

#include <vector>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>
#include <type_traits>
#include <cstdlib>
#include <utility>

#include "graph_context_base.hpp"
#include "node_base.hpp"
#include "port_def.hpp"

namespace liblombardi {

/// GraphContext - manages nodes and datum pool for a data-flow graph
/// Inherits from GraphContextBase which provides datum pool management
class GraphContext : public GraphContextBase {
public:
    GraphContext() = default;
    ~GraphContext() = default;

    /// Create a node (user-owned)
    template <typename TNode, typename... Args>
    std::shared_ptr<TNode> create_node(Args&&... args) {
        // Create the node - TNode should be a Node specialization with specific port types
        auto node = std::make_shared<TNode>(std::forward<Args>(args)...);

        // Set the context pointer
        node->set_context(this);

        // Assign a stable id (its index in the node list) so the scheduler can
        // map raw pointers back to slots without searching.
        node->set_id(_nodes.size());

        // Add to nodes list using NodeBase pointer
        _nodes.push_back(node);

        // Topology changed: the memoized schedule is stale.
        _dirty = true;

        return node;
    }

    /// Link two ports together via statically-typed PortRef handles.
    /// Allocates a shared datum if the producer port is not yet bound.
    template <typename FromRef, typename ToRef>
    void link(FromRef from, ToRef to) {
        static_assert(std::is_same_v<typename FromRef::port_def::buffer_type,
                                     typename ToRef::port_def::buffer_type>,
                      "Cannot link ports of different buffer types");

        using BufferType = typename FromRef::port_def::buffer_type;
        const uint from_port = port_index_from_def_v<typename FromRef::port_def>;
        const uint to_port = port_index_from_def_v<typename ToRef::port_def>;

        datum_index_t global_id;
        if (auto existing = binding_for(from.node.id(), from_port)) {
            global_id = *existing;
        } else {
            global_id = allocate_datum<BufferType>();
            bind_port(from.node.id(), from_port, global_id);
        }

        bind_port(to.node.id(), to_port, global_id);
        record_edge(from.node, to.node);
    }

    /// Record a producer-consumer edge for scheduling (no port wiring).
    void record_edge(NodeBase& from, NodeBase& to) {
        const auto edge = std::make_pair(from.id(), to.id());
        if (std::find(_edges.begin(), _edges.end(), edge) != _edges.end()) {
            return;
        }
        _edges.push_back(edge);
        _dirty = true;
    }

    /// Clear all recorded edges (e.g. before rewiring a render graph).
    void clear_edges() {
        if (_edges.empty()) {
            return;
        }
        _edges.clear();
        _dirty = true;
    }

    /// Force the schedule to be rebuilt on the next run().
    void mark_dirty() { _dirty = true; }

    const std::vector<std::pair<size_t, size_t>>& edges() const { return _edges; }

    /// Get a datum by global ID (returns shared_ptr to any Datum subclass)
    std::shared_ptr<Datum> get_datum(datum_index_t global_id) const {
        return GraphContextBase::get_datum(global_id);
    }

    /// Get all nodes in the graph
    const std::vector<std::shared_ptr<NodeBase>>& nodes() const {
        return _nodes;
    }

    /// Run the graph - executes nodes in topological order.
    /// The topological sort is memoized: it is only recomputed when the graph
    /// is marked dirty (a node is added or a link is made). The common case of
    /// running the same graph repeatedly (e.g. per frame) is just a walk over
    /// the cached schedule plus one virtual compute() call per node.
    void run() {
        if (_dirty) {
            rebuild_schedule();
            _dirty = false;
        }

        for (NodeBase* node : _schedule) {
            node->compute();
        }
    }

    /// Recompute and cache the execution order via topological sort.
    /// Uses explicit edges recorded by link() / record_edge().
    virtual void rebuild_schedule() {
        const size_t n = _nodes.size();
        _schedule.clear();
        if (n == 0) {
            return;
        }

        std::vector<int> indegree(n, 0);
        std::vector<std::vector<size_t>> adjacency(n);

        for (const auto& edge : _edges) {
            const size_t from = edge.first;
            const size_t to = edge.second;
            if (from >= n || to >= n) {
                continue;
            }
            adjacency[from].push_back(to);
            indegree[to]++;
        }

        std::vector<size_t> ready;
        for (size_t i = 0; i < n; ++i) {
            if (indegree[i] == 0) {
                ready.push_back(i);
            }
        }

        while (!ready.empty()) {
            const size_t current = ready.back();
            ready.pop_back();
            _schedule.push_back(_nodes[current].get());

            for (size_t next : adjacency[current]) {
                if (--indegree[next] == 0) {
                    ready.push_back(next);
                }
            }
        }

        for (size_t i = 0; i < n; ++i) {
            if (indegree[i] > 0) {
                continue;
            }
            NodeBase* node = _nodes[i].get();
            if (std::find(_schedule.begin(), _schedule.end(), node) == _schedule.end()) {
                _schedule.push_back(node);
            }
        }
    }

    /// Get number of datums
    size_t datum_count() const {
        return GraphContextBase::datum_count();
    }

protected:
    const std::vector<NodeBase*>& execution_schedule() const { return _schedule; }
    bool execution_dirty() const { return _dirty; }
    void clear_execution_dirty() { _dirty = false; }

    /// Nodes in the graph
    std::vector<std::shared_ptr<NodeBase>> _nodes;

    /// Cached topological execution order, rebuilt only when _dirty.
    std::vector<NodeBase*> _schedule;

    /// Explicit producer-consumer edges (node id pairs).
    std::vector<std::pair<size_t, size_t>> _edges;

    /// Set when topology changes; triggers a schedule rebuild on next run().
    bool _dirty = true;

private:
};

} // namespace liblombardi
