#pragma once

#include <vector>
#include <memory>
#include <unordered_map>
#include <utility>
#include <optional>
#include <functional>
#include <stdexcept>

#include "datum_pool.hpp"
#include "port_def.hpp"

namespace liblombardi {

struct PortBindingHash {
    size_t operator()(const std::pair<size_t, uint>& key) const {
        return std::hash<size_t>{}(key.first) ^ (std::hash<uint>{}(key.second) << 1);
    }
};

/// Base class for all graph contexts
/// Provides datum pool management without node-specific details
class GraphContextBase {
public:
    GraphContextBase() = default;
    virtual ~GraphContextBase() = default;

    /// Allocate a datum of type T (templated for type-safety)
    template <typename T>
    datum_index_t allocate_datum(size_t size = 0) {
        return _global_data.allocate_datum<T>(size);
    }

    /// Get a datum by global ID (returns shared_ptr to any Datum subclass)
    std::shared_ptr<Datum> get_datum(datum_index_t global_id) const {
        return _global_data.get_datum(global_id);
    }

    /// Get data from a datum of type T (throws if type mismatch)
    template <typename T>
    std::vector<T>* get_data(datum_index_t global_id) {
        return _global_data.get_data<T>(global_id);
    }

    template <typename T>
    const std::vector<T>* get_data(datum_index_t global_id) const {
        return _global_data.get_data<T>(global_id);
    }

    /// Get number of datums in pool
    size_t datum_count() const {
        return _global_data.size();
    }

    /// Look up an existing port binding without allocating.
    std::optional<datum_index_t> binding_for(size_t node_id, uint port_id) const {
        const auto it = _port_bindings.find(std::make_pair(node_id, port_id));
        if (it == _port_bindings.end()) {
            return std::nullopt;
        }
        return it->second;
    }

    /// Bind a port to a global datum id.
    void bind_port(size_t node_id, uint port_id, datum_index_t global_id) {
        _port_bindings[std::make_pair(node_id, port_id)] = global_id;
    }

    /// Ensure a port has a datum, allocating one of type T if needed.
    template <typename T>
    datum_index_t ensure_datum_for(size_t node_id, uint port_id) {
        if (auto existing = binding_for(node_id, port_id)) {
            return *existing;
        }
        const datum_index_t global_id = allocate_datum<T>();
        bind_port(node_id, port_id, global_id);
        return global_id;
    }

    /// Clear all port bindings (e.g. before rewiring a graph).
    void clear_bindings() {
        _port_bindings.clear();
    }

protected:
    /// Global pool of datums
    DatumPool _global_data;

    /// Maps (node_id, port_id) -> global datum id.
    std::unordered_map<std::pair<size_t, uint>, datum_index_t, PortBindingHash> _port_bindings;
};

} // namespace liblombardi
