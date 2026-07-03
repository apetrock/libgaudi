#pragma once

#include <vector>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

namespace liblombardi {

/// Slab allocator for fixed-size memory blocks
/// O(1) allocate and deallocate
template <typename T>
class SlabAllocator {
public:
    SlabAllocator(size_t block_size = 64) : _block_size(block_size) {}

    /// Allocate a single element
    size_t allocate() {
        if (_free_slots.empty()) {
            // Allocate a new block
            _blocks.push_back(std::vector<T>(_block_size));
            for (size_t i = 0; i < _block_size; ++i) {
                _free_slots.push_back(_blocks.size() - 1);
            }
        }

        // Get a free slot from the most recent block
        size_t block_id = _free_slots.back();
        size_t local_idx = _blocks[block_id].size();
        _blocks[block_id].push_back(T{});
        _free_slots.pop_back();

        return block_id * _block_size + local_idx;
    }

    /// Deallocate a previously allocated element
    void deallocate(size_t idx) {
        if (idx >= total_allocated()) {
            throw std::out_of_range("Invalid slab index");
        }

        size_t block_id = idx / _block_size;
        size_t local_idx = idx % _block_size;

        _blocks[block_id].resize(local_idx);
        _free_slots.push_back(block_id);
    }

    /// Get pointer to element
    T* data(size_t idx) {
        size_t block_id = idx / _block_size;
        size_t local_idx = idx % _block_size;
        return &_blocks[block_id][local_idx];
    }

    const T* data(size_t idx) const {
        size_t block_id = idx / _block_size;
        size_t local_idx = idx % _block_size;
        return &_blocks[block_id][local_idx];
    }

    /// Get total allocated elements
    size_t total_allocated() const {
        size_t total = 0;
        for (const auto& block : _blocks) {
            total += block.size();
        }
        return total;
    }

    /// Get number of blocks
    size_t num_blocks() const {
        return _blocks.size();
    }

    /// Get block size
    size_t block_size() const {
        return _block_size;
    }

    /// Clear all allocations
    void clear() {
        _blocks.clear();
        _free_slots.clear();
    }

private:
    std::vector<std::vector<T>> _blocks;
    std::vector<size_t> _free_slots;
    size_t _block_size;
};

} // namespace liblombardi
