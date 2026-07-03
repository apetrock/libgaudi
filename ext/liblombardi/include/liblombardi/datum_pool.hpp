#pragma once

#include <memory>
#include <vector>
#include <type_traits>

namespace liblombardi {

/// Index into the datum pool
using datum_index_t = size_t;

/// Base datum interface
/// Implement custom datums by inheriting from this class
/// Example:
///   struct MyDatum : public Datum {
///       std::vector<int> _data;
///       MyDatum() { _data.reserve(100); }
///       
///       virtual void resize(size_t size) override {
///           _data.resize(size);
///       }
///       
///       virtual size_t size() const override {
///           return _data.size();
///       }
///       
///       virtual void clear() override {
///           _data.clear();
///       }
///       
///       virtual void* get_data() override {
///           return _data.data();
///       }
///       
///       virtual const void* get_data() const override {
///           return _data.data();
///       }
///   };
class Datum {
public:
    Datum() = default;
    virtual ~Datum() = default;

    virtual void resize(size_t size) = 0;
    virtual size_t size() const = 0;
    virtual void clear() = 0;

    /// Get raw pointer to data
    virtual void* get_data() = 0;
    virtual const void* get_data() const = 0;

private:
    datum_index_t _index;
};

/// Implementation of datum for a specific type
template <typename T>
class DatumImpl : public Datum {
public:
    DatumImpl() = default;
    DatumImpl(const std::vector<T>& data) : _data(data) {}

    virtual void resize(size_t size) override {
        _data.resize(size);
    }

    virtual size_t size() const override {
        return _data.size();
    }

    virtual void clear() override {
        _data.clear();
    }

    virtual void* get_data() override {
        return _data.data();
    }

    virtual const void* get_data() const override {
        return _data.data();
    }

    std::vector<T>& data() { return _data; }
    const std::vector<T>& data() const { return _data; }

private:
    std::vector<T> _data;
};

/// Pool for managing arbitrary datum types
class DatumPool {
public:
    DatumPool() = default;

    /// Insert a datum into the pool
    template <typename T>
    datum_index_t insert_datum(const std::vector<T>& data) {
        auto datum = std::make_shared<T>(data);
        datum_index_t index = _datums.size();
        _datums.push_back(datum);
        return index;
    }

    /// Allocate a new empty datum
    template <typename T>
    datum_index_t allocate_datum(size_t size = 0) {
        auto datum = std::make_shared<T>();
        datum->resize(size);
        datum_index_t index = _datums.size();
        _datums.push_back(datum);
        return index;
    }

    /// Allocate a datum from data
    template <typename T>
    datum_index_t allocate_datum_from_data(const std::vector<T>& data) {
        return insert_datum(data);
    }

    /// Get datum by index
    const std::shared_ptr<Datum>& get_datum(datum_index_t index) const {
        if (index >= _datums.size()) {
            throw std::out_of_range("Invalid datum index");
        }
        return _datums[index];
    }

    std::shared_ptr<Datum>& get_datum(datum_index_t index) {
        if (index >= _datums.size()) {
            throw std::out_of_range("Invalid datum index");
        }
        return _datums[index];
    }

    /// Get data of a specific type
    template <typename T>
    std::vector<T>* get_data(datum_index_t index) {
        auto datum = get_datum(index);
        if constexpr (std::is_same_v<DatumImpl<T>, std::decay_t<decltype(*datum)>>) {
            return &static_cast<DatumImpl<T>*>(datum.get())->data();
        } else {
            throw std::bad_cast();
        }
    }

    template <typename T>
    const std::vector<T>* get_data(datum_index_t index) const {
        auto datum = get_datum(index);
        if constexpr (std::is_same_v<DatumImpl<T>, std::decay_t<decltype(*datum)>>) {
            return &static_cast<const DatumImpl<T>*>(datum.get())->data();
        } else {
            throw std::bad_cast();
        }
    }

    /// Get number of datums in pool
    size_t size() const {
        return _datums.size();
    }

    /// Clear all datums
    void clear() {
        _datums.clear();
    }

private:
    std::vector<std::shared_ptr<Datum>> _datums;
};

} // namespace liblombardi
