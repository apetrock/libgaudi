#ifndef GAUDI_SMALLVEC_HPP
#define GAUDI_SMALLVEC_HPP

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace gaudi {

/// Stack-first vector (Rust `SmallVec` style): inline storage, spills to heap when full.
template <typename T, std::size_t INLINE = 16> class smallvec {
public:
  using value_type = T;
  using size_type = std::size_t;
  using reference = T &;
  using const_reference = const T &;
  using pointer = T *;
  using const_pointer = const T *;
  using iterator = T *;
  using const_iterator = const T *;

  static constexpr size_type inline_capacity = INLINE;

  smallvec() = default;

  smallvec(std::initializer_list<T> init) {
    reserve(init.size());
    for (const T &v : init)
      push_back(v);
  }

  explicit smallvec(size_type count, const T &value = T()) {
    reserve(count);
    for (size_type i = 0; i < count; ++i)
      push_back(value);
  }

  template <typename InputIt,
            typename = typename std::enable_if<
                !std::is_integral<InputIt>::value>::type>
  smallvec(InputIt first, InputIt last) {
    reserve(static_cast<size_type>(std::distance(first, last)));
    for (; first != last; ++first)
      push_back(*first);
  }

  smallvec(const smallvec &other) { assign_from(other); }

  smallvec(smallvec &&other) noexcept { move_from(std::move(other)); }

  smallvec &operator=(const smallvec &other) {
    if (this != &other) {
      clear();
      assign_from(other);
    }
    return *this;
  }

  smallvec &operator=(smallvec &&other) noexcept {
    if (this != &other) {
      clear();
      move_from(std::move(other));
    }
    return *this;
  }

  ~smallvec() { clear(); }

  void clear() {
    destroy_range(0, size_);
    heap_.clear();
    size_ = 0;
    spilled_ = false;
  }

  bool empty() const { return size_ == 0; }
  size_type size() const { return size_; }
  size_type capacity() const {
    return spilled_ ? heap_.capacity() : inline_capacity;
  }

  bool is_spilled() const { return spilled_; }

  reference operator[](size_type i) { return data()[i]; }
  const_reference operator[](size_type i) const { return data()[i]; }

  reference front() { return (*this)[0]; }
  const_reference front() const { return (*this)[0]; }

  reference back() { return (*this)[size_ - 1]; }
  const_reference back() const { return (*this)[size_ - 1]; }

  pointer data() { return spilled_ ? heap_.data() : inline_ptr(); }
  const_pointer data() const {
    return spilled_ ? heap_.data() : inline_ptr();
  }

  iterator begin() { return data(); }
  iterator end() { return data() + size_; }
  const_iterator begin() const { return data(); }
  const_iterator end() const { return data() + size_; }
  const_iterator cbegin() const { return begin(); }
  const_iterator cend() const { return end(); }

  void reserve(size_type n) {
    if (n <= inline_capacity)
      return;
    if (!spilled_)
      spill_to_heap(n);
    else if (heap_.capacity() < n)
      heap_.reserve(n);
  }

  void push_back(const T &value) {
    if (!spilled_ && size_ < inline_capacity) {
      new (inline_ptr() + size_) T(value);
      ++size_;
      return;
    }
    if (!spilled_)
      spill_to_heap(size_ + 1);
    heap_.push_back(value);
    size_ = heap_.size();
  }

  void push_back(T &&value) {
    if (!spilled_ && size_ < inline_capacity) {
      new (inline_ptr() + size_) T(std::move(value));
      ++size_;
      return;
    }
    if (!spilled_)
      spill_to_heap(size_ + 1);
    heap_.push_back(std::move(value));
    size_ = heap_.size();
  }

  template <typename... Args> reference emplace_back(Args &&...args) {
    if (!spilled_ && size_ < inline_capacity) {
      new (inline_ptr() + size_) T(std::forward<Args>(args)...);
      ++size_;
      return inline_ptr()[size_ - 1];
    }
    if (!spilled_)
      spill_to_heap(size_ + 1);
    heap_.emplace_back(std::forward<Args>(args)...);
    size_ = heap_.size();
    return heap_.back();
  }

  void pop_back() {
    if (size_ == 0)
      return;
    if (spilled_) {
      heap_.pop_back();
      --size_;
      if (size_ <= inline_capacity)
        compact_to_inline();
    } else {
      inline_ptr()[size_ - 1].~T();
      --size_;
    }
  }

  void resize(size_type n, const T &value = T()) {
    if (n < size_) {
      while (size_ > n)
        pop_back();
      return;
    }
    reserve(n);
    while (size_ < n)
      push_back(value);
  }

private:
  size_type size_ = 0;
  bool spilled_ = false;
  alignas(T) unsigned char inline_bytes_[sizeof(T) * INLINE];
  std::vector<T> heap_;

  T *inline_ptr() { return reinterpret_cast<T *>(inline_bytes_); }
  const T *inline_ptr() const {
    return reinterpret_cast<const T *>(inline_bytes_);
  }

  void destroy_range(size_type from, size_type to) {
    if (spilled_) {
      heap_.erase(heap_.begin() + static_cast<std::ptrdiff_t>(from),
                  heap_.begin() + static_cast<std::ptrdiff_t>(to));
      return;
    }
    for (size_type i = from; i < to; ++i)
      inline_ptr()[i].~T();
  }

  void spill_to_heap(size_type want_cap) {
    std::vector<T> tmp;
    tmp.reserve(std::max(want_cap, size_));
    for (size_type i = 0; i < size_; ++i) {
      tmp.push_back(std::move(inline_ptr()[i]));
      inline_ptr()[i].~T();
    }
    heap_ = std::move(tmp);
    spilled_ = true;
    size_ = heap_.size();
  }

  void compact_to_inline() {
    if (!spilled_ || size_ > inline_capacity)
      return;
    for (size_type i = 0; i < size_; ++i)
      new (inline_ptr() + i) T(std::move(heap_[i]));
    heap_.clear();
    heap_.shrink_to_fit();
    spilled_ = false;
  }

  void assign_from(const smallvec &other) {
    clear();
    reserve(other.size());
    for (size_type i = 0; i < other.size(); ++i)
      push_back(other[i]);
  }

  void move_from(smallvec &&other) noexcept {
    size_ = other.size_;
    spilled_ = other.spilled_;
    if (other.spilled_) {
      heap_ = std::move(other.heap_);
    } else {
      for (size_type i = 0; i < size_; ++i) {
        new (inline_ptr() + i) T(std::move(other.inline_ptr()[i]));
        other.inline_ptr()[i].~T();
      }
    }
    other.size_ = 0;
    other.spilled_ = false;
  }
};

} // namespace gaudi

#endif // GAUDI_SMALLVEC_HPP
