#pragma once

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "node_base.hpp"

namespace liblombardi {

template <typename T>
struct add_op {
  T operator()(const T &a, const T &b) const { return a + b; }
};

template <int N, typename BufferType, typename CombineOp>
class junction_node : public Node {
public:
  static_assert(N >= 2, "junction_node requires at least 2 inputs");
  static_assert(N <= 16, "junction_node currently supports up to 16 inputs");
  using ptr = std::shared_ptr<junction_node>;

  using value_type = typename BufferType::value_type;

  enum class PortId {
    Input0 = 0,
    Input1,
    Input2,
    Input3,
    Input4,
    Input5,
    Input6,
    Input7,
    Input8,
    Input9,
    Input10,
    Input11,
    Input12,
    Input13,
    Input14,
    Input15,
    Output
  };

  template <size_t I>
  using InputPortDef = PortDef<BufferType, static_cast<PortId>(I)>;
  using OutputPortDef = PortDef<BufferType, PortId::Output>;

  uint port_count() const override { return static_cast<uint>(N + 1); }

  void compute() override { compute_impl(std::make_index_sequence<N>{}); }

  template <size_t I>
  PortRef<junction_node, InputPortDef<I>> input() {
    static_assert(I < static_cast<size_t>(N), "input index out of range");
    return {*this};
  }

  template <size_t I>
  PortRef<junction_node, InputPortDef<I>> input_port() {
    static_assert(I < static_cast<size_t>(N), "input index out of range");
    return {*this};
  }

  PortRef<junction_node, OutputPortDef> output() { return {*this}; }

private:
  template <typename DatumPtr>
  static const std::vector<value_type> &values_from(const DatumPtr &datum) {
    return datum->values();
  }

  template <size_t... Is>
  void compute_impl(std::index_sequence<Is...>) {
    auto out = get_datum<OutputPortDef>();
    std::array<std::shared_ptr<BufferType>, N> inputs{get_datum<InputPortDef<Is>>()...};

    size_t out_size = 0;
    for (const auto &input : inputs) {
      if (input) {
        out_size = std::max(out_size, input->size());
      }
    }

    auto &out_values = out->values();
    out_values.assign(out_size, value_type{});

#ifdef _OPENMP
#pragma omp parallel for if(out_size > 1024)
#endif
    for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(out_size); ++i) {
      bool seeded = false;
      value_type acc{};
      for (const auto &input : inputs) {
        if (!input || static_cast<size_t>(i) >= input->size()) {
          continue;
        }
        const auto &values = values_from(input);
        if (!seeded) {
          acc = values[static_cast<size_t>(i)];
          seeded = true;
        } else {
          acc = _combine(acc, values[static_cast<size_t>(i)]);
        }
      }
      out_values[static_cast<size_t>(i)] = seeded ? acc : value_type{};
    }
  }

  CombineOp _combine{};
};

} // namespace liblombardi
