#ifndef __GAUDI_HEPWORTH_BLOCK_SOLVER_NODE__
#define __GAUDI_HEPWORTH_BLOCK_SOLVER_NODE__

#include <cstddef>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include "gaudi/duchamp/field_nodes.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

#include "liblombardi/node_base.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

namespace detail {

template <size_t I, typename... DofBlocks>
using slot_block_t = std::tuple_element_t<I, std::tuple<DofBlocks...>>;

template <size_t I, typename... DofBlocks>
using slot_input_datum_t = typename slot_block_t<I, DofBlocks...>::input_port::datum_type;

} // namespace detail

template <typename... DofBlocks>
class block_solver_node : public liblombardi::Node {
public:
  static constexpr size_t N = sizeof...(DofBlocks);
  static_assert(N >= 1, "block_solver_node requires at least one DOF block");
  static_assert(N <= 16, "block_solver_node currently supports up to 16 DOF blocks");

  using ptr = std::shared_ptr<block_solver_node>;
  using config_type = block_solver_config<DofBlocks...>;

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
    Input15
  };

  template <size_t I>
  using InputPortDef =
      liblombardi::PortDef<detail::slot_input_datum_t<I, DofBlocks...>,
                           static_cast<PortId>(I)>;

  block_solver_node() = default;
  explicit block_solver_node(config_type config) : _config(std::move(config)) {}

  void set_config(config_type config) { _config = std::move(config); }
  const config_type &config() const { return _config; }

  /// Mutate outer step without rebuilding the constraint graph.
  void set_dt(real h) { _config.dt = h; }
  real dt() const { return _config.dt; }
  void set_damping(real d) { _config.damping = d; }
  real damping() const { return _config.damping; }
  void set_iterations(int n) { _config.iterations = n; }
  int iterations() const { return _config.iterations; }

  void compute() override { compute_impl(std::make_index_sequence<N>{}); }

  unsigned int port_count() const override { return static_cast<unsigned int>(N); }

  template <size_t I>
  liblombardi::PortRef<block_solver_node, InputPortDef<I>> input_at() {
    static_assert(I < N, "input index out of range");
    return {*this};
  }

private:
  template <size_t... Is>
  void compute_impl(std::index_sequence<Is...>) {
    solver_context ctx;
    prepare_solver_step(_config, ctx);
    (apply_slot<Is>(), ...);
    solve_solver_step(_config, _solver, ctx);
    flush_external_inputs(_config);
  }

  template <size_t I>
  void apply_slot() {
    auto datum = get_datum<InputPortDef<I>>();
    if (!datum) {
      return;
    }
    auto block = std::get<I>(_config.dof_blocks);
    if (block) {
      block->apply_external_input(datum->data());
    }
  }

  config_type _config;
  projection_solver _solver;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_BLOCK_SOLVER_NODE__
