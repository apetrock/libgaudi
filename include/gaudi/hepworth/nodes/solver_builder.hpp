#ifndef __GAUDI_HEPWORTH_SOLVER_BUILDER__
#define __GAUDI_HEPWORTH_SOLVER_BUILDER__

#include "gaudi/hepworth/block/solver_composition.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

template <typename... DofBlocks>
class block_solver_builder {
public:
  static block_solver_builder create() { return block_solver_builder{}; }

  block_solver_builder &with_blocks(std::shared_ptr<DofBlocks>... blocks) {
    _config.dof_blocks = std::make_tuple(blocks...);
    return *this;
  }

  block_solver_builder &with_presolve(presolve_fn fn) {
    if (fn) {
      _config.presolve.push_back(std::move(fn));
    }
    return *this;
  }

  block_solver_builder &with_recompute(constraint_recompute_fn fn) {
    if (fn) {
      _config.recompute.push_back(std::move(fn));
    }
    return *this;
  }

  block_solver_builder &with_bundle(const constraint_bundle &bundle) {
    for (const auto &fn : bundle) {
      with_recompute(fn);
    }
    return *this;
  }

  block_solver_builder &dt(real h) {
    _config.dt = h;
    return *this;
  }

  block_solver_builder &damping(real d) {
    _config.damping = d;
    return *this;
  }

  block_solver_builder &iterations(int n) {
    _config.iterations = n;
    return *this;
  }

  block_solver_config<DofBlocks...> build() const { return _config; }

private:
  block_solver_config<DofBlocks...> _config;
};

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_SOLVER_BUILDER__
