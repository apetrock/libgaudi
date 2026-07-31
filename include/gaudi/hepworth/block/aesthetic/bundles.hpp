#ifndef __HEP_ROD_AESTHETIC_BUNDLES__
#define __HEP_ROD_AESTHETIC_BUNDLES__

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/hepworth/block/aesthetic/rod_aesthetic_constraints_init.hpp"
#include "gaudi/hepworth/block/aesthetic/symmetric_bezier_constraint.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

template <size_t... Is>
inline constraint_recompute_fn make_rod_symmetric_bezier_recompute(
    asawa::rod::rod::ptr rod, real w) {
  return [rod, w](solver_context &ctx) {
    init_symmetric_bezier(*rod, ctx.constraints, w, select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_uniform_step_recompute(
    asawa::rod::rod::ptr rod, real w, real induced_twist = 0.0) {
  return [rod, w, induced_twist](solver_context &ctx) {
    init_uniform_step(*rod, ctx.constraints, w, select_blocks<Is...>(ctx),
                      induced_twist);
  };
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_helicity_recompute(asawa::rod::rod::ptr rod,
                                                           real w) {
  return [rod, w](solver_context &ctx) {
    if (w <= 0.0)
      return;
    init_helicity(*rod, ctx.constraints, w, select_blocks<Is...>(ctx));
  };
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __HEP_ROD_AESTHETIC_BUNDLES__
