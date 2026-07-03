#ifndef __GAUDI_HEPWORTH_CONSTRAINT_BUNDLES__
#define __GAUDI_HEPWORTH_CONSTRAINT_BUNDLES__

#include <memory>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

inline constraint_recompute_fn make_shell_stretch_recompute(shell_position_block::ptr shell,
                                                            real w) {
  return [shell, w](solver_context &ctx) {
    init_triangle_strain(*shell->mesh, ctx.constraints, shell->xs->get(), w,
                         ctx.blocks);
  };
}

inline constraint_recompute_fn make_shell_bending_recompute(shell_position_block::ptr shell,
                                                            real w) {
  return [shell, w](solver_context &ctx) {
    init_bending(*shell->mesh, ctx.constraints, shell->xs->get(), w, ctx.blocks);
  };
}

inline constraint_bundle make_shell_physics_bundle(shell_position_block::ptr shell,
                                                   real stretch_w,
                                                   real bending_w) {
  return {make_shell_stretch_recompute(shell, stretch_w),
          make_shell_bending_recompute(shell, bending_w)};
}

template <size_t IPos, size_t IQuat>
inline constraint_recompute_fn make_rod_stretch_shear_recompute(asawa::rod::rod::ptr rod,
                                                                real w) {
  return [rod, w](solver_context &ctx) {
    init_stretch_shear(*rod, ctx.constraints, rod->l0(), w,
                       {ctx.blocks[IPos], ctx.blocks[IQuat]});
  };
}

template <size_t IQuat>
inline constraint_recompute_fn make_rod_bend_twist_recompute(asawa::rod::rod::ptr rod,
                                                             real w) {
  return [rod, w](solver_context &ctx) {
    init_bend_twist(*rod, ctx.constraints, w, {ctx.blocks[IQuat]}, false);
  };
}

template <size_t IPos>
inline constraint_recompute_fn make_rod_collisions_recompute(
    asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic, real w) {
  return [rod, dynamic, w](solver_context &ctx) {
    auto &x = ctx.blocks[IPos];
    init_collisions(*rod, *dynamic, ctx.constraints, w, {x, x});
  };
}

template <size_t IPos, size_t IQuat>
inline constraint_bundle make_rod_physics_bundle(asawa::rod::rod::ptr rod,
                                                 asawa::rod::dynamic::ptr dynamic,
                                                 real stretch_w, real bend_w,
                                                 real collision_w) {
  return {make_rod_stretch_shear_recompute<IPos, IQuat>(rod, stretch_w),
          make_rod_bend_twist_recompute<IQuat>(rod, bend_w),
          make_rod_collisions_recompute<IPos>(rod, dynamic, collision_w)};
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __GAUDI_HEPWORTH_CONSTRAINT_BUNDLES__
