#ifndef __GAUDI_HEPWORTH_CONSTRAINT_BUNDLES__
#define __GAUDI_HEPWORTH_CONSTRAINT_BUNDLES__

#include <memory>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/hepworth/block/generic_constraints_init.hpp"
#include "gaudi/hepworth/block/rod_constraints_init.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/rod_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

// Solver-step wrappers around init_* in shell/rod constraints_init.
// Keep init there for direct use (knotted_surface, tests); bundles route
// joint-matrix blocks via select_blocks<Is...>(ctx).

template <size_t... Is>
inline constraint_recompute_fn make_shell_stretch_recompute(shell_position_block::ptr shell,
                                                            real w) {
  return [shell, w](solver_context &ctx) {
    init_triangle_strain(*shell->mesh, ctx.constraints, shell->xs->get(), w,
                         select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_recompute_fn make_shell_bending_recompute(shell_position_block::ptr shell,
                                                            real w) {
  return [shell, w](solver_context &ctx) {
    init_bending(*shell->mesh, ctx.constraints, shell->xs->get(), w,
                 select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_bundle make_shell_physics_bundle(shell_position_block::ptr shell,
                                                   real stretch_w,
                                                   real bending_w) {
  return {make_shell_stretch_recompute<Is...>(shell, stretch_w),
          make_shell_bending_recompute<Is...>(shell, bending_w)};
}

template <size_t... Is>
inline constraint_recompute_fn make_shell_laplacian_recompute(
    shell_position_block::ptr shell, real w, laplacian_mode mode,
    laplacian_stencil stencil) {
  return [shell, w, mode, stencil](solver_context &ctx) {
    init_laplacian(*shell->mesh, ctx.constraints, shell->xs->get(), mode, stencil, w,
                   select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_bundle make_shell_laplacian_bundle(
    shell_position_block::ptr shell, real w, laplacian_mode mode,
    laplacian_stencil stencil = laplacian_stencil::cotan) {
  return {make_shell_laplacian_recompute<Is...>(shell, w, mode, stencil)};
}

template <size_t... Is>
inline constraint_recompute_fn make_shell_area_recompute(shell_position_block::ptr shell,
                                                         real w, area_mode mode) {
  return [shell, w, mode](solver_context &ctx) {
    init_area(*shell->mesh, ctx.constraints, shell->xs->get(), w,
              select_blocks<Is...>(ctx), mode);
  };
}

template <size_t... Is>
inline constraint_bundle make_shell_area_bundle(shell_position_block::ptr shell, real w,
                                                area_mode mode = area_mode::rest) {
  return {make_shell_area_recompute<Is...>(shell, w, mode)};
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_pin_recompute(rod_position_block::ptr rod,
                                                      real w) {
  return [rod, w](solver_context &ctx) {
    init_pinned(*rod->rod, ctx.constraints, rod->rod->x(), w,
                select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_bundle make_rod_pin_bundle(rod_position_block::ptr rod, real w) {
  return {make_rod_pin_recompute<Is...>(rod, w)};
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_stretch_shear_recompute(asawa::rod::rod::ptr rod,
                                                                real w) {
  return [rod, w](solver_context &ctx) {
    init_stretch_shear(*rod, ctx.constraints, rod->l0(), w,
                       select_blocks<Is...>(ctx));
  };
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_bend_twist_recompute(asawa::rod::rod::ptr rod,
                                                             real w) {
  return [rod, w](solver_context &ctx) {
    init_bend_twist(*rod, ctx.constraints, w, select_blocks<Is...>(ctx), false);
  };
}

template <size_t... Is>
inline constraint_recompute_fn make_rod_collisions_recompute(
    asawa::rod::rod::ptr rod, asawa::rod::dynamic::ptr dynamic, real w) {
  static_assert(sizeof...(Is) == 1, "rod collisions use one position block twice");
  return [rod, dynamic, w](solver_context &ctx) {
    auto blocks = select_blocks<Is...>(ctx);
    init_collisions(*rod, *dynamic, ctx.constraints, w, {blocks[0], blocks[0]});
  };
}

// IPos then IQuat: fans into the right sub-packs for stretch / bend / collision.
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
