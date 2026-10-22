# ARP Procrustes Pyramid + BVH MG TDD

## Status

**Phase B (pyramid):** implemented.  
**Phase C (BVH hierarchical TP):** experiments-only — no line connectivity, cannot
run the same rod TP + line search at every level. Prefer topology
[`rod_residual_forest_tdd.md`](rod_residual_forest_tdd.md) for uniform V-cycles.

- Pyramid: [`procrustes_pyramid.hpp`](../include/gaudi/arp/procrustes_pyramid.hpp)
- Points TP: [`point_tangent_point.hpp`](../include/gaudi/calder/point_tangent_point.hpp)
- MG (legacy): [`mg_tangent_point.hpp`](../include/gaudi/duchamp/modules/mg_tangent_point.hpp)

Topology residual forest: **Phase 1 implemented** (`residual_forest.hpp`).

## MG control flow

1. **Ascent** `j = 0 .. max_sim_L`: `extract_at_depth(j)` → Calder `tangent_point_gradient_points_cauchy` → \(R_{\mathrm{sim}}=\mathrm{Kabsch}(x_0\to x_1)\) → \(T\leftarrow R_{\mathrm{sim}}T\) → prolong_x into hierarchy  
2. **Descent**: `prolong_from` at each depth; Kabsch on children before/after prolong to update leaf \(T\)  
3. Optional leaf `tp_propose_curvature_step` refine  

`max_sim_L = max(0, max_L - max_sim_offset)` (default offset 4).

## Inter-level vs sim rotation

- **Inter-level \(R\)**: pyramid residual frames (`δ = Rᵀ(x−c)`).  
- **\(R_{\mathrm{sim}}\)**: computed each level solve from observed \(x_0\to x_1\); transports tangents.

## Tests

See `procrustes_pyramid_tests.hpp`, `mg_tangent_point_tests.hpp`.
