# Construct / Init Split — Full Migration Plan

**Status:** Phase 0–4 complete (2026-08-31). Phase 5–6 pending.
**Stashed (orthogonal):** [`docs/rbf-mode-basis-redesign.md`](rbf-mode-basis-redesign.md) — active RBF plan; prerequisite split done.  
**Defer:** persistent constraint pools (perf-only), RBF mode-basis redesign

---

## Problem

Projection constraints naturally have **indices + blocks + live `q`**. Rest data (`d_rest`, `Q_r`, `_vg`, `_O`, pierce spokes, dipole normals) is baked in constructors fed by parallel geometry:

| Init path | Projection path |
|-----------|-----------------|
| `const vector<vec3>& x` from `shell_xs->get()` / `asawa::get_vec_data` | `_blocks[0]->get_vec3(vi, q)` |
| `R.__u` / `R.__x` read directly from asawa rod | `_blocks[k]->get_quat(ii, q)` |

Even when numbers match at t=0, **code paths diverge** — different helpers, no guarantee of shared evaluation (slab field `Q_r` vs `Q_k` is the worst case).

## Goal

1. **Constructor** — topology + config only (stencil ids, weights, `w`, `blocks`, scalar params like `l0` where not geometric).
2. **`init_rest(q_rest, ctx)`** — bake rest pose via the **same block accessors** as `project(q)`.
3. **Solver hook** — after constraint *create* (recompute), pack `q_rest` from blocks once, call `init_rest` on every constraint before `solver.step`.

Pierce / slab topology: **frozen at first successful `init_rest`** (current semantics).

---

## Is this worthwhile?

**Yes.** This is the correct structural fix, not a slab-prism band-aid.

| Benefit | Why it matters |
|---------|----------------|
| Single geometry API | Rest init and projection share `get_vec3` / `get_quat` |
| Fixes field identity class | `Q_r` centers/coeffs from same stencil gather as live refit |
| Safer multi-block coupling | Shell+rod demos (knotted surface) stop mixing raw `__x` and blocks |
| Testable | Parity: `init_rest(pack_q(blocks))` ≡ old vec path at t=0 |
| Low runtime cost | One `pack_q` + one `init_rest` pass per step (same order as today’s full rebuild) |

**Risk:** Medium surface area, but **mechanical** — same pattern repeated per constraint family. Pilot on slab prism, then batch shell/rod. Topology-only constraints need no change.

**Not in scope:** changing variational energies, persistent constraint caching, RBF basis redesign.

---

## Demo-driven scope

Constraints touched by production demos (migrate these):

### `knotted_surface.hpp` (shell + rod + coupling)

| Recompute | Constraint / init | Rest baked in ctor? |
|-----------|-------------------|---------------------|
| weld | `edge_edge_weld` | No (live-only) |
| dipole | `dipole_weld`, `triangle_dipole_tunneling` | **Yes** (`Nr`, `Nm`) |
| rod strain | `stretch_shear` | Scalar `l0` only |
| rod kink | `min_kink` | No |
| rod bend/twist | `bend`, `twist` | **Yes** (`R.__u` → `_O`) |
| helicity | removed | Was optional knotted-surface feature (never enabled) |
| angle | `angle` | Axis `z`, `phi` only |
| pin | `pinned` | **Yes** (`x[iv]`) |
| rod collision | `rod_collision` | No (pair topology) |
| shell strain | `triangle_strain` | **Yes** (`_rest`) |
| shell bend | `bending` | **Yes** (`_vg`) |
| willmore | `edge_willmore` | Partial (uses mesh const x) |
| area | `area` | **Yes** (rest area) |
| shell collision | `pnt_tri_collision` | **Yes** (`_Nv`) |

### `braid_planar_demo.hpp` / rod graph demos

| Bundle / recompute | Constraints | Rest bake? |
|--------------------|-------------|------------|
| `make_rod_edge_stretch_recompute` | `edge_stretch` | `l0` scalar |
| `make_rod_smooth_recompute` | `smooth` | No |
| `make_rod_min_kink_recompute` | `min_kink` | No |
| `make_rod_collisions_recompute` | `rod_collision` | No |

### `rod_constraints_solver.hpp`, `torque_loop.hpp`

| Bundle / recompute | Constraints | Rest bake? |
|--------------------|-------------|------------|
| `make_rod_physics_bundle` | stretch_shear, edge_stretch, bend, twist, collision | bend/twist via `init_rest` |
| smooth / min_kink / squad_smooth | squad_smooth: quat fairing on S³ (not min-jerk) |

### `cubic_stylization.hpp`

| Bundle | Constraints |
|--------|-------------|
| stretch | `triangle_strain` |
| cubification | `arap_normal_align` (+ external `N_rest`, `N_goal`) |
| pin | `pinned` |

### `slab_prism_bending.hpp`

| Bundle | Constraints |
|--------|-------------|
| `make_shell_slab_prism_bundle` | pierce stencil, slab ARAP, twist, trivector/field divergence |

### `rod_forces.hpp`

**No block solver** — force field only. No migration. Demos that *combine* TP forces + block solver (`braid_planar`) are covered above.

---

## Constraint taxonomy (implementation batches)

### Batch A — Already aligned (no `init_rest` work)

`smooth`, `min_kink`, `straight`, `edge_edge_weld`, `pnt_pnt_weld`, `rod_collision`, `angle` (axis/φ are config, not sampled geometry).

### Batch B — Shell rest from positions

`edge_strain`, `triangle_strain`, `area`, `bending`, `willmore`, `arap`, `laplacian`, `pinned`, `face_normal_flow`, `edge_normal_flow`.

**Pattern:**

```cpp
// ctor: ids, edge_weights, w, blocks
void init_rest(const vecX &q_rest, solver_context &) override {
  one_ring_gather(_blocks[0], _ids, q_rest, _edge_weights, [&](...) {
    _vg = ...;  // or _rest, _A0, _p, etc.
  });
}
```

### Batch C — Rod Cosserat rest from quat block

`bend`, `twist` — today `init_bend` passes full `R.__u`; must use `_blocks[0]->get_quat(i, q_rest)`.

`stretch_shear`, `edge_stretch` — keep `l0` in ctor (material rest length, not a position sample). Optional: recompute `l0` from `q_rest` if policy is geometric rest length.

### Batch D — Slab prism + field

`slab_*`, `shape_field_divergence_strain`, `trivector_divergence_strain`, `relative_plane_twist`.

Remove `x_rest` from `init_slab_prism_constraints` / bundles; pierce + `Q_r` in `init_rest`.

### Batch E — Coupling / collision with stored normals

`dipole_weld`, `triangle_dipole_tunneling`, `pnt_tri_collision`, `point_edge_creep`.

Normals and offset geometry sampled in `init_rest` from shell/rod blocks at `q_rest`.

### Batch F — Removed (2026-08-31)

`helicitiy`, `uniform_step`, `symmetric_bezier` block constraints deleted (unused).
`squad_smooth` retained (quat SQUAD fairing; distinct from `min_kink` position min-jerk).

---

## Infrastructure (Phase 0)

### 1. Base API — `block_constraint.hpp`

```cpp
virtual void init_rest(const vecX &q_rest, solver_context &ctx) {}
virtual bool rest_initialized() const { return _rest_initialized; }

protected:
  bool _rest_initialized = false;
```

Constraints that bake rest set `_rest_initialized = true` at end of `init_rest`. Slab pierce may skip re-init if already frozen.

### 2. Geometry helpers — `block_geometry.hpp`

```cpp
vec3 block_vec3(sim_block &b, index_t vi, const vecX &q);
quat block_quat(sim_block &b, index_t vi, const vecX &q);

void gather_shell_stencil(const sim_block &shell_block,
                          const std::vector<index_t> &ids,
                          const vecX &q, std::vector<vec3> &out);

vecX pack_block_positions(const std::vector<sim_block::ptr> &blocks);
// Full q for multi-block: map_to_x on scratch vecX (or reuse pre-step q before integrate)
```

Shared stencil walkers (`one_ring_iter`) should accept `(block, ids, q)` — refactor shell_constraints.hpp helpers to dual-path: legacy `vector<vec3>& x` delegates to block path for tests.

### 3. Solver hook — `solver_composition.hpp`

After all `recompute` fns populate `ctx.constraints`:

```cpp
vecX q_rest;
pack_q_rest(ctx.blocks, q_rest);
for (auto &c : ctx.constraints)
  c->init_rest(q_rest, ctx);
```

**Rest policy enum** (optional, per-constraint or global):

- `rest_at_step_start` — default for knotted surface / braid (rebuild each step, rest = pre-step q)
- `rest_frozen_once` — slab pierce + `Q_r` after first init
- `rest_external` — cubification `N_goal` stays user-provided datum

---

## Migration phases

### Phase 0 — Infrastructure + tests skeleton

- [x] `init_rest` on `block_constraint`
- [x] `block_geometry.hpp`
- [x] `init_all_rest(q_rest, ctx)` in solver_composition
- [x] Test helper: `shape_field_init_rest_matches_legacy_ctor`

### Phase 1 — Slab prism pilot (highest ROI)

- [x] Split slab constraint ctors (topology only)
- [x] `init_slab_prism_constraints` → create + `init_rest`
- [x] Update `make_shell_slab_prism_bundle` — drop `x_rest` arg
- [x] `slab_prism_smoke` (25 tests) passes
- [ ] Remove RBF fixed-basis band-aids (optional follow-up)

### Phase 2 — Shell Batch B (unblocks knotted surface + cubic)

- [x] `triangle_strain`, `bending`, `area`, `edge_strain`
- [x] `pinned`
- [x] `edge_willmore`, `willmore`, `arap`, `laplacian`
- [x] `arap_normal_align` (positions from block; `N_rest`/`N_goal` from datums at init time)
- [x] `shell_constraints_init.hpp` — create without `x`; solver hook calls `init_rest`
- [x] `bundles.hpp` — `shell_positions_for_ctx` + no `shell->xs->get()` in recompute lambdas

### Phase 3 — Rod Batch C (unblocks knotted surface + all rod demos)

- [x] `bend`, `twist` — lightweight ctor + `init_rest` from quat block
- [x] `init_bend` / `init_twist` / `init_bend_twist` — no longer pass `R.__u`
- [x] `make_rod_physics_bundle` path compiles (torque_loop, braid, rod_constraints)

### Phase 4 — Coupling Batch E (knotted surface completion)

- [x] `dipole_weld`, `triangle_dipole_tunneling` — `Nr`/`Nm` in `init_rest`
- [x] `pnt_tri_collision` — `_Nv` in `init_rest` (uses `ctx.shell_mesh`)
- [x] `knotted_surface.hpp` — sets `ctx.shell_mesh`, `ctx.rod_normals` in recompute lambdas
- [x] `init_rod_shell_weld` — live-only (no change)

### Phase 5 — Cleanup (was aesthetic audit)

- [x] Removed dead aesthetic constraints (`helicitiy`, `uniform_step`, `symmetric_bezier` block)
- [x] Removed legacy `hepworth/rod` `cylinder` helicity init
- [ ] `rod_guided_deformation_with_morph.hpp` direct `init_*` calls audit (if any remain)

### Phase 6 — Cleanup

- [ ] Deprecate `x` parameters on `init_*` (thin wrappers calling block path for legacy tests)
- [ ] Document rest policy per demo in this file
- [ ] Optional: `init_rest` skip if constraint pointer persisted (future perf work)

---

## Call-site checklist

| File | Action |
|------|--------|
| `bundles.hpp` | Remove `shell->xs->get()` / `x_rest` from all recompute lambdas |
| `knotted_surface.hpp` | No direct change if init_* updated; verify dipole cache vs `q_rest` |
| `braid_planar_demo.hpp` | Auto via bundle changes |
| `cubic_stylization.hpp` | Cubification: pass datum ids, sample in init_rest |
| `slab_prism_bending.hpp` | Drop `_x_rest` member / bundle arg |
| `rod_*` demos | Auto via Phase 3 |
| `slab_prism_tests.hpp` | Parity + identity tests |

---

## Safety strategy

1. **Dual-path period** — keep old `create(ids, x, ...)` as deprecated inline forwarding to `create(ids,...)` + manual `init_rest` for tests until parity suite passes.
2. **One batch per PR** — Phase 1 alone is mergeable; Phase 2 can split by constraint family.
3. **Parity tests** — for each migrated constraint: at fixed mesh/rod, assert `project(q)` and rest targets unchanged vs baseline at t=0.
4. **No solver algorithm change** — still clear + rebuild constraints each step; only *where* rest is sampled moves.
5. **CI targets** — `slab_prism_smoke`, knotted-surface smoke (if exists), braid demo compile; run slab tests every batch.

---

## Success criteria

- [x] Zero `vector<vec3>` rest args on constraint ctors in slab + shell + rod Batch B/C (legacy forwarding ctors kept for tests)
- [x] `init_bend_twist` never reads `R.__u`; uses quat block via `q_rest`
- [x] `make_shell_slab_prism_bundle` has no `x_rest` parameter
- [ ] RBF low-order mode-basis redesign — see [`docs/rbf-mode-basis-redesign.md`](rbf-mode-basis-redesign.md)
- [x] All listed demos compile; `slab_prism_smoke` + `gaudi_tests` pass at t=0

---

## Estimated effort

| Phase | Constraints | Demos unblocked |
|-------|-------------|-----------------|
| 0 | infra | — |
| 1 | ~6 slab/field | slab_prism |
| 2 | ~10 shell | knotted surface (shell), cubic |
| 3 | 2 rod (+ init glue) | knotted surface (rod), braid+, torque_loop |
| 4 | ~4 coupling | knotted surface (full) |
| 5 | ~4 aesthetic | rod_constraints module |
| 6 | cleanup | all |

Total: large but ** repetitive** — suitable for 4–6 focused PRs after Phase 0–1 prove the hook.
