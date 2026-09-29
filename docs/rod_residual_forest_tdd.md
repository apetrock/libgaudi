# Rod Procrustes Residual Forest TDD

## Status

**Phase 1 (forest on fat rod, x-channel):** implemented —
`include/gaudi/asawa/rod/residual_forest.hpp` and
`include/gaudi/test/residual_forest_tests.hpp` (five tests green).

Preferred hierarchy for **uniform V-cycles** (same Calder TP + CFL/line search
at every level via `pickle_rod` → relax → `set`). ARP BVH Procrustes pyramid MG
remains experiments-only (no line connectivity).

**Data-less rod** (old Phase A) deferred until pickle must project arbitrary
datums beyond `x`.

**Related code (today):**
[`rod.hpp`](../include/gaudi/asawa/rod/rod.hpp) (fat attribute member arrays),
[`dynamic.hpp`](../include/gaudi/asawa/rod/dynamic.hpp) (BVH / split-collapse adapter),
[`shell.hpp`](../include/gaudi/asawa/shell/shell.hpp) (topology + datum pool — target parity),
[`slab_prism_operators.hpp`](../include/gaudi/hepworth/constraints/slab/slab_prism_operators.hpp)
(`fit_procrustes_rotation` — Kabsch SO(3) pattern),
[`procrustes_pyramid.hpp`](../include/gaudi/arp/procrustes_pyramid.hpp) (`procrustes_rotation`).

**Header:**
`include/gaudi/asawa/rod/residual_forest.hpp`

---

## Concept

### One-line pitch

Hold a **frozen hierarchy of residuals** on top of an L0 rod while connectivity is
stable. **V-cycle** by projecting absolute state to coarser levels (`pickle`),
relaxing there, then writing back (`set`) which **re-encodes residuals** and
reconstructs finer levels through **k-window Procrustes** frames. On L0 topology
change, **invalidate** and rebuild — do not incrementally repair the forest.

### Why this is uncommon

| Existing tool | What it does | Gap |
|---------------|--------------|-----|
| El Topo / Los Topos | Surface tracking + topology ops | No residual pyramid / restrict–prolongate for MG |
| Progressive meshes / ZBrush multires | Hierarchy + local-frame detail | Typically static connectivity |
| Classical geometric MG | Nested complexes + transfers | Assumes fixed topology |

We want **tracking-style topology changes** interleaved with **MG transfers** on a
frozen connectivity window:

```
topology change → create(forest) → V-cycle (pickle / relax / set)* → topology change → invalidate → …
```

Vertex **positions** may mutate freely while the forest is valid. Resampling,
insert/delete, strand split/merge, or any L0 index reshuffle **invalidates**.

### Adapter pattern (like `rod::dynamic`)

`residual_forest` does **not** own the rod. It holds `rod::ptr`, owns hierarchy side
data (parent maps, residuals, rest snapshots), and answers project/write-back
queries — same relationship `dynamic` has to the rod for BVH / remesh ops.

```cpp
auto F = asawa::rod::residual_forest::create(R);
```

### No freelist / slab node pool

Hierarchy is **SoA + `index_t`**: `parent`, `keep`, `rest`, `delta`. Topology change
triggers full rebuild. A freelist node pool only pays off for incremental
alloc/free of hierarchy nodes — out of scope.

---

## Prephase: data-less rod (deferred)

Was required so `pickle(level, datum_id)` can project arbitrary channels.
**Phase 1** ships x-only on the fat rod (`pickle` / `pickle_rod` / `set` on positions).
Revisit data-less rod when non-x datums need residual transfer.

---

## Residual math

### Reconstruction

At each fine level \(\ell\) with parent level \(\ell+1\):

\[
x_\ell = c_{\ell+1} + R_{\ell+1}\,\delta_\ell
\]

- \(c_{\ell+1}\): parent anchor (parent sample position, or window centroid —
  **lock: use parent live sample** for \(\delta = x - x_p\) consistency).
- \(R_{\ell+1}\): live Kabsch / Procrustes from **rest parent window → live parent
  window** (ARAP local step; not a global Poisson solve in Phase B).
- \(\delta_\ell\): stored **residual** (detail), not an absolute position.

**Procrustes rotates residuals.** Writing full absolute coordinates into the
residual slots is incorrect.

### Encoding

Given absolute fine \(x_\ell\) and parent \(x_{\ell+1}\):

\[
\delta_\ell = R^\top (x_\ell - c_{\ell+1})
\]

At `create` (rest): \(R = I\), so \(\delta_\ell = x_\ell^{\mathrm{rest}} - x_{p(\ell)}^{\mathrm{rest}}\).

On `set(ℓ, abs)`: recompute \(\delta_\ell\) from the new absolute vs the **current
parent** level (live Procrustes inverse), then reconstruct levels \(\ell-1,\ldots,0\)
with existing finer residuals unchanged (unless those levels are also `set`).

### Window (3D)

Two samples fix only an edge → twist about \(T\) is free. Use a **centered odd
window of \(k=5\)** coarse samples about the parent (`k=3` near open ends).

Near-degenerate (collinear / tiny singular values): fall back to \(R = I\) and
parent translation. Primary tests use a **non-planar helix**.

Kabsch: same SVD SO(3) pattern as `fit_procrustes_rotation` (det \(+1\)).

### Collapse / parents (pointerless)

- Every-other keep; dropped fine index → nearest kept coarse index (closed loop wraps).
- Coarse rest sample = kept fine rest sample (**no averaging**) so restrict/inject
  of survivors is exact.
- `parent[ℓ][i]` = `index_t` into level \(\ell+1\) only — no pointers.

---

## API (normative)

Header: `include/gaudi/asawa/rod/residual_forest.hpp`  
Namespace: `gaudi::asawa::rod`

### Lifecycle

```cpp
auto F = residual_forest::create(R);   // snapshot rest x, build levels, encode δ
// … while F->valid() …
F->invalidate();                       // L0 topology changed
F = residual_forest::create(R);         // rebuild
```

Caller-driven invalidate in Phase B; optional hooks from `dynamic` ops in Phase C.

### Pickle (read / project)

```cpp
rod::ptr Rℓ = F->pickle(ℓ);                 // coarsened topology + reconstructed abs datums
std::vector<vec3> xℓ = F->pickle(ℓ, x_id);  // one channel at level ℓ
```

- `pickle` is **read-only**. It reconstructs **absolute** state at \(\ell\) from
  coarser anchors + residuals (`parent abs + R δ` for x).
- `pickle(ℓ)` returns a temporary data-less `rod` suitable for coarse relax
  (corner count = keep count; single open/closed strand linking).

### Set (write-back)

```cpp
F->set(ℓ, x_id, xℓ);   // after simulating at level ℓ
```

Semantics of `set(ℓ, datum_id, absolute_data)`:

1. Treat `absolute_data` as the new absolute state **at level ℓ**.
2. **Re-encode** \(\delta_\ell\) vs the **current parent** level (for x).
3. **Reconstruct** levels \(\ell-1,\ldots,0\) by applying stored finer residuals
   through live Procrustes; write into L0 rod datums.

| Call | Effect on residuals |
|------|---------------------|
| `set(5, x, x5)` after coarse solve | Updates L5 anchors; **keeps** \(\delta_4,\delta_3,\ldots\) so detail follows parent |
| `set(4, x, x4')` after L4 sim | **Re-encodes** \(\delta_4\) from `x4'` vs current L5; finer \(\delta_3,\ldots\) unchanged; L0 rebuilt |

Optional later: `set_residual(ℓ, datum, δ)` to write \(\delta\) directly. Phase B
implements `set` with internal encode; add `set_residual` only if a call site needs it.

**Phase B residual channel:** encode/recompute \(\delta\) for the **x** datum only.
Other datums: inject survivors only (no Procrustes) until Phase C.

### Sketch: V-cycle (caller owns relax)

```text
F = create(R)
relax(R)                              // L0
for ℓ = 1 .. N:
  Rℓ = F->pickle(ℓ)
  relax(Rℓ)                           // or many iters at ℓ == N
  F->set(ℓ, x_id, positions_of(Rℓ))
relax(R)                              // L0 again
```

Phase B ships forest + transfers; **not** a full V-cycle driver or TP/hepworth wiring.

---

## Storage sketch (pointerless)

```text
residual_forest {
  rod::ptr R;                      // L0, not owned
  bool valid;
  bool loop;
  int k;                           // window, default 5
  int n_levels;

  vector<vector<index_t>> parent;  // parent[ℓ][i] -> index in ℓ+1
  vector<vector<index_t>> keep;    // fine indices that become coarse verts
  vector<vector<vec3>> rest_x;     // rest positions per level at create
  vector<vector<vec3>> delta;      // residuals at ℓ vs parent (x-channel)
}
```

Live absolute buffers at intermediate levels may be cached for Procrustes windows
or recomputed on demand from `delta` + coarsest live — implementation choice;
API does not require exposing them.

---

## Phases (implementation checklist — later)

### Phase A — Data-less rod

- [ ] `rod`: topology + datum pool; strip baked `__x,__v,…`
- [ ] `x()` / friends as shims to datums
- [ ] `pack()` remaps all datums uniformly
- [ ] `dynamic` interp on split/collapse generic over datums
- [ ] Braid / TP / knotted / rod tests compile and pass

### Phase B — `residual_forest` adapter

- [ ] `create` / `invalidate` / `valid`
- [ ] Collapse pyramid + parent/keep maps
- [ ] Encode \(\delta\) at create (x datum)
- [ ] `point_window_procrustes` (k-window Kabsch)
- [ ] `pickle(ℓ)` / `pickle(ℓ, datum_id)`
- [ ] `set(ℓ, datum_id, abs)` with residual re-encode + reconstruct below
- [ ] Tests listed below green

### Phase C — Later

- [ ] V-cycle driver / relax integration (TP, hepworth, …)
- [ ] Auto-invalidate inside `dynamic` topology ops
- [ ] Per-datum residual encodings beyond x
- [ ] Mesh residual forest (same policy)
- [ ] Taylor / parallel level approximations
- [ ] Optional `set_residual`

---

## Tests (specify now; implement with Phase A/B)

Fixture: closed **helix** (~64 samples), `k=5`, several collapses (e.g. 3–4 levels).

| Test | Assert |
|------|--------|
| `rod_dataless_split_resizes_datums` | After dynamic edge split, all datum sizes match corner count; `x()` shim works |
| `rod_forest_identity_recover` | `pickle(N,x)` then `set(N,x)` with same data → L0 matches create snapshot (\(\infty\)-norm tiny) |
| `rod_forest_rigid_coarse` | `set(N, x, R xN + t)` → L0 ≈ `R x0 + t`; finer \(\delta\) magnitudes unchanged |
| `rod_forest_midlevel_reencode` | `set(N,…)`, edit `pickle(N-1)`, `set(N-1, x')` → \(\delta_{N-1}\) changes; move N again → L0 tracks **new** residual |
| `rod_forest_invalidate` | `invalidate` ⇒ `!valid`; operations refuse or no-op until `create` |
| `rod_forest_pickle_topology` | `pickle(ℓ)->` active corner/vert count equals keep count; one strand |

---

## Design tie-breakers

- Forest is an **adapter** (`create(rod::ptr)`), not a free-floating point pyramid.
- **Data-less rod before forest** so pickle can project any datum.
- Pointerless index maps; **rebuild** on topology change (no incremental repair).
- `pickle` reads reconstructed **abs**; `set` writes abs → **re-encode \(\delta\)** → reconstruct below.
- Position transfer = Procrustes on residuals; other datums inject-only in Phase B.
- No freelist / slab node pool for hierarchy nodes.
- Thin `vector<vec3>`-only prototype without rod adapter: **rejected / stashed**.

---

## Out of scope (this TDD pass)

- Implementing Phase A/B/C
- Wiring forest into braid / TP solvers
- Mesh forest

---

## Acceptance (for this document)

- API (`create` / `pickle` / `set` / `invalidate`) and residual re-encode semantics
  are unambiguous enough to implement Phase A then B without re-litigating design.
- Mid-level L5→L4 residual update is explicit.
- Data-less rod is called out as a required prephase, not an afterthought.
