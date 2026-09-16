# Divergence Shape-Field Bending Energy TDD

## Status

Draft — **volume-style shell bending** via Hepworth projection, using **paired Darboux
fields** and a **Green / divergence field energy** on a **dual one-ring slab stencil**.
Replaces:

- [`medial_axis_projected_bending_energy_tdd.md`](medial_axis_projected_bending_energy_tdd.md) (deleted)
- [`osculating_torus_medial_bending_tdd.md`](osculating_torus_medial_bending_tdd.md) (deleted)

**Not this doc:** Halloween cubification (ARAP + ℓ¹ normals) —
[`reference/MANIFEST.md`](../reference/MANIFEST.md) (Liu & Jacobson 2019).

**Related code:**
[`darboux_cyclide.hpp`](../include/gaudi/albers/darboux_cyclide.hpp),
[`shape_operator_constrained_fit.hpp`](../include/gaudi/albers/shape_operator_constrained_fit.hpp),
[`shell_constraints.hpp`](../include/gaudi/hepworth/block/shell_constraints.hpp),
[`shell.hpp`](../include/gaudi/asawa/shell/shell.hpp),
[`smallvec.hpp`](../include/gaudi/smallvec.hpp),
[`ray_simplex.hpp`](../include/gaudi/arp/ray_simplex.hpp),
[`ray_traversal.hpp`](../include/gaudi/arp/ray_traversal.hpp),
[`brute_force_tree.hpp`](../include/gaudi/arp/brute_force_tree.hpp),
[`scripts/sympy/darboux.py`](../scripts/sympy/darboux.py).

---

## Concept

### One-line pitch

At each foot vertex, build a **slab stencil** (top one-ring + inward pierce + bottom
one-ring). Fit a **rest Darboux field `Q_r`**, embed deformed samples in a **live field
`Q_k`**, and penalize **field mismatch** `ΔD = D(Q_r,·) − D(Q_k,·)` (and its gradient)
over stencil samples — **not** hard projection onto `{D = 0}`. Local resolve → Hepworth
**difference rows** (ring spokes + spine), same currency as ARAP.

### What we are not doing

| Rejected | Why |
|----------|-----|
| Per-face macro-element + xyz ridge march | Heavy; medial hinge decoupled from field energy |
| Osculating-torus pivot + `{D=0}` Newton project | Level-set projection is stiff; hinge is secondary |
| Spine spring `‖d − d₀‖` | Spine row pins **gradient-resolved** difference, not rest length |
| Raw `‖Q_r − Q_k‖` in ℝ¹⁴ | Gauge-heavy; use **evaluated field** `Q·v` at samples |

---

## Stencil topology

### Construction (per foot `v₀`)

1. **Top ring:** `shell.get_one_ring(v₀)` + foot index.
2. **Inward pierce:** ray from foot along **inward** mesh normal; first triangle hit;
   closest vert on hit face → pierce foot `v'₀`. Use
   [`arp::ray_pierce_tree`](../include/gaudi/arp/ray_traversal.hpp) on a
   `bvh_tree<3>` / `aabb_tree<3>` / `brute_force_tree<3>` over face simplices:

   ```cpp
   arp::Ray ray(foot, -n_in, t_min, t_max);
   auto hit = arp::ray_pierce_tree(*tri_tree, ray);  // id, t_ray, tri_local{u,v}
   ```
3. **Bottom ring:** `get_one_ring(v'₀)` + pierce foot.
4. **Centers:** `c_top = x[v₀]`, `c_bot = x[v'₀]` (foot verts, not centroids v1).
5. **Spine:** vector `d = c_bot − c_top` (difference DOF, not a rest spring).

Cache pierce pairing at rest; invalidate on remesh. Reuse tri/edge intersection
patterns from shell collision / `ray_triangle_intersect`.

### Stencil storage

Use [`smallvec`](../include/gaudi/smallvec.hpp) for ring ids, constraint `ids`, and
deduped slab vert lists — avoid heap on typical valence-6/7 rings.

---

## Two fields, two gauges

| Symbol | Role |
|--------|------|
| **`Q_r`** | Rest Darboux fit on dual ring at **`t=0`**; cached for sim |
| **`Q_k`** | Live field — “rest implicit for current geometry” (transport and/or infrequent refit) |
| **`x_i`** | Current stencil vertex positions |

**Embedding:** `D(Q_k, x_i) ≈ 0` at ring samples (value rows; full normal rows only at
foot + pierce if available).

**Shape memory:** at the same samples,

```text
ΔD(x) = D(Q_r, x) − D(Q_k, x)
```

For fixed **`Q`**, **`D(Q,x)`** is a **polynomial in `x`** (Darboux). Energy is a
coupled polynomial over all **`x_i`** in the stencil.

### Rest gauge vs advected gauge

| Gauge | Field | Sample | Meaning |
|-------|-------|--------|---------|
| **Rest (frozen)** | `Q_r` fixed | current `x_i` | `f_i^rest = D(Q_r, x_i)` — motion off rest surface |
| **Advected** | `Transport(Q_r)` | current `x_i` | `f_i^adv = D(Q_adv, x_i)` — residual after comoving field |
| **Live embed** | `Q_k` | current `x_i` | `g_i = D(Q_k, x_i)` — current geometry owns live implicit |

**Distortion beyond advection:** `g_i − f_i^adv`. **Field drift:** `ΔD(x_i)`.

`Transport(Q)` is approximate (local motor / SE(3) jet on coeffs, or refit warm-start);
exact polynomial pushforward is an open item.

Orient **`Q_r`** once: flip so `∇D(Q_r, foot) · n_in > 0`.

---

## Variational formulation (authoritative)

**This section is the derivation chain.** Discrete weights, SymPy codegen, and Hepworth
`project()` all follow from here — not ad hoc term picking.

### Domain and field mismatch

Slab domain **`Ω`** (between complementary pierce rings: top sheet, bottom sheet, and
interior between spine samples). Coefficients **`Q_r`** (rest) and **`Q_k`** (live) are
**fixed during the local vertex solve**; optional outer loop updates **`Q_k`**.

```text
ΔD(x) = D(x; Q_r) − D(x; Q_k)
```

For Darboux, **`D(·;Q)`** is a polynomial in **`x`**; **`ΔD`** is therefore a fixed
polynomial in **`x`** once **`Q_r, Q_k`** are fixed.

### Primary bulk energy (what we argmin)

Local stencil solve:

```text
x* = argmin_{x ∈ stencil samples}  E_bulk[x]
```

with

```text
E_bulk[x] = ∫_Ω |∇(ΔD)|² dV  +  λ_gap ∫_Ω (ΔD)² dV
```

Optional gauge / embedding terms (same variational step, separate weights):

```text
E_embed[x] = μ ∫_Ω (D(Q_k, x))² dV
E_adv[x]   = ν ∫_Ω (D(Q_k, x) − D(Q_adv, x))² dV
```

**`E_bulk`** is the divergence-derived shape-field energy between the two Darboux states.
**Not** level-set projection **`{D=0}`**, **not** a rest spine spring **`‖d−d₀‖`**.

### Divergence / Green identity (exact rewrite of the gradient term)

For smooth **`ΔD`** on **`Ω`** with boundary **`∂Ω`** (both ring sheets):

```text
∫_Ω |∇(ΔD)|² dV
  = ∫_∂Ω ΔD (n · ∇(ΔD)) dA  −  ∫_Ω ΔD · Δ(ΔD) dV
```

(Green’s first identity with **`φ = ψ = ΔD`**.)

So the **same** **`argmin E_bulk`** can be analyzed or discretized as:

| Piece | Integral | Role |
|-------|----------|------|
| Boundary flux | **`∫_∂Ω ΔD (n·∇ΔD) dA`** | Top + bottom ring samples |
| Bulk correction | **`− ∫_Ω ΔD · Δ(ΔD) dV`** | Few interior slab quadrature points |
| Gap term | **`λ_gap ∫_Ω (ΔD)² dV`** | Added directly (co-area / thin-shell limit ≈ **`h ∫ (ΔD)² dA`**) |

Implementation may quadrature **either** the left side **`∫ |∇(ΔD)|²`** **or** the
Green split — they are **equivalent in the continuum limit** when discretized consistently.
SymPy derives **`∂E_bulk/∂x`** from the **bulk integrand** (left side + gap); the Green
split is used to **validate** boundary vs interior quadrature and to choose sample placement.

### First variation → local gradient descent

With **`Q_r, Q_k`** fixed, **`ΔD`** depends on stencil vertex positions **`x_i`**. For a
sample point **`x`** tied to vert **`i`**:

```text
∂(ΔD)/∂x_i = (∂/∂x) D(Q_r, x) − (∂/∂x) D(Q_k, x)     // chain through Darboux eval
```

Standard calculus of variations on **`E_bulk`** (bulk form):

```text
δE_bulk / δx  ∋  −2 Δ(ΔD) ∇(ΔD)  −  2 λ_gap ΔD ∇(ΔD)     // from |∇(ΔD)|² and (ΔD)² terms
```

(Plus **`E_embed` / `E_adv`** terms if enabled.) **SymPy generates this exactly** from
the summed discrete **`E`** — no hand-wavy gradients.

**Local algorithm (per stencil, per Hepworth outer iter):**

1. Fix **`Q_r`**, **`Q_k`** (and optional **`Q_adv`**).
2. Quadrature: ring verts + optional interior **`p_k`** along spine → discrete **`E`**.
3. **`g = ∇_x E`** via generated code (or SymPy lambdify in tests).
4. Gradient descent / few Newton steps on **stencil verts only** → **`x*`**.
5. **`project(q → p)`** from **`x*`** (below).

### Discrete quadrature (implementation of the integral)

Replace **`∫_Ω f(ΔD, ∇ΔD) dV`** by:

```text
Σ_i w_i f(ΔD(x_i), ∇(ΔD(x_i)))     // ring verts (boundary contribution)
+ Σ_k w_k f(ΔD(p_k), ∇(ΔD(p_k)))   // interior slab (bulk correction path)
```

Weights **`w_i`** from cotan / valence / slab thickness estimate **`h`**. Foot + pierce
may carry extra flux weight for the **`n·∇(ΔD)`** boundary term when using the Green
split explicitly.

Scale: fix **`‖Q_r‖`** (homogeneous **`D`**) once per stencil.

---

## Local resolve → Hepworth projection

Same pattern as ARAP: **differences**, not absolute pins. **Targets come from the
variational argmin above**, not guessed springs.

**Casting split (important):**

| Role | Geometry |
|------|----------|
| Fit **`Q_r` / `Q_k`** | Full pierce stencil (foot + pierce one-rings) |
| Energy quadrature + GD DOFs | **Foot one-ring only** (`top_ring_positions`) |
| Hepworth `fill_A` / `p` | **Top hub-and-spokes only** |

Do **not** project pierce-side verts from the field energy — dual-ring field rows
couple distant mesh DOFs and fight relative-plane twist. Outer loop may still
refit **`Q_k`** on the full stencil; local GD freezes both fields and re-places
foot samples inside that gauge (divergence / Green story).

### Local step

Minimize discrete **`E_bulk`** on **foot-ring** DOFs (SymPy **`∇E`**, few GD steps)
with **`Q_r, Q_k` fixed** → **`x*`**.

### `project(q → p)`

**Not** `p = w · x*`. Emit **relational targets** on the foot ring only:

```text
Top ring:    p_j = w_t · (x_j* − x_ctop*)
```

**`fill_A`:** rows **`x_j − x_hub`** for the foot one-ring. Physics lives in **`p`**.

Field forces enter the local resolve; optional **normal-only** wiring of **`∂E/∂x`**
into spoke targets (tangential freedom like ARAP).

### ADMM (optional)

Auxiliaries **`z_i ≈ ΔD(x_i)`** or **`z_i ≈ flux_i = ΔD · (n·∇ΔD)`** with dual **`u`**,
as in Liu–Jacobson cubify — on **field samples**, not **`R n`**.

---

## SymPy derivation

Location: [`scripts/sympy/slab_field_energy.py`](../scripts/sympy/slab_field_energy.py) —
extends [`darboux.py`](../scripts/sympy/darboux.py).

Implements the **same** **`E_bulk`** density as the variational section:

```python
delta = D(Q_r, x) - D(Q_k, x)
E = |grad(delta)|² + lam_gap * delta**2
```

Also emits **Green split** densities **`ΔD (n·∇ΔD)`** and **`−ΔD Δ(ΔD)`** for quadrature
cross-checks against **`|∇(ΔD)|²`**.

Stencil: sum **`E`** over samples; **`sp.diff(E, x_i)`** per vert. Emit via
[`albers_printer.py`](../scripts/sympy/albers_printer.py) / [`generate.py`](../scripts/sympy/generate.py):

- `slab_field_energy(Qr, Qk, x, lam_gap, out)`
- `slab_field_energy_grad(Qr, Qk, x, lam_gap, out_g)`

Validate: `python3 scripts/sympy/slab_field_energy.py --fd-check` (FD vs SymPy grad).

---

## `smallvec` (C++ adjacency storage)

Rust [`SmallVec`](https://docs.rs/smallvec) pattern for stack-first, heap-spill collections.

**Header:** [`include/gaudi/smallvec.hpp`](../include/gaudi/smallvec.hpp)

```cpp
template <typename T, std::size_t INLINE = 16>
class smallvec { /* inline [INLINE] + spill to std::vector */ };
```

| Use | Typical `INLINE` |
|-----|------------------|
| `get_one_ring` copy | 12–16 |
| Stencil vert ids for constraint | 24 |
| Pierce hit face verts | 8 |

**API (v1):** `push_back`, `size`, `data`, `begin`/`end`, `clear`, `reserve`, move-only
or rule-of-five as needed. **No** full Rust feature parity (insert_remove, drain, etc.)
until needed.

**Migration:** replace hot-path `std::vector<VertId>` temporaries in
`build_slab_stencil`, `shell_constraints_init`, constraint `ids` assembly.

---

## Contrast table

| | Shell `bending` | ARAP | Divergence shape-field bend |
|---|-----------------|------|----------------------------|
| Stencil | Vertex one-ring | Spokes | **Dual ring + spine** |
| Rest data | `_vg` | `_d_rest`, optional `R` | **`Q_r`**, pierce cache |
| Local target | Scaled Laplacian | **`R d_rest`** | **Field-resolved differences** |
| Volume | No | No | **`ΔD` slab quadrature** |
| Projection | Vector at hub | Edge diffs | **Spokes + spine + field grad** |

---

## Implementation phases

| Phase | Deliverable | Status |
|-------|-------------|--------|
| 1 | [`smallvec.hpp`](../include/gaudi/smallvec.hpp) | done (header) |
| 2 | `build_slab_stencil` (pierce + dual ring, `smallvec`) | pending |
| 3 | SymPy [`slab_field_energy.py`](../scripts/sympy/slab_field_energy.py) + generated C++ | sympy draft |
| 4 | `darboux_slab_field_bend` `block_constraint` (local resolve + `project`) | pending |
| 5 | Rest **`Q_r`** init bundle; transport for **`Q_k`** | pending |
| 6 | Tests: flat slab, cylinder shell, pierce failure skip | pending |
| 7 | Vermeer / medial demo wire-up | pending |

---

## Open questions

1. **Pierce failure:** skip stencil vs fallback to single-ring + weak bending?
2. **`Q_k` schedule:** transport-only vs refit every `N` steps?
3. **Overlap:** valence-weighted constraint merge on shared verts?
4. **Interior slab samples:** count and placement along spine?
5. **Concave / no pierce:** MAT-gated skip (inherited policy from older medial drafts)?

---

## Test plan

| Test | Expect |
|------|--------|
| SymPy grad vs FD on single sample | `< 1e-6` rel err |
| Flat plate, small bend | Mostly tangential motion; low **`E_gap`** |
| Thick cylinder, inward pierce | Stable thickness; **`ΔD`** small under moderate bend |
| Pierce miss | Constraint skipped, no NaN |
| **`smallvec`**: push past inline cap | Correct size, no leak |

---

## See also

- Dipole tunneling (different energy, same projection idiom):
  [`dipole_rod_tunneling_projection_tdd.md`](dipole_rod_tunneling_projection_tdd.md)
- Darboux fit: [`darboux_cyclide.hpp`](../include/gaudi/albers/darboux_cyclide.hpp)
- Medial axis **demo** code (xyz ridge — reference only, not this energy):
  [`darboux_cyclide_medial.hpp`](../include/gaudi/duchamp/darboux_cyclide_medial.hpp)
