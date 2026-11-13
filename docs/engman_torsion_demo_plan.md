# Engman torsion demo: torsional buckling and Hopf-fiber shapes

A new, small demo on the Engman solver (`include/gaudi/engman/`, `MATH.md`
section 17) for torsional buckling of tubes and helical ridges along the
Hopf fibers (Villarceau circles) of a torus. It is also the place where TODO
item 18 phases 4 and 5 happen: a row gets in only once it does what it says
under the hard coupling. `monge_fairing_demo` stays as it is, as the
projection-solver reference.

Background: `docs/kkt_admm_curvature_solver_plan.md` (solver),
`docs/twist_and_tension_constraints.md` (twist readout, tension rows).

---

## As built: composable Engman scene (2026-10-27)

| Piece | File |
|---|---|
| Energies: `as_energy` (any constraint, soft), `store_bending` (`plate`, `helfrich`), `torsion_target`, `deviator_hold`, `twisted_spin2`, `twist_rate`; readouts `read_twist`, `read_torsion` | `engman/energies.hpp` |
| Hard per-face area `face_area` | `engman/constraint.hpp` |
| Terms: `coupling`, `memory`, `store_metric`, `bending`, `torsion`, `twist_rate`, `twisted_dirichlet`, `anisotropy_hold`, `face_area`, `conformal`; `term` variant, `find<T>` | `engman/terms.hpp` |
| Scene: asawa shell, per-frame remesh (`asawa::shell::dynamic`), grouped config (`mesh`, `step`, `rest`, `remesh`, `terms`, `readout`, `draw`), presets `torsion` (default), `torsion_hopf`, `torsion_elastic`, `twist_rate`, `relax` | `engman/scene.hpp` |
| Demo (`python gaudi.py run vermeer_engman_demo`, preset by `ENGMAN_PRESET`) | `duchamp/engman_demo.hpp`, `..._adapter.hpp`, `projects/vermeer_engman_demo/` |
| Tests | `test/engman_tests.hpp`; opt-in `engman_sweep` (`GAUDI_SWEEP=1`) and `engman_remesh_probe` |

A scene is a mesh plus a list of terms. Each term is a small struct with its
own knobs and a `build`; constraint terms carry `mode off | soft | hard`. A
preset is a designated-initializer term list:

```cpp
scene_config cfg{.terms = {terms::coupling{}, terms::memory{}, terms::bending{.w = 0.1},
                           terms::torsion{.w = 10, .tau = 4}}};
```

Findings:

- **The twist energy acts only on torsion.** `torsion_target` is
  `w ∫ (W(a, N×a) − τ)²`, one linear row per edge on the store: mean
  curvature, anisotropy magnitude and principal angle are free, so curvature
  changes, and under the hard coupling positions follow. Its linear part is
  the Helfrich-Prost chiral term. A first version completed the chiral term
  into a pin on the whole tensor, which prescribed k0 and k1; it was replaced.
- **The twist is an instability.** A round torus is stationary for any chiral
  term by mirror symmetry; torsion grows from the startup noise. At `τ = 4`,
  `w = 10` the torsion rises steadily (0.002 to 0.09 over 18 frames,
  displacement about 0.3 of the tube radius) with the coupling exact to
  `1e-14`. At `τ = 8` it is stronger but the coupling breaks by frame 12.
- **Frame step.** `h = 0.05` lets a strong drive move the shape too far per
  frame and the coupling stops converging; `h = 0.02` keeps it exact. More
  iterates or a lower `ρ` do not fix `h = 0.05`.
- **The twist-rate energy breaks the coupling solve** under strong drive (the
  store departs from the mesh by 40%), and readouts from the store overstated
  the twist. Readouts now come from the mesh's own read-out.
- **Weights scale with the frame step.** Inertia is `M / h²`, so going from
  `h = 0.05` to `0.02` stiffens it about six times; the twist-rate preset's
  weight went from 1 to 10 with it.
- **Remesh thresholds.** `asawa::shell::avg_length` returns about a quarter
  of the mean edge length (it divides the edge sum by twice the corner
  count), which `monge_fairing_demo`'s `C = 2` absorbs. The scene uses the
  true mean and a band of 0.5 to 1.6 mean edges; remeshing costs about 3 ms a
  frame and the vertex count settles near 600 to 700.
- **Solver fix.** With no constraint rows the automatic `ρ` was infinite and
  the merit became `∞ · 0`; `ρ` is 1 then.

## 0. What we want to see, and what physics produces it

| Shape | Mechanism | Regime |
|---|---|---|
| Helical diamond buckles on a twisted tube (Kresling / Yoshimura) | metric held, ends twisted, bending cheap | solid shell: hard area, high shear modulus `μ`, small bending `κ` |
| Helical ridges along the Hopf fibers of a torus | material shear along the Villarceau winding, metric otherwise held | same, closed surface, no pins |
| Piecewise developable creased tubes | as the first row with `κ` very small: bending concentrates in creases | large Föppl-von Kármán number `γ = μ R² / κ` |
| Pearling, budding (sanity regime) | spontaneous curvature at fixed area, no shear stiffness | fluid membrane: hard area, `μ = 0` |
| Ruffled growth | curvature-driven area targets held hard, `μ > 0` | growth-induced buckling |

Two facts shape the presets:

- **Gaussian curvature is intrinsic.** If the metric is held (hard area plus
  stiff shear), `K` stays where it started. A cylinder starts developable and
  stays near-developable as it buckles, so torsion on a cylinder gives
  developable facets for free. A torus cannot become developable while
  smooth: smooth flat tori in R³ exist only in `C¹` (Nash-Kuiper, Borrelli's
  corrugated torus), and Hopf tori are flat in S³, not R³. "Developable Hopf
  shapes" therefore means piecewise developable: creases along the fibers,
  which needs the metric soft enough to give up `K` and `κ` small enough for
  creases.
- **A symmetric start is a saddle.** Every preset seeds noise
  (`perturb_vertices_noise`, as in `monge_fairing_demo`) or an azimuthal
  bump. Multipliers do not remove saddle points.

## 1. Inventory

### Built and used as is

| Piece | Where | Role here |
|---|---|---|
| Engman step, KKT, report | `engman/solver.hpp`, `engman/kkt.hpp` | the integrator |
| Hard coupling `R̄(x) = S` | `engman::coupling` | the store is the shape operator |
| Hard position pins | `engman::position_pin` | twisting the caps (2.4) |
| Store metric (`λ_H L_τ + λ_z L_z`) | `engman/store_metric.hpp` | store proximal term |
| Memory | `hb::dirichlet_memory_gn` | position proximal term |
| Soft shear modulus `μ` | `hb::triangle_conformal` | traceless strain energy |
| Store pins and goals | `hb::edge_shape_fair` (`pin`, `tensor_metric`), `willmore_target` | bending (2.2) |
| Villarceau winding | `stretch_winding {.ring = 1, .tube = 1}` (the angle only) | Hopf-fiber shear axis |
| Material directions | `kusama/material_direction.hpp` | carried fiber, twist readout |
| CR block and store | `hb::shell_cr_block` | store DOFs, stencil table, vertex stars |
| Remesher | `asawa::shell::dynamic` | closed meshes only (no boundary handling found), so every preset is closed |

### Not used: the Monge local fits

The demo does not use the per-vertex Monge-patch machinery from
`monge_fairing_demo`: `monge_axes`, `fan_chart`, `ring_frame`, `frame_fair`,
`albers_ring_fair`, `frame_minimizer`, the `axis_ratio` and
`curvature_shear_target` goals, the Hopf jets in `hopf_sweep.hpp`, or
`vertex_frame_cache_compute`. The shape operator is the store, read through
the CR coupling, so no local fit is needed. Only geometric helpers carry over:
`stretch_winding` for the Villarceau angle, `perturb_vertices_noise`, and mesh
construction. Copy those, or move them to `duchamp/utils/`. Do not extract
anything else from `monge_fairing_demo`.

The datum transfer of 2.6 needs a tangent frame per vertex. Define it fresh in
`engman/multiplier_datums.hpp`: the area-weighted vertex normal and the edge
to the lowest-id neighbor, projected and normalized, with tensors moved through
`kusama::weingarten::vertex_stars`. Do not use the Monge-era `vertex_chart`.
The chart only has to be deterministic and identical between scatter and
gather.

`hopf_sweep.hpp` can optionally build a starting mesh later: a Pinkall Hopf
torus from an elastic curve on S², as a `torus_hopf` initial shape.

### New

| Piece | Proposed file | Section |
|---|---|---|
| Hard per-face area constraint | `engman/constraint.hpp` (`face_area`) | 2.1 |
| Edge-length shear energy with an axis source | `hepworth/block/edge_metric_target.hpp` | 2.3 |
| Twist readout `Ω` and Uzawa multiplier | `engman/twist.hpp` | 2.5 |
| Region pins re-selected each frame | in the demo (`twist_caps`) | 2.4 |
| Multiplier transfer across remesh | `engman/multiplier_datums.hpp` | 2.6 |
| Developability energy on the store (optional) | `edge_shape_fair` row kind `det` | 2.7 |
| Curvature-driven area targets | `engman/area_targets.hpp` | 2.8 |
| The demo, adapter, vermeer project | `duchamp/engman_torsion_demo.hpp`, `..._adapter.hpp`, `projects/vermeer_engman_torsion_demo/` | 3 |

## 2. Components

### 2.1 `engman::face_area` (hard)

Per face `C_f = (A_f(x) − a*_f) / a_ref`, `a_ref` the mean face area at the
frame start (so `ρ_rel` means the same at any resolution). Jacobian
`∂A_f/∂x_i = ½ n_f × (x_{i+2} − x_{i+1})`. Targets: the frame-start areas
(incompressible), or the curvature rule of 2.8.

- **Count.** About `2V` rows on `3V` positions. It does not lock (full edge
  isometry would: about `3V` rows on a closed mesh).
- **Redundancy.** Faces whose three vertices are all pinned get no row. Any
  remaining redundancy with the pins is absorbed by the `−I/ρ` block of the
  quasi-definite system.
- **Multiplier.** `λ_f` is the local surface tension field. Draw it signed.

Tests: areas met to `1e-10` after one step from a perturbed state; zero
residual and no motion at rest; rigid motions in the Jacobian's null space;
an end twist through `position_pin` still moves the interior (no locking).

### 2.2 Bending: two knobs on the store

Under the hard coupling there is no `bend` row. Bending is a store energy that
reaches positions through the coupling's tangent space.

- **Plate:** `edge_shape_fair` `pin` toward zero on the full tensor with
  `tensor_metric = true`, weight `κ`. Energy `κ ∫ |W|²`, the thin-shell
  bending energy. Default for the solid presets.
- **Helfrich:** the trace pin toward `c₀`, weight `κ`. Default for the fluid
  preset.

Both rows already exist. The new code is only the preset wiring.

### 2.3 `edge_metric_target`: shear and scale as edge-length targets

One energy for every in-plane metric target built from a strain tensor `ε` in
a tangent frame `(a, b)` (`docs/twist_and_tension_constraints.md` 1.3, 2.6):

```
ℓ*_e² = ℓ_e² (1 + 2 ε(d̂_e, d̂_e))         shear:  ε = γ/2 (a⊗b + b⊗a)
                                          scale:  ε = ½ diag(s_a² − 1, s_b² − 1)
```

Gauss-Newton residual `r_e = w sqrt(ℓ_e) (|d_e|² − ℓ*_e²) / ℓ_e²`, rest at the
frame start, so `γ` is a per-frame increment. Axis sources:

- `winding`: the Villarceau direction from `stretch_winding` (Hopf fibers);
- `material`: `kusama::material_direction`, carried with the surface;
- `store`: the smoothed `κ₂` line field of the frame-start store
  (twist doc 1.2), sign-free.

This is a port of `shear_target` from the one-ring projection to an edge-length
energy, so it composes with `face_area` without double-counting the area. Keep
`shear_target` in `monge_fairing_demo` untouched.

Tests: a straight cylinder under a `winding` shear slides without shape change
to `1e-10`; flipping the sign of `a` on random edges leaves the rows unchanged.

### 2.4 Twist by rotating caps

`dynamic` has no boundary handling, so tubes are closed capsules. Each frame,
select the vertices in the two cap regions (`|z − z_end| < r_cap`, re-selected
after every remesh) and pin them with `engman::position_pin` at their current
positions rotated by `±dθ/2` about the axis, optionally moved together by
`dz` (axial compression; Kresling patterns need some). Hard pins mean the
imposed twist is exact and the cap torque is the pin multiplier: plot torque
against total twist to see the buckling point as a drop in stiffness.

### 2.5 Twist readout and Uzawa (torus preset)

`engman/twist.hpp` computes the axis field, `Δ_e`, `ω_e` and the global
`Ω = Σ A_e |z_e| ω_e / Σ A_e |z_e|` from the frame-start store (twist doc
1.2). The shear amplitude is a multiplier: `γ_e = λ_twist R_e`, with
`λ_twist += ρ_twist (Ω* − Ω)` once per frame and clamped to
`±γ_max / R_max`. No KKT row: a dense row would fill the factor, and on a
surface of revolution `∂Ω = 0`, which would send a hard multiplier to
infinity. The capped Uzawa update only grows until symmetry breaks.

Variant B (later, behind a flag): the energy `−λ_twist Ω(S)` on the store as a
fresh-gradient linear term each iterate, no axis field, same Uzawa update.

### 2.6 Multiplier transfer across remesh (TODO 18, phase 4)

At the end of a frame, scatter the intensive coupling multiplier
`μ_e = λ_e / c_e` (`MATH.md` 17.5) and the face-area multiplier
`λ_f / A_f` to vertex datums (weights `A_e / 3`, `A_f / 3`; tensors through
`vertex_stars` charts). After the remesh, gather and rescale. Also carry the
store offset `δ_e = S_e − R̄_e(x)` the same way. `λ_twist` is a scalar and
persists as is.

Tests: on a rigidly turned torus a transferred tensor round-trips to `1e-12`;
over 20 frames with a remesh each frame, the first-iterate `|C|` after a
remesh is within 2× of the last before it.

### 2.7 Developability on the store (optional)

Per edge `r_e = w sqrt(c_e) det(S_e) / κ_ref²`, `det S = k11 k22 − k12²`.
Gauss-Newton, Jacobian linear in `S`. It acts on positions through the coupling
and pushes `K` toward zero. It fights a held metric (Gauss), so use it only
with a soft metric (`face_area` off or curvature-driven, `μ` moderate). This
is the knob for creased, piecewise developable tori. Expect crease formation
to depend on the remesher keeping dihedrals (`max_flip_dihedral`).

### 2.8 Curvature-driven area targets

`a*_f = a_f (1 − α + β (|W_f|² − κ_c²) / κ_c²)` from the frame-start store,
clamped to `[1 − α_max, 1 + β_max]`, fed to `face_area`. With `μ = 0` this is
the surface tension we settled on. With `μ > 0` it drives growth-induced
buckling.

### 2.9 Growth as the drive: what buckles and what does not

`monge_fairing_demo`'s position `shear` is `axis_ratio` at ratio
`(1 + g) / (1 − g)` with an `area_band` (default `[1, 1]`), rest reset every
frame. It is traceless growth: each frame the ring's metric is stretched
along the axis and compressed across it at fixed area. The compression is
relieved by bending, and with bending cheap against the membrane the
relief is ruled folds. That is why it looks developable.

Under Engman the same drive splits by trace, like the curvature store:

| Growth part | Carried by | Hard or soft |
|---|---|---|
| trace (area) | `face_area` target `a*_f` | hard by default, toggle (3.1) |
| traceless (ratio, shear) | `edge_metric_target` (2.3) | soft, weight `w_growth` |

"Increase the area and let the constraints do the work" works only if the
growth is incompatible with the current metric:

- **Uniform area growth on a free closed surface does nothing.** It is
  realized exactly by a dilation, and `∫|W|²` is scale-invariant, so nothing
  buckles.
- **Uniform growth with confinement buckles.** Pinned caps (2.4) or a held
  volume stop the dilation, and the surface buckles like a rod growing
  between walls.
- **Nonuniform area growth buckles.** A gradient in `a*_f` changes the
  target `K` (`Δ log` area acts like `−2K`): edge growth ruffles, center
  growth domes. The curvature rule of 2.8 is one such field.
- **Traceless growth buckles,** which is the monge mechanism, now as
  `edge_metric_target` with `ratio_per_frame` and the area held by
  `face_area` separately.

Preset `growth` gets a `growth.mode`: `uniform_confined`, `gradient`,
`curvature` (2.8), `traceless`.

### 2.10 Bending and membrane energies: which to use

- **Bending is a store energy.** The `pin` toward zero on the full tensor
  (`tensor_metric`) is `κ ∫|W|²`, and the trace pin is Helfrich. Both are
  quadratic in `S`, so `H_SS` is constant and block-diagonal, and all the
  nonlinearity sits in the coupling's exact symbolic Jacobian. The CR
  mid-edge read-out is the consistent discretization (Wardetzky et al. 2007).
- **Legacy `edge_bend` / `willmore_bend` (hinge) and the block `bending`
  row stay out of the energy list.** They discretize the same energy in
  position space a second time. The two disagree at mesh resolution, and
  under the hard coupling the multiplier absorbs the difference as noise.
  The hinge Hessian is 12 × 12 and nonconvex per edge. Keep the hinge as a
  diagnostic overlay: its energy against `κ ∫|W|²` on a refinement sequence.
- **ARAP is a membrane energy, not bending.** It holds the full metric (trace
  and traceless together) against its rest and is rotation-invariant, so it
  does not resist bending. Offer it as `membrane.mode = arap` beside
  `area_conformal` (hard `face_area` plus `μ` `triangle_conformal`).

**Rest persistence: elastic versus viscous.** Every membrane rest above is
rebuilt at the frame start, so a membrane resists only the strain of one frame.
Shear accumulates freely across frames: the sheet behaves as a viscous fluid,
not an elastic solid. Viscous sheets buckle too (rope coiling, folding of
viscous sheets, Ribe; Teichman and Mahadevan). By the Stokes-Rayleigh analogy
a viscous sheet is an elastic one with strain replaced by strain rate. Here
the per-frame drive increment plays the strain rate, so the `η = 1` regime
buckles when the drive per frame passes a threshold, with a wavelength set by
the same `sqrt(κ / μ)` scale. That is how the monge shear row buckles today.
The regimes differ in what holds the pattern: an elastic buckle is an
equilibrium that stores energy and pushes back (cap torque persists), while a
viscous buckle is frozen flow with no restoring force. The threshold depends on
the amount of strain for an elastic sheet and on its rate for a viscous one.
The ratio of relaxation time to drive time (the Deborah number,
about `(1 − η) / η` frames over frames per drive cycle) moves between them.
Add a persistent rest metric `ḡ_v`, a symmetric 2 × 2 per vertex
in its chart, carried on vertex datums like the multipliers (2.6), with a
Maxwell relaxation knob:

```
ḡ_v ← ḡ_v + η (g_v − ḡ_v) + growth increment       η = 1 viscous (today), η = 0 elastic
```

ARAP and `triangle_conformal` read their rest from `ḡ_v`. `face_area` holds
`det ḡ_v` per face. This is the length store of `MATH.md` 17.7, sited on
vertices so it survives the remesh.

### 2.11 Willmore track: minimal `∫H²` shapes (Soliman et al. 2021)

The energy is the Helfrich trace pin toward `c₀ = 0` (2.2), so `κ Σ_e A_e τ_e² / 4`,
reaching positions through the hard coupling. Three levels of constraint,
each a preset:

1. **Free Willmore flow.** Plus total area held by an Uzawa multiplier
   (`∫H²` is scale-invariant, so this only fixes the null direction). A torus
   should flow to the Clifford torus: `R / r → √2`, `W → 2π²`. A sphere-like
   start goes to the round sphere, `W → 4π`. This is the validation preset.
2. **Reduced volume.** Plus enclosed volume by Uzawa (2.3 of the twist and
   tension doc), at reduced volume `v = 6 √π V / A^{3/2} < 1`. This gives the
   Helfrich vesicle phase diagram: prolates, discocytes, stomatocytes.
3. **Constrained Willmore (Soliman).** Willmore at a fixed conformal class.
   Per interior edge, a hard length-cross-ratio constraint
   `C_e = log lcr_e(x) − log lcr̄_e`, `lcr_e = (l_il l_jk) / (l_lj l_ki)` over
   the edge's diamond, an `engman::constraint` with a sparse Jacobian through
   four lengths. Rank is about `E − V` on `3V` positions (lengths are fixed
   only up to a per-vertex scale), so it does not lock. Far from the Clifford
   class, tori cannot reach `2π²` and buckle into constrained Willmore tori.
   The conformal class is the knob: start from a torus and set `lcr̄` from a
   rectangular or skewed flat torus to choose it.

**Conformal class across a remesh.** Discrete cross ratios depend on the
triangulation. Carry a per-edge conformal rest length `ℓ̄_e`, from which
`lcr̄` is computed:

- **Flip:** the Ptolemy relation `ℓ̄_new = (ℓ̄_a ℓ̄_c + ℓ̄_b ℓ̄_d) / ℓ̄_old`
  keeps the discrete conformal class exactly (Bobenko, Pinkall, Springborn;
  Gillespie, Springborn, Crane). The `datum::flip` hook exists.
- **Split:** compute the new lengths in the rest triangles by the law of
  cosines.
- **Collapse:** breaks the class, so disable collapses in this preset.

Edge datums have no split or collapse transfer in `datums.hpp` today. Run the
conformal preset with `remesh_every = 0` first, and add Ptolemy flips and
rest-metric splits after that works.

Readouts: `W / 4π` (sphere) or `W / 2π²` (torus), `R / r`, reduced volume,
cross-ratio residual, the area and volume multipliers.

### 2.12 Hard constraints in shape space

Instead of fixing the conformal class with metric cross ratios (the legacy
`hepworth/legacy` constraint stays where it is), constrain the store. Under
the hard coupling, a constraint on `S` is a constraint on the realized shape
`R̄(x)`. All of these rows are linear or algebraic in `S`, so the store block
of their Jacobian is cheap and the positions see them only through the
coupling.

**Counting rule.** The coupling's `3E` rows pair one-for-one with the `3E`
store unknowns and cost the positions nothing. Store constraints are charged
against the `3V` positions:

| Family | Rows | Verdict |
|---|---|---|
| per edge scalar (or per-edge phase) | about `3V` | locks, or is infeasible under Codazzi; soft only |
| per face scalar | about `2V` | borderline; never together with `face_area` |
| per vertex scalar (through `vertex_stars`) | `V` | leaves about `2V` freedom; the useful family |
| global scalar | 1 each | free; Uzawa, never a KKT row |

Do not stack families past `3V` (for example `face_area` plus a per-vertex
constraint). A hard constraint must also be feasible on the topology:
`K = 0` everywhere cannot hold on a smooth closed surface (Gauss-Bonnet), and
all-umbilic holds only on spheres, so those stay soft.

**Candidates, in order of interest:**

1. **Bilayer couple, global.** `∫H dA = M*` with area and volume (Svetina and
   Žekš; Seifert's phase diagram). The one well-studied shape-space
   constraint: budding, pears, stomatocytes. Validates the global path.
2. **Pipe surface, per vertex.** `κ_hi(S_v) = 1/r`: one principal curvature
   constant, so the surface is a tube of constant radius about a free
   centerline. With Willmore and held total area, the tube's bending energy
   becomes the centerline's rod bending. Driven by twist (2.4, 2.5) with the
   ends held, this should exchange twist for writhe (Călugăreanu,
   `Lk = Tw + Wr`) and coil into plectonemes. Torsional buckling of an
   elastic rod produced from a shape constraint alone. The eigenvalue rows
   `dκ/dS` are already derived (`MATH.md` 16, Eigensatz); guard near umbilics,
   which a tube does not have.
3. **Canal surface, per vertex.** `e₁ · ∇κ₁ = 0` for one family only (one
   of the two `dupin` rows of `edge_shape_fair`, made hard): envelopes of a
   one-parameter family of spheres, which are tubes of varying radius. This is
   the family of tubular undulations. Both families hard is Dupin cyclides,
   which is about `4V` rows and too much.
4. **Linear Weingarten, per vertex.** `a H + b K = c`, with `c` either fixed or
   a single global unknown. `b = 0` is constant mean curvature (unduloids,
   pearls as exact surfaces). The pipe surface is the principal-curvature
   version. Slide `(a, b)` live.
5. **Curvature caps, per vertex, inequality.** `|κ| ≤ κ_max`: a minimum
   feature radius, which also keeps the mesh locally injective. Needs
   complementarity (`λ ≥ 0`, active set per iterate). Engman handles
   equalities today, so this is the last one.

`shape_det_conserve.hpp` already conserves and redistributes the store's
determinant as soft rows. Its per-vertex Gauss sum is the natural first test of
the per-vertex wiring made hard.

## 3. The demo

`duchamp/engman_torsion_demo.hpp`, modeled on `monge_fairing_demo` but small.
Shared plumbing is limited to mesh setup, noise, the remesh call and the
snapshot (see "Not used" in section 1).

### 3.1 Config

```cpp
enum class torsion_preset { torus_hopf, capsule_twist, crease_tube, fluid_pearl, growth };

struct engman_torsion_config {
  torsion_preset preset = torsion_preset::capsule_twist;
  mesh      { torus R, r, res | capsule L, R, res; noise amplitude; seed bump m, amplitude }
  solver    { engman::solver_config; iterations; dt; damping; remesh_every }
  membrane  { hard_area (bool); mu (triangle_conformal weight) }
  bending   { mode plate | helfrich; kappa; c0 }
  twist     { caps: dtheta, dz, r_cap  |  shear: axis source, gamma or Uzawa (Omega*, rho, gamma_max) }
  develop   { w_det }               // 2.7, 0 = off
  growth    { alpha, beta, kappa_c } // 2.8, off unless preset growth
  draw      { overlay; scale }
};
```

Every weight is dimensionless: `μ` against inertia at the mean face, `κ`
through `γ = μ R² / κ` (shown in the UI).

**Tunable at runtime.** Every field above is a live UI control, as in
`monge_fairing_demo`, and constraints and energies are rebuilt each frame,
so a change takes effect on the next frame. What changes is what a number
means:

- **Energy weights** (`μ`, `κ`, `c₀`, `w_growth`, `w_det`, `w_shear`) are
  tuned exactly as now. Under the hard coupling they trade only against each
  other and inertia, never against a constraint, so a weight change moves
  the shape without loosening the coupling.
- **Drive rates** (`dθ`, `dz`, `ratio_per_frame`, `α`, `β`, `Ω*`) set how
  fast the design moves.
- **Hard constraints have no accuracy weight.** Their `ρ` is a conditioning
  knob only (`solver_config::rho_rel`).
- **Hard or soft per family.** Each constraint family (coupling,
  `face_area`, cap pins) gets a switch `hard | soft(w)`. Soft wraps the
  constraint's residual and Jacobian as a Gauss-Newton energy with weight
  `w` (`engman::as_energy`, a few lines on `constraint::evaluate`). This
  restores the old penalty behavior one family at a time for comparison, and
  is the only way to make, for example, area a tunable compromise.
  `solver_config::dual_updates = false` already does this globally.
- **Presets are starting points.** Switching preset reloads the defaults.
  Editing a field afterwards does not.

### 3.2 Frame

1. Prepare `shell_cr_block`; gather multipliers and store offset (2.6).
2. Build energies (inertia, memory, store metric, `triangle_conformal`,
   bending, shear, `det`) and constraints (coupling, `face_area`, cap pins).
3. Update `λ_twist` (torus) or the cap goals (capsule).
4. `engman::solver::step`.
5. Scatter multipliers; remesh every `remesh_every` frames.

### 3.3 Overlays

- coupling multiplier `|μ_e|`: where the design fights realizability;
- face-area multiplier `λ_f`, signed: the tension field;
- angle defect per vertex: developability, near zero off the creases;
- principal lines (existing drawing);
- readouts: `Ω`, `λ_twist`, cap torque against twist, `γ`, solver report
  (iterates, restorations, final `|C|`, factor time).

### 3.4 Presets

| Preset | Mesh | Constraints | Energies | Drive | Expect |
|---|---|---|---|---|---|
| `capsule_twist` | capsule `L/R = 4` | coupling, area, cap pins | plate `κ`, `μ` high (`γ ≈ 1e3`) | `dθ` per frame, small `dz` | helical diamond buckles past a critical twist; torque drops |
| `crease_tube` | same | same | `γ ≈ 1e5`, small `κ` | same | sharp creases, angle defect concentrated on them |
| `torus_hopf` | torus `R/r = 2`, `cos 3θ` bump | coupling, area | plate, `μ` high | `winding` shear, Uzawa on `Ω*` | ridges spiral along the Villarceau family; sign of `Ω*` picks the handedness |
| `fluid_pearl` | torus | coupling, area | Helfrich `c₀ > 1/r`, `μ = 0` | none | tube pearls after noise |
| `willmore_free` | torus `R/r = 3`, or sphere | coupling, total area (Uzawa) | Helfrich `c₀ = 0` | none | Clifford torus `R/r = √2`, or round sphere |
| `willmore_vesicle` | sphere | coupling, area and volume (Uzawa) | Helfrich `c₀ = 0` | reduced volume `v` | prolate, discocyte, stomatocyte by `v` |
| `willmore_conformal` | torus | coupling, length cross ratios (the legacy constraint, if wired), total area | Helfrich `c₀ = 0` | conformal class from a chosen flat torus | constrained Willmore tori; optional |
| `bilayer` | sphere | coupling, area, volume, `∫H dA` (Uzawa) | Helfrich `c₀ = 0` | `M*` | buds, pears, stomatocytes |
| `pipe_coil` | capsule | coupling, pipe `κ_hi = 1/r` per vertex, total area, cap pins | Willmore | cap twist | twist to writhe, plectonemes |
| `canal` | torus with noise | coupling, canal rows per vertex | Willmore, `μ` | none or twist | tubular undulations as exact canal surfaces |
| `growth` | sphere or torus | coupling, curvature-driven area | plate, `μ` moderate | 2.8 | ruffles at crests |

The wrinkle wavelength must cover about 6 edges or the mesh will alias it.
Set the remesher's target edge length from `sqrt(κ / μ)`, the shell thickness
scale, and show the ratio in the UI.

## 4. Milestones

Each one lands with its tests in `test/engman_tests.hpp` (built through
`tests/weingarten_tests.cpp`), then a demo checkpoint.

### Sprint 0: helical pipe first (do this before the milestones below)

The quickest route to a helical shape, built from a shape-space constraint
alone: a constant-radius pipe surface (2.12, candidate 2) twisted at its caps
buckles into a helix (Greenhill's problem), then into plectonemes. No metric
constraints, no remesh, no twist readout.

1. **Capsule mesh.** A closed capsule of length `L` and radius `r`, `L/r`
   about 10, 2 to 4k vertices, small noise.
2. **Skeleton.** The new demo with inertia, memory, store metric, the hard
   coupling, and plate bending (2.2). No remesh.
3. **Cap pins.** 2.4: rotate the caps by `±dθ/2` per frame about the axis,
   separation held. Readouts: cap torque (pin multipliers), total twist.
4. **Pipe constraint.** Per vertex `C_v = κ_hi(S_v) − 1/r`, `S_v` gathered
   from the incident edges through `vertex_stars`, `dκ/dS` from the Eigensatz
   rows (`MATH.md` 16). Exclude the cap vertices (the caps are spheres,
   umbilic). Start it soft through `engman::as_energy` at a large weight, then
   switch it hard.
5. **Expect.** Straight and twisting below a critical twist, then a helical
   centerline, then plectonemes. The torque drops at the buckle.

Sobolev preconditioning (Soliman) is not required here: the Engman step is an
implicit Newton-KKT step whose Hessian carries the bending pin through the
coupling, `J_xᵀ c J_x`, already a fourth-order operator, and the store metric
adds the `H¹` term on `S`. Watch iterate counts and restorations on the first
runs; if the step stalls, raise `refresh_every` to 1 before adding a
preconditioner.

Optional second toy on the same skeleton: the torus with the Villarceau shear
(2.3, `winding` axis) and no area constraint. That is the viscous helical
ridge the monge demo already shows, as a check that the drive carries over.

1. **Skeleton.** Torus, coupling, memory, inertia, store metric, no drive,
   no remesh. Idle at rest: drift below `1e-12` over 20 frames. Overlays
   for the coupling multiplier and the solver report.
2. **`face_area` and the hard/soft switch.** 2.1 with its tests, plus
   `engman::as_energy` with a test that soft at large `w` matches hard to
   `1/w²`. `fluid_pearl` preset reaches the pearling instability.
3. **Membrane and bending knobs.** `triangle_conformal` `μ`, plate and
   Helfrich bending, `γ` readout.
4. **Capsule and cap pins.** `capsule_twist` buckles; torque against twist
   recorded. This is the first "cool shape" checkpoint.
5. **Remesh with multiplier transfer and the persistent rest.** 2.6 with its
   tests, and the rest metric `ḡ_v` with `η` (2.10). With `η = 0` a twisted
   capsule held for 20 frames keeps its cap torque, and with `η = 1` the
   torque decays. At `η = 1` a fast enough twist rate still buckles, and the
   buckle stays when the drive stops. Record the threshold twist per frame
   against `η`. All presets run with `remesh_every = 1`.
6. **Hopf shear.** `edge_metric_target` with the `winding` axis, then the
   Uzawa twist controller. `torus_hopf` preset.
7. **Creases and growth.** `crease_tube`, the `det` energy, the
   curvature-driven area targets, `growth` preset.
8. **Willmore track** (2.11; independent of milestones 4 to 7, can follow
   milestone 3). Global area and volume Uzawa, then `willmore_free` with the
   Clifford-torus test (`R / r` within 1% of `√2` at 4k vertices), then
   `willmore_vesicle`. Then shape-space constraints (2.12): `bilayer`
   (global), the per-vertex wiring through `vertex_stars` tested by making the
   vertex Gauss sum hard, then `pipe_coil`, `canal`, linear Weingarten, and
   curvature caps last. The conformal preset and Ptolemy flips are optional
   and use the legacy cross-ratio constraint, not a new one.
9. **Audit.** For every row ported from `monge_fairing_demo`, record in
   `MATH.md` section 17 whether it does what it says under the hard coupling.
   That record closes TODO 18 phase 5.

## 5. Risks

- **Remesher against creases.** Flips and collapses smooth creases and can
  re-introduce `K`. Check `max_flip_dihedral` and the collapse threshold on
  `crease_tube`. Fall back to `remesh_every > 1` and a dihedral-aware split
  rule.
- **Hard area at high `μ` near locking.** Area plus stiff shear approaches
  isometry. If the capsule refuses to buckle, lower `μ` before suspecting the
  solver, and compare against the CR-midpoint isometry route
  (English and Bridson 2008) as phase 6 of TODO 18.
- **Quasi-definite conditioning.** Three hard constraint families at once
  (coupling, area, pins) with `ρ_rel = 1e6`. Watch restorations and the final
  `|C|`. If they rise, scale each family's rows to unit mean (as `a_ref`
  does for area) before touching `ρ`.
- **Twist readout near umbilics.** `Ω` is weighted by `|z|`. On a nearly round
  tube start, the readout is noise until the bump grows. The `γ_max` clamp
  keeps the controller from winding up.
- **Factor time.** `kkt_mode::full` costs 3.9 s at 8192 vertices (solver
  header). Interactive play wants a 2 to 4 k vertex default, with
  `refresh_every > 1` once milestone 4 is stable.
