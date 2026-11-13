# Curvature solver on multipliers: ADMM, then a KKT block solve

A plan to take the Weingarten coupling (`docs/cr_curvature_hepworth_plan.md`
section 0, `MATH.md` sections 9 to 16) off penalty rows and onto multipliers,
in two stages: first ADMM inside the existing `projection_solver` (small,
reversible, measurable in a day), then a saddle-point block solver that owns
the constraints and uses the spin-2 CR Dirichlet operator as the metric for
the store step and the dual update. The second question, whether the energy
should be written as the split `∫|∇W|² = ∫|∇H|² + ∫|∇W₀|²`, is answered in
section 4: the split is already built (`edge_shape_fair` rows `dirichlet` and
`dirichlet_trace`); what changes is how it is used.

**Revisions 2026-10-26** (the formulation is `MATH.md` section 17; the build
is the Engman plan, `include/gaudi/engman/`):

- Stage A (duals inside `projection_solver`) is skipped; the same experiment
  is Engman's first phase, so the `projection` path stays the reference.
- Every term declares a role: soft energies keep their weights (ADMM on them
  converges to the same weighted compromise, only faster), constraints get
  multipliers. Only the coupling is a constraint at first; shape edits,
  memory and inertia stay energies. No `bend` row in Engman.
- Exit criteria (1.5, 2.6) are restated: under a hard coupling an edit
  competes with inertia and memory, so "store on the pin to 1e-9" cannot
  hold. Check feasibility `|R(x) - S| -> 0`, `rho` independence and
  realizable targets instead.
- 1.4 is withdrawn: the memory is rotation invariant (normalized frames); a
  fixed `M + alpha L_vec` in `x` would penalize rotation, i.e. bending. Its
  runaway was step control, now bounded by backtracks and the ascent caps.
- The Newton-KKT step is the main inner step, not an escape: one
  `SimplicialLDLT` of the quasi-definite `[H J^T; J -I/rho]`, which makes the
  energies orthogonal to the constraints. Block Gauss-Seidel (2.2 a) is the
  fallback; the Schur PCG (2.3) is for scale.
- What persists is `lambda`, transferred as the intensive `lambda / c`
  (MATH 17.5); `S - R(x)` vanishes at feasibility.

References: Overby, Brown, Li, Narain 2017 (ADMM ⊇ projective dynamics);
Nocedal and Wright chapter 17 (augmented Lagrangian) and 16.2 (Schur
complement, null-space method); Stein et al. 2020 (CR operators); Wang, Liu,
Tong 2012 (reconstruction from fundamental forms); Soliman et al. 2021
(cross-ratio conformal constraint).

---

## 0. Diagnosis, in one paragraph

Every row in the solve is a penalty: projection rows `w A q = w p` and
Gauss-Newton rows `w r(q)`. With no multiplier anywhere, the solver can only
land on a weighted compromise, and each listed limitation follows from that:
the store settles at the blend `S = (w² R̄ + Σ w_k² g_k) / (w² + Σ w_k²)` (an
edit lands halfway); `k_coupling` has to be calibrated against `bend` so the
compromise is even; the coupling diverges above about `1e4 ×` inertia because
a penalty needs its weight large for accuracy and small for conditioning;
and the per-frame restart zeroes the only state that could accumulate. An
augmented Lagrangian keeps a moderate weight and gets exactness from the
multiplier; the multiplier is the state that persists.

## 1. Stage A: ADMM inside `projection_solver` (question 1, part 1)

Minimal change, same loop, same factor. Stage A is the thread adapter: it
answers whether the compromise was the problem before any new solver is
written.

### 1.1 Projection rows get a dual

Per Overby et al.: keep a dual vector `u` the size of the row space (`id0`),
project on the shifted argument, and feed the global solve the shifted
target.

```
// projection_solver::step, inner loop
project_all(q):   p_i = Π_i(A_i q + u_i)          // was Π_i(A_i q)
b = Mh s + Aᵀ (p − u)                              // was Aᵀ p
q = S->solve(b)
u += A q − p                                       // after the solve
```

`project` takes `q` today and computes its own `A_i q` internally; add an
optional per-row offset so `Π_i` sees `A_i q + u_i`. The simplest route is a
`vecX u` member of the solver and a `project(const vecX &q, const vecX &u,
vecX &p)` overload with a default that ignores `u` (every existing row keeps
working; rows opt in one at a time). The weight `w_i²` becomes the ADMM `ρ_i`
for that row: it sets convergence speed, not where the solution lands.

### 1.2 Gauss-Newton rows get a multiplier offset

The coupling `weingarten_gn` is the constraint `C_e(x, S) = R̄_e(x) − S_e = 0`.
Its augmented Lagrangian is `ρ/2 |C + λ/ρ|²`, so the row only needs its
residual shifted:

```
r_e  = w sqrt(c_e) (R̄_e(x) − S_e) + λ_e        // λ_e in the same scaled units
fill_J, add_gradient, residual_sq unchanged otherwise
after the inner loop:   λ_e += r_e (q_final)       // i.e. λ += ρ C
```

Implementation: `block_constraint` gains `std::vector<real> dual` sized to
its rows, applied inside `linearize` as `r += dual`, and a
`update_dual(const vecX &q)` hook the solver calls once after the loop.
`dirichlet_memory_gn` does not get one (section 1.4). `exponent_coupling`
gets the same offset on its curvature rows if it stays in play.

### 1.3 Persistence across the remesh

`shell_cr_block::prepare` resets `shape` to the live reading every frame and
the coupling starts idle. Keep the restart (every frame retopologizes) but
restart from transferred state. Write the store as `S_e = R̄_e(x) + δ_e` and
carry `δ` and `λ` as **vertex** datums: `datum_t<vec3>` with `prim_type::VERTEX`
already interpolates on split (`subdivide`) and merges on collapse
(`collapse`) by `va::mix`, while the edge branches are `/*not implemented*/`,
and flips leave vertices alone. Edge quantities are the least stable thing
under remeshing; vertex quantities ride the velocity code path.

- Representation: a symmetric 2×2 in each vertex's chart frame, stored as
  `vec3 (k11, k12, k22)`. `vertex_stars` and `vertex_chart` already define the
  edge-to-vertex polar; reuse them for the scatter and gather.
- Scatter (end of frame, before remesh): `δ_e = S_e − R̄_e(x)`, rotated into
  each endpoint's chart and accumulated with weight `A_e / 3`; same for
  `λ_e / c_e` (transfer the density, not the integrated value, or a split
  halves a multiplier and a collapse doubles it).
- Gather (`prepare`, after remesh): `S_e = R̄_e(x) + ½ (δ_v0 + δ_v1)` rotated
  into the edge frame; `λ_e = c_e · ½ (λ_v0 + λ_v1)`. Degenerate-face edges
  (`find_degenerate`) get `δ = λ = 0` as now.
- Damping: on the first `gn_refresh_every` iterates after a remesh, scale the
  gathered `λ` by `0.5` if an oscillation is seen at remesh events; do not
  zero it.
- Global multipliers (section 3.4) are scalars and persist as is.

The memory rest (`build_memory_rest` at the frame start) stays rebuilt from
the live shape each frame: it is a proximal anchor, nothing to carry.

### 1.4 The memory becomes the metric of the position step

`dirichlet_memory_gn` is a Dirichlet energy with its rest at the previous
frame, i.e. a proximal term `½ ‖x − x_prev‖²` in the norm `M + αL_vec`
pretending to be a row, which is why it can run away on large remeshed
meshes. Make it the metric of the step: add its (PSD, constant per frame)
matrix to `Mh` in the global system rather than to `AᵀA`, with the right-hand
side `(Mh + αL_vec) s` instead of `Mh s`. It then scales with the step by
construction, keeps "rough metric change costs more than smooth", and loses
its weight knob to `α` (units: `h²` times a stiffness, so it is
resolution-independent in the same way as inertia).

### 1.5 Exit criteria for stage A

| Check | Now | Target |
|---|---|---|
| `weingarten_store_is_the_blend_of_its_rows` | store at the blend | store on the pin to `1e-9` at any `w_coupling ≥ 1` (rewrite the test) |
| `weingarten_target_reaches_the_mesh_at_the_series_weight` | ratio `1.8` | mesh reaches the realizable target, independent of `w` |
| `weingarten_gn_torus_tube_over_frames`, 20 frames | `2e-12` | same or better in fewer frames; no `k_coupling` dependence |
| `weingarten_gn_cube_to_sphere_over_frames` | `6e-4` | `< 1e-6` |
| Stiffness ceiling | diverges above `1e4 ×` inertia | no divergence at `w = 1e2`, exact to `1e-9` |
| New: remesh persistence | state zeroed | torus tube target reached across 20 frames with a flip or split every frame, no visible flicker at remesh events |

If the first four hold, the compromise was the problem and stage B is worth
building. If they do not, stop and look at the coupling rows before writing
any solver.

## 2. Stage B: the saddle-point block solver (question 1, part 2)

A new class beside `projection_solver`, not a rewrite of it:
`hepworth/block/saddle_solver.hpp`. It keeps the `sim_block` DOF layout and
the `projection_constraint` row interface, and adds one new interface:
a constraint block.

### 2.1 Variables and blocks

| Block | DOFs | Hessian contributions |
|---|---|---|
| `x` positions | `3V` | inertia `Mh`, proximal metric `αL_vec` (1.4), `bend`, membrane, collisions, external forces |
| `S` curvature store | `3E` | shape rows (`edge_shape_fair` pins and goals: block-diagonal `3×3` per edge), fairing `L_S` (spin-2 element, sparse SPD) |
| `ℓ` metric store (section 3.1) | `E` | length goals; absent until buckling work starts |
| multipliers `λ` | one per constraint row | none |

Constraint blocks (each a `constraint_block` with `residual(q)`,
`jacobian(q)`, `dual`, `rho`):

- coupling `R̄_e(x) − S_e = 0`, `3E` rows, Jacobian `[J_x, −I]` (the symbolic
  `cr_curvature_edge_readings` Jacobian, as in `weingarten_gn::linearize`);
- metric `|x_j − x_i|² − ℓ_e² = 0`, `E` rows, or Soliman's cross ratio on each
  edge for a conformal class;
- global scalars (total Willmore, total `|W₀|²`, writhe): one row each, dense.

Shape edits stay soft energies on `S`. That is the design part and should
remain a compromise by intent, between the request and realizability.

### 2.2 The step

One frame is one augmented-Lagrangian solve with warm-started multipliers:

```
for outer k = 1..K_outer (2 to 4):
  for inner j = 1..K_inner:
    // primal step on (x, S) at fixed λ:  min E_x + E_S + Σ ρ/2 |C + λ/ρ|²
    option (a) block Gauss-Seidel: x-solve with cached factor, then S-solve
    option (b) one Newton-KKT step (2.3)
  λ += ρ C(x, S)          // preconditioned, 2.4
  ρ ← min(ρ_max, τ ρ) if |C| did not fall by a factor 4
```

Option (a) is stage A's loop with the store split out. Option (b) is what
makes the blocks stop fighting; build (a) first, keep it as the fallback.

### 2.3 The Newton-KKT step and its Schur complement

With `H = diag(H_xx, H_SS)` and `J = [J_x, −I]` for the coupling, the KKT
system is

```
[ H_xx    0     J_xᵀ ] [dx]   [−g_x]
[ 0     H_SS   −I   ] [dS] = [−g_S]
[ J_x   −I      0   ] [dλ]   [−C  ]
```

Eliminate `dS` and `dx`:

```
(J_x H_xx⁻¹ J_xᵀ + H_SS⁻¹) dλ = J_x H_xx⁻¹ g_x − H_SS⁻¹ g_S − C
```

- `H_xx⁻¹` is the frozen Cholesky `projection_solver` already keeps
  (`_cache_solver`): two triangular solves per matvec.
- `H_SS⁻¹` is a `3×3` block solve per edge when fairing is off and a sparse
  Cholesky of `3E × 3E` when it is on (`build_spin2` plus the pin diagonal).
- The dual system is SPD and lives on edge rows. Solve it with PCG. Its
  natural preconditioner is `M_S + βL_S`, the spin-2 CR Dirichlet operator
  (section 4.3), prefactored once per frame.
- Metric rows add a second block in the same pattern. Global rows (a handful
  of dense columns) go through the null-space trick: `Z` is one column per
  global constraint, so project the step onto `null(J_global)` explicitly.

Cost per Newton step: one PCG on `3E` unknowns with a prefactored
preconditioner and `O(10)` matvecs, each matvec two triangular solves on
`3V`. That is the same order as one iterate of the present loop; expect
interactive rates at a few hundred thousand edges without multigrid.
Multigrid on tensor fields over unstructured, retopologizing meshes is its
own project; keep it out until the refactor cap bites.

### 2.4 The dual update, preconditioned

`λ += ρ C` is gradient ascent on the dual with conditioning
`J H⁻¹ Jᵀ`. Preconditioned: `λ += ρ (M_S + βL_S)⁻¹ C`. One triangular solve
pair on `3E` per outer iterate. This is the piece that lets option (a) of
2.2 compete with the Newton-KKT step on most frames, so the expensive step
is only taken when the preconditioned ascent stalls (measured as `|C|` not
falling by `4×` in an outer iterate).

### 2.5 What the integrator looks like

Unchanged in spirit: implicit Euler with the momentum prediction `s` as the
proximal anchor in the `x` block (`Mh (q − s)`), velocities read back by
`map_from_x`. What changes is that the constraint rows are solved to
feasibility within the frame instead of to a weighted blend, the store and
duals are warm-started from the vertex datums (1.3), and `ρ` is a conditioning
parameter with no accuracy meaning.

### 2.6 Exit criteria for stage B

- Stage A checks pass under `saddle_solver` with option (a); option (b) then
  halves the inner iterates on cube to sphere at equal accuracy.
- Realizable targets (torus tube, torus `R`, cube to sphere, cube to
  ellipsoid) reach `1e-12` in one frame, no restarts, from `λ = 0`.
- Frozen-Hessian warm restart: a second frame on an unchanged target takes
  `< 5` inner iterates.
- A frustrated target on a fixed metric buckles (3.1): flat plate with a
  sphere-cap read-out and isometry rows gives the cylinder ratio
  `q_min / q_max < 0.1` as the exponent solver did at `km = 1e2`, without
  that solver.

## 3. Modeling additions the solver makes possible

### 3.1 The metric store (needed for buckling)

With a free metric, more curvature is cheapest as uniform shrink, so there is
no curvature-domain buckling. Add a per-edge length store `ℓ_e` and the
isometry constraint `|x_j − x_i|² = ℓ_e²` as a constraint block (or the
cross ratio for a conformal class, which is the knob Soliman turns). A store
edit that changes `det W` is then unrealizable (Gauss) and lands on a
compromise; one that changes anisotropy at fixed `det` must be realized by
bending, which is where the modes appear. The symmetry-breaking noise is
still needed; a multiplier does not remove a saddle point. The Gauss check is
cheap and worth logging: angle defect of the live mesh against `det S` per
vertex.

### 3.2 Codazzi as a soft term on `S`

At convergence `S = R̄(x)` so Gauss-Codazzi hold implicitly. A soft Codazzi
residual on `S` (`div W₀ − ½∇H` per face, CR gradient of the store) keeps the
store subproblem near integrable fields during the inner iterates and is the
right regularizer for an unrealizable edit. Low priority; add after 3.1.

### 3.3 Dupin and directional smoothness

`edge_shape_fair` rows `dupin` already hold each principal curvature's
derivative along its own line with a frozen principal frame. Under the split
of section 4 that is a directional restriction of `∫|∇W₀|²`; nothing changes
except that the same rows can be weighted against the full traceless Dirichlet
term rather than against a pin.

### 3.4 Global constraints

Total Willmore, total `|W₀|²`, writhe: one scalar multiplier each, dense
gradient on `S` or `x`, matvec only. These persist across remeshes
unchanged, so push design intent into them where it is global by nature
(an excess Willmore energy above the minimum at fixed conformal class is the
cleanest buckling driver).

## 4. The energy split (question 2)

### 4.1 What exists

`MATH.md` section 14 already writes the smoothness energy as
`∫|∇z|² + ½ ∫|∇τ|²` with `τ = k11 + k22` (spin 0) and `z = ((k11 − k22)/2,
k12)` (spin 2), and `edge_shape_fair` builds both: rows `dirichlet` are the
full spin-2 element `E_f` (`spin2_face_element`), rows `dirichlet_trace` are
`½ τᵀ K τ` with `K = cr_face_stiffness`, which equals the trace part of `E_f`.
So the answer to "should we modify the energy to follow the split" is: it
already does, implicitly. What is missing is the traceless part on its own
and separate weights.

### 4.2 Change: expose `dirichlet_traceless` and two weights

Add `edge_fair_rows::dirichlet_traceless` built as `E_f − ½ Kᵀ⊗(trace)` (or
directly from the spin-2 frame maps with the trace projected out), and let
callers weight the two parts independently:

```
E_fair(S) = λ_H ∫|∇H|²  +  λ_z ∫|∇W₀|²
```

Readable meaning: `λ_H` smooths mean curvature; `λ_z` smooths the principal
direction field weighted by anisotropy `(κ₁ − κ₂)/2`, so near umbilics the
directions are free to turn. Dupin rows are the directional restriction of the
`λ_z` term. Test: on the analytic torus, `dirichlet_trace + dirichlet_traceless`
equals `dirichlet` to round-off; the traceless element is zero on any
constant-direction, constant-anisotropy field and on the umbilic sphere.

### 4.3 Use it as the metric, not mainly as the energy

Two cautions on the energy use. First, Simons' identity ties `∫|∇W|²` on an
integrable field to the Hessian of `H` plus Willmore-like terms, so once the
coupling holds the Dirichlet term mostly flattens `∇²H`; its unconstrained
minimizers (parallel second fundamental form: planes, spheres, cylinders) are
much narrower than the surfaces that appear once `S = R̄(x)`. Keep it as a
weak prior. Second, the fairing numbers in `MATH.md` section 14 show
over-smoothing on cube to sphere. The heavy use is as a metric:

- the store's proximal step in 2.2 (a) and the `H_SS` block in 2.3 use
  `M_S + λ_H L_τ + λ_z L_z`;
- the dual preconditioner in 2.4 is `M_S + βL_S` with the same split, because
  `λ` is an edge curvature quantity and this is its natural inner product;
- a Sobolev step on `S` (solve `(M_S + αL_S) g = ∇_S E`) gives
  mesh-independent step sizes for Willmore-type edits, which is the
  `H²`-flow argument one derivative down since `S` is already two derivatives
  below `x`. (The `sobolev` step that was tried and removed acted on `x`
  with an `H¹` operator, which is the wrong order for a fourth-order energy;
  on `S` it is the right order.)

### 4.4 Decision

Modify the energy only by exposing the split with two weights (one new row
kind, no new element). Modify the solver to use the same operator as the
metric for the store step and the dual update. Do not raise the energy weight
to get smoothness; raise `λ_z` in the metric and let the coupling decide the
shape.

## 5. Files

| Piece | File | Change |
|---|---|---|
| ADMM duals on projection rows | `hepworth/block/solver.hpp` | `vecX u`, shifted `project`, `u` update after solve |
| multiplier offset on GN rows | `hepworth/block/block_constraint.hpp`, `weingarten_gn.hpp` | `dual` vector, `r += dual`, `update_dual` |
| memory as proximal metric | `solver.hpp`, `weingarten_gn.hpp` | memory matrix added to `Mh`, rhs `(Mh + αL) s` |
| persistence | `hepworth/blocks/shell_cr_block.hpp` | vertex datums `shape_offset`, `shape_dual`; scatter in `commit_edge_tangents`'s sibling, gather in `prepare` |
| traceless rows | `hepworth/block/edge_shape_fair.hpp`, `kusama/weingarten.hpp` | `dirichlet_traceless`, `spin2_traceless_face_element` |
| saddle solver | `hepworth/block/saddle_solver.hpp` (new), `constraint_block.hpp` (new) | stage B |
| metric store | `shell_cr_block.hpp`, `hepworth/block/edge_length_constraint.hpp` (new) | 3.1 |
| demo wiring | `duchamp/monge_fairing_demo.hpp` | solver choice, `ρ`, `λ_H`, `λ_z` knobs |
| tests | `test/weingarten_tests.hpp` | section 1.5, 2.6, 4.2 |

## 6. Order of work

1. Stage A, 1.1 and 1.2: duals on projection rows, offset on the coupling
   row. Rewrite `weingarten_store_is_the_blend_of_its_rows` as the pin test.
   Measure 1.5.
2. 1.4: memory to proximal metric. Confirm the runaway on large remeshed
   meshes is gone.
3. 1.3: vertex datums for `δ` and `λ`; remesh persistence test.
4. 4.2: traceless rows and the two-weight fairing; parity test against the
   full element.
5. Stage B option (a): `saddle_solver` with block Gauss-Seidel, preconditioned
   dual update (2.4) using 4.3. Stage A tests pass under it.
6. 3.1: metric store and isometry rows; buckling exit criterion.
7. Stage B option (b): Newton-KKT step with the Schur PCG. Compare inner
   iterate counts against (a) and keep (b) as the stall escape.
8. 3.4 global multipliers; 3.2 Codazzi prior; multigrid only if the refactor
   cap is the bottleneck.

## 7. Risks and open questions

- **Duals on nonconvex projections.** The SVD projectors are prox operators of
  nonconvex sets; ADMM has no convergence guarantee there, only the empirical
  record of Overby et al. If a projector row oscillates with a dual, cap its
  `u` or leave that row without one; the coupling row is smooth and is the one
  that matters.
- **`ρ` versus `bend`.** `k_coupling` was calibrated so the coupling and
  `bend` have equal position stiffness. Under a multiplier the coupling's
  stiffness is `ρ`, a conditioning choice; start it at the calibrated value
  so the first inner iterate matches today's behavior, then let the adaptive
  rule move it.
- **Remesh warm start quality.** A transferred `λ` is the dual of a slightly
  different problem. Expect a few re-convergence iterates per remesh event;
  the damping in 1.3 is the fallback. Measure `|C|` at the first iterate after
  a remesh against the last before it.
- **The exponent solver.** It already reaches `1e-14` on realizable targets
  and buckles under a stiff `km`. It is a position-side subproblem solver
  (rotations through the exponential), and could serve as the `x` block in
  2.2 (a) where Gauss-Newton on positions stalls. Keep it as an option, not a
  dependency.
- **Store frame on vertex datums.** Vertex charts rotate with the mesh;
  scatter and gather must use the same chart convention as
  `vertex_frame_cache_compute`, or a transferred tensor picks up a spurious
  rotation. Test on a rigidly turned torus: transferred `δ` round-trips to
  `1e-12`.
