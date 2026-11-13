# Two constraints: material twist from the curvature store, and surface tension with a volume multiplier

Companion to `docs/kkt_admm_curvature_solver_plan.md`. Both constraints fit
the existing `projection_constraint` / `block_constraint` interface and the
multiplier offsets of stage A. Both read the frame-start store `h0`
(`shell_cr_block::shape`) and act on positions. Neither puts a row on the
store's phase, which Codazzi rejects (`cr_curvature_hepworth_plan.md`
section 0, "twist turns the curvature lines only slightly").

Notation per edge `e`, frame `(t_e, b_e, N_e)` from `build_edge_frames`:
store `h_e = (k11, k12, k22)`, trace `τ = k11 + k22`, traceless part as a
complex number `z = (k11 − k22)/2 + i k12`, so `|z| = (κ₁ − κ₂)/2` and
`½ arg z` is the angle of the `κ₁` direction in the edge frame. Doubled-angle
complex numbers represent line fields: `u = e^{2iφ}` is the line at angle `φ`.

---

## 1. `material_twist`: read the phase in the curvature frame, shear in the material frame

### 1.1 What it does

Each frame it reads the axial line field and the twist rate from the store,
compares the measured twist rate against a target, and applies a small
material shear about the axis as an edge-length target relative to the
current shape. The twist accumulates in the embedding; nothing is wrapped
or stored as an angle. On a surface of revolution the shear is absorbed by
sliding (cross-sections rotate, shape unchanged); with azimuthal asymmetry
it spirals the ridges.

### 1.2 Frame-start build

1. **Axis line field.** Raw doubled-angle axis per edge `w_raw,e = −z_e / |z_e|`
   (the `κ₂` direction is at 90° to the `κ₁` direction, which is `−z` in
   doubled angle), weight `|z_e|`; edges with `|z_e| < z_floor` get weight 0.
   Smooth with the spin-2 CR Dirichlet operator (`build_spin2` applied to the
   2-vector `(Re w, Im w)`, mass `A_e / 3`):
   `(M + λ_a L₂) w = M w_raw`, then normalize `w_e`. `a_e` is the unit tangent
   at angle `½ arg w_e` in `(t_e, b_e)`, `b'_e = N_e × a_e`. The sign of `a_e`
   is arbitrary and nothing below depends on it.
2. **Radius proxy.** `R_e = 1 / max(|κ_max,e|, κ_floor)`, `κ_max = τ/2 + |z|`
   (the circumferential curvature on a tube). Clamp `R_e ≤ R_max`.
3. **Twist readout.** Local phase deviation
   `Δ_e = ½ atan2(Im(z̄_e w_e), Re(z̄_e w_e))` (degree 0 in `|z|`). Twist rate
   along the axis: with `∇_a z` the CR directional derivative of the store's
   traceless part along `a` in the face (the construction of the `dupin`
   rows in `edge_shape_fair`, applied to `z` with the spin-2 frame maps),
   `ω_e = Im(z̄_e ∇_a z_e) / |z_e|²`. Global:
   `Ω = Σ_e A_e |z_e| ω_e / Σ_e A_e |z_e|` (anisotropy-weighted so umbilics
   do not vote).
4. **Controller.** Shear increment per unit length for this frame:
   `τ_step = clamp(k_p (Ω* − Ω) + τ_ff, −τ_max, τ_max)`. `Ω*` carries the
   handedness. Per-edge shear strain `γ_e = τ_step R_e`.

### 1.3 Rows

A pure shear `γ` in the `(a, b')` frame is the metric `ḡ = I + γ (a⊗b' + b'⊗a)`.
Its exact effect on an edge of length `ℓ_e` at angle `ψ_e` to `a` is

```
ℓ_target,e² = ℓ_e² (1 + γ_e sin 2ψ_e)
```

so edges at +45° to the axis stretch, at −45° shrink, along and across the
axis unchanged. Rows are the existing edge-length projection: `A` the edge
difference, `p = ℓ_target,e d̂_e` after a rotation-free projection, weight
`w_twist sqrt(ℓ_e)` so the row is weight times geometry like every other
row. The rest is the current shape (frames are restarts), so `γ_e` is an
increment, not an accumulated strain.

Gradient/Hessian version (`gauss_newton() = true`): residual
`r_e = w sqrt(c_e) (|d_e|² − ℓ_target,e²) / ℓ_e`, Jacobian `2 w sqrt(c_e) d_eᵀ / ℓ_e`
on the two endpoints; multiplier offset from stage A applies unchanged.

Alternative row form, when the one-ring Procrustes machinery is preferred:
per vertex `G_v = I + γ_v (a⊗b' + b'⊗a)` in the vertex tangent plane,
identity along `N`, ring rest at the frame start; local step is the
Procrustes fit onto `G_v d̄_j`. Same energy to first order, more rows.

### 1.4 Knobs

| Knob | Meaning | Units |
|---|---|---|
| `Ω*` | target twist rate, sign is handedness | rad / length |
| `k_p`, `τ_ff` | controller gain and feed-forward | 1, rad / length |
| `τ_max` | cap on the per-frame increment | rad / length |
| `λ_a` | axis smoothing | length² |
| `z_floor`, `κ_floor`, `R_max` | umbilic and flat-region guards | curvature, curvature, length |
| `w_twist` | row weight, i.e. `ρ` under multipliers | 1 |

### 1.5 Tests

- **Invariance.** Scale every `|z_e|` by 2: `Δ`, `ω`, `Ω` unchanged. Flip
  `a_e` on random edges: rows unchanged to round-off.
- **Straight cylinder.** `n` frames of `τ_step`: shape unchanged to `1e-10`
  (pure sliding); a material fiber carried as a vertex datum turns by
  `n τ_step` times the axial length, read as the doubled angle to `a`.
- **Seeded tube.** Cylinder with an azimuthal bump (`cos 3θ` radial
  perturbation) and `Ω* ≠ 0`: ridges spiral; `Ω` tracks `Ω*` within 10%
  after 20 frames; reversing the sign of `Ω*` reverses the handedness.
- **Remesh.** Same with a flip or split every frame: `Ω` continuous across
  remesh events (the axis field is re-smoothed each frame, nothing carried).

### 1.6 Pitfalls

- A row on `h_e = k12 cos 2θ − ((k11 − k22)/2) sin 2θ` constrains
  `|z| sin 2Δ`; the solver meets it by scaling `|z|`, and it is non-monotone
  past 45°. Do not use it as a constraint; `Δ` and `ω` are the readouts.
- Near umbilics the axis is undefined; the `|z|` weighting and `z_floor`
  handle it, and the smoothing fills the gap from neighbors.
- The store's phase is never a row. The material shear is the only lever.

---

## 2. `surface_tension` and `enclosed_volume`: tension in the global step, volume as a global multiplier

### 2.1 What it does

Tension minimizes area; alone it collapses everything. Against a fixed
enclosed volume, area can only leave necks for bulges, which is the
pearling instability of a tube. With a spontaneous-curvature pin on the
store's trace (`edge_shape_fair` rows `pin` on `τ` to `c₀`, already built)
this is the Helfrich energy with tension and pressure:

```
E = σ Area(x) + (κ/2) Σ_e A_e (τ_e − c₀)² − p (V(x) − V*)
```

Pearling wavelength scales as `sqrt(κ / σ)`.

### 2.2 Surface tension rows (no local step)

Area equals the Dirichlet energy of the embedding when the operator is
built from the surface itself, and the Dirichlet energy of a map between
surfaces bounds area from above with equality when the map is conformal
(Pinkall and Polthier). So the frame-start Dirichlet energy is a majorizer of
area that touches it at the frame start, and minimizing it lowers area
every frame. Per face `f` with frame-start area `A_f⁰` and gradient basis
`∇φ_i` (the P1 hat gradients at the frame start):

```
E_D(x) = Σ_f (A_f⁰ / 2) ‖ Σ_i x_i ⊗ ∇φ_i ‖²_F   =  ½ xᵀ L_cot(x⁰) x
```

Rows: per face, `sqrt(σ A_f⁰ / 2) [∇φ_0 | ∇φ_1 | ∇φ_2]ᵀ` applied to each
coordinate (6 rows per face, 3 coordinates), `p = 0`. Written per face as
the Jacobian of the embedding the rows are PSD for any triangle, so obtuse
angles need no clamping (a per-edge `sqrt(cot)` form would). In
`projection_solver` this is a projection constraint whose `project` writes
nothing, so the term lives entirely in `AᵀA`; it is exactly one implicit
mean curvature flow step per frame (Desbrun's implicit fairing).

Weight: `σ` enters as `w² = σ`. Relative to inertia the dimensionless knob
is `σ h² / m`, so hold that fixed across resolutions.

### 2.3 Enclosed volume as a global multiplier

`V(x) = (1/6) Σ_f (x_0 × x_1) · x_2`, gradient per vertex
`∇_v V = (1/3) Σ_{f ∋ v} A_f n_f` (the area-weighted vertex normal). One
Gauss-Newton row:

```
r = sqrt(ρ_V) (V(x) − V*) + λ_V       λ_V += r after the step
```

Dense in one row (`1 × 3V`), handled by the stage A multiplier offset. The
multiplier `λ_V` is the pressure. It is a scalar, so it persists across
remeshes unchanged. Open-loop alternative: a prescribed pressure `p` as the
force `p ∇_v V` in `integrate_inertia`, no row; use this for inflation,
the multiplier for pearling.

Under the stage B solver the volume row is one dense column; use the
null-space trick (`Z` is that one column) rather than putting it through
the Schur PCG.

### 2.4 Knobs

| Knob | Meaning | Units |
|---|---|---|
| `σ` | tension | force / length, dimensionless as `σ h² / m` |
| `V*` or `p` | volume target, or open-loop pressure | length³, force / length² |
| `c₀`, `κ` | spontaneous curvature pin and its weight (existing `pin` rows on `τ`) | curvature, force · length |
| `ρ_V` | volume row weight | 1 |

### 2.5 Tests

- **Area descent.** Tension plus inertia only: `Area` decreases every frame
  on a noisy sphere; the frame-start Dirichlet energy equals the area to
  round-off at the frame start.
- **Sphere.** Fixed `V*`, tension on: a sphere is stationary to `1e-10`;
  `λ_V` converges to `2σ / r` (Laplace).
- **Rayleigh-Plateau.** Cylinder of radius `r`, radial perturbation of
  wavelength `λ`, fixed `V*`: `λ > 2πr` grows, `λ < 2πr` decays. With the
  `c₀` pin at weight `κ`, the fastest-growing wavelength moves with
  `sqrt(κ / σ)`.
- **Torus.** Fixed `V*`, tension on: the tube pearls after symmetry-breaking
  noise; without noise it stays a torus (a saddle).

### 2.6 Relation to the one-ring scale constraint

The "flat regions shrink their one-ring" half of the curvature-driven scale
rule is tension to first order (`σ ≈ k (1 − s)` for a uniform ring scale
`s`). Use the rows above for that half; keep the ring or edge-length
machinery only for the anisotropic part `G_v = diag(s(κ₁), s(κ₂))`, written
as the same edge-length target form as section 1.3 with
`ℓ_target,e² = ℓ_e² (s_a² cos²ψ_e + s_b² sin²ψ_e)`. Section 1 and this
anisotropic scale are one constraint type, a metric target from the store,
differing only in the strain tensor they build: shear `γ (a⊗b' + b'⊗a)` or
scale `diag(s_a² − 1, s_b² − 1)`.
