# Local incompressible material (φ on slab stencil)

**Status:** Active design (2026-08-31).  
**Patch:** top one-ring + bottom one-ring + pierce spine + hub normals `n_top`, `n_bot`.  
**Valence:** **n_top** and **n_bot** come from the mesh (`get_one_ring`); they need not match (e.g. triangle foot, hex pierce, irregular rings). All sums below are over actual stencil size, not fixed 6+6.  
**Supersedes** RBF live-refit and “frozen centers” sketches for shape-field bending.

---

## Model (one sentence)

Fit **Q_r** once (LS quadric + normal gauge); live state is **incompressible** deformation **φ = id + u**, **div u = 0**, on the slab **Ω** between the two rings; relax control **q** with NR using analytic **∇f**, **H** and wedge volume; emit ARAP spoke targets.

---

## Rest (init once per foot)

| Quantity | How |
|----------|-----|
| **Q_r** | Normal-constrained LS quadric on ring samples + foot row (`∇f·n_top = 1`) |
| **x^rest** | Ring + hub positions from block `q` at `init_rest` |
| **V₀** | Sum of fan-wedge oriented volumes **τ^rest** (same topology as trivector path) |
| **A₀^top, A₀^bot** | Cap areas from top/bottom ring polygons (optional area term) |
| **φ** | Identity |

---

## Live parameterization: incompressible φ

```text
φ_q(x) = x + u_q(x)
div u_q = 0     (exact linear incompressibility on Ω)
```

**Parameterize u** (pick one for v1):

1. **Solenoidal:** `u = curl Ψ`, **Ψ** linear in foot-local coords, coeffs = **q** (3–6 DOF after div-free).
2. **Div-free basis:** precomputed nullspace of discrete div on wedge interior + boundary; **u = B q**.
3. **Micro-FFD projected:** trilinear lattice velocity, project each control increment onto div-free subspace.

**Do not** refit **Q** on live geometry. Live field:

```text
f(x) = f_r( φ_q(x) ; Q_r )
S_live = { x : f(x) = 0 } = φ_q⁻¹( S_rest )
```

---

## Volume conservation

**Primary (v1):** linear **div u = 0** ⇒ first-order volume preservation on Ω.

**Consistency check:** wedge volumes **V_k(q) = (det Dφ_q)|_{x_k} · V_k^rest** (or compare **τ** to **τ^rest** with deformed geometry) — same samples as trivector wedges.

**Later (optional):** penalty or multiplier on **det Dφ_q − 1** for finite strain.

---

## Local energy (NR on q)

Sum over ring samples **i** and wedge quadrature points; SymPy → **∂E/∂q**, **∂²E/∂q²**.

```text
E(q) = μ E_gauge(q) + ν E_area(q) + η E_strain(q)
```

| Term | Definition |
|------|------------|
| **E_gauge** | Σ w_i ( \|∇ΔD_i\|² + λ_gap ΔD_i² ), **ΔD = f_r(φ_q(x))** (target 0) |
| **E_area** | (A_top(q) + A_bot(q) − A₀)² **or** co-area on samples **Σ \|∇f\| σ(f/ε)** |
| **E_strain** | Σ \|sym(Dφ_q − I)\|² on wedge points (optional; small-strain regularizer) |

**Constraints:** **C q = 0** enforces **div u = 0** (if not built into basis).

**Solve:** Newton / Gauss–Newton on **q** → **q*** → boundary **x*** → Hepworth **p = w·(x* − hub)**.

No marching cubes in the solver loop (MC optional for viz only).

---

## Relation to other paths

| Path | Role |
|------|------|
| **Trivector wedge spring** (`τ − τ^rest`) | Discrete volume memory; can be **consistency check** or light penalty alongside **div u = 0** |
| **RBF / live Q_k refit** | Deprioritized — gauge mismatch, not incompressible material |
| **Shape-field ΔD with two fits** | Replace live refit with **f_r ∘ φ_q** |

---

## Implementation order

1. **`scripts/sympy/slab_phi_incompressible.py`** — quadric **f_r**, linear **u_q**, **ΔD**, **E_gauge**, **∂E/∂q**; FD check.
2. **Div-free basis** on the actual stencil (**n_top**, **n_bot** wedges; template or per-valence precompute).
3. **C++ bundle** `incompressible_phi_field` next to `darboux_field_bundle` (eval **f**, **∇f** via chain rule).
4. **Constraint** replacing or augmenting `shape_field_divergence_strain` inner solve: NR on **q**, then spoke targets.
5. **Three.js probe** — drag foot, step NR, draw **f_r∘φ_q = 0** (variable **n_top** / **n_bot** UI).

---

## Appendix: wedge div-nullspace basis (v1)

At `init_rest`, for each foot:

1. Collect ring verts (top + bottom hubs and one-rings); **n_dof = 3 × n_verts**.
2. Build fan wedges via `build_slab_fan_wedges`; one discrete **div u = 0** row per wedge at the barycenter (P1 barycentric grads).
3. **B = ker(D)** from SVD of **D**; **dim(q) = n_dof − rank(D)** (typically ~6–15).
4. **Q_r** and **B** share the same fixed foot origin **`foot_rest`** (world position at init). All gauge samples use **local = x − foot_rest**, not the live hub.

Runtime: **u(vi) = B vi × q** at ring verts; **x* = x + u** after GD on **q**.

---

## References

- Hirota et al., volume-preserving FFD (integral VP context): `10.1145/304012.304036`
- nLab: volume-preserving diffeomorphism (pointwise **det Dφ = 1** — nonlinear follow-up)
- Existing: `quadric.py`, `slab_field_energy.py`, `slab_wedge_strain.hpp` (wedge **V₀**)
