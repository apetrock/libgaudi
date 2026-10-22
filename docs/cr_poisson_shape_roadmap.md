# Shape-operator targets on L_cot and L_dir, as projective dynamics

## Premise

Prescribe the shape operator `S` and recover positions with constant, prefactored
Laplacians, keeping every nonlinearity in local steps, the way projective dynamics
(Bouaziz et al. 2014) and ARAP (Sorkine and Alexa 2007) do.

The Crouzeix–Raviart (CR) derivative of the embedding, `G = M^{-1} D X`
(Stein et al. 2020, eq. 7 and App. A.2, `reference/stein2020-vector-dirichlet-energy.pdf`),
gives two vectors per edge: the derivative of `X` along the edge and across it.
Their tangential half is what the vector Dirichlet energy measures; their normal
half is the shape operator. Per face, with frozen rest coefficients and normals,

`S_f = -sum_{e in f} (s_e l_e / A_f) [ (N_f . g_par,e) vperp vpar^T + (N_f . g_perp,e) vperp vperp^T ]`

is linear in `G`. So `S = f(D X)` with `f` linear once the normals are frozen,
and the inverse `g(S_goal)` is a linear least-squares solve.

Positions follow from rotated rest tangents by a Poisson solve. This is the loop of
Crane, Pinkall, Schröder 2013 (`reference/crane2013-conformal-willmore-flow.pdf`):
curvature, then tangents, then `Delta f = div T`. Wang, Liu, Tong 2012 reconstruct
frames from the second fundamental form and positions by a Poisson solve in the same
spirit.

## Operators

- `L_cot`: cotan Laplacian on positions; ARAP spoke rows (`A^T A = L_cot`).
- `L_dir`: vector Dirichlet on CR edge DOFs (`kusama::laplace::vector_dirichlet`),
  `2E x 2E`; regularizes the per-edge curvature components.
- `D0`, `M0`: the CR derivative (`2E x n` per coordinate) and edge mass `A_e / 3`.
  `D0^T M0 D0` is the CR Dirichlet operator on vertex functions.

## Formulation as PD constraints

Each term is a constraint with constant rows `A` (built from the step-start rest)
and a local projection `p`; the global step is
`(M/h^2 + sum w A^T A) q = (M/h^2) s + sum w A^T p(q)`, factored once per step.

| Constraint | Acts on | Rows `A` | Projection `p` |
|---|---|---|---|
| ARAP | `X` | spoke differences (`L_cot`) | Procrustes rotation times rest spokes |
| CR derivative | `X` | `sqrt(M0) D0` | rotated rest tangents plus normal curvature part `N c^T` |
| Curvature, linear | `c` | frozen map `c -> S` | `S_goal` |
| Curvature, Newton | local derivative | local derivative rows | closest point on `{g = S_goal}` (Lagrange–Newton with `g'`, `g''`) |
| Curvature smoothness | `c` | `L_dir` | 0 (on the change) |

`c_e` are the normal components of the two edge derivative vectors, an
intrinsic 2-vector per edge in the Stein edge frame, so `L_dir` applies to it.

`L0` is never inverted explicitly; the per-step factor is back-substituted.
Freezing operators within a step is first-order accurate in metric strain.

## Phases

### A. Roadmap and TODO — done
This document; `docs/TODO.md` item.

### B. Scratch base system in Python — in progress
`scripts/sympy/cr_poisson/`. One outer iteration:
1. read the current edge derivatives `G = D0 X` and their normal parts `c`;
2. c-step: `(M0 + alpha L_dir + lambda B^T B) dc = lambda B^T (S_goal - S(c))`, prefactored;
3. R-step: per-vertex Procrustes on cotan-weighted spokes against the ARAP rest;
4. X-step: `(w_a L_cot + w_d D0^T M0 D0 + eps M) X = w_a b(R) + w_d D0^T M0 F(R, c) + eps M X`,
   prefactored, with `F_e = R_e G0_tan,e + N_e c_e^T`.

Exit criterion: on isometric known answers (a flat strip bent to a cylinder of
radius r, an open cylindrical patch rebent from radius 1 to 1.5) the loop reaches
the analytic shape, agrees with a `least_squares` reference, and is idle when
`S_goal` equals the rest curvature. Non-isometric targets (a scaled torus) get
the ARAP / curvature compromise and are reported, not asserted.

### C. Composable with hepworth PD — built (part 2)

What exists and what it found: `docs/cr_curvature_hepworth_plan.md` section 0.
The edge curvature store, `weingarten_gn`, `dirichlet_memory_gn`, unified
`edge_shape_fair` and the uncoupled `kusama::weingarten_solver` are in; the
demo's default curvature solver is Weingarten. Open: stability under
remeshing (vertex growth and spin, as with the principal solver).

Implementation plan: `docs/cr_curvature_hepworth_plan.md` (live CR read-out
as Gauss–Newton rows coupled to the curvature block, CR-ARAP edge memory,
frames as restarts; standalone `kusama::cr_curvature_solver`). The original
sketch below predates the findings in `scripts/sympy/cr_poisson/MATH.md` §9–12.

The constraints above as hepworth constraints over the existing blocks, with
curvature and positions updated every PD iterate; the CR-derivative rows on the
positions block; `c` either a tracked block in `q` (with `L_dir` rows) or
presolved per step; optionally the Newton projection on `g`.

## Open decisions

- Read-out: CR derivative of the embedding (`D X`, gives full `S`) or of a tangent
  field (`N . D V = S V`, gives `S` applied to the field).
- Curvature components `c` as a tracked block in `q` or presolved per step.
- Whether the Newton projection on `g` replaces the Eigensatz Gauss–Newton rows.

## References

Stein, Wardetzky, Jacobson, Grinspun 2020; Crane, Pinkall, Schröder 2013;
Sorkine and Alexa 2007; Wang, Liu, Tong 2012; Bouaziz et al. 2014.
