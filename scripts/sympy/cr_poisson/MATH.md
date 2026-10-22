# X = L^{-1} D^T M F(J^{-1}(S_g)): shape-operator targets through the CR derivative

Source: Stein, Wardetzky, Jacobson, Grinspun 2020, *A Simple Discretization of
the Vector Dirichlet Energy*, eq. (9) and App. A.2
(`reference/stein2020-vector-dirichlet-energy.pdf`). Code it mirrors:
`hepworth::block::cr_coordinate_gradients` and `cr_face_operators`
(`include/gaudi/hepworth/block/cr_shape.hpp`), ported in `cr_ops.py`.

## 1. The shape operator in the paper

Continuous (eq. 9): the shape operator is the normal part of the covariant
Hessian of the embedding coordinates,

    S(v, w) = -N . (grad grad X)(v, w).

With `T = dX` the tangent frame field, `S = -N . grad T`.

Discrete (App. A.2), two parts:

1. Gradient of each coordinate from hat functions to CR tangent vectors,
   `G = M^{-1} D(X)`. Per edge `e`, stacked over x, y, z,
   `G_e = [g_par,e, g_perp,e]` (3 x 2) is the tangent frame at the edge
   midpoint: the derivative of `X` along the edge and across it. `M` is the
   CR mass, `M_e = A_e / 3` (`A_e` the area of the incident faces).
2. The CR covariant derivative of that frame field, normal part only. Per face,
   over its edges `e_i` (local orientation, `s_i = +1` when it agrees with the
   global orientation),

       S_f += (s_i |e_i| / |f|) [ (N_f . g_par,i) v_perp,i v_par,i^T
                                + (N_f . g_perp,i) v_perp,i v_perp,i^T ],

   with `v_par,i` the unit edge vector, `v_perp,i = N_f x v_par,i`. The code
   (and this note) uses `K_CR_SHAPE_SIGN = -1` so that `S = -dN`.

Two facts from the formula:

- `g_par,e` is the edge's own unit direction (both incident faces give
  `(x_head - x_tail) / |e|` up to `s`), so `N_f . g_par,e = 0` and the
  `v_par` term vanishes for every edge, which is the paper's remark that it
  is zero away from the boundary. Verified symbolically in `derive.py`.
- Write `n_f in R^6` for the face's normal components
  `(N_f . g_par,i, N_f . g_perp,i)` over its three edges and `s_f in R^4` for
  the coefficients `(T S T, T S B, B S T, B S B)` in the face frame `(T, B)`.
  Then `s_f = K_f n_f` with `K_f` (4 x 6) built from `s_i |e_i| / |f|` and the
  outer products in `(T, B)`.

So the read-out is

    S = J(X) D X,     J(X) = blockdiag_f( K_f(X) . (N_f(X)^T on each edge column) ),

linear in the frame `G = D X`, with `J` depending on `X` through the face
normals, edge directions, and `s_i |e_i| / |f|`.

## 2. What J can and cannot see

`J` reads only the normal components `N_f . G_e`. The tangential part of the
frame is invisible to `S`. A target `S_g` therefore fixes at most the normal
components; the tangential part has to come from somewhere else. Here it comes
from the rest frames rotated with the mesh (a continual ARAP of the rest
tangents):

    R_e = argmin_{R in SO(3)} |G_e - R G0_e|_F = polar(G_e G0_e^T),
    T_e = R_e G0_e,

with `G0 = D X0` the rest frames. A rigid rotation of a neighborhood rotates
`G` and `N` together, so `S` is unchanged; only the relative rotation between
neighboring frames (the connection) changes `S`, which is what `R_e` tracks.

## 3. Energy

    E(X) = 1/2 |D X - F(X)|^2_M,

where `F(X)` is the closest frame field to `D X` in the set

    C(X) = { F : tangential part of F_e equals T_e,  K_f n_f(F) = s_g,f for every f },

with `n_f(F)` the normal components of `F` against the live face normals.

### Local step: the inverse of J

Because `N_f . g_par = 0` (section 1), only the three perpendicular normal
components `c_f = (N_f . g_perp,i)_{i=1..3}` enter, each through
`v_perp,i v_perp,i^T`. Those three rank-one matrices span the symmetric
2 x 2 matrices for any non-degenerate triangle, so the read-out restricted to
a face is a 3 x 3 invertible map `K3_f` from `c_f` to `(k11, k12, k22)`:

    c_f* = K3_f^{-1} s_g,f.

The local inverse of `J` is exact and unique per face; there is no
least-norm freedom left. (The general closed form
`n* = n + K^+ (s - K n)`, verified in `derive.py`, reduces to this.) It is
the Newton step for `J`: one step is exact because the map is linear once
`K3_f` is fixed; `K3_f` itself is re-linearized by the outer iteration. The
correction to edge `e` from face `f` lies along `N_f`:
`dF_e^(f) = N_f (c_f* - c_f(T))_e` on the perpendicular column.

### Recombination

Each interior edge frame is shared by two faces with different normals. Two
options:

- Local (default): `F_e = T_e + sum_{f ni e} w_fe N_f (n_f* - n_f(T))_e^T`,
  `w = 1/2` for interior edges, 1 on the boundary. Face `f` then reads
  `N_f . F_e = n_f(T)_e + 1/2 dn_fe + 1/2 (N_f . N_g) dn_ge`; the error
  against `n_f*` is `1/2 ((N_f . N_g) dn_ge - dn_fe)`, small where
  neighboring faces ask for similar corrections and the dihedral is small.
- Exact (global): per-edge corrections `d_e` (3 x 2) jointly,
  `min |d|^2_M s.t. K_f n_f(T + d) = s_g,f` for all faces, a sparse
  least-norm problem with a live matrix (not prefactored).

### Global step

With `F` fixed, `E` is quadratic in `X`:

    L X = D^T M F,     L = D^T M D.

`D` (hence `L`) uses only edge lengths and areas, so it is intrinsic: frozen
at the rest, it stays exact under isometric bending. `L` has constants in its
null space (translation), fixed by a small `eps I` with `eps X_k` on the right.

### Derivative and descent

    dE/dX = D^T M (D X - F) - (dF/dX)^T M (D X - F).

If `C` were a fixed set, the second term would vanish (the residual
`D X - F` is normal to `C` at its closest point, envelope theorem),
stationarity would be `L X = D^T M F(X)`, and the local/global iteration
would decrease `E` monotonically (PD's argument). Here `C(X)` moves with `X`:
the ARAP rotations are themselves closest points and are covered, but the live
normals and `K_f` are not. So the iteration is a fixed-point scheme whose
descent is checked, not guaranteed. That is the frozen-Hessian situation:
keep `L` (and `G0`) factored, evaluate `E` every iterate, and rebuild at the
current `X` when `E` does not decrease.

    X = L^{-1} D^T M F( J^{-1}(S_g) ),   iterated.

## 4. Exact versus first order

| Piece | Status |
|---|---|
| Read-out `S = J(X) D X` | exact (the paper's discretization) |
| `L` frozen at rest | exact under isometry; first order in metric strain otherwise |
| ARAP tangential target | closest point, exact per edge |
| Normal correction per face | exact inverse of the live `K3_f` (3 x 3) |
| Local recombination | first order in neighboring-correction mismatch |
| `K_f`, normals within an iterate | frozen, re-linearized every iterate |

## 5. Integrability

A frame field is the derivative of an immersion only if it is integrable
(Gauss–Codazzi). `L X = D^T M F` returns the least-squares immersion for any
`F`; an `S_g` without a realizing surface (with the rest metric) gets the
compromise. Isometric targets (a flat strip bent to a cylinder, an open patch
rebent to another radius) are reachable; a torus scaled or with a new tube
radius changes the metric and is not, unless `G0` is reseeded.

## 6. Rebuild levels

- `L` only: current metric, rest tangents `G0` kept (ARAP memory survives).
- `L` and `G0`: reseed, as the demo reseeds per frame.

## 7. Measured: frame targets with frozen normals cannot reach `S_g`

Strip (17 x 6, 2 x 1) flat to cylinder r = 1.5, `curvature_residual` 0.865
flat, 0.00025 on the exact cylinder.

| Step from flat | curvature residual |
|---|---|
| local recombination, `F` read out | no change (E ~ 1e-11) |
| global least-norm `d`, `lam` 1e6, `F` read out | 0.825 |
| global, then `X = L^{-1} D^T M F` | 0.840 |
| one Gauss–Newton step on `X` (frozen `B`) | 0.133 |

Why. `c_fe = N_f . g_perp,e` is a per-(face, edge) quantity, two values per
interior edge, while `g_perp,e` is one vector per edge. A face reads its
edges with the orientation sign `s_i`, so at a flat start the two faces of
an edge ask for opposite corrections along the same normal: the local average
cancels exactly (section 3's mismatch term is the whole correction), and the
global solve has one effective unknown per edge against three equations per
face (`nE ~ 1.5 m` against `3 m`). Curvature in this read-out lives in the
face normals tilting against the shared edge frame, i.e. in how `N_f` moves
with `X`, which a frame-only target with frozen normals discards.

Consequence: the local inverse of `J` has to move the faces, not the edge
frames. Two forms that keep the PD shape:

- Per-face stencil projection: face `f` and its three opposite vertices
  (the support of `S_f`), project onto `S_f(p) = s_g,f` by local
  Newton–Raphson (3 equations, min-norm over the stencil), global step over
  the stencil differences. This is PD with a nonlinear local projection.
- Per-face rotations: unknown `R_f` per face, edge targets from the rotated
  rest edges, global step `L_cot X = div(R T0)` (Poisson on rotated
  frames). The curvature is the relative rotation of neighboring faces.

## 8. Both face-moving forms (`face_loops.py`)

### Rotations

The per-face inverse `c*_f = K3_f^{-1} s_g,f` is a pair of tilts per interior
edge, one from each face, against the shared `g_perp`. Their difference is
the dihedral, and the dihedral is the relative rotation of the two faces.
Both hold exactly on the cylinder (residual 1e-16 and 5e-6):

    theta_e = asin(c_+ / |g_perp,e|) - asin(c_- / |g_perp,e|),
    R_- = R_+ Rot(a0_e, theta_e),

with `+` the face whose winding runs along `E[e, 0] -> E[e, 1]`. One iterate:

1. `theta*_e` from `c*` (exact local `J^{-1}`, live `K3`, live `|g_perp|`).
2. Rotation sync, linear least squares then polar per face:
   `min sum_e |R_- - R_+ P_e|^2 + mu sum_f A_f |R_f - Rproc_f|^2`,
   `P_e = Rot(a0_e, theta*_e - theta0_e)`. The matrix depends on `P_e` only,
   so it is fixed while the target is.
3. Per-face ARAP Poisson `L_cot X = div(R T0)`, prefactored at rest.

Step 3 can use the CR Dirichlet matrix instead (`global_step="cr"`):
`D = sum_f D_f` splits by face, each `D_f` is rigid-equivariant, so
`D^T M D X = D^T M sum_f R_f D_f X0` has the same fixed point when the faces
are rigidly rotated. Measured: identical to `L_cot` (rigid error 2.561e-4, the
same ~0.5 contraction per iterate from flat and from the inverted strip), so
the rate is set by the local step, not the global matrix. The 2.5e-4 floor is
the reference cylinder itself: it wraps arc length, so its chords are shorter
than the flat edges and it is not an exact isometric image of the strip
(`sum_f R_f D_f X0` against `D X_cyl` differs by 2.9e-4).

This is surface reconstruction from discrete fundamental forms (lengths
through `T0`, dihedrals through `P_e`), with the dihedrals read through the
CR shape operator instead of prescribed.

### Stencil projection

Projection by Newton–Raphson converges quadratically per face (3 steps to
1e-15). Boundary faces have `N_f . g_perp = 0` on their boundary edge
identically (that `g_perp` comes from the face alone), so their `J` has rank 2
and the step must be a rank-truncated pseudo-inverse; the plain min-norm step
blows up as soon as the stencil leaves the exact flat. With that, PD is
monotone, but the projected stencils keep their current tilt, so bending
travels by diffusion through the overlaps.

### Strip flat to cylinder r = 1.5 (17 x 6)

| Form | Iterates | Rigid error | Curvature residual |
|---|---|---|---|
| reference least squares (`reference.py`) | — (98 s) | 2.7e-4 | |
| first loop, ARAP competing (`base_loop.py`) | converged | 4.4e-3 | |
| frame targets, frozen normals (`cr_loop.py`) | stuck / diverges | 0.21 | 0.86 |
| rotations, `mu` 1e-2 | 12 to 3e-4 curvature, ~0.5 per iterate | 2.56e-4 | 7.8e-5 |
| stencil projection, `eps` 1e-3 | 3000 (53 s) | 6.1e-3 | 2.8e-3, flattening |

Rotations reach the same fixed point (2.56e-4) from flat, from r = 1, and
from the strip bent the other way (r = -1.5); it is unchanged for `mu` from
1e-4 to 0.1 and biased at `mu = 1` (4e-2), where the Procrustes tie competes
with the target. The exact cylinder reads 2.5e-4 against its own target
because the cylinder is not exactly the read-out's fixed point.

## 9. The goal energy on tangents (`cr_goal.py`)

Section 7's roadblock came from freezing the normals at `X`. Reading the
whole operator from the tangent field removes it: per face, the oriented edge
vectors are `v_i = s_i l0_i g_par,i`, the normal and frame come from them, and

    S_f(T) = K sum_i s_i |v_i| / A(v) (N(v) . g_perp,i) v_perp,i v_perp,i^T.

With `T = D X` this is the embedding's read-out (3e-15). With `T` free, a face
tilts by changing its `g_par`, so `S(T) = s_g` is solvable locally.

ARAP advection must be per face. The CR `g_perp` of a bent edge averages
the two faces' in-plane perpendiculars and shortens by about cos(theta / 2);
the per-edge `R_e G0_e` keeps it unit length, so even the exact answer is not a
fixed point (bias 4e-4). `T0 = sum_f R_f D_f X0` (face Procrustes) equals
`D X` whenever faces move rigidly.

Two solvers, with test targets read from an exactly isometric cylinder (columns
on the circle at chord = flat spacing) and from a sphere-cap patch (radius 2,
exp map, 3.4% strain from flat):

- `tangent_loop` (PD): `T* = argmin |T - T0|_M s.t. S(T) = s_g` by
  Newton–Raphson (3 steps to 1e-14), then `X = L^{-1} D^T M T*`; `inner`
  iterates per frozen system, then rebuild (`L`, or `L+G0` reseeding the ARAP
  rest).
- `goal_gn_loop`: Gauss–Newton on
  `1/2 |S(D X) - s_g|^2_W + w/2 |D X - T0(X)|^2_M`, system
  `(J D)^T W (J D) + w L` factored per outer step and kept for `inner`
  iterates with a fresh gradient, refactored early on ascent.

| Case | Solver | Result |
|---|---|---|
| isometric cylinder, flat start | `tangent_loop`, frozen | err 8.8e-9, E 1e-17 |
| isometric cylinder, inverted start | `tangent_loop`, frozen | err 3.2e-7 |
| isometric cylinder | `tangent_loop`, reseed every 30 | err 1.7e-4 (flat), 4.8e-3 (inverted) |
| isometric cylinder, flat start | `goal_gn_loop`, w 1e-2 | err 4e-13 in 22 iterates |
| isometric cylinder, inverted start | `goal_gn_loop`, w 1e-2 | E 1e-27, err 0.21 (another zero) |
| isometric cylinder, flat start | GN on positions, no ARAP | E 1e-29, err 0.059 (another zero) |
| sphere cap | `tangent_loop`, any schedule | stalls, E 5e-5, err 9e-3 |
| sphere cap | `goal_gn_loop`, w 1e-4 | E 2e-11, err 7e-3, strain 0.027 (target 0.034) |

Findings:

- With per-face advection the formulation is exact on isometric targets and
  needs no rebuild: `D` is intrinsic, so the rest system stays correct.
  Reseeding the ARAP rest mid-solve bakes an unconverged shape into the
  memory and costs accuracy.
- `S(D X) = s_g` does not pin the shape: Gauss–Newton without ARAP finds a
  different zero, and the cap solve meets `s_g` to 2e-11 on a surface 7e-3
  from the source with less strain. The read-out has more solutions than the
  rigid motions; the ARAP memory is what selects among them (the least-strain
  one).
- The PD loop stalls on non-isometric targets: its min-norm tangent
  correction is mostly non-integrable and the Poisson step drops it. The
  goal Gauss–Newton does not stall there.
- `goal_gn_loop` at w 1e-4 from the inverted strip diverges (full GN steps
  without a line search) — needs damping before it goes near hepworth.
  (Since fixed: Armijo backtracking, refactor on failure.)

## 10. Where the ARAP rotation lives

The CR basis of edge `e` is `1 - 2 lambda_k` (`k` opposite `e`), so per
triangle the CR stiffness is 4 x P1 re-indexed to edges: two edges sharing
vertex `k` couple by `2 cot theta_k`. The vector Dirichlet matrix is that
scalar stiffness times the Levi-Civita transport: every off-diagonal 2 x 2
block of `vector_dirichlet` is `2 |cot theta|` times a rotation by `theta`
(checked numerically). An edge's stencil is its diamond (both faces).

Memories in `cr_goal.py` (all rest-invariant and rigid-invariant, gradients
checked against finite differences; Hessians trace-normalized to `D^T M D` so
`w` is comparable):

| Memory | Rotation | Energy / Hessian |
|---|---|---|
| `vertex_tangent` | vertex ARAP, `polar((R_i + R_j) / 2)` to edges | `|D X - R_e G0|^2_M`, `D^T M D` |
| `vertex_arap` | vertex ARAP | vertex ARAP, `L_cot` |
| `face_tangent` | face Procrustes, `sum_f R_f D_f X0` | `|D X - T0|^2_M`, `D^T M D` |
| `edge_arap` | diamond Procrustes (cotan spokes of both faces) | CR-ARAP `sum_e sum_{f ni e} A_f/3 |grad_f X - R_e grad_f X0|^2`, `L_cot`, exact gradient |
| `edge_tangent` | diamond Procrustes | `|D X - R_e G0|^2_M` |
| `dir_tangent` | `polar(sum w_ee' G_e' G0_e'^T)`, vector Dirichlet weights | `|D X - R_e G0|^2_M` |
| `dir_arap` | Procrustes on frame spokes `G_e' - G_e` | ARAP on frame spokes, `D^T L_w D`, exact gradient |

Goal Gauss–Newton, frozen system, `w = 1e-2` (`edge_compare.py`, plot
`edge_compare.png`); shape = mean distance to the target surface after
rigid alignment / size (rigid vertex error / size for the open cases):

| Case | vertex_tangent | face_tangent | edge_arap | edge_tangent | dir_tangent | dir_arap |
|---|---|---|---|---|---|---|
| isometric cylinder | 1e-4 | 0 | 1e-4 | 1e-4 | 1e-4 | 7e-4 |
| sphere cap | 6.0e-3 | 5.9e-3 | 6.1e-3 | 6.0e-3 | 6.0e-3 | 5.4e-3 |
| torus tube 0.7 -> 0.9 | 0.032 | 0.033 | 0.033 | 0.032 | 0.032 | 0.036 (seam) |
| torus R 2 -> 2.5 | 0.030 | 0.030 | 0.029 | 0.030 | 0.030 | 0.030 |
| cube -> sphere | 0.031 | 0.022 | **0.016** (E 0.053) | 0.024 | 0.024 | **0.012** (E 0.088) |
| cube -> ellipsoid | 0.014 | 0.017 | 0.013 | 0.016 | 0.016 | 0.015 |

PD loop on the isometric cylinder (tangent base, frozen): face 3e-9, CR-ARAP
edge 1.4e-4, vector-Dirichlet edge 3.1e-4, vertex 1.05e-3.

Findings:

- Edge rotations beat vertex rotations wherever the placement matters: 8x
  smaller bias in the PD loop, and on cube -> sphere the edge energies
  (`edge_arap`, `dir_arap`) halve the shape error and the goal residual of
  the vertex memory, with visibly fewer corner dents. The gain comes from the
  energy forms (exact gradient, rotation fitted on the same stencil the energy
  uses); the tangent forms (`edge_tangent`, `dir_tangent`) behave like
  `vertex_tangent`.
- On both tori and on cube -> ellipsoid every memory lands on the same
  stationary point: there the limit is the metric (section 9), not where the
  rotation lives. The ellipsoid keeps the same crease under every memory.
- `dir_arap` is a memory of frame variation (bending-like, `D^T L_w D`), so it
  partly resists curvature change; it leaves a seam on the torus tube.
- Remaining edge bias on isometric targets is the `g_perp` shortening of a
  bent edge (`R_e G0_e` keeps it unit length).

## 11. The read-out must follow the metric

Sections 9 and 10 read `S` through `D` frozen at the rest metric. That is
not a frozen-Hessian approximation: it changes the residual. Under stretch the
rest `D` gives the wrong `g_perp`, and the Jacobian never sees `dS/dX` through
the lengths and areas in `D`. The remaining residual on the non-isometric
targets is exactly that part, so restarting from the converged shape did not
help: the torus tube restart measured 2.32 at the end of solve 0 and 0.76 for
the same shape with `D` rebuilt, then stalled at tube 0.77 for 10 restarts.

`S_f` depends only on its 6-vertex stencil (face and the three opposite
vertices), including `D`'s dependence on the metric (each `g_perp` row uses
the two faces of its edge). `face_loops.stencil_s3_live` recomputes `D` from
the stencil (matches `cr_derivative` at any metric, 1e-15), and
`goal_gn_loop(readout="live")` takes the Jacobian by per-face central
differences over the 18 stencil coordinates (gradient checked against finite
differences).

Restart protocol (`live_compare.py`, plot `live_compare.png`): solve the
frozen system to convergence (CR-ARAP memory, `w = 1e-2`), make the result
the new rest (metric and ARAP memory), repeat.

| Case | start | after solve 1 | after solve 5 |
|---|---|---|---|
| isometric cylinder | 0.114 | 2.3e-4 | 2.2e-4 |
| sphere cap (open) | 0.098 | 4.5e-3 | 4.2e-3 |
| torus tube 0.7 -> 0.9 | 0.035 | 4.0e-4 | 0 (E 6e-16) |
| torus R 2 -> 2.5 | 0.049 | 1.6e-3 | 0 (E 8e-14) |
| cube -> sphere | 0.061 | 7.9e-4 | 0 (E 3e-9) |
| cube -> ellipsoid | 0.052 | 9.4e-4 | 0 (E 3e-11) |

(Shape: mean distance to the target surface after rigid alignment / size;
rigid vertex error / size for the open cases.) With the frozen read-out the
same cases stopped at 1.2–3.5%. The closed targets converge to the exact
shape; each restart cuts the error about 10x. The open patches meet `S_g` (E
1e-9) on a shape slightly off the source: an open boundary does not pin the
shape. The earlier collapse under reseeding was the frozen read-out
mismeasuring scale (`S ~ 1/size`), not missing an area anchor.

Caveat: these targets are read-outs of meshes with the same connectivity, so
they are exactly realizable. A prescribed `S_g` field (conformal edits) need
not be, and will land on a compromise.

## 11b. Curvature on the right-hand side only

`goal_gn_loop(lhs="memory")`: `H = w L` (CR-ARAP, constant), the curvature
enters only through the gradient `-J^T W r` (a Sobolev-preconditioned
descent, PD-shaped: the matrix never sees `J`). Frames of 30 iterates, rest
reseeded each frame, shape after frames 1, 2, 3 and 10:

| Case | `H = J^T W J + w L`, w 1e-2 | `H = w L`, w 0.1 | `H = w L`, w 1 |
|---|---|---|---|
| torus tube 0.7 -> 0.9 | 4e-4, 0, 0 ... 0 | 3.5e-3, 1.6e-3, 9e-4 ... 4e-5 | 1.4e-2, 7e-3, 4e-3 ... 1.8e-4 |
| cube -> sphere | 0.29, 0.11, 3e-3 ... 0 | 0.26, 0.20, 0.13 ... 3.6e-3 | 3.5e-2, 2e-2, 1.1e-2 ... 2.4e-3 |

It converges, about 10x slower per frame, and the weight matters (at w 1e-2
the torus stalls at 8e-4; a mass term `mu M` slows it further).

## 12. Analytic Jacobian of the live read-out (`cr_curvature_jacobian.py`)

Per face, symbolic in the 18 stencil coordinates and the three boundary
flags `h_i`: `s_f = (k11, k12, k22)` and `J_f = d s_f / d P` (3 x 18). The
edge orientation sign enters `g_perp` and the read-out once each and cancels,
so the face winding is used directly. On a boundary edge the caller puts the
reflected dummy `P_{3+i} = P_i + P_{i+1} - P_{i+2}` in the empty slot: the
masked neighbor area then stays nonzero and its derivative finite (a
degenerate dummy gives `0 * inf`).

Checks (random perturbations of torus, cube, open strip): value against
`stencil_s3_live` 7e-15; `J` against central differences 9e-10; translation
and rotation in the null space of `J`; `J p = -s` (`S` is homogeneous of
degree -1 in the positions) to 6e-13.

`--emit` writes `scripts/sympy/generated/cr_curvature_stencil_generated.hpp`
(`hepworth::block::cr_generated::cr_curvature_stencil(P, h, s, J)`): 1217
CSE temporaries, about 5.7k operations for `s` and `J` together. Compiled
check: `J` against central differences 4.8e-10 over 200 random stencils (a
third with a boundary edge), `J p + s` 3e-15.

## 13. The memory, compared under the working solver

Section 10 compared memories with the read-out frozen at the rest metric
and without restarts (the setup section 11 showed to be broken; the torus
seam came from it), and its `dir_arap` was not Stein's operator: it kept the
block weights `2 |cot theta|`, dropped the transport `Rot(theta)` and the
orientation signs, and compared extrinsic 3D frame differences, which made
it a bending memory. This section redoes the comparison.

### The intrinsic vector Dirichlet memory (`memory="dir_intrinsic"`)

Per edge, from the current CR frame `G_e = [g_par, g_perp]`:

    t = g_par / |g_par|,  n = (g_par x g_perp) / |.|,  b = n x t
    U_e = [t, b]^T G_e        (2 x 2; b . g_par = 0)

`t, b` is exactly the Stein edge frame (`g_perp` lies along `N x t`), so
`U_e` is in the coordinates `ops.vector_dirichlet` expects. The energy is

    E = 1/2 sum_{j = par, perp} c_j^T (L_dir + mu M2) c_j,   c_j = column j of U(D X) - U0

with `L_dir = ops.vector_dirichlet(m, X_rest)` as written (weights, block
rotations, orientation signs, diagonal) and `mu = 1e-3 tr(L_dir) / tr(M2)`.
`U_e` depends only on `G_e`, so `dc / dT` is one 4 x 6 block per edge; the
Gauss–Newton memory Hessian is `J_c^T K J_c` (`J_c = dc/dT D`), frozen per
restart, gradient fresh every iterate.

Checks (strip, 17 x 6): rest 0; rigid motion 2e-27; exactly isometric
cylinder 1.8e-4 (the `cos(theta / 2)` shortening of `g_perp`); 10% stretch
2.8; 10% uniform scale 3.5; gradient against central differences to 1e-9.
A smooth in-plane stretch costs 1.6 and rough in-plane noise of similar size
10.4: unlike CR-ARAP, which pays every strain alike, this memory charges
rough metric change much more than smooth change.

### Comparison (`memory_restart_compare.py`)

Live read-out, `J^T W J` on the left, 30 iterates per frame, rest reseeded
every frame, 8 frames; each memory at `w` in `{1e-3, 1e-2, 1e-1}` (Hessians
trace-normalized). Shape error (surface distance after rigid alignment /
size; rigid vertex error / size for the open cases), best `w` per memory:

| Case | CR-ARAP | `dir_intrinsic` | vertex tangent | `dir_arap` |
|---|---|---|---|---|
| isometric cylinder (open) | 6.6e-3 (w 0.1) | **3.4e-4** (w 0.1) | 8.6e-3 | 5.5e-3 |
| sphere cap (open) | 5.5e-3 | 4.6e-3 | 4.5e-3 | **3.6e-3** |
| torus tube 0.7 -> 0.9 | 1e-16 | 1e-15 | 1e-16 | 3e-16 |
| torus R 2 -> 2.5 | 9e-17 | 2e-15 | 1e-16 | 1e-16 |
| cube -> sphere | 3.5e-7 | **1.6e-10** | 7.0e-7 | 3.1e-7 |
| cube -> ellipsoid | 4.2e-10 | **4.6e-13** | 8.9e-9 | 9.1e-10 |

Per frame (`memory_restart_compare_curves.png`): on the cube cases
`dir_intrinsic` falls about 10x per frame from the first frame (cube ->
sphere 7.5e-4 at frame 2, against 3.6e-2 to 4.2e-2 for the others, which
plateau for three or four frames first). On the tori every memory reaches
machine precision; `dir_intrinsic` is fastest on the tube at frame 1
(1.2e-4 against 1.3e-3). On the open strip it is 20x better than CR-ARAP at
the same budget, with dihedral error 2.8e-5 against 4.9e-4. Final shapes
show no artifacts (`memory_restart_compare.png`); `w = 1e-3` is best on the
closed targets, `1e-1` on the open strip.

### Decision

`dir_intrinsic` is the memory for the hepworth build, with CR-ARAP kept as
an option. Its Jacobian depends on the shape, so in hepworth it is a
Gauss–Newton constraint per edge (4 rows over the edge's diamond, 12
columns, frozen per frame) rather than a constant projection row; the demo
rebuilds every frame anyway. (The design doc refines this to one constraint
per face, `r_f = K_f^{1/2} c_f`, since `L_dir` assembles from positive
semidefinite per-face elements.)

## 14. Where the curvature tensor lives: edges

A symmetric 2-tensor `(k11, k12, k22)` splits into its trace
`tau = k11 + k22` (spin 0) and trace-free part `z = ((k11 - k22) / 2, k12)`
(spin 2); the smoothness energy is `int |grad z|^2 + 1/2 |grad tau|^2`.

### Elements (`shape_dirichlet.py`)

- **Edge (CR).** Each edge stores its tensor in its frame `(t_e, N x t_e)`
  (unfolded per face; flipping `t_e` flips `b_e` too, so the orientation does
  not matter). Per face, map each edge's tensor into the face frame (`M_a`,
  the spin-2 frame map for the edge's angle) and interpolate with the CR
  basis: `E_f = sum_ab K_ab q(M_a h_a, M_b h_b)`, `K` the CR scalar stiffness
  (4 x the P1 stiffness of the opposite vertices), `q` the inner product
  above. The same construction at spin 1 (`M_a` a rotation, `q` the dot
  product) reproduces `ops.vector_dirichlet` exactly (scale 1.000000,
  difference 3e-16 on a jittered planar mesh and a torus), so the spin-2
  element is Stein's construction applied to tensors.
- **Face (dual graph).** Per interior edge, the two faces' tensors read in
  the unfolded edge frame, `w_e q(d, d)`, with circumcentric
  `w = 2 / (cot alpha + cot beta)` or barycentric `w = |e| / |c_f - c_g|`.

All three: zero on a constant tensor field, invariant under rotation; the
edge element's trace part equals the scalar CR Laplacian.

### Correctness (`shape_dirichlet_compare.py`)

Planar field with an exact reference (`z = a(y) e^{2 i pi x}`,
`a = 1 + 0.5 sin(2 pi y)`, `tau = 1 + 0.3 x`), five families, five levels:

| | Edge | Face circumcentric | Face barycentric |
|---|---|---|---|
| order (right grid, lattice, jittered, 4:1, random Delaunay) | 2.00, 2.01, 1.99, 1.84, 2.00 | NaN, 2.11, 1.02, -31.6, NaN | -0.04, -0.17, 1.03, -0.10, -0.24 |
| finest relative error | -1.9e-3 to -1.2e-4 | NaN (cocircular quads), up to 1.5e10 | -7.9e-3 to +1.86 |
| cylinder (parallel `S`, exact 0) | 3.7e-4, 4.6e-5, 5.7e-6 | NaN | 0 (trivially, piecewise constant) |
| umbilic sphere (exact 0) | 7e-14 | 9e-14 | 2e-14 |
| torus, regular / flipped diagonals | 29.43, 29.87, 29.98 (both) | NaN | 33.98 / 32.36 at 96 x 48 |

The circumcentric face weight is singular on cocircular quads (every
diagonal of a right-triangle grid) and huge near them; the barycentric weight
is finite but inconsistent (no convergence, 186% error at 4:1 aspect, 13%
high and triangulation-dependent on the torus). This is the face-butterfly
distortion seen in fairing. The edge element converges at second order on
every family and passes the curved checks. (The torus regular and flipped
triangulations agree exactly for the edge element by symmetry: the grid's
quads are isosceles trapezoids; the planar jittered and random families are
the sensitivity test.)

### End to end (`shape_site_e2e.py`)

Clean target tensor read at the site, 10% noise, fairing
`(M + 2 lambda L) h = M h_noisy` (`lambda` in units of the mean edge length
squared), recovery with the live read-out, `dir_intrinsic` memory and
restarts (8 frames). Coupling for edges: the mean of an edge's two faces'
readings against `h_e` ("edge_avg", exactly idle at the target); the
per-(face, edge) variant ("edge_pair") has a floor at the exact target (the
two faces' readings differ by 8.6% rms on the torus).

| Case | Site | clean | noisy | best faired (lambda) | target error, best |
|---|---|---|---|---|---|
| torus tube | face (barycentric) | 0 | 0.0023 | 0.0042 (0.03) | 0.073 |
| torus tube | edge_avg | 0 | **0.0016** | **0.0018** (0.03) | **0.051** |
| torus tube | edge_pair | 0.0029 | 0.0025 | 0.0270 (0.3) | 0.087 |
| cube -> sphere | face | 0 | **0.0084** | **0.0021** (1.0) | 0.183 |
| cube -> sphere | edge_avg | 0 | 0.0166 | 0.0093 (0.3) | **0.071** |
| cube -> sphere | edge_pair | 0.0108 | 0.0107 | 0.0042 (0.3) | 0.081 (smoothness 28x clean) |

(Shape: mean surface distance / size.) On the torus the edge site wins on
every measure, and fairing with the face energy makes the shape worse. On
cube -> sphere the face site's lower shape error comes from over-smoothing
(target error 0.18, smoothness 0.09 of the clean field): the field is pushed
toward a constant, right for a sphere but not for varying curvature. The edge
energy is the better denoiser of the tensor itself on both cases. Open: the
edge site's unfaired cube -> sphere run is worse than the face site's (0.0166
against 0.0084, a low-frequency lopsided error), to revisit with fairing
weights in part 2.

### Decision

The curvature tensor lives on **edges** (`shape_dof` in `shell_cr_block`,
in each edge's unfolded frame), coupled to positions by the **mean of the
edge's two faces' live readings** (3 rows per edge over the union of the two
face stencils, exactly zero at a realizable target), and faired by the
**CR spin-2 element** (`edge_dirichlet(m, X, 2)`: CR stiffness x spin-2 frame
maps), which converges at second order on every triangulation tested. The
face site is ruled out by its smoothness energy: neither dual weight is
correct on general triangulations.

## 15. C++ (part 2)

Sections 11, 13 and 14 are implemented in `include/gaudi/kusama/weingarten.hpp`
(builders, `weingarten_solver`) and `include/gaudi/hepworth/block/weingarten_gn.hpp`,
`edge_shape_fair.hpp`, `edge_shape_goals.hpp`; the map and the findings are in
`docs/cr_curvature_hepworth_plan.md` section 0. Parity with this study
(`include/gaudi/test/weingarten_tests.hpp`):

| Check | Python | C++ |
|---|---|---|
| spin-1 element vs `vector_dirichlet` | equal | 6e-16 |
| spin-2 planar field, right grid, rel. error n 11 / 21 | -8.19e-3 / -2.05e-3 | -8.190e-3 / -2.053e-3 |
| torus tube 0.7 -> 0.9, restarts | exact | rigid error 7e-14 |
| cube -> sphere (clean), 8 frames of 30 | 0 (E 3e-9) | max radial error / size 1.2e-7 |
| same through `projection_solver`, 20 frames | n/a | torus 2e-12, cube 6e-4 (inertia anchors each frame) |
