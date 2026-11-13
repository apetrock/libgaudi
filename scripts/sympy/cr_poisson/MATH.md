# X = L^{-1} D^T M F(J^{-1}(S_g)): shape-operator targets through the CR derivative

Source: Stein, Wardetzky, Jacobson, Grinspun 2020, *A Simple Discretization of
the Vector Dirichlet Energy*, eq. (9) and App. A.2
(`reference/stein2020-vector-dirichlet-energy.pdf`). Code it mirrors:
`hepworth::block::cr_coordinate_gradients` and `cr_face_operators`
(`include/gaudi/hepworth/block/cr_shape.hpp`), ported in `cr_ops.py`.

**Current formulation.** This file is a lab notebook; section numbers are
cited from code, so sections keep their numbers. The live pieces: the CR
read-out and its Jacobian (sections 1, 11, 12), the vector Dirichlet memory
(13), the edge store and its spin-2 fairing (14), the `weingarten` coupling in
C++ (15), the exponent solver (16), and the constrained position + shape
formulation, Engman (17), shape as a Lie sphere coordinate (18), and a
Poisson solver in Dupin space (19).

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

## 16. The exponent solver: curvature drives position (`exponent_solver.py`)

Lineage: Stein et al. 2020 supply the operators (CR derivative `D`, mass
`M`, the Dirichlet operator `D^T M D`, the face read-out of App. A.2); Crane
et al. 2013 (`reference/crane2013-conformal-willmore-flow.pdf`) supply the
loop shape (edit curvature, turn it into new tangents by per-element
similarities `T~ = lambda-bar T lambda`, recover positions by one Poisson /
Dirichlet solve); Wang, Liu, Tong 2012 supply the reconstruction from
discrete fundamental forms (lengths, dihedrals) and its integrability
condition. Motivation: the Gauss-Newton coupling (sections 11 to 15) reaches
any realizable target to machine precision, but shape-space curvature edits
did not buckle the mesh (TODO 2): its memory (`dir_intrinsic`) picks the
smoothest least-strain shape and its convex model gives positions no outlet.

### Rotation form with fresh targets (section 8 revisited)

Sync of the per-face rotations: Jacobi polar sweeps converge far too slowly
(strip flat -> cylinder 1.5, 200 iterates: rigid error 0.11 / 0.06 / 0.015 at
1 / 2 / 4 sweeps); the exact factorized sync reproduces 2.561e-4; dihedral
targets read once per rest (they depend on positions only through the
metric, exact under isometry) and the sync factored once give 2.5 to 3.6e-4
at half the cost. With fresh targets every iterate (live `K3`), frozen
factors per frame and one defect correction of the sync per iterate, the
isometric cylinder reaches 8e-12 by frame 2 (Gauss-Newton best 3.4e-4).
Rotations never change the metric: torus tube / R stall at 3.5e-2 / 5e-2.

Buckling under frustrated curvature (flat 13 x 13 sheet, target a sphere cap
read-out, frame-0 rest kept): the rotation form snaps toward the bimetal
cylinder (shape ratio q_min / q_max of the fitted height quadric 0.16 at 3.5%
strain, sync anchor `mu = 3`); Gauss-Newton with a CR-ARAP memory dishes
(soft, ratio 0.99, 18% strain) or freezes (stiff, curvature 0.07 against 1).

### Stretch: what failed and why

Every curvature mismatch can be met by bending (dihedral) or stretching
(length), `k ~ theta / l`; which one is right is global (an open strip
bends, a closed tube must stretch, a plate asked for a sphere must stretch or
buckle). Variants tried in the prototype:

| Variant | Result | Cause |
|---|---|---|
| per-face split with a stretch prior | faces shrink instead of bending from flat | local choice cannot see topology |
| stretch what bending left, per frame (gain `beta`) | closed band 1.298 of 1.3, tube 9e-4, oscillates | a fixed-point guess, not a Newton step |
| per-face Newton / naive global block on stretch | diverge / worse than no stretch | double counting: bend and stretch both answer the full residual |
| rest-length relaxation between frames | tori 1e-3, cubes diverge | fixed-dihedral step against a reconstruction that undoes it |
| joint row in (theta, log l) | dihedral part not realized on closed surfaces | the step is not integrable |
| joint row + Wang's vertex closure | tube 1.3e-3, band and torus R stall | rotations rebuilt through the mu-tied sync |
| body-frame preconditioner `P` (sync + gather) | `<J P r, r> / |r|^2` +0.34 (cylinder), -0.02 (tube) | `(omega, L)` from dihedrals alone is not a displacement gradient |

### The factorization (checked)

With each face's first-order deformation split into rotation `omega_f` and
in-plane log-stretch `L_f` (torus, random perturbation 1e-6):
`dtheta_e = a_e . (omega_- - omega_+)` to 3.9e-6 relative,
`dlog l_e = t^T L t` to 1.8e-5, and `ds = J_theta C omega + J_ell E Lambda`
to 1.6e-5, so `J_X = [J_theta C, J_ell E] [Omega; Lambda]`. In the body frame
the fundamental-form Jacobian is rigid-invariant (3e-10 against 73% for the
world-frame `J_X`) and changes 3x less under isometric bending (4 / 9 / 17%
against 15 / 29 / 55% at r = 6 / 3 / 1.5). The Gauss-Newton Hessian is a
metric-only core conjugated by the per-face rotations,
`J_X^T W J_X = Omega_0^T R (C^T B^T K^T W K B C) R^T Omega_0`; the rotations
around the core are why Gauss-Newton refactors.

### Read-out in the fundamental forms (`cr_ff_jacobian.py`)

Each face's stencil is rebuilt from 9 lengths (its edges, its wings' edges)
and 3 dihedrals: the face in z = 0 by the law of cosines, each wing opposite
the third vertex turned about its edge by `-(theta + pi)` (theta in the sync
convention of `face_loops.dihedrals`, the same from both faces). The read-out
is rigid-invariant, so `J_ff = J_P dP/d(l, theta) diag(l, 1)` with `J_P` the
generated stencil Jacobian (section 12). Rebuilt read-out against live 1e-14
(torus, cube, strip with boundary dummies); `J_ff` against central
differences 1.4e-8 / 3e-10 / 3e-10. `--emit` writes
`scripts/sympy/generated/cr_ff_rebuild_generated.hpp` (187 CSE temporaries,
about 730 operations).

### The exponent formulation

Unknowns per face: `xi_f = [omega_f]x + P_f L_f P_f^T`, the face's
displacement gradient (rotation 3, in-plane log-stretch 3). One constant
system per restart at the current `X`:

- curvature rows `sqrt(W) ([J_theta C, J_ell E] xi - r)`, `r = s_g - s(X)`;
- compatibility `xi_f d_e = xi_g d_e` on every interior edge (stiff, `1e8`):
  the field is the gradient of a displacement, which also gives Wang's closure
  and edge-stretch agreement;
- priors `kb |omega|^2` (tiny; removes the rigid null space) and, for an
  elastic material, `km |L_cur + dL|^2` with `L_cur` each face's current
  log-stretch against a persistent rest.

Each iterate: fresh residual, `xi` by one backsolve, `F_f = expm(a xi_f)` per
face (rotations leave the solve through the exponential), advected frames
`sum_f F_f g_f(X) / M_e` gathered into positions by `D^T M D` (constant per
restart), line search on `|r|_W^2 + km tau sum A_f |L_cur|^2`. To first order
this is Gauss-Newton on positions with restarts; the exponential carries large
rotations exactly. All matrices are symmetric positive definite and constant
per restart (LDL^T in C++; LU in the prototype).

Parity (restart every 30 iterates, `w_c = 1e8`, `km = 0`):

| Case | Iterates | Residual | Shape | Gauss-Newton best (section 13) |
|---|---|---|---|---|
| isometric cylinder | 5 | 8.6e-1 -> 2.7e-9 | 2.0e-6 | 3.4e-4 |
| sphere cap (open) | 90 | 9.2e-1 -> 1.3e-7 | 2.9e-4 | 3.6e-3 |
| torus tube 0.7 -> 0.9 | 4 | 2.6 -> 9.8e-9 | 3.4e-10 | 1e-15 |
| torus R 2 -> 2.5 | 4 | 1.3 -> 1.4e-8 | 6.3e-10 | 1e-16 |
| cube -> sphere | 30 | 7.4 -> 4.9e-9 | 4.5e-10 | 1.6e-10 |
| cube -> ellipsoid | 10 | 8.1 -> 5.1e-9 | 2.2e-10 | 4.6e-13 |
| closed band r 1 -> 1.3 | 6 | 5.8e-1 -> 6.6e-9 | 1.5e-9 | n/a |

The residual floor near 1e-9 is the same with finite-difference and symbolic
`J_ff`, so it is not the Jacobian. The C++ port reaches 1e-14 on the same
cases, so the floor is the prototype's numerics.

Buckling (plate as above, elastic against the flat rest; ratio 1 dish, 0
cylinder):

| Target curvature | km | Ratio | Principal curvatures | Max strain |
|---|---|---|---|---|
| 0.667 | 1 | 0.996 | 0.663 / 0.661 | 6.8e-2 |
| 0.667 | 1e2 | 0.085 | 0.800 / 0.068 | 7.2e-3 |
| 0.667 | 1e4 | 0.007 | 0.868 / 0.006 | 1.6e-4 |
| 1.0 | 1 | 0.991 | 1.005 / 0.996 | 1.2e-1 |
| 1.0 | 1e2 | 0.067 | 1.464 / 0.098 | 8.3e-3 |
| 1.0 | 1e4 | 0.008 | 1.560 / 0.013 | 2.0e-4 |

`km` is the material: soft stretches into the dish, stiff holds the metric
and the frustrated isotropic curvature snaps to the cylinder, the only
developable shape that meets one principal curvature. The cylinder's
curvature overshoots the target (0.87 for 0.667, 1.56 for 1.0); not yet
understood.

### Eigensatz on the same system

Principal-curvature rows (Eigensatz 2010 eq. 3.12) are the eigenvalues of each
face's `s_f`, `kappa = m +- d`, `m = (k11 + k22) / 2`,
`d = sqrt(((k11 - k22) / 2)^2 + k12^2)`; their rows are
`dkappa/ds . [J_theta C, J_ell E]` with

    dkappa_1,2 / ds = (1/2, 0, 1/2) +- ((k11 - k22) / (4 d), k12 / d, -(k11 - k22) / (4 d)),

and within `umbilic (|m| + d)` of an umbilic both rows are `dH`, the paper's
A.6. Its metric terms (eq. 3.11) are rows on the stretch: conformal on the
trace-free part `((l11 - l22) / 2, l12)`, areal on the trace `l11 + l22`. With
`km` they form one 3 x 3 stiffness `M` on `(l11, l12, l22)`, in `H` as
`tau A_f / mean A M` and on the right-hand side as `-tau A_f / mean A M L_cur`,
the merit carrying the same `L_cur^T M L_cur`. Everything else is shared.

| Check | Result |
|---|---|
| `dkappa/ds` against central differences | 4.4e-10; umbilic rows equal `dH` exactly |
| tube 0.7 -> 0.9 (read-out eigenvalues of the target mesh) | 14 iterates, 2.3 -> 2.9e-14, shape 1.3e-14 |
| metric rows, same targets, elastic | trace-free stretch 0.138 -> 0.0021 at `w_conformal` 10; trace 0.228 -> 0.016 at `w_area` 10 |
| torus x1.2 against the kusama solver (analytic targets) | kusama R 1.189, tube 0.396; exponent R 1.205, tube 0.393 |

The analytic targets are not exactly realizable by the discrete read-out, so
the exponent merit stops at 0.41 there.

### C++ port (`hepworth/block/exponent_coupling.hpp`)

- `E` reads each length from the stencil face that owns it (the face for its
  own edges, the wing face for the wing edges), not the mean of the edge's two
  faces: the compatibility rows make them agree, and `H` halves (118 to 60
  nonzeros per row; 1.9 to 0.64 s to factor on 2304 faces). Same results.
- The gather and the projection rows are in increment form,
  `x + (D^T M D)^-1 D^T M (G - G(x))` and `D x + G_step - G(x)`, with `D` at
  the step start. The absolute form reproduces `x` only where `D` was built;
  once the iterate moved, every line search was offset and stalled.
- `factor` (curvature rows, compatibility, priors, `H` factored) is reused
  across frames, its steps turned from each face's build frame `[T B N]_0` to
  the live one by `R_f = [T B N] [T B N]_0^T`. On a rigidly turned torus the
  reused step matches a fresh one to 9e-7 of the step, which is the floor
  between two fresh factors: the rigid part is held only by the tiny `kb`. A
  step that stalls on a reused factor refactors at once. If that also stalls
  (the large-residual balance of an elastic target, where Gauss-Newton has no
  descent direction left), the cache is marked futile until its age refresh.
  The standalone solver restarts when an iterate cuts the merit by less than
  10x (a frozen factor contracts linearly, a fresh one quadratically).
- The stretch prior is area-weighted as the merit; unweighted, a fresh step
  was not a descent direction on the torus (faces differ 2x in area).
- Weights: the gather rows are `w sqrt(k_gather M_e) D`, `k_gather = 34.5`
  calibrated against the cotan `bend` row like `k_coupling` (unscaled the
  rows are resolution-independent at 0.029 of `bend`); the store is anchored
  at its frame-start reading at `w_coupling`. One weight then means the same
  under `weingarten`, `exponent` and `eigensatz`.
- A reused factor keeps its rows and `H` together. Rebuilding only the rows
  from the live mesh each frame (fresh gradient, frozen `H`) is unstable:
  `H` holds the rigid turn only by `kb = 1e-6`, rows and `H` from the same
  build leave the right-hand side no component there, and live rows leak one
  that `H^-1` amplifies about `1e6` (every face hits `max_step`; roughness
  0.19 to 4.5 in 10 frames). The demo default is a fresh factor every frame,
  so all tangent geometry is reseeded; reuse (10 frames at 5% pressure: tube
  0.1991 against 0.1992 fresh) is opt-in.
- Demo frames (48 x 24 torus, 45 iterates): 0.8 s refactoring every frame,
  0.11 s at `refresh_every` 8 or 16 with the same shapes to 3 or 4 digits; idle
  0.08 s and exact to 1e-15. A mass on the exponents is stable but only slows
  the response.

## 17. Constrained position + shape formulation (Engman)

Plan: `docs/kkt_admm_curvature_solver_plan.md` (revised) and the Engman plan.
Code: `include/gaudi/engman/`. Every row in `projection_solver` is a penalty
(projection rows `w (A q - p)`, Gauss-Newton rows `w r(q)`), so a solve lands
on a weighted compromise: the store sits at the blend
`S = (w^2 R + sum w_k^2 g_k) / (w^2 + sum w_k^2)`, `k_coupling` has to be
calibrated against `bend`, and a stiff coupling trades accuracy against
conditioning. Engman keeps the same variables and rows but solves the
coupling as a constraint.

### 17.1 Variables and roles

`q = (x, S)`: positions `x` (`3V`) and the edge store `S` (`3E`, section 14).
Every term declares a role.

- **Energies** `E_k(q) = 1/2 w_k^2 |r_k(q)|^2`, a compromise by intent:
  inertia `1/2 |x - s|^2_{M_h}`, the intrinsic vector Dirichlet memory
  (section 13, Gauss-Newton; rotation invariant, so it charges metric change
  and not bending), the store rows (`edge_shape_fair` pins and goals,
  fairing), and optional position rows (stretch, area, ...). No cotan `bend`
  row: under a hard coupling the store and its fairing carry all curvature.
- **Constraints** `C_j(q) = 0` with multipliers. First the coupling,

      C_e(x, S) = c_e (R_e(x) - S_e),    c_e = k_coupling A_e,

  scaled like the `weingarten_gn` row and the store pins (both `c_e` times a
  curvature), so `C` is an integrated curvature;
  later isometry (17.7) and global scalars.

The problem is `min_q E(q)` subject to `C(q) = 0`, with Lagrangian
`L = E + lambda^T C`.

The roles matter. An augmented Lagrangian or ADMM on a soft term converges to
the same weighted minimizer as projective dynamics (Overby et al. 2017), only
faster: the weight still decides where it lands. Only a constraint becomes
exact, and its `rho` is then a conditioning parameter with no accuracy
meaning. Making every row a constraint is infeasible (stretch, area and Monge
landings conflict).

### 17.2 Why the store stays a variable

The coupling is explicit in `S` (`S = R(x)` at feasibility), so `S` and
`lambda` could be eliminated: `min_x E_x(x) + E_S(R(x))`, with Gauss-Newton
Hessian `H_xx + J_x^T H_SS J_x`. `J_x` reaches each edge's two face stencils
(a 2-ring) and the fairing in `H_SS` a ring of edges, so the eliminated
system fills to a 3-ring and is fourth order. Keeping `S` keeps
`H = diag(H_xx, H_SS)` with the coupling Jacobian
`J = c [J_x, -I]` sparse, and the position block equal to today's
factor. At feasibility `S - R(x) = 0`: the store is a splitting variable, not
state. The persistent state is `lambda` (17.5).

### 17.3 The Newton-KKT step

At `(q, lambda)` with `H` the Gauss-Newton Hessian of the energies, `g` their
gradient and `J = dC/dq`, one step of the proximal method of multipliers is

    [ H    J^T     ] [ dq      ]     [ g + J^T lambda ]
    [ J   -I / rho ] [ dlambda ] = - [ C              ]

and `lambda <- lambda + dlambda`. The second row reads
`C + J dq = dlambda / rho`: the linearized constraint, relaxed by `1 / rho`.

**Orthogonality.** Split `dq = dq_n + dq_t` with `dq_n` in `range(J^T)` and
`dq_t` in `null(J)` (columns `Z`). For `rho -> inf`, `J dq_n = -C` restores
feasibility and

    Z^T H Z dq_t = -Z^T (g + J^T lambda + H dq_n):

the energies move the solution only along the constraint manifold. They are
still least squares (Gauss-Newton), but orthogonal to the constraints: an
energy cannot buy a constraint violation, and several constraints stacked in
`J` are met together instead of each pulling with its own weight on the
whole space, as penalties do. The coupling has full row rank always (the
`-c I` in the store columns), so its normal space is never degenerate;
`-I / rho` regularizes any other rank deficiency.

**Factorization.** With `H` positive definite the matrix is symmetric
quasi-definite (Vanderbei 1995): an `LDL^T` exists for every symmetric
permutation, with `D` having `dim q` positive and `dim C` negative entries,
so `SimplicialLDLT` factors it with its fill-reducing ordering and no
pivoting. `H` positive definite needs every variable held by some energy:
inertia holds `x`; the store needs a floor where no row touches an edge, a
proximal `mu M_S |S - S_prev|^2` (the store damping, 17.6).

**Iterates.** As in `projection_solver`: the factor is kept for several
iterates with a fresh right-hand side (`g`, `J^T lambda`, `C`), refactored on
a schedule and on ascent. Each step takes a line search on the augmented
Lagrangian merit `phi = E + lambda^T C + rho / 2 |C|^2`, with the same
backtrack, fresh-ascent stop and refactor cap. The step ends with a
feasibility restoration: KKT solves with no energy gradient (rhs `[0; -C]`,
the smallest correction in the `H` metric), Newton on `C = 0`, so every step
ends on the constraints even when its iterates ran out descending (cube to
sphere: frames 0 and 1 end at `|C|` 1.4e-2 and 5.7e-3 after 30 iterates,
1e-14 after restoration). Usable `rho_rel` is up to about 1e6; at 1e9 the
constraint block of the factor is lost and steps crawl.

Measured (`test/engman_tests.hpp`, torus 24 x 12 and cube 6): with the duals
off and `rho = w^2` the solve lands on `projection_solver`'s `weingarten_gn`
result (positions 1.6e-10, store 8.6e-12); with them on the coupling is met to
1e-14 and the result does not depend on `rho_rel` in 1e2 to 1e6 (spread
1e-15); a ring of vertices pinned exactly beside the coupling leaves both
exact (6e-18, 5e-15) where as penalties they trade (0.11, 0.074); over 20
frames the torus tube reaches 3.1e-15 (`projection_solver` 5.8e-13) and cube
to sphere 3.2e-13 (8.6e-12), every frame feasible to 1e-14.

**Cost, measured** (`engman_kkt_cost`; memory, store pins and the coupling
on a torus). One factor of the full KKT matrix has the fill of
`projection_solver`'s matrix with the `weingarten_gn` rows and is a little
faster: 48 x 24, `L` 1.51M either way, 0.13 s; 128 x 64 (8192 vertices), `L`
18.2M against 19.0M, 3.9 against 4.4 s. The direct Schur reduction onto
positions (`kkt_mode::schur`: per edge `D_e = c_e^2 H_SS,e^-1 + I / rho`,
then `H_xx + J_x^T c D^-1 c J_x`; exact, matches the full solve to 8e-16,
and applies only while every store row is per-edge) is 7x smaller (24576
unknowns) but fills more (`L` 21.4M, 5.9 s): eliminating the store couples
each vertex to its 2-ring, the fill 17.2 keeps the split to avoid. The full
factor is the default.

**Scale.** For large meshes, eliminate `dq` instead: with
`J = c [J_x, -I]`,

    (c (J_x H_xx^-1 J_x^T + H_SS^-1) c + I / rho) dlambda = rhs,

SPD on edge rows, solved by PCG with the position factor for `H_xx^-1` and
the preconditioner `M_S + beta L_S` (17.6). Not needed until the direct
factor is the bottleneck.

### 17.4 Block Gauss-Seidel (fallback)

Position step on `H_xx + rho J_x^T c^2 J_x` (cached), then the store step
(per edge `3 x 3`, or the fairing system), then `lambda += rho C`. It is
ADMM on the split and reaches the same point; kept as the diagnostic for the
KKT step and its fallback.

### 17.5 The multiplier as persistent pressure

At a solution `grad E + J^T lambda = 0`. Its store block gives
`lambda_e = grad_{S_e} E / c_e`: the force the shape energies exert
through the coupling, the curvature pressure the mesh still owes the design.
For a pin `1/2 w^2 c_e^2 |S_e - g_e|^2`, `lambda_e = w^2 c_e (S_e - g_e)`,
which scales like the area. Transfer the intensive part,
`mu_e = lambda_e / c_e` (a curvature): at the end of a frame scatter
`mu` to vertex datums (each vertex's chart frame through `vertex_stars`,
weight `A_e / 3`), and after the remesh gather `lambda_e = c_e mu_e`
with `mu_e` the mean of its ends rotated into the edge frame. Vertex datums
interpolate on split and merge on collapse; edge datums do not
(`datums.hpp`). A split then neither halves a multiplier nor a collapse
doubles it.

### 17.6 The store's metric

The smoothness energy splits by spin (section 14),

    E_fair = lambda_H int |grad tau|^2 + lambda_z int |grad z|^2,

`tau = k11 + k22`, `z = ((k11 - k22) / 2, k12)`. `edge_shape_fair` has the
full spin-2 element (`dirichlet`) and the trace part (`dirichlet_trace`);
the traceless part on its own (`dirichlet_traceless`) gives the two weights.
`lambda_z` smooths the principal direction field weighted by anisotropy, so
directions turn freely near umbilics; the Dupin rows are its directional
restriction. Use `M_S + lambda_H L_tau + lambda_z L_z` as the metric more
than as an energy: the store's proximal term in `H_SS` (17.3), and the dual
preconditioner (17.3, Scale), since `lambda` is an edge curvature density
and this is its natural inner product. Raising the energy weight
over-smooths (section 14, cube to sphere).

Built (`test/engman_tests.hpp`). The split is exact: with `t = (1, 0, 1)`,
`tensor_inner = q_z + 1/2 t t^T` and every frame map keeps the trace
(`t^T M = t^T`), so `E_f = K_ab M_a^T q_z M_b + 1/2 K_ab t t^T`
(`spin2_traceless_face_element`; element 5e-16, row energies 2e-16 on a noisy
torus field). The traceless rows (`edge_fair_rows::dirichlet_traceless`,
demo `shape_distortion::w_grad_traceless`) vanish on a constant planar field
(4e-16) and on any umbilic field (1e-14 against 18 for the trace rows on a
varying mean curvature). The metric is `engman::store_metric`: the two row
kinds toward the live store (`goal_live`), zero residual and gradient, so
they add `lambda_H L_tau + lambda_z L_z` to `H_SS` and leave the solution
(2.8e-15). Under the direct KKT step its smoothing is small (first-step store
roughness 0.471 to 0.462 at `lambda` 1e2): the coupling ties `S` to `R(x)`,
so the metric only reshapes steps the positions take too. Its larger role is
the dual preconditioner of the Schur path.

### 17.7 The buckling outlet

With a free metric, more curvature is cheapest as uniform shrink (the pressure
runs, section 16 and the TODO). A per-edge length store `l_e` and the
isometry constraint `|x_j - x_i|^2 - l_e^2 = 0` (or Soliman's cross ratio for
a conformal class) close that outlet: an edit that changes `det S` at fixed
metric is unrealizable (Gauss) and lands on a compromise; one that changes
anisotropy at fixed `det` has to bend, which is where modes appear. A
multiplier does not remove a saddle point: symmetry-breaking noise is still
needed.

### 17.8 What it does not change

The memory stays the rotation-invariant Gauss-Newton energy: a fixed
quadratic metric `M + alpha L_vec` in `x` would penalize rotation, which is
the bending curvature edits ask for. Its runaway on large remeshed meshes was
step control on a nonlinear row (bounded by `gn_backtracks`, the fresh-ascent
stop and the refactor cap). The `projection` path keeps every row, so the two
solvers compare like for like until Engman is the default; a later audit
(the Engman plan, phase 5) keeps only the rows that do what they say under a
hard coupling.

## 18. Shape as a coordinate: Lie sphere geometry and the Dupin cyclide

Code: `include/gaudi/albers/lie_sphere.hpp`, `include/gaudi/albers/lie_dupin.hpp`,
tests `include/gaudi/test/lie_sphere_tests.hpp` (`tests/lie_sphere_tests`).
Reference: T. Cecil, *Lie Sphere Geometry* (2008).

The ring rows (`osculating_torus_ring`, `osculating_dupin_ring`) fit a model
shape per vertex and land the one-ring on it. The bare torus row works; the
volume modes and the cyclide row do not hold up, and the corotational volume
mode helps only a little. This section records why we think the one-ring fit
is missing structure, and the coordinate that supplies it: a torus or
cyclide as a single element of a symmetric space, with exp/log, closed-form
fits and transport, the way a unit quaternion is a rotation.

### 18.1 Lie sphere coordinates

R^{4,2} with `<x, y> = -x0 y0 + x1 y1 + x2 y2 + x3 y3 + x4 y4 - x5 y5`.

    sphere (c, r):  [ (1 + |c|^2 - r^2)/2, (1 - |c|^2 + r^2)/2, c, r ]
    point x:        sphere (x, 0)
    plane n.x = h:  [ h, -h, n, 1 ]

All are null, coordinates are projective, and oriented contact is
`<s, t> = 0` (`<s, t> = ((r1 - r2)^2 - |c1 - c2|^2) / 2`). Orientation: a
sphere touches the plane with normal n at x when its centre is x + r n, so the
sphere of normal curvature k at (p, n) has centre p + n/k, radius 1/k (k > 0
bends towards n, the monge sag's sign).

A contact element (p, n) is the null 2-plane L spanned by its point sphere P
and tangent plane T, and the sphere of curvature k is

    k P + T,

with no division: k = 0 is the tangent plane, k = inf the point. **The shape
operator at a point is two marked points on this null line (the two curvature
spheres) plus the frame.** Flat points and umbilics are not singular here; the
`1 / k_hi` spine points and their `apex_cap` clamp are an artifact of
Euclidean coordinates.

### 18.2 The group: transports g

Lie sphere transformations are `A^T eta A = eta` (O(4,2), 15-dimensional),
acting on sphere vectors by multiplication. Built: rotation, translation,
rigid, uniform scale, Laguerre normal offset `(p, n) -> (p + d n, n)`,
orientation flip, reflection in a non-null vector (mirror planes, Möbius
inversion in a sphere), the so(4,2) generators, Cayley transform (exact group
element, cheap), exp and log. `A^-1 = eta A^T eta`.

The edge transport `g_ij` between two vertex shape coordinates is an element
of this group; there is no separate "g coordinate". Volume is not invariant
under it (only the rigid subgroup preserves volume), so a volume hold is a
constraint on the update (18.7), not a property of the coordinate.

### 18.3 The cyclide sigma

A Dupin cyclide (torus, cylinder and cone of revolution, and their Möbius and
Laguerre images) is a split `R^{4,2} = U + U^perp`, both of signature (2,1):
one curvature-sphere family is the conic of null vectors in U, the other in
U^perp (Cecil section 5.4). Its coordinate is the involution

    sigma = +1 on U, -1 on U^perp,   sigma^2 = I,  eta sigma symmetric,  tr sigma = 0.

The space of these is the symmetric space `O(4,2) / (O(2,1) x O(2,1))`,
dimension 15 - 3 - 3 = 9: the classical 9-parameter family. `-sigma` is the
same surface with the families swapped. All Dupin cyclides are Lie equivalent,
so one group parametrizes all of them, the torus included.

Canonical torus (axis z, spine radius R, tube r, outward normals): U is
spanned by `[(1 + R^2 - r^2)/2, (1 - R^2 + r^2)/2, 0, 0, 0, -r]`, `e2`, `e3`
(the tube spheres, centres on the spine circle, radius -r), and
`sigma = 2 Pi_U - I` with `Pi_U = B (B^T eta B)^-1 B^T eta`.

- **Membership is linear in sigma.** L lies on the cyclide iff `sigma L = L`,
  i.e. `<sigma a, w> = 0` for a in L and w in `L^perp` (4-dimensional,
  contains L): eight rows per contact element.
- **Curvatures.** The members of L in U and in U^perp: solve
  `(I - sigma)(alpha P + beta T) = 0` (and with `I + sigma`) projectively,
  `k = alpha / beta`. On the torus this gives `-1/r` and
  `-cos(phi) / (R + r cos(phi))` to 1e-10, including k_v = 0 on top of the
  tube.
- **Advection is conjugation**, `sigma <- A sigma A^-1`: a rigid motion moves
  the torus exactly, a Laguerre offset by d is the torus with tube r + d, a
  scale scales both radii, an inversion gives a non-torus cyclide (tube
  curvature varies around it). Conjugation never leaves the manifold.
- **Normalization is the matrix sign function**: eta-symmetrize, then Newton
  `X <- (X + X^-1)/2` (stays eta-self-adjoint, quadratic). The analogue of
  normalizing a quaternion, or of the SVD in Procrustes.

### 18.4 Closed-form fit

With `sigma = eta M`, M symmetric (21 unknowns), each membership row is
`a^T M w`: a homogeneous linear system. The identity is always in its null
space (it fixes every contact element); a trace-0 row removes it, and the
smallest eigenvector of the normal equations, normalized, is the fit.
Measured null space of the rows on a torus, by number of generic contact
elements: 1 -> 14, 2 -> 7, 3 -> 3, 4 -> 2 (= span{I, sigma}). **Four contact
elements determine the cyclide linearly; three leave a pencil.** Four is
poorly conditioned (eigenvalue gap ~1e-7 on random torus points); six or more
is comfortable, and a one-ring patch of seven points recovers the torus to
1e-4. Exact recovery of a torus and a Möbius-image cyclide from 12 points.
Generic contact elements lie on no common cyclide: three or four random ones
give a null space of 1 (the identity only). "Four determine it" holds for
contact elements already on a cyclide. On an edge's diamond (two triangles,
four contact elements, the quad face of a cyclidic net) the fit is a least
squares one and its defect is meaningful: zero on a torus, and on a generic
surface the smallest normal-equation eigenvalue ratio falls from 1.5e-10
(diamond size 0.2) to 3.6e-17 (0.025). The next eigenvalue falls as fast
(1.3e-9 to 5.5e-16), so at mesh scale a diamond pins only the second-order
part and leaves the third-order part nearly free: per-edge sigma needs the
advected state and the Dirichlet coupling as regularizers.

Related closed forms: a Darboux cyclide (`albers/darboux_cyclide.hpp`) is
linear in 13 coefficients; cyclidic nets (18.6) fix a cyclide patch from four
concyclic points and one frame.

### 18.5 Geodesics, transport, holonomy, tangent space

- **Transvection** `g = (sigma_b sigma_a)^(1/2)` satisfies
  `g sigma_a g^-1 = sigma_b` and lies in the group; `log g` is in so(4,2) and
  anticommutes with `sigma_a` (a tangent vector). This is the edge transport.
- **Geodesic** `exp(t xi) sigma_a exp(-t xi)`, `xi = log g`: the cyclide's
  slerp. The midpoint m is equidistant: `m sigma_a m = sigma_b`.
- **Holonomy** of a face, `g_ki g_jk g_ij`, fixes `sigma_i` (it lies in the
  stabilizer O(2,1) x O(2,1)); identity for coincident cyclides, nonzero for a
  generic triangle of tori.
- **Tangent space** at sigma: the algebra elements anticommuting with sigma,
  `(X - sigma X sigma)/2`, dimension 9. To first order
  `sigma(xi) = sigma + 2 xi sigma`.
- **Corrector** (`lie_dupin_correct`): damped Gauss-Newton in that tangent
  space, `sum_k w_k |rows_k|^2 + prior |xi|^2`, updated by Cayley. Recovers a
  perturbed torus to 1e-7 from ten contact elements; a heavy prior holds the
  start.
- Measured: the Laguerre offset is a group element that maps tori to tori but
  **not** a geodesic (its log does not anticommute with sigma), so the
  geodesic midpoint of tubes r and r + d is not the tube r + d/2.

The metric on the symmetric space is indefinite (the invariant form is
positive on rotations, negative on boosts). For penalties use a Frobenius norm
in a gauge, or the Euclidean torus subcase (dual quaternion axis frame plus
(R, r), 7 parameters).

### 18.6 Why the one-ring fit misses something

1. **The osculating shape is not determined by the 2-jet.** A torus has 7
   parameters against 6 second-order conditions (1 free: where the axis sits);
   a cyclide 9 against 6 (3 free: exactly the `b` the Dupin row fits from the
   store's cubic slopes). The rows re-decide that freedom every step, so it is
   a convention, not state.
2. **Neighbouring fits are not coupled by integrability.** Two fitted models
   are compatible only if the field of shape operators satisfies Codazzi,
   `(nabla_X S) Y = (nabla_Y S) X`. Position rows, ARAP and a Dirichlet energy
   on S do not see it. The exact gluing condition is flatness: per edge a
   transition `g_ij` in the group, per face holonomy `g_ki g_jk g_ij = I`
   (Cartan: a surface modelled on a homogeneous space; the Cartan curvature
   measures the departure from one global model). A Dirichlet energy on S is
   never zero on a torus (S varies along it); the holonomy is.
3. **The glue penalizes values, not fluxes.** The numerical analogue is a
   Trefftz / discontinuous Galerkin method: local exact solutions (here,
   tori), error only in the interface jumps. Penalizing only the value jump
   (projective dynamics' global step reconciles positions) is consistent but
   blurry; convergence needs the flux jump (normal and frame) with the
   Nitsche / interior-penalty consistency terms. Two neighbouring tori can
   agree on every shared vertex and disagree on the tangent plane and the
   principal frame there.

The natural energy is a **harmonic map into the cyclide space**:
`sum_edges d(sigma_i, sigma_j)^2`, zero iff the patch is a single cyclide.
Precedent one level down: Bryant's conformal Gauss map (each point to its mean
curvature sphere, in de Sitter space) has the Willmore energy as its area, so
"make neighbouring spheres agree" is already a conformally invariant fairing
energy; this moves it from spheres to cyclides.

The discrete theory that builds the gluing in is the **cyclidic net**
(Bobenko, Huhnen-Venedey 2012): a circular quad (four concyclic points) with a
curvature-line frame at one vertex determines a unique Dupin cyclide patch;
the frame at the next vertex is the reflection of the current one in the
edge's perpendicular bisector plane (exact discrete parallel transport of the
shape frame), and adjacent patches meet C^1. The circularity condition is a
discrete Codazzi.

Edge elements in the same direction: Grinspun, Gingold, Reisman, Zorin 2006
(*Computing discrete shape operators on general meshes*: normals at edge
midpoints, CR degrees of freedom, a shape operator per face); Wardetzky,
Bergou, Harmon, Zorin, Grinspun 2007 (*Discrete quadratic curvature energies*:
hinge bending is a CR Laplacian).

### 18.7 Advect, then correct

Carry sigma_i per vertex (and later g_ij per edge) as state:

1. **Advect** (predict): `sigma <- A sigma A^-1` with A the lift of the ring's
   best-fit motion since the last step. The rigid part is the weighted
   Procrustes `(R, t)` the corotational mode already computes (`lie_rigid`);
   a conformal part could come from motor estimation in CGA (Valkenburg,
   Dorst 2011).
2. **Correct** (minimal move): `lie_dupin_correct` toward the ring's contact
   elements, the prior holding the advected state, plus the edge coupling
   (jump or holonomy) and the volume constraint as one linear condition
   `dV . xi = 0` on the 9-dimensional step (Steiner's formula gives dV
   analytically). Volume then no longer fights the landing the way the
   positional corrections did.
3. **Land** the positions on their cyclides as now.

What this fixes: the free axis and `b` parameters persist instead of being
refit; near an umbilic the advected state holds where the data stops
determining the cyclide (the `b` blow-up the ridge `b_ridge` papers over);
flat regions need no apex cap. Remeshing: a split takes the geodesic midpoint
of its two vertices, a collapse the Karcher mean.

Caution from 17.1: ADMM on a soft term converges to the same weighted
minimizer as projective dynamics. The edge jumps only become exact as
constraints, which is Engman's machinery, not a dual variable on a penalty.

### 18.8 Status and plan

Built and tested (16 tests): coordinates, the group, sigma for the torus and
its images, membership, curvatures, normalization, the linear fit,
transvection, log/exp/geodesic, holonomy, tangent basis, corrector. Not wired
into any row. Plan in TODO item 19. A solver built on it: section 19.

## 19. A Poisson solver in Dupin space (`dupin_poisson.py`)

Code: `lie_ops.py` (numpy port of `lie_sphere.hpp` / `lie_dupin.hpp`, checked
against every value in `lie_sphere_tests.hpp` by `check_lie_ops.py`, plus the
pieces below), `dupin_poisson.py` (stage 1 field tests, the solver),
`dupin_compare.py` (the shared suite: `edge_compare.cases()` and the closed
band, same metrics as section 16).

Idea: the shape state is one cyclide `sigma_i` per vertex (a local shape
function, the Trefftz view of 18.6), coupled across the one ring; positions
chase points and normals only. No CR operators in the solver and no
derivative of curvature with respect to positions: the existing solvers
differentiate the read-out (second order in `x`, so a fourth-order gradient),
this one solves second-order Poisson problems on constant factors.

### 19.1 The coupling is a Poisson problem, and needs no gauge

The invariant form on the tangent space `p = {X : X sigma = -sigma X}`
(maps `U -> U^perp`) is `<,>_U (x) <,>_{U^perp}` = (2,1) x (2,1), signature
**(5,4)**: a geodesic harmonic map `sum d(sigma_i, sigma_j)^2` has four
negative stiffness directions and is a saddle, not a Poisson problem.

The ambient form avoids it. `sigma = eta M` (M symmetric) is a linear space,
and `sum_ij w_ij |sigma_i - sigma_j|^2 + mu sum_i M_i |sigma_i - sigma_i*|^2`
is the scalar Laplacian on each entry: one constant SPD factor
`(w L + mu M)`, 36 right-hand sides, then the matrix sign function per vertex
(`normalize`; eta-self-adjointness and trace 0 are linear and survive the
solve, only `sigma^2 = I` is restored). An extrinsic harmonic map, the
ARAP / Alouges pattern.

**No gauge is needed.** A quadratic of Kronecker form `L (x) G` has the same
minimizer for every constant metric `G`, so the Frobenius choice does not
matter, and the sign function commutes with conjugation. The step is
therefore exactly equivariant under any global Lie sphere transformation
(conjugate the field, solve, conjugate back): measured on a noisy torus
field translated by 0.5 / 1 / 3 times the mesh size, change 8e-14 / 8e-11 /
5e-7 (round-off: entries grow to 1.4e4). The (5,4) indefiniteness never
enters. Couplings stay nonnegative with clipped cotan weights.

Stage 1 (torus 2 / 0.7, 24 x 12): an exact field is a fixed point (3e-15); a
field with 1e-2 noise recovers 2.6e-2 -> 3.3e-3 / 1.1e-3 / 6.3e-4 at
`w` 1 / 10 / 100, about 15 ms per solve.

### 19.2 Reading a cyclide off a mesh

- **The linear fit (18.4) is unusable on mesh data.** With analytic normals
  the one-ring fit is exact (3.6e-8). It tolerates normal error up to about
  1e-4 (one ring) / 1e-3 (two ring), linear amplification about 500x / 7x;
  vertex normals of the 24 x 12 torus are off by about 1e-2 rad (0.7
  degrees), and the fit returns non-cyclides at every vertex (error 50).
  Contact elements read off a cubic jet are consistent but the jet is worse
  on a coarse tube (3.6e-2 rad) and the fit still fails. The linear fit
  searches the 20-dimensional projective space for a 9-dimensional manifold:
  data on no cyclide pulls it off the manifold.
- **On the manifold it works.** `osculating_torus(p, n, k1, k2, e1)` builds
  the second-order osculating torus in Lie coordinates as the span of three
  tube spheres (centres on a circle, so linear in (1, cos t, sin t); the
  orientation follows the signs of k1, k2; one free parameter, the axis, by
  convention along e1). Correcting it on the two ring's contact elements
  (`lo.correct`, 9-dim Gauss-Newton, tiny prior) gives 2.5e-2 (24 x 12) and
  7.7e-3 (48 x 24) against the analytic torus with mesh normals, exact
  (1e-15) with analytic normals. The floor is the normals' third-order error.
- The osculating torus alone (no correction) is second order only: the ring
  sits 6.6e-2 off it on the coarse torus, so landing targets are wrong.

### 19.3 Positions

Targets are pose-free: the target mesh's read-out conjugated into each
vertex's body frame `B_i = lie_rigid([t1 t2 n], x_i)` and re-posed on the
live mesh every iterate (so no target positions leak in). Each iterate:
pose, couple (19.1), land every ring neighbour on its vertex's cyclide, and
gather the spokes `p_ij - p_ii` with one constant factor
`A^T W A + eps M` (A the directed spoke incidence).

Landing: a cyclide's contact elements are the pairs (u, v) of null vectors,
u on the conic of U and v on the conic of U^perp (always orthogonal, so
span{u, v} is a contact element): the surface is parametrized by two angles,
and the closest point is a two-parameter Gauss-Newton (`lo.dupin_land`,
exact to 1e-9 on the torus). No warm start across iterates: the conic bases
come from an eigendecomposition and are not continuous in sigma.

What failed, in order (torus tube 0.7 -> 0.9, 24 x 12, shape error 3.4e-2 at
the start):

| Variant | Result | Cause |
|---|---|---|
| landing + ARAP-predicted neighbours | no growth (3.4e-2 flat) | the prediction re-imposes rest lengths; a closed tube must stretch (curvature x circumference = 2 pi) |
| ring pushed by the transvection `(sigma_t sigma_c)^(1/2)` | first iterate 9.8e-3, then diverges | the read-out's third-order noise (19.2) moves `g` by about 0.09 per iterate against a 0.15 signal; and the geodesic of an indefinite metric minimizes nothing Euclidean (large boosts) |
| same, osculating read-out, `g^alpha` | diverges | as above |
| pure landing (projective dynamics on `sum dist(x_j, sigma_i)^2`) | 3.4e-2 -> 2.2e-2 in 60, then drifts back up | slow and unstable |
| pure landing + shape coupling (`w` 10) | 4.2e-3 after one iterate, 2.5e-3 by 10, then a growing mode (about 1.15x per iterate) | the body frame's t1 from the first ring edge twists with the mesh and re-poses the cyclide (feedback); a weak ARAP energy only delays it |
| pure landing + coupling + principal-direction frames | **2.5e-3 by iterate 6, stable over 80** | |

So the shape coupling is what makes it work: one iterate removes most of the
error, which the uncoupled landing never does. That is the role 18.4 gave the
Dirichlet coupling (the third-order part of a per-vertex read-out is noise
the neighbours average out).

### 19.4 The shared suite (first cut)

`dupin_compare.py`, defaults (`corrected` read-out, principal frames, cotan
spokes, `w_shape` 10), 40 iterates:

| Case | Shape | Fixed point | Iterates | Note |
|---|---|---|---|---|
| torus tube 0.7 -> 0.9 | 3.4e-2 -> 2.5e-3 | 2.5e-3 | 6 | at the read-out floor |
| closed band r 1 -> 1.3 | 1.2e-1 -> 4.9e-3 | 4.4e-3 | 6 | at the read-out floor |
| torus R 2 -> 2.5 | 4.9e-2 -> 1.4e-3 (it 10, relax 0.25), then overshoots | 6e-4 | | no merit picks it (landing residual is not monotone in shape) |
| isometric cylinder | fails (fixed point 1.4e-1) | | | flat start, k2 = 0 everywhere, open boundary |
| sphere cap, cube -> sphere / ellipsoid | every cyclide invalid | | | umbilic targets: a sphere is not a cyclide |

Against section 16 (exponent solver: residual to 1e-8 and shape 1e-10 in 4
to 30 iterates on every case): where it works it converges in a similar
number of iterates, but to the read-out floor (2.5e-3 at 24 x 12), not to
machine precision, because the target is a per-vertex cyclide read-out of
the target mesh and the solver reproduces that read-out, not the CR one.
Python cost is about 0.9 s per iterate, nearly all landing (two-angle
Gauss-Newton per spoke with a finite-difference Jacobian); the factors are
constant and the solves are milliseconds.

### 19.5 Open

- **Degenerate targets.** Umbilics (spheres) and flat / cylindrical regions
  have no cyclide coordinate (U degenerates). Sphere elements (the curvature
  sphere pair `s_a = k_a P + T`, or the mean curvature sphere) are defined
  there; `|grad_{e_a} s_a|^2 = |e_a . grad k_a|^2` in that representation,
  which is the existing Dupin row, so sphere elements reproduce today's
  rows and cyclide elements add the third order. A hybrid by vertex
  (cyclide where anisotropic, sphere where umbilic) is the next experiment.
- **A fixed energy or merit.** Posing targets on the live mesh makes the
  iteration a fixed-point map, not descent; torus R overshoots and the
  landing residual does not select the best iterate. Candidates: freeze the
  pose for a few iterates (a fixed energy, plain projective dynamics, then
  re-pose), or the merit on the coupled field's Dirichlet energy plus
  landing.
- **The read-out floor.** 2.5e-3 at 24 x 12 from mesh-normal error at third
  order; better normals (from the landed cyclides themselves, iterated) or
  the coupling inside the read-out.
- **Speed.** Vectorized landing (all spokes at once, analytic Jacobian of
  the decode) before any timing comparison.
- **Lie Procrustes (route 2).** Built and tested in `lie_ops`
  (`lie_procrustes`, generalized polar `group_polar`): on contact elements of
  one cyclide the rows leave the pencil {A, A sigma_0} (sigma_0 fixes each
  of them), both transporting sigma_0 the same; a prior projected on the two
  smallest eigenvectors picks one. Not yet used in the solver.

### 19.6 The monge-demo port (design only)

- **Roles.** Hard: the coupling only (store = mesh read-out; store rows then
  act on the mesh exactly, no `k_coupling` calibration). Soft (material, so
  the shape can grow and stretch): volume about the Dupin spine (per ring,
  spine = curvature-sphere centre from the vertex cyclide, no `1/k_hi` /
  `apex_cap`), position shear, area, area ratio. Soft fairing: `w_dupin`,
  `w_minimizer`, `w_grad`.
- **All-soft is not orthogonal.** With no hard rows the KKT step is
  Gauss-Newton on the weighted sum and lands on the same compromise as
  `projection_solver` (17.1); orthogonality exists only against hard rows.
  Counting: each per-ring hard scalar costs about V on 3V positions (shear
  plus area is about Soliman's 2V conformal class, leaving V of bending for
  the fairing energies; a third locks), the Dupin rows are about 4V.
  Lexicographic (priority) solving of the fairing energies is the route if
  they should resolve perpendicular to each other.
- **Surface tension.** Soft `gamma sum A_f`, gradient the mean-curvature
  normal: second-order smoothing, faster than the fourth-order fairing rows.
  In this solver it is `gamma L_cot` in the constant position matrix
  (`SolverConfig::tension`, implicit mean curvature flow); it balances
  against soft area as a membrane at a preferred area. In Engman,
  `face_area`'s multiplier is the local tension.
- **Store fork.** The edge store, its coupling and the rows `w_dupin` /
  `w_minimizer` / `w_grad` are CR. Either bridge (edge tensors read from the
  vertex cyclides by `lie_dupin_curvatures` at edge midpoints) or replace
  (the vertex sigma field is the store, Lie-space versions of the rows; the
  Dupin rows are then redundant since each element is a cyclide). Decided
  after the degenerate cases work.

### 19.7 Edge spheres: a finite-element cyclidic net (`lie_edges.py`)

Framing: Bobenko and Huhnen-Venedey's cyclidic nets
(`reference/1101.5955v2.pdf`) as finite elements instead of surface patches.
In a cyclidic net, vertices are contact elements with curvature-line frames,
each edge's endpoints share a sphere exactly (principal contact element net),
each circular quad carries a cyclide patch, patches meet C^1, and the mesh
must follow curvature lines. Here, vertices carry cyclides (shape functions,
solved for), edges carry spheres, faces carry planes, and any triangle mesh
works; the net's exact conditions become measured residuals.

**What each element holds.** A face holds a plane `T_f` (null) and its
points, never a sphere. An edge holds what its two faces give,
`W_e = span{T_f, T_g}` (signature (1,1)): one dihedral. The spheres tangent
to both planes are the null vectors of `W_e^perp` (signature (3,1)), a
2-parameter family; a length scale picks one. So an edge holds one sphere,
one directional curvature. Twist, the full tensor and third order are
composed at vertices from several edges. Oriented spheres carry saddles
(edges of opposite sign).

**Two edge spheres.**

- `hinge` (primal, across, from the face planes):
  `k_perp = 2 tan(theta / 2) / l*`, `l* = (|e| / 2)(cot alpha + cot beta)`
  the circumcentric dual length (a sphere through a face's vertices meets
  the face plane's parallel at the circumcentre; exact for a polygon's
  incircle). `l* = 0` on the diagonal of a cyclic quad (four concyclic
  points are coplanar): such a hinge measures nothing across itself, so
  hinges are weighted by `max(l*, 0) |e|`. A first version with
  `d = (A_f + A_g) / (2 |e|)` read 2 to 2.7 relative error on the cube-sphere
  from exactly those diagonals. Converges O(h^2) on reliable hinges
  (5.2e-2, 1.5e-2, 3.9e-3 at cube 4, 8, 16).
- `contact` (dual, along, from the endpoint contact elements): the sphere in
  `L_i cap L_j`. It exists iff the normal lines meet, and then
  `x_j - x_i = (n_i - n_j) / k` (the discrete Rodrigues relation). Closed form
  `k_par = -(dn . dx) / |dx|^2`, exact on circles (1e-15 on a sphere and on
  the torus's curvature-line edges with exact normals). The Lie least-squares
  common sphere (`common_sphere`: `s = a P_i + b T_i`, rows against
  `L_j^perp`, 4 x 2) agrees on curvature-line edges, and its residual is
  0 there (9e-15) and 0.3 to 0.4 on the torus's diagonals: the twist the
  hinge cannot see. With mesh vertex normals it reads 0.27 off on the
  cube-sphere (valence-3 corners); the targets use the same read-out, so it
  is consistent.
- Euler is linear on the contact line:
  `s(phi) = cos^2 phi s_1 + sin^2 phi s_2` (1e-16), so a vertex cyclide
  predicts every incident edge's sphere linearly. (Not yet used: 19.8.)

Checks: `check_edge_spheres.py` (signs, exactness, twist residual, saddle
signs +0.99 / -0.99, hinge projection to 3e-11).

**Solver.** Edge curvatures are Euclidean invariants, so targets read off
the target mesh are pose free (no body frames) and the merit
`sum w (k(X) - k*)^2` is a fixed energy. Projective dynamics with one
constant factor: the hinge turns its wings about the edge to the target
dihedral (closed form; `l*` unchanged); the contact edge takes a length
blended from the Rodrigues length `turning / k*` (trusted as
`(|k*| / mean |k*|)^2`, capped at 1) and a length memory (rest re-seeded
every 20 iterates, plastic), capped at 25% against the memory. 1 ms per
iterate in Python (against 0.9 s for the vertex-cyclide landing loop), and
every fixed point is exact (2e-17 to 2e-16).

Shared suite (`lie_compare.py`, 400 iterates, keep-best by merit):

| Case | Shape | Note |
|---|---|---|
| torus tube 0.7 -> 0.9 | 3.4e-2 -> **1.1e-4** | the Rodrigues lengths are the stretch outlet |
| closed band r 1 -> 1.3 | 1.2e-1 -> 6.4e-3 | still descending |
| cube -> sphere | 6.1e-2 -> 1.0e-2 | the first Lie variant to move it (umbilic target) |
| torus R 2 -> 2.5 | 4.9e-2 -> 2.1e-2 | |
| cube -> ellipsoid | 5.2e-2 -> 4.4e-2 | |
| isometric cylinder, sphere cap, saddle (flat starts) | fail (2e-1 to 5e-1) | below |

Without the length memory every open case blew up (strain to 8): lengths
came only from `turning / k*`, ill-conditioned where `k*` is small and wrong
before the normals have turned, and the per-iterate cap compounded.

**Two limits, both seen before.**

1. **Stretch or bend is a global choice** (section 16, "faces shrink instead
   of bending from flat"). Along an edge, `k = turning / L` is met by turning
   (bend) or shortening (stretch). On a closed tube topology fixes the
   turning, so the Rodrigues length is right; on an open sheet from flat the
   turning is still 0 and the local projection shrinks.
2. **Hinge bending from flat is diffusive** (section 8, "stencil projection
   is monotone but diffusive"). Every vertex is a hinge endpoint and a wing,
   so wing turns cancel in the gather and only the boundary moves: hinges
   alone with elastic memory take the cylinder's merit 0.63 -> 0.22 in 600
   iterates, and the sphere cap does not move.

Both point to the same missing global unknown: **the face planes**. Faces
hold planes, and an edge sphere prescribes the relative turn of its two
planes (the dihedral) and, through Rodrigues, how far apart their contact
elements sit. Section 8's rotation form is that lift: per-face rotations
synchronized by relative rotations across edges (one constant factor),
gathered by `L_cot X = div(R T0)`, strip to cylinder in about 12 iterates
from flat. It answers the earlier question of how to lift local sphere
projections globally: through the face planes (rotations), with the
Rodrigues lengths as the stretch outlet, not through per-edge frames
(section 7).

### 19.8 Next

- Face-plane sync: per-face rotations with the hinge spheres' dihedrals as
  relative rotations (section 8, factored once), positions by
  `L_cot X = div(R T0)` with `T0` carrying the Rodrigues / memory lengths.
- Vertex cyclides as the composer: Euler rows from each vertex's curvature
  spheres to its edges' spheres (anisotropy-weighted), for edits and the
  third-order coupling (tried in 19.9).
- The suite with both, plus the stretch-or-bend choice made globally
  (lengths as unknowns in the sync, or continuation: bend first, then
  admit stretch).

On the torus R stall (edges only), measured per parallel: every parallel
reaches about 60% of its needed stretch (outer 1.107 of 1.185, top 1.151 of
1.25, inner 1.250 of 1.385; tube edges 1.011 of 1.0), so the zero normal
curvature on top (the parallel's curvature sphere is the tangent plane) is
not the bottleneck; the Rodrigues / memory blend is. Structurally, though,
edge spheres carry normal curvature only: the top parallel is a circle of
radius R with purely geodesic curvature, which no edge sphere sees.
Cyclidic-net edges are circles (curvature-line arcs), not spheres.

Why one sphere per edge (the `contact` sphere) cannot act alone here: the
common-sphere condition is the vector Rodrigues relation
`k (x_j - x_i) = n_i - n_j`, linear in positions and normals together, but
the normals are derived from positions, so only its length along the edge is
used: it sets lengths and never turns a normal. Both failures of 19.7 follow
(stretch chosen locally, bending only by cancelling wing turns). With the
vertex contact elements as state (a principal contact element net proper:
`x` and `n` unknown, unit normals restored in the local step), one sphere
per edge is a linear joint system on a constant factor that makes the
stretch-or-bend choice globally; exact on curvature-line edges, the twist on
others composed at the vertex. Not built.

### 19.9 The coupled net (`lie_net.py`): vertex cyclides + edge spheres

The original sketch: posed vertex cyclides (principal frames, corrected
read-out, shape Poisson) predict their edges' along / across spheres by
Euler (`e_U = dp / dth_v` from the conic parametrization: the U sphere held
fixed traces the U curvature line); each edge sphere is glued as the blend
of the edge's target and its endpoints' predictions, weighted by anisotropy
`|k_u - k_v| / (|k_u| + |k_v|)` (0 at umbilics); positions from one constant
factor with three terms (ring landing blended by anisotropy, hinge wing
turns, Rodrigues + plastic memory lengths).

| Case | edges only, 400 it | net, 150 it | net fixed point (last iterate) |
|---|---|---|---|
| torus tube | **1.1e-4** | 3.6e-3 | 6.7e-4 |
| closed band | **6.4e-3** | 1.6e-2 | 4.6e-3 |
| cube -> sphere | **1.0e-2** | 1.2e-2 | 6e-16 |
| torus R | **2.1e-2** | 3.7e-2 | 3.1e-4 |
| cube -> ellipsoid | **4.4e-2** | 4.7e-2 | 5.6e-3 |
| cylinder / cap / saddle | fail | fail | 9.4e-2 / 1.7e-2 / 4.9e-1 |

Edges-only fixed points are exact (1e-14); the net's drift. The cyclides'
Euler predictions disagree with the edge targets on the target mesh itself:
each cyclide is an independent read-out from mesh normals at third order
(about 2.5e-2 at 24 x 12, 19.2), while the edge read-out is consistent by
construction. Blended in, the cyclides add bias, not information, and cost
about 1000x per iterate (ring landing in Python). Cube -> sphere is
unaffected because umbilic vertices get zero weight.

What this does not test: edits (no target mesh, where a composer is needed)
and third-order consistency. What it says about the design: the vertex
cyclide must not be a second, independent read-out. It should be composed
from its edges (fitted to its incident edge spheres by Euler rows plus its
contact element), so that at a target mesh its predictions equal the edge
targets by construction and the fixed points stay exact; its job is then to
carry what edges cannot (twist, the geodesic curvature of curvature lines,
third order), not to re-measure what they can.

### 19.10 The contact element net (`contact_net.py`)

The principal contact element net (Bobenko and Suris) as a solver: each
vertex holds a contact element with **both** position and normal unknown,
and each edge holds one sphere, the common sphere of its endpoints.

**Read-out (exact, pose free).** Per edge, with `nbar = normalize(n_i + n_j)`
(for unit normals `dn . nbar = 0` exactly) and `dx_t` the part of `dx`
normal to `nbar`,

    dn = -k dx_t - tau (nbar x dx).

`k` is the curvature along the edge (the common sphere; exists exactly when
`tau = 0`), `tau` the twist (geodesic torsion), and the two are all of `dn`
(reconstruction 4e-16; `tau` 3e-15 on the torus's curvature-line edges,
median 0.67 on its diagonals). Per directed spoke the normal offset
`c_ij = n_i . d_ij / |d_ij|` (0 for a smooth surface, `O(h)` on a mesh).
All three are read off the target.

**Rows.** None mixes coordinates (cross terms lagged into the right-hand
side), so the global step is one constant factor on `(x_c, n_c)`, solved
three times per iterate:

- Rodrigues, per edge: `k (x_j - x_i) + (n_j - n_i) = k (dx . nbar) nbar
  - tau (nbar x dx)`. Positions and normals in the same row: the global
  solve splits each curvature into turning and length.
- Tangency, per spoke (positions follow normals): `x_j - x_i` toward its
  part tangent to `n_i` plus `c_ij |d| n_i`.
- Normal fit, per vertex (**normals follow positions**): the unit normal that
  best meets the ring's offsets for the current spokes,
  `min |D n - c|^2 + mu |n - n_i|^2`, normalized.
- Length memory (plastic, `w` 0.1), anchor.

Every target is a fixed point (to 1.6e-10, the factor's conditioning).

**Two findings on the way.**

1. Without the twist the read-out is incomplete and even fixed points drift
   (7e-2 on the tube): on a triangle mesh most edges are not curvature-line
   edges, and their sphere is only part of `dn`.
2. With one-way tangency (positions follow normals, normals held only by
   unit length and the Rodrigues differences), the tori and the ellipsoid
   turn around after 12 to 87 iterates: the state normals detach from the
   surface (2.5 degrees at iterate 30, 19 at 100, flipped by 200) and drag
   their spokes. That is the director-shear mode of TODO item 16. The normal
   fit (two-way tangency) removes it; every case then descends
   monotonically.

**Shared suite** (`lie_compare.py --solver=cnet`, 10000 iterates, about 2 to
4 ms per iterate in Python):

| Case | Shape | Iterates to merit / start < 1e-3 / 1e-6 / 1e-9 | Exponent solver (16) |
|---|---|---|---|
| torus R 2 -> 2.5 | 4.9e-2 -> **3.9e-15** | 900 / 1719 / 2536 | 6.3e-10 |
| torus tube 0.7 -> 0.9 | 3.4e-2 -> **4.1e-14** | 1102 / 2908 / 4713 | 3.4e-10 |
| cube -> sphere | 6.1e-2 -> **1.6e-13** | 427 / 1847 / 3341 | 4.5e-10 |
| cube -> ellipsoid | 5.2e-2 -> **3.5e-12** | 865 / 4151 / 8880 | 2.2e-10 |
| closed band r 1 -> 1.3 | 1.2e-1 -> 1.2e-9 | 2909 / 7319 / - | 1.5e-9 |
| isometric cylinder | 1.1e-1 -> 3.2e-4 | 3957 / - / - | 2.0e-6 |
| sphere cap | 9.8e-2 -> 6.3e-4 | 2835 / - / - | 2.9e-4 |
| saddle a 0.8 | 1.6e-1 -> 3.1e-3 | - / - / - | |

Against every earlier Lie variant (19.3 to 19.9) this is the first that
reaches machine precision, and it does so on the umbilic (cube -> sphere)
and non-Dupin (ellipsoid) targets too: the edge spheres plus normals as
state carry what neither the vertex cyclides nor the face-normal edge
spheres could. It confirms 19.8's argument: one sphere per edge suffices
once the contact elements are the state, because the Rodrigues relation is
then a joint linear row in positions and normals and the stretch-or-bend
choice is global.

**Cost and what is left.**

- Linear convergence of projective dynamics: hundreds of iterates to 1e-3,
  thousands to 1e-9. Each iterate is a backsolve on a constant
  `2V x 2V` factor plus per-vertex 3 x 3 normal fits, so per-iterate cost in
  C++ is a backsolve; the exponent solver needs 4 to 90 iterates but a fresh
  factor per restart. An Anderson corrector (TODO item 12) is the obvious
  next step on the iterate count.
- Open sheets converge slowest (cylinder, cap, saddle stall near 1e-3
  relative merit): the boundary has incomplete rings and no hinge memory;
  to look at.
- Targets here are read off a target mesh (all three read-outs), so this is
  the reconstruction problem from discrete contact-element data, Wang et
  al.'s fundamental forms in Lie form. For edits, the targets `k`, `tau`, `c`
  come from the shape layer: the vertex cyclides (fitted to the net's own
  contact elements, the Bobenko step, not to mesh normals) composing the
  edge spheres by Euler and supplying `tau` (`k_u - k_v`) `sin cos`.
- Not yet done: stage 2 (cyclides on the converged net as a diagnostic is
  uninformative on this suite: the net reproduces the target's own
  contact elements, so the fit repeats 19.2), timing against C++, remeshing
  (contact elements are vertex datums).

### 19.11 Accelerating the net

Two causes of the slow rate: lagged coupling (twist, tangency directions and
the normal fit frozen per iterate) and low-frequency modes relaxing through
local projections (worst on the open sheets). Three accelerators, same rows
and targets:

- **Anderson** (Peng et al. 2018, window 5, their merit safeguard) on the
  fixed-point map `z = (x, n) -> G(z)`, one constant factor. It needs the map
  to be a function of `z`: the plastic memory reseeded every 20 iterates
  changes the map and clears the history, so the memory is reseeded every
  iterate instead (a proximal damping). With that it is 3 to 20 times faster
  than plain and reaches 1e-10 to 1e-14 on every case including the open
  sheets.
- **Chebyshev** (Wang 2015, `rho` from the plain iterates): 2 to 5 times
  faster than plain on closed shapes, little on the open sheets.
- **Gauss-Newton / LM on the un-lagged rows** (`ContactNetGN`): Rodrigues with
  its `nbar` and twist, scalar tangency `n_i . d / |d| - c` (both ways at
  once, replacing the normal fit), unit length, length memory reseeded every
  iterate; vectorized central-difference Jacobian over each row's local
  unknowns, factor of `J^T J + lambda diag` kept 5 iterates (the hepworth
  pattern). With the memory at 0.1 it crawls on the flat starts (the
  per-iterate proximal term limits the step); at 0.001 it converges in
  Weingarten-like iterate counts on everything.

Iterates to merit / start below 1e-9 (shape error at the end; 24 x 12 meshes):

| Case | plain | Chebyshev | Anderson | GN, refresh 5 (factors / backsolves) |
|---|---|---|---|---|
| torus tube | 4713 | 1378 | 334 | **8** (2 / 10), 7.6e-16 |
| torus R | 2536 | 215 | 147 | **8** (2 / 10), 1.1e-15 |
| cube -> sphere | 3341 | 760 | 171 | **9** (5 / 16), 4.7e-16 |
| cube -> ellipsoid | 8880 | 1917 | 380 | **20** (15 / 35), 2.9e-16 |
| closed band | - | 1510 | 1e-6 at 136 | **11** (3 / 14), 3.4e-14 |
| isometric cylinder | - | - | 1e-6 at 1413 | **44** (18 / 72), 1.1e-12 |
| sphere cap | - | - | 1503 | **48** (19 / 83), 2.5e-12 |
| saddle | - | - | 1831 | **34** (12 / 51), 9.6e-13 |

(plain and Chebyshev with the 0.1 plastic memory; Anderson and GN with the
0.001 per-iterate memory, now the defaults.)

Cost, measured in scipy on tori (factor includes assembly for PD):

| Vertices | PD factor (2V) | PD backsolve | GN factor (6V) | GN backsolve |
|---|---|---|---|---|
| 288 | 7 ms | 0.06 ms | 27 ms | 0.31 ms |
| 1152 | 32 ms | 0.34 ms | 0.22 s | 2.5 ms |
| 4608 | 0.23 s | 2.1 ms | 3.4 s | 17 ms |

At 1152 vertices that puts a torus-tube solve at about 0.5 s either way
(GN: 2 factors and 10 backsolves; Anderson: about 1500 backsolves on one
factor), the same order as the exponent solver's 0.8 s refactoring frame
(C++, 2304 faces) and the Engman KKT factor (0.13 s, C++). The scipy numbers
are not CHOLMOD; a C++ timing decides. GN wins on iterates and on the open
sheets; Anderson keeps one factor for the whole solve (the Rodrigues `k` sit
in the matrix, so either refactors when the targets change).

Every fixed point stays exact (4e-11 against the factor's conditioning).
The Gauss-Newton rows are also the right object for the Engman route: a
constraint `C(x, n)` there is one of these residual blocks, and the
multipliers then report which targets an edit cannot realize.

**Correction (fading memory).** The open-sheet rows of the table above were
an artifact: the length memory, reseeded every iterate, is a proximal term
on every step, and on the soft bending modes of an open sheet even 0.001
of it turns quadratic convergence linear (boundary and interior errors were
equal; only 3 to 10 LM steps were rejected). Without it Gauss-Newton takes
the cylinder, cap and saddle to 1e-9 in 6 or 7 iterates; with it the
closed cube -> ellipsoid is faster (14 against 27). The memory now fades
with the merit (`w_m * merit / merit_start`, updated at each refactor), at
weight 0.01 and refresh 5: every case reaches 1e-9 in 7 to 15 iterates
with 2 to 6 factorizations (`bench_contact_net.py cold`).

### 19.12 Warm frames, the C++ port and the demo pathway

**Warm frames** (`bench_contact_net.py warm`: the target interpolated over
20 frames, each warm started): the projection solve does not get cheaper,
700 to 1500 iterates per frame to its own 1e-9 (linear rate, the same
whatever the start), and the reference-`k0` factor (target change lagged as
`(k - k0) dx` so the factor survives) loses accuracy as the targets drift
(cube -> ellipsoid 3.5e-3); Gauss-Newton per frame costs about 15
factorizations against 2 to 6 for the whole deformation cold. The
backsolve-only per-frame loop does not materialize at that tolerance.

**C++** (`kusama/contact_net.hpp`, rows generated by
`contact_net_jacobian.py` into `generated/contact_net_generated.hpp`,
`test/contact_net_tests.hpp`, all passing): read-out exact, twist 0 on
curvature lines, generated Jacobians against differences, every suite target
a fixed point under both solvers. Gauss-Newton reaches 1e-9 in 7 to 15
iterates (2 to 6 factorizations, 7 to 29 ms per solve at 100 to 300
vertices, shapes 1e-16) on torus tube / R, cube -> sphere / ellipsoid, cap,
saddle, cylinder; the projection solve with Anderson in 123 to 430 iterates
(45 to 92 ms), Chebyshev slower.

**Demo pathway** (`curvature_solver::contact_net`,
`hepworth/block/contact_net_coupling.hpp`): the exponent pattern, the CR
derivative gather rows with `k_gather` (same calibration against `bend`),
targets in increment form from the edge store (`k += dh_0`, `tau += dh_1`,
`c += dk |d| / 2`; the store's frame is the net's convention, so no edit is
an exact fixed point), the store anchored at `w_coupling`. Probe, torus test
target (R 0.75 / tube 0.25 -> 1 / 0.4), 6 frames of 45 iterates at
`w_coupling` 1:

| Solver | R | tube | s / frame |
|---|---|---|---|
| contact_net, inner solve every iterate | 1.026 | 0.404 | 47 to 79 |
| contact_net, 10 solves per frame | 1.029 | 0.390 | 12.5 |
| contact_net, 6 | 1.018 | 0.372 | 9 |
| contact_net, 3 (default) | 0.974 | 0.336 | 5 |
| eigensatz | 0.874 | 0.446 (overshoot) | 3 |
| exponent | 0.881 | 0.440 | 3.2 |
| weingarten | 0.751 | 0.346 | 6 |

The net tracks the target best (R especially); its cost is the inner solves
per frame (each a few Gauss-Newton factorizations of the `6V` system on the
1152-vertex torus), set by `solves_per_frame` (the store settles over the
outer iterates and each move re-solves, warm started from the last net).

**Global accelerators** (`projection_solver::accel`, demo
`fairing_step::accel`): Anderson and Chebyshev as a corrector after the
global step, Peng's safeguard, history cleared on refactor. In the demo they
changed energies and shapes only in the 3rd or 4th digit and doubled the
frame time of `eigensatz` / `exponent` (3.1 -> 7.2 s): the safeguard's extra
projection runs the stateful couplings' inner steps, and the map is not a
pure function of the iterate there (2 to 26 of 44 corrections taken). Off
by default; they pay where the map is fixed (the contact net's `pd`).

### 19.13 Native rows in the projection solver

The inner pathway solves the net on top of the projection solver; the native
one (`hepworth/block/contact_net_gn.hpp`, demo `contact_net_mode::native`,
the default) puts the net's rows in it as Gauss-Newton rows, as
`weingarten_gn` is: positions, per-vertex normals (a `vec3_block` on
`shell_cr_block`, `normal_site`, reseeded from area-weighted face normals each
frame, mass 1e-8) and the edge store in one factor.

**Sphere rows fail.** Rows in Lie form (each neighbour on its edge's
curvature sphere, `ContactNetSphereGN` in Python) are exact on curvature
lines only: off them the neighbour slides on the sphere, which the rows do
not see. The native rows keep the Rodrigues form per edge
(`contact_rodrigues_kt`, `k = k_0 + dS_0`, `tau = tau_0 + dS_1`, store
columns), tangency per spoke as point on sphere (`n_i . d / |d| = |d| kappa /
2`, live `|d|`), unit rows on the normals; row scale `k A_e / |e|`,
`k_contact_gn = 14.2` from the `bend` calibration (0.97 / 1.09 of `bend`'s
response; at `k = 1` it was 0.005).

**The rows are right.** With the target's own tangency offsets, rows rebuilt
every frame at the current state, the closed shapes reach machine precision in
12 frames: torus tube 8.6e-14 (24 iterates, 3 refactors), R 1.0e-13,
cube -> sphere 1.4e-11 (77, 13), cube -> ellipsoid 8.9e-14.

**The tangency target.** The frame start's offset (`c - k |d| / 2` held,
`k` the edge's) holds the old shape: the one-sided `n_i . d / |d|` is the
vertex's normal curvature toward `j`, not the edge's `k`. Tori survive it
(9.9e-8, 2.1e-8), cube -> sphere ends at 1.2 and the cap at 3.5. The target
now comes from the store's vertex tensor (`tangency::vertex`, default):
`kappa_ij = kappa_0 + sum_e a_e . dS_e`, the spoke's normal curvature in
vertex `i`'s tensor (`kusama::weingarten::vertex_tensor`, the weighted mean
of its edge tensors in its plane), coefficients frozen at the frame start, so
linear in the store and with all three store components (the edge rows see
only `k11`, `k12`; `k22` enters only here).

**What else the rows need.** On their own the rows have a collapsed zero:
all points and normals equal satisfy Rodrigues, tangency and unit rows, and
without a memory the cubes shrink to a point (shape error 0.5, independent of
mass). Positions also need damping (`projection_solver` has no LM): mass
1e-3 or the demo's Dirichlet memory turns the open sheets from stalls into
convergence. Vertex target, store moved by CR edge-reading increments (as
the demo's is), mass 1e-3, 12 frames, shape error:

| case | memory at the demo's ratio (= `w`) | no memory |
|---|---|---|
| torus tube | 1.9e-3 | 2.1e-6 |
| torus R | 6.3e-3 | 8.3e-5 |
| cube -> sphere | 6.0e-2 | collapse |
| cube -> ellipsoid | 2.4e-1 (6e-2 at frame 3) | collapse |
| sphere cap | 4.2e-3 | 1.8e-2 |
| saddle | 4.3e-2 | 3.9e-2 |
| cylinder | 4.8e-3 | 2.2e-2 |

The floors of a few percent are the read-out mismatch: the store moves by CR
read-out increments, and the parts the net cannot express are a bias
(Weingarten's rows read the CR quantity itself, so it has none). The memory
slows the closed shapes (it is a proximal term per frame).

**The demo.** Probe, torus target, 6 frames: native (vertex or edge target)
and `weingarten` give the same shape (frame 5 R 0.687 / 0.688 / 0.698, tube
0.392 / 0.390 / 0.393), native 6 to 7 s per frame against 5 at 45 iterates.
The frames were converged: at 150 iterates the energies agree with 45 to
five digits and the shapes to three (`|dq|` 1e-4 against 0.03, motion
the energy does not see), and 10 or 20 iterates give R 0.69 to 0.70, tube
0.390 to 0.393 at 1.4 to 3.8 s per frame. So the demo's shape is the
minimum of its energy, not of the solver: R moves away from the target
(0.75 -> 0.69 against 1) because the other rows hold it, and no curvature
pathway or accelerator changes that.

### 19.14 Lie edge net: two spheres per edge, composed at the vertices (`lie_edge_net.py`)

Why the earlier Lie rows failed, and what this changes. Pure sphere sharing
(Bobenko-Suris principal nets) holds only on curvature-line edges; a triangle
has three edge directions and a surface two principal ones, so most edges
carry twist, and without it even the target drifts 7e-2 (19.10). Incidence
rows alone are scalars: they do not say which way an edge leaves its vertex,
so neighbours slide on their spheres. Vertex cyclides were a second,
independent read-out (bias, 19.9) and do not exist at umbilics (19.4).

State: contact elements (`x`, `n`) and per edge `H_e = (a, tau, b)`, the
along sphere `a P + T`, the twist, the across sphere `b P + T` (the CR edge
store in Lie form, so a C++ port keeps its edge store). Vertices compose: an
auxiliary tangent tensor `S_i` whose Euler samples must match every incident
edge (`s(phi) = cos^2 phi s_1 + sin^2 phi s_2`, linear in R^{4,2}; checked
to 2e-16). Rows: Rodrigues with `(a, tau)` (direction, the sliding fix), Lie
incidence `<a P(x_i) + T(x_i, n_i), P(x_j)> = 0` from both ends, hinge
`2 tan(theta / 2) = b l*`, vertex Euler consistency, ambient Poisson coupling
of `S` (no frame gauge), unit normals, targets on `H`, length memory.
`check_lie_edge_net.py`: all pass (Lie rows against `lie_ops`, Euler,
Rodrigues and hinge exact at the targets, Jacobians, rigid invariance).

The read-outs are not mutually consistent at a target (rms bias: incidence
1e-2 to 2.6e-2, vertex composition 1e-2 to 0.11), so targets are not fixed
points and the solve lands on a compromise. That is the accepted lower
precision. Consequences: a faded memory lets open sheets stretch to absorb
the leftover (strain 0.55 to 0.69, shapes 0.11 to 0.49), so the memory is
kept (`w_memory` 0.1, not faded); the open-sheet floor sits on the boundary
(drift / size 1e-2 at the boundary, 1.3e-3 inside on the cylinder and cap),
helped a little by scaling boundary vertices' incidence and vertex rows
(`w_boundary` 0.1). No single row sets the floor (ablations of incidence,
vertex, coupling, hinge, target weight all within a factor 2).

Suite (`lie_compare.py --solver=lsn`, 60 iterates, defaults):

| case | shape | fixed point | iterates | `cnet_gn` (19.11) |
|---|---|---|---|---|
| torus tube | 3.4e-2 -> 1.1e-3 | 1.1e-3 | 28 | machine precision |
| torus R | 4.9e-2 -> 2.4e-3 | 2.3e-3 | 36 | machine precision |
| cube -> sphere | 6.1e-2 -> 3.0e-3 | 4.6e-3 | 60 | machine precision |
| cube -> ellipsoid | 5.2e-2 -> 4.6e-3 | 2.7e-3 | 60 | machine precision |
| closed band | 1.2e-1 -> 1.4e-2 | 1.4e-2 | 60 | 1.2e-9 |
| isometric cylinder | 1.1e-1 -> 2.6e-2 | 1.4e-2 | 60 | 3.2e-4 (pd) |
| sphere cap | 9.8e-2 -> 3.4e-2 | 7.5e-3 | 60 | 6.3e-4 (pd) |
| saddle | 1.6e-1 -> 1.9e-2 | 2.7e-2 | 60 | 3.1e-3 (pd) |

Every case converges, umbilic (cube -> sphere) included, where every cyclide
variant failed; closed shapes at 1e-3 to 5e-3, open sheets and the band at
1.4e-2 to 3.4e-2 (the pass bar was about 1e-2). Python cost 25 to 500 ms per
iterate (finite-difference Jacobians on `9V + 3E` unknowns).

True edits (`edits_lie_edge_net.py`, no target mesh): cube -> sphere's
sphere with every curvature times 1.2 reaches radius 1.162 (want 1.152, the
memory resists the shrink) and stays round to 7.6e-3 (radius spread over
mean), 27 factorizations; `ContactNetGN` on the same edit 1.155, 1.7e-3, 13.
The shape error rises (0.10 -> 0.13) only because vertices slide on the
sphere. An incompatible edit (the torus told it is umbilic everywhere)
stays finite but crumples (mean |tau| 5.3): nothing reports that the edit
cannot be met (Engman multipliers would).

So: a working Lie-sphere solver at about 1e-3 to 3e-2, robust where the
cyclides were not, but not better than the contact net on any case, and
slower in Python. What would change that: consistent read-outs (the across
sphere and the vertex composition from the same contact elements, so targets
are fixed points), boundary rows, and edits from the shape layer rather than
a scaled read-out. Stage 2 (vertex cyclides where anisotropic) not started.

### 19.15 Lie cyclide net: edge spheres, vertex cyclides, Poisson coupling (`lie_cyclide_net.py`)

The Lie edge net (19.14) plus a Dupin cyclide `sigma_i` per vertex, composed
from the net rather than read off the mesh (the 19.9 lesson): fitted (9-dim
Gauss-Newton in the tangent space, local coordinates) to its own contact
element, its ring's contact elements (the net's normals) and the curvatures
of its edge-composed tensor `S_i`; then the ambient Poisson coupling of the
field (19.1, one factor, matrix sign). Block Gauss-Seidel: (B) fit, couple,
read each cyclide's curvature spheres and principal direction at its vertex;
(A) the edge net's Gauss-Newton with two extra rows per cyclide vertex,
weighted by anisotropy: `S_i` toward the cyclide's tensor, and the ring's
contact elements in contact with the cyclide (third order; a vectorized
smooth residual, the contact line's distance from holding a curvature
sphere of each family). Vertices with anisotropy under 0.05 keep `S_i` and
stay out of the coupling; flat surfaces are out of scope.

Checks (`check_lie_cyclide_net.py`, all pass): the fit recovers
`lo.torus(2, 0.7)` from its contact elements and curvatures (1.7e-8, the
finite-difference curvature rows); the contact residual is 0 on the cyclide
(3e-16) and not off it; the coupling commutes with a global Lie
transformation including scale (6.6e-15) and keeps a constant field; the
joint Jacobian matches a directional difference.

**The ring contact rows destabilize.** With them (weight 1) the second
sweep's geometry block raises the merit five-fold on the tori (3.7e1 ->
2.1e2), strain 0.24, after which no cyclide fit is accepted; torus tube
6.9e-2, R 4.4e-2, band 0.30. Alone they do the same (merit 4e2 in sweep 2,
fixed point 1.4e-2). At the torus target they sit at 5.9e-2 rms and the
tensor rows at 8e-2: the cyclides composed from the net's state do not hold
the ring's contact elements, so the third-order rows pull against the edge
rows. Off by default (`w_ring` 0).

**Tensor rows only** (6 sweeps of 8 Gauss-Newton iterates):

| case | 19.14 edge net | cyclide net | fixed point (19.14 -> here) |
|---|---|---|---|
| torus tube | 1.1e-3 | 1.6e-3 | 1.1e-3 -> 8.5e-4 |
| torus R | 2.4e-3 | **1.1e-3** | 2.3e-3 -> 1.0e-3 |
| cube -> ellipsoid | 4.6e-3 | 4.7e-3 | 2.7e-3 -> 2.8e-3 |
| closed band | 1.4e-2 | 1.4e-2 | 1.4e-2 -> 1.6e-2 |
| cube -> sphere | 3.0e-3 | 3.2e-3 | 4.6e-3 -> 4.0e-3 |

Stable and nearly monotone (merit rises by under 1% in two sweeps on the
tori). It halves torus R's error and lowers both tori's fixed points; the
ellipsoid, band and sphere are unchanged. The pass bar (a lower floor on
every closed anisotropic case) is met on the tori only. On the ellipsoid
and band the accepted cyclide count flickers between sweeps (218 / 0 / 210 /
0 / 98 on the ellipsoid), so the fit is often rejected there (defect or
anisotropy test); those cases mostly run as the edge net. Cost: 12
factorizations over 6 sweeps on the tori, 11 s (the edge net 8 to 19 s);
20 to 55 on the cubes.

What is left: why the composed cyclides do not hold their rings (the fit
weights the vertex's own contact element 10 to 1 and the curvatures, so
the ring is a compromise, and the coupling then smooths it further); a
fit that is accepted consistently on the ellipsoid and band; then the ring
rows again, weak or ramped in once the field is settled.

### 19.16 Lie shape projection and cheap position routes (`lie_shape_projection.py`)

Reframing: keep the edge targets, project them each frame onto what Lie
elements express (central spheres everywhere, the conformal Gauss map whose
Dirichlet energy is the Willmore energy, Bryant 1984; implicit vertex
cyclides where anisotropic), couple both fields with one Poisson factor
(the fidelity trade `w`), then derive positions cheaply. Losses outside the
Willmore / Dupin class accepted if the surface's intention is kept.

**Implicit vertex cyclides.** A Dupin cyclide's U-family curvature spheres
span one 3-dim subspace (signature (2, 1)), so `U_i` = the top three
eigenvectors of the ring's sphere moments (minus the U^perp family's), in a
local frame, then `from_subspace`. Exact on exact data (5e-14 on the
analytic torus). Per-edge principal decompositions fail (each edge mixes
the contact and hinge read-outs; families swap on diagonals: an edge read
(-0.71, +0.17) where the torus has (-1.11, -0.3); cyclides 0.42 off even
after smoothing). From each vertex's ring-fitted tensor instead, with the
vertex's own spheres weighted 10 to 1: the read-back at the vertex is
faithful (curvatures within 0.004, directions within 5 degrees of the
tensor), but the third-order part is not (cyclide 0.4 to 0.5 off the
torus): one ring's spheres are nearly collinear, so the subspace's third
direction comes from tiny variations and amplifies the read-outs' ~10%
curvature bias.

**The projection is lossy beyond intention.** At the exact torus target
(representable: a torus is a cyclide) the projected edge targets differ
from the input by 15 to 48% (along / twist / across), with or without
cyclides, at any `w` from 0.1 to 10: the edge -> vertex tensor -> edge round
trip itself does it (the read-outs disagree edge to edge, the across
hinge values most).

**Positions** (closed cases, 8 projection sweeps, shape error; start in
brackets):

| route | torus tube (3.4e-2) | torus R (4.9e-2) | cube -> sphere (6.1e-2) | cube -> ellipsoid (5.2e-2) | band (1.2e-1) |
|---|---|---|---|---|---|
| contact net, raw targets | **2.0e-3** | **4.9e-3** | **1.0e-2** | **4.2e-3** | **1.7e-2** |
| contact net, central sphere | 1.6e-2 | 5.6e-2 | 2.4e-2 | 5.1e-2 | 4.6e-2 |
| contact net, cyclides | 1.5e-1 | 1.2e-1 | 4.8e-2 | 4.6e-2 | 1.2e-1 |
| spoke gather, raw | 6.6e-2 | 2.9e-2 | 3.2e-2 | 2.1e-1 | 1.3 |
| mean-curvature Poisson, raw | 1.1e-1 | 1.9e-2 | 3.9e-1 | 1.7e-1 | 2.4e-1 |
| envelope landing, raw | 6.2e-1 | 1.0 | 1.9 | 1.5e-1 | 2.7e-1 |

Only the contact net reconstructs, and the projection degrades it on every
case. The cheap routes (one prefactored factor, or no factor) are fast
(0.1 to 0.2 s for 8 sweeps against 1.3 to 2.7 s) but move away from the
target or diverge. The spoke gather holds an exact target (drift 1.7e-3)
but cannot reach one: with rest lengths reseeded every 10 iterates and a
Rodrigues blend, the tube does not grow and the band blows up (the
stretch-or-bend choice is local, as in 19.7). The mean-curvature Poisson
carries `H` only and drifts even at the target (2.2e-2: the edge read-out's
`H` is not the cotan Laplacian's). Envelope landing is Jacobi-local: even
relaxed and step-capped it does not settle (moves stay at 3e-2 to 6e-2).

So positions from curvature need the global stretch-or-bend choice the
contact net makes (positions and normals solved together). Not built yet:
the sphere-block variant and the Willmore / Dupin speed test (the
hypothesis that these objectives converge faster here than in the forward
solvers). Both would need a position route that works; the contact net is
the only one so far.

### 19.17 Lie group synchronization: global positions in one solve (`lie_sync.py`)

Each vertex holds a transformation `F_i` placing a canonical contact
element (origin, normal `e_z`, tangent `e_x` toward its first fan
neighbour); each edge prescribes `G_ij = F_i^-1 F_j`; one sparse least
squares `min sum |F_j - F_i G_ij|^2` synchronizes them, and positions and
normals are read off `F_i` (translations are group elements, so positions
come out of the same solve). Row form: each row of `[R_i | x_i]` obeys
`f_j = f_i G_ij`, so the three rows share one factor. Stage 1 is the
Euclidean subgroup (SE(3), 4 x 4; the 6 x 6 Lie form is only needed once
inversions enter). `G_ij` in the edge's midpoint frame: chord of length
`l`, endpoint frames turned by `-+ omega / 2` with `omega = tau l t - a l b`
(the Rodrigues relation as a rotation of the contact element, exact
magnitude), tangent axes at the spokes' polar angles (corner angles scaled
to `2 pi`). So `G` uses only the metric (lengths, angles) and each edge's
along sphere and twist; the across sphere enters through neighbouring
edges.

Two fixes on the way: the residual is homogeneous, so a weak gauge pull
lets every frame shrink (6e-2 even from the target); one frame is pinned
exactly. And the linear relaxation still shrinks the rotation blocks a
little, so the translations are re-solved with the projected rotations
(a Laplacian; two factorizations in all).

Checks (`check_lie_sync.py`, all pass): the model's relative rotations
match the target mesh's own frames (median 2.8e-5 on the cylinder, 8.6e-3
on the tori: the polar angles' `2 pi` scaling where Gauss curvature is
not zero); one sync from the start with the target metric reproduces the
target; a rigid motion of the start moves the result rigidly (1e-13).

One sync from the start shape, target edge spheres (shape error):

| case | start | target metric | start metric |
|---|---|---|---|
| isometric cylinder | 1.1e-1 | 3.6e-3 | 3.6e-3 |
| sphere cap | 9.8e-2 | 5.2e-3 | 9.0e-3 |
| torus tube | 3.4e-2 | **7.4e-4** | 3.3e-2 |
| torus R | 4.9e-2 | **6.5e-4** | 4.8e-2 |
| cube -> sphere | 6.1e-2 | 3.2e-3 | 1.2e-2 |
| cube -> ellipsoid | 5.2e-2 | 7.4e-3 | 3.8e-2 |
| closed band | 1.2e-1 | 4.0e-3 | 1.2e-1 |
| saddle | 1.6e-1 | 2.7e-3 | 4.6e-2 |

0.01 to 0.04 s in Python (the contact net: 1.3 to 2.7 s for its 8
sweeps), two factorizations, no iteration. With the metric given, a global
reconstruction from the contact-sphere network in one linear solve, close
to the contact net's accuracy (2e-3 to 1.7e-2) and better on the tori. With
the start's metric, isometric and nearly isometric changes work (cylinder,
cap, cube -> sphere at 1.2e-2); changes that need stretch (tori 15% length
change, band 17%) stay where they started: the curvature data carry no
metric, unlike the contact net, whose Rodrigues lengths choose stretch
globally.

**Conformal stretch, tried and off.** Scale factors `u` from the
linearized discrete-conformal equation (`Omega + L u = K_t A`, angle
defects growing with `u`, Springborn, Schroeder, Pinkall 2008), lengths
times `exp((u_i + u_j) / 2)`, re-synchronized: worse on every case,
including the cylinder that needs none (3.6e-3 -> 0.13). The target
Gauss curvature from the target spheres' vertex tensors is too noisy (a
product of two ~10%-biased curvatures), and the Poisson turns the noise
into stretch. `stage = "similarity"` kept, not the default.

Where this fits: in the demo the metric is not unknown, it is what the
material rows (area, axis ratio, shear, growth) prescribe in position
space, while the store carries curvature. A sync from the store's spheres
and the material rows' metric is one solve per frame. Next: the metric from
those rows' targets (or edge-length targets), then the Moebius stage
(relative transformations that are conformal; the basis for the Moebius
constraints).

### 19.18 The sync as a warm start for the contact net (`warm_start_sync.py`)

One sync (start metric, target edge spheres; positions and normals) as the
first iterate of `ContactNetGN`, against the start shape. Merit levels are
relative to the cold start's merit.

| case | sync start merit (of cold) | cold: its to 1e-3 / 1e-6 / 1e-9, factorizations | warm |
|---|---|---|---|
| isometric cylinder | 1.1e-2 | 8 / 11 / 14, 6 | **2 / 7 / 9, 3** |
| sphere cap | 2.8e-2 | 10 / 12 / 15, 6 | 7 / 14 / 15, 5 |
| cube -> sphere | 8.6e-2 | 5 / 10 / 10, 5 | 4 / 10 / 12, 4 |
| cube -> ellipsoid | 1.5e-1 | 5 / 9 / 10, 5 | 6 / 10 / 11, 4 |
| saddle | 7.8e-2 | 8 / 13 / 14, 6 | 13 / 24 / 25, 9 |
| torus tube / R | 0.93 / 0.92 | 3 to 5 / 6 / 7, 2 to 3 | 5 / 9 / 9, 3 |
| closed band | 0.97 | 6 / 8 / 10, 3 | 9 / 10 / 14, 13 |

Final shapes identical (machine precision either way). The sync lowers the
starting merit 10 to 90 times on the near-isometric cases but the
Gauss-Newton needs about as many iterates: it already converges in 7 to 15
from the start shape, and its last digits take the same few iterates from
anywhere. Only the cylinder gains clearly (factorizations 6 -> 3); the
saddle and band lose (the fading memory and the damping start from the
sync's state, and the band's sync, which cannot stretch, is not closer).
The sphere edit (curvatures times 1.2) is not reached either way.

So for the standalone contact net the warm start does not pay. Where it
might: in the demo, whose frames carry many competing rows and whose cost is
the joint factor's refactors; to be measured there before building it into
C++.

### 19.19 AMAP conformal rows and surface tension in the monge demo

**AMAP** (Vaxman, Mueller, Weber 2015, `reference/vaxman2015-conformal-mobius-deformations.pdf`,
section 6.2, the 3D edge-deviation form): per vertex a quaternion `Y_i` (not
unit: its modulus is the local scale), per edge the residual
`w_ik - conj(Y_i) q_ik Y_k` (`w` the live edge, `q` the frame-start edge,
pure quaternions; 4 rows) plus `alpha (Y_i - Y_k)` (the inversion damping).
`hepworth/block/amap_gn.hpp` (Gauss-Newton rows, analytic Jacobian, scaled
`w sqrt(A_e) / |q_ik|`); the quaternions are a `vec4_block` site on
`shell_cr_block` (`amap_site`, identity each frame). Demo
`position_distortion::w_amap`, `amap_inversion`; probe `wamap=`, `amapinv=`.
`triangle_conformal` (`w_conformal`) is a different notion: each face's
deformation gradient projected onto its closest similarity (first order,
per face, Euclidean); AMAP is vertex-based and Moebius-exact.

Tests (`amap_tests`, all pass): the frame start is a fixed point (0); the
gradient matches differences (2e-10); a mesh moved by a Moebius map (two
inversions) admits quaternions with AMAP energy 3.4e-29, an anisotropic
stretch 1.3 leaves 0.76 (the defining property).

**Surface tension.** Position space: per-triangle area rows with a zero
target (`area_mode::zero`, the area gradient being the mean-curvature
normal), demo `position_distortion::w_tension`, probe `wtension=`. Curvature
space: `goal_mean_tension` (`edge_shape_goals.hpp`), each edge's mean
curvature moved to `h0` (`constant`, Young-Laplace) or to the store's mean
(`mean`, constant-mean-curvature fairing), trace-free part kept; linear in
the store. Demo `shape_distortion::w_tension`, `tension`; probe `wcmc=`,
`cmch0=`. Tests: the goal sets the mean and keeps the trace-free part
(6e-17); zero-area rows shrink a torus (54.3 -> 18.8 in 3 steps).

Probe, torus target, 6 frames of 20 iterates, every knob explicit (baseline:
AMAP and tension off, `w_conformal` 1), frame 5:

| setting | R | tube | area | rough | s / frame |
|---|---|---|---|---|---|
| baseline | 0.774 | 0.295 | 9.00 | 0.182 | 3.9 |
| AMAP 1 | 0.775 | 0.295 | 8.99 | 0.183 | 8.0 |
| AMAP 10 | 0.807 | 0.276 | 8.72 | 0.182 | 7.8 |
| AMAP 1, `w_conformal` 0 | 0.773 | 0.296 | 9.00 | 0.181 | 7.3 |
| AMAP 10, inversion 1 | 0.807 | 0.276 | 8.72 | 0.177 | 7.2 |
| position tension 0.5 | 0.740 | 0.272 | 7.93 | 0.168 | |
| curvature tension 10 (mean) | 0.773 | 0.296 | 9.02 | 0.160 | |
| curvature tension 10, `h0` -1.5 | 0.770 | 0.294 | 8.92 | 0.170 | |

AMAP at 1 barely moves the shape and doubles the frame time (4 more
unknowns per vertex in the joint factor). At 10 it resists the tube's
growth (0.276 against 0.295) and lets R move further (0.807): a torus with
a larger tube is not conformally equivalent to the start, so a conformal
row holds the tube and gives way elsewhere. The inversion damping barely
matters: per frame, neighbouring quaternions stay nearly equal (each frame's
change is locally almost a similarity). Position tension shrinks area and
the tube, as surface tension should. Curvature tension evens the
curvature (roughness 0.182 -> 0.160 in mean mode) at almost no change of
size. The tension-only runs' times were not consistent (machine load), so
they are left out.

A harness note, again: under zsh an unquoted `$v` holding two flags
passes them as one argument, so only the first takes effect. The first
probe matrix was wrong for that reason; use `${=v}`.

### 19.20 Position-space Dupin circle rows (`dupin_circle.hpp`)

Dupin cyclides are the surfaces whose curvature lines in both families are
circles (`reference/2312.04341v1.pdf`, Bobenko and Fairley, section 4). Any
three points lie on a circle, so a row needs four; four points concyclic is
Moebius-invariant (real cross ratio, Vaxman et al. 2015). Groups: for an
edge within 25 degrees of one family's principal direction at both ends
(the store's vertex charts), the neighbours of its endpoints that best
continue the line, `(a, i, j, b)`; weight alignment times anisotropy.
Row `circle4` (projective, 12 rows: the four positions minus their
centroid): each iterate fits the live circle (least-squares plane, the
algebraic circle fit in it, a line when nearly straight) and a mode picks
the targets on it. Drive modes seed their targets once from the frame-start
circle and each iterate takes the seeds' closest points on the live circle
(no drift within a frame; the refit keeps the line circular); modes with a
direction along the line take it from the shear axis (oriented) or a
reference vector.

Modes: `closest` (the plain row); `slide` (rotation of the circle, `delta`
edge lengths of arc per frame); `even` (ends kept, inner points at thirds);
`graded` (arcs `1 : g : g^2`); `hyperbolic` and `parabolic` (Moebius maps
of the circle fixing the group's middle and its antipode, or the antipode
alone; cross ratio kept); `loxodromic` (the circle moved on its curvature
sphere by the map fixing `x_i`: turns the line about `x_i` by `psi`).

Tests (`amap_tests`, all pass): points on a circle fixed (5e-16),
perturbed points land on a circle (3e-16), near-collinear groups fit a
line; slide, even, graded angles exact (9e-16); hyperbolic and parabolic
keep the cross ratio (6e-16) while moving the points (0.61); loxodromic
keeps `x_i` fixed, stays on the sphere and on one circle, turns the
tangent by exactly `psi`; no drift (two local steps agree exactly). The
canonical torus (grid on curvature lines): 1024 groups, residual 5.5e-16 (a
Dupin cyclide); a noisy torus's residual falls 1.73e-2 -> 1.70e-3 under the
rows.

Probe, torus target, 6 frames of 20 iterates (shear on, so the drive modes
take its axis), frame 5; `dcres` the mean circle residual over the mean
edge length (frame 0: about 1e-4):

| setting | R | tube | area | rough | dcres | s / frame |
|---|---|---|---|---|---|---|
| off | 0.772 | 0.296 | 9.00 | 0.179 | 3.29e-3 | 3.8 |
| closest 1 | 0.772 | 0.296 | 8.98 | 0.179 | 3.22e-3 | 4.9 |
| closest 10 | 0.774 | 0.293 | 8.92 | 0.188 | 2.77e-3 | 4.8 |
| slide 0.2 | 0.780 | 0.288 | 8.81 | 0.181 | 2.31e-3 | 4.8 |
| even | 0.780 | 0.287 | 8.78 | 0.179 | 2.21e-3 | 4.8 |
| graded 1.5 | 0.774 | 0.285 | 8.67 | 0.178 | 2.13e-3 | 4.7 |
| hyperbolic 1.3 | **1.126** | 0.222 | 9.88 | 0.180 | 3.97e-3 | 5.0 |
| parabolic 0.2 | 0.738 | 0.287 | 8.29 | 0.279 | 2.51e-3 | 4.8 |
| loxodromic 0.1 | 0.856 | 0.299 | 10.06 | 0.165 | 2.02e-3 | 5.6 |

The plain row is weak against the demo's other rows (residual down 2% at
1, 16% at 10) and costs about a second per frame. The drive modes do
reshape: even, graded and slide lower the circle residual a third and pull
the tube in; loxodromic turns the lines and grows the area (10.06) while
smoothing (roughness 0.165) and keeps the lines most circular; hyperbolic
regrading at 1.3 per frame is strong (R 0.77 -> 1.13, tube 0.30 -> 0.22,
anisotropy up); parabolic roughens (0.279). The amounts are per frame and
uncalibrated; these are first readings of each mode's character, not tuned
settings.

**Fix: the row keeps each group's size.** The plain row shrank the surface.
The circle condition is scale-free, so the row has no restoring force on
size, and the algebraic circle fit underestimates the radius on short arcs:
each iterate pulled the four points slightly inward and the bias
accumulated (circle rows alone on a noisy torus: volume -0.53% in 20
frames, against -0.15% from removing the noise). The targets are now
rescaled about their centroid to the live points' RMS spread, so the row
changes only the group's shape (-0.37% in the same test). In the demo
(later defaults, so compare within this pair), frame 5 with the plain row
at 10 against off: area 8.90 against 8.80 (no shrink), circle residual
4.8e-4 against 3.9e-3 (eight times rounder lines; before the fix the row
fought the others and moved it 16%), roughness 0.215 against 0.243.

**Why it shrinks, isolated** (`dupin_circle_study.py`: a replica of
`circle4` on a torus whose grid lines are curvature lines, every group at
weight 1, projective dynamics with a weak mass). The rows constrain only
the group's shape (positions minus centroid), so a uniform offset of the
fitted circle is invisible to them. Alone they do not shrink: the clean
torus is an exact fixed point, and from normal or tangential noise of 0.1
edge the volume ends at or slightly above its start (+0.2% from the
spread rescale counting noise as size); the algebraic and geometric
circle fits agree. One failure mode exists: when the four points zig-zag
sideways by more than the arc's sagitta (fine meshes, groups following
edges rather than lines), the least-squares plane follows the zig-zag and
the row flattens the line: 32 x 16 at 0.25 edge, volume -17% and area
+85%; 64 x 32 at 0.4 edge, -24% and +150%. The plane of the chord and the
sagitta (inner pair's mean off the ends' mean), in which an alternating
zig-zag cancels, removes it (within 0.2%, and 3.5% / 11% in the worst
case) and matches the least-squares plane under random noise; the plane
through the three chord midpoints also cancels it but collapses under
random noise (volume 0.47 in 20 frames). In C++ the zig-zag does not
reproduce (test `dupin_circle_rows_keep_volume_under_zigzag`: 1.0008 either
way): the 25 degree alignment test drops the groups where the zig-zag
dwarfs the sagitta, and the sagitta plane was slightly worse under random
noise (0.9941 against 0.9963), so the C++ fit is unchanged. The demo's
shrink is therefore not the row alone; the interaction with the other
rows (coupling, memory, remeshing) is the next suspect, by probe ablation.

**Ratio and shear as settings** (same study). Ratio: each family's
targets scaled from the frame-start spread, the tight (tube) family by
`sqrt(rho)` and the loose by `1 / sqrt(rho)` per frame. Well posed and it
sets size, so it cannot shrink; the torus stays a torus while its aspect
changes (`rho` 1.02 per frame: `R / r` 2.86 -> 2.45 in 10 frames, about 80%
of the drive; 1.05: 1.96). Scaling from the live spread instead compounds
over iterates and diverges. Shear: one family's group slid along the other
family's direction by `sigma` times arc position turns that family's lines
by only 0.2 to 0.35 degrees (`sigma` 0.05 to 0.1): the slid points leave the
circle, the next fit undoes them, and curvature lines must stay orthogonal,
so a differential turn is not representable. Turning both families
together is the loxodromic mode, which already exists.

**Variational rows and the Moebius shift** (`dupin_mobius_shift.py`). The
current drive modes re-project frame-start seeds onto the live circle and
rescale to the live spread, so their targets follow the iterate and they
are not variational. Variational form: at the frame start each group's
target shape `T` is its points on their (geometric) circle mapped by the
mode, sized to the frame-start group; each iterate the live points are
matched to rigid copies of `T` (Kabsch). The local-global energy never
rose in any run. The shift: rays from the circle's centre reprojected onto
a circle shifted by `delta` keep the points on a circle but change their
cross ratio (second order, about `0.4 |delta / r|^2`), so they are not
Moebius; the conformal form is the unit-disk automorphism
`z -> (z - a) / (1 - conj(a) z)`, `a = delta / (2 r)`, which agrees with
the reprojection to first order (angle gap `0.25 |delta / r|^2`), keeps
cross ratios exactly and composes as a group. On a square cell, `a` on a
diagonal gives a right kite (sides 0.78 / 1.18 / 1.18 / 0.78, 90.00 degrees
at the side corners), along an edge direction an isosceles trapezoid.

Size: a Moebius map also scales the arc between its fixed points (by
`(1 - |a|^2) / (1 + |a|^2)` midway), so the image's own size shrank the
torus (line groups, amount 0.1: volume -12%); sized to the frame-start
spread, line groups move volume only -0.06 to -0.24% (they regrade points
along their circles, tangentially). Cells: one kite asked of every cell
is not realizable (identical kites do not tile a grid by translation): the
angles did not change (3.6 degrees off square throughout) and volume fell
(-6% at 0.2). Alternating the axis in a checkerboard (any quad tiles by
half turns about its edge midpoints) realizes them: angles 6 to 23 degrees
off square at amounts 0.02 to 0.1. Sized by spread that still lost volume
(-0.5% to -7.6%: a kite of equal spread encloses less area); sized by the
cell's area it is volume neutral on the diagonal axis (+0.01%, +0.3%,
+2.6%), while the tube axis (alternating trapezoids) changes the torus's
proportions (`R` 2.00 -> 2.13, `r` 0.70 -> 0.64 at 0.1). Lessons for a C++
row: targets fixed per frame and Procrustes-matched (variational); size
from the frame start, by area for cells; per-element kites only realize as
a compatible pattern (a checkerboard, or an integrable Moebius field), so
the axis field has to be chosen with that in mind.

**Kites from one integrable field** (`dupin_kite_field.py`; torus with a
bump, the direction from the vertex store's principal axis sign-aligned
frame to frame and turned 45 degrees, standing in for the Dirichlet edge
tangents; `g = alpha |z| t`, one Poisson solve `L u = div g`, each cell's
shift from its `grad u`; kite measured against the same run at `alpha = 0`,
6 frames):

| mode | realized | correlation | roughness | volume |
|---|---|---|---|---|
| checkerboard (alpha 0.5) | 84 to 92% (7 to 13 degrees) | 0.94 to 0.97 | 1.8 to 1.9 | +0.2 / 0.0% |
| integrable, area held | 1 to 2% | 0.3 to 0.5 | 0.5 | 0 |
| per-cell, no Poisson, area held | 1 to 2% | 0.2 to 0.35 | 0.7 | 0 |
| integrable with its scale `e^{2u}` (alpha 0.02) | 29 to 50% | 0.3 to 0.6 | 0.04 to 0.11 | +0.6 to +1.5% |

A Moebius field is a conformal growth: its kites are compatible only with
the cell sizes `e^{2u}` that come with them (held area removes exactly the
part that makes neighbours agree). Smooth kites are therefore small: a
cell's kite is about `|grad u| h`, and any visible angle compounds into
exponential growth across the surface (alpha 0.5 with the scale diverged).
Visible kites (degrees per cell) exist only as a mesh-scale alternating
texture, which the checkerboard realizes at about 90%, volume neutral.
Two features, one direction field: geometry (integrable conformal growth,
`u` into edge or area growth targets, kites as its trace) and texture
(checkerboard kites oriented by the field).

**C++: `kite_growth.hpp`** (demo `position_distortion::kite`, `w_kite`; probe
`wkite= kalpha= kangle= kmaxu=`, readout `kiteu`; tests in `amap_tests`: the
Poisson solve recovers a linear function to 1e-7 and drops a divergence-free
ring field to 3e-6, zero rate is a fixed point, the area rows follow `2 u`
at correlation 0.98). The Dirichlet edge tangents (sign-aligned to the
carried datum, turned by `angle`) times `alpha |z|`, one Poisson solve,
per-triangle scale-only `area` rows at `e^{2u}` (rotation and conformal
class free). Probe, torus target, frame 5: `u` spans only 0.002 to 0.004 per
frame at `alpha` 0.02 to 0.05, whatever the angle: on a torus the
curvature-line field circulates around closed loops, and its loop-integrated
part has no potential (a scale gradient that never closes up), so the
Poisson solve keeps only the small variation along the loops. With `u` near
zero the rows act as per-triangle area locks at weight 1 (area 9.98 -> 9.86
against off): the Moebius growth shares the scale channel with the demo's
growth, as predicted.

Rows replaced by the edge form: per edge an `edge_strain` row at
`e^{(u_i + u_j) / 2}` times its frame-start length (discrete conformal
scaling: a triangle whose vertices carry different `u` changes shape, the
kite at the triangle level), and an optional direction `source`: the
tangents, or the gradient of the store's `H` or `|z|` (`mean_gradient`,
`aniso_gradient`), which have sources and sinks. Tests: edge lengths follow
`(u_i + u_j) / 2` (correlation 0.99, gain 0.90, total area +1%); on an exact
torus the ring tangents give a `u` range of 3e-7, the `|z|` and `H`
gradients 0.079. Probe (the demo's later defaults, `alpha` 0.2, `w_kite`
1), frame 5:

| setting | R | tube | area | roughness | `u` range per frame |
|---|---|---|---|---|---|
| off | 0.816 | 0.217 | 6.85 | 0.68 | 0 |
| tangent | 0.759 | 0.245 | 7.29 | 0.22 | 0.009 |
| `|z|` gradient | 0.764 | 0.234 | 7.10 | 0.51 | 0.17 |
| `H` gradient | 0.776 | 0.245 | 7.31 | 0.54 | 0.16 |
| `|z|` gradient, alpha 0.1 | 0.767 | 0.242 | 7.30 | 0.37 | 0.08 |

The gradient sources drive (`u` near the 0.1 cap each way); the tangents do
not. At weight 1 the edge rows also hold every edge's frame-start length
up to `e^u`, a metric memory: with tangents (no drive) they only stabilize
the frame (roughness 0.22 against 0.68 off, which roughens in this
configuration), and with a drive they raise roughness to 0.37 to 0.54.
(These gradient runs had the 45 degree turn applied; a gradient turned by
`theta` is `cos theta` gradient plus `sin theta` circulation, and the
Poisson solve keeps only the first, so the gradient sources now ignore
`angle` and set their own direction.)

**Scalar sources replace the gradient sources.** The gradient sources
shrank the torus: the Poisson solve integrates `alpha grad s` back to
`alpha s` plus a constant, and the constant was set so `u` had zero mean,
so growth where `H` or `|z|` peaks (the small inner side) was paid for by
shrinking the large outer side (`R` and volume down). The sources now set
`u` directly: `mean_curvature`, `u = alpha (-H - level)` (each element's
Moebius map anchored on its mean-curvature sphere: convex grows, concave
shrinks, saddles stay; outward normals make convex `H` negative), and
`anisotropy`, `u = alpha (|z| - level)`; `balanced` restores the zero mean.
The tangent source keeps the Poisson route (always mean zero). Tests: on
the exact torus `u` is 0.045 on the outer equator (`0.05 x 0.90`) and 0.017
on the inner, positive everywhere; the edge rows alone inflate the volume
12% in 3 frames at `alpha` 0.02. Probe (demo defaults, `alpha` 0.1,
`w_kite` 0.05), frame 5: off R 0.709, tube 0.391, area 11.03, roughness
0.66; mean R 0.732, tube 0.400, area 11.52, roughness 0.48; mean balanced
R 0.735, tube 0.382, area 11.08 (no net growth), roughness 0.49. At
`alpha` 0.1 on this torus (`-H` about 2.6) `u` sits at the `max_u` cap
over most of the surface, so the growth is nearly uniform there (the kites
come from the variation of `u`; keep `alpha |H|` below the cap to keep it).

### 19.21 Complex curvature store and curvature advection (`hz_store.py`, `curvature_advection.py`)

**Complex form.** In an orthonormal tangent frame `S = H I + [[a, b], [b, -a]]`;
with `z = a + i b = (k11 - k22) / 2 + i k12` and `w` a tangent vector as a
complex number, `S w = H w + z conj(w)`. Turning the frame by `phi` sends
`z -> z e^{-2 i phi}` (spin 2), `H` is a scalar. A vertex store `(H, z)` (3
reals per vertex against 3 per edge: 864 against 2592 unknowns on the
48 x 24 torus) maps to edges by `B`: endpoint mean, `z` turned into the edge
frame by `e^{-2 i phi}`; `P` is the least-squares projection of edge
tensors onto the vertex store. Checks (`hz_store.py`): the identity 1e-16,
frame rotation 2.5e-16, frame invariance 5e-16. Representational loss
`|h - B P h| / |h|` on the suite targets: 0.002 to 0.087 on the CR read-out,
0.13 to 0.65 on the contact read-out (its edge values carry twist and
along-edge detail the vertex store cannot hold).

**Weingarten solver** (`hz_weingarten.py`, `shape_site_e2e` protocol,
targets `h` against `B P h`):

| case | edge clean | vertex clean | edge noisy | vertex noisy | best CR-faired noisy |
|---|---|---|---|---|---|
| torus tube 0.7 -> 0.9 | 0 | 1.9e-3 | 1.6e-3 | 2.3e-3 | 2.6e-3 (w 0.03) |
| cube -> sphere | 0 | 7.8e-3 | 2.28e-2 | 1.27e-2 | 4.2e-3 (w 1.0) |

The vertex store is a projection with a fixed, mild smoothing: it costs
0.2 to 0.8% shape on clean targets and denoises about like a light spin-2
CR fairing, but no better than a tuned one, and the tuning that wins on
the sphere (w 1.0) ruins the torus (8e-2).

**Contact net** (`hz_contact_net.py`). Projected targets into `ContactNetGN`:
the raw edge read-out recovers every case to 1e-16 in 9 to 18 iterates;
`B P`-projected targets land at 6e-3 to 9e-2 in 14 to 44 iterates with 3 to
10 times the factorizations (the targets are no longer realizable). Joint
solve with the store as unknowns: the vertex store has 35 to 40% fewer
unknowns but more fill (each vertex value enters every edge row of its star
through the phase factors: L+U 1.09e6 against 6.8e5 on the torus, 8.9e5
against 4.6e5 on cube -> sphere), 4 to 7 times the factor time, twice the
iterates, and a 2e-3 to 5e-2 shape floor. **Decision:** the complex form is
a notation (and the right variable for transport), not a cheaper store;
keep the edge store in the solvers.

**Complex on the edges** (`hz_edge_solvers.py`). The store kept on edges,
3 reals each, as `(k11, k12, k22)` or `(H, Re z, Im z)` (`k = C s`), in both
joint solvers (Weingarten: positions + store, read-out rows `R(x) - C s`;
contact net: positions + normals + store), times the store rows' metric:
component `|dk|^2 = 2 dH^2 + 2 dRe(z)^2 + dIm(z)^2` (the C++ default: the
deviator along the edge counts twice the diagonal one, so it depends on
how the edge sits against the axes) or invariant `2 dH^2 + 2 |dz|^2`
(`tensor_metric`). Results on 2 Weingarten and 8 contact cases, clean and
with 10% component noise:

- Parametrization: no change. Clean targets reach 1e-16 in the same
  iterates in every variant; noisy, the shape errors agree to three digits
  in all 10 case-solver pairs (positions to 1e-6 to 1e-4 on the contact
  net; up to 5e-3 on Weingarten, a slide along a near-free direction under
  differently scaled damping; cube -> sphere noisy hits the 40-iterate cap
  either way). Iterates and factorizations differ by a few because the
  Marquardt damping `diag(J^T J)` is not invariant. Factor fill moves with
  the ordering: complex 15 to 30% more on the contact net (`k = H + Re z`
  couples two unknowns per row), 7 to 15% less on Weingarten; factor times
  follow.
- Metric: the invariant metric changes the noisy answer by 1 to 18%
  (here slightly worse, but the noise is independent per component, which
  is the component metric's own noise model).

So acting on complex numbers directly changes nothing in the solve; only
the store rows' metric is a real choice, and it is available today as
`tensor_metric`.

**Quaternion frames as the store** (`quat_store.py`, 19.22 in spirit). Each
vertex carries a quaternion `Y` (frame and normal `n = Y e3 conj(Y)`; scale
`|Y|^2` unused here); across an edge the normals give the edge's bending
and twist by the Rodrigues relation, `k11 = -dn . t / |dx_t|`,
`k12 = -dn . b / |dx_t|` (three edge directions per vertex fix the full
tensor; the rotation about the normal is a gauge). A frame field holds only
curvature some frame field realizes: no slack between store and mesh.

- Contact net (store: positions, normals, edge store, goal rows on the
  store; frames: positions and `Y`, goals on the Rodrigues and tangency
  rows directly): clean targets exact on all 8 suite cases in about the
  same iterates with half the unknowns (e.g. 2016 against 4320 on the
  torus); fill lower on the open sheets and the band (band 204k against
  413k), 24 to 32% higher on the tori and cubes; noisy shapes within about
  10% of the store's on 7 cases (torus R stalled after one iterate, 1.1e-2
  against 4.8e-3). The edge store is redundant next to the normals.
- Weingarten (frames coupled to the CR read-out's `(k11, k12)`): a floor
  even on clean targets (7.6e-3 on the torus tube; 1.6e-3 with the
  frame-start offset between the two read-outs held, the
  `contact_net_coupling` increment form), and failure on the open sheets
  from flat (shape errors of order 1). With frames, the normals and the
  tangency rows already couple positions and curvature (the contact net);
  the CR read-out is a second discretization that fights it.
- Frustrated targets (a uniform twist added to the torus target): both
  shrink (the contact rows' collapse mode once the length memory fades); at
  twist 0.8 the store crumples harder (mean edge 0.18 against 0.26, realized
  twist roughness 1.0 against 0.18); with the memory held at 0.01 the frames
  solve stalls after 2 to 3 iterates, at 0.1 both diverge. No clean
  curvature-space buckling yet: it needs a firm metric anchor (hard edge
  lengths or areas) so frustration goes into bending rather than shrinking.

**Per-edge rotations as the store (E2)** (`darboux_store.py`). Each edge
carries `q_e`, the rotation taking vertex `i`'s frame to vertex `j`'s. Exact
read-out: bending and twist from its action on the normal (the Rodrigues
chord form on `n_i` and `q_e n_i conj(q_e)`, gauge independent; matches the
contact read-out to 1e-15); `kappa_g l` the remaining turn about the image
normal relative to the shortest arc (intrinsic, gauge-fixed by the vertex
frames' tangent field). Contact net with positions, normals and `q_e`:
rotation rows `n_j - q_e n_i conj(q_e)`, goals on `(k, tau)`, a weak gauge
row holding `kappa_g` at its frame-start value, tangency, unit rows, length
memory. Holonomy (the rotations composed around a face) is the identity
for rotations that come from frames, in any gauge; its off-normal part
measures bending/twist frustration, its turn about the normal intrinsic
frustration.

- Suite: clean exact on all 8 cases; noisy best or tied on 7 (torus tube
  3.26e-3 against store 3.85e-3 and E3 3.52e-3; torus R 3.67e-3 against
  4.56e-3 / 4.29e-3; cube -> sphere 5.23e-3 against 5.81e-3 / 5.64e-3; the
  ellipsoid 0.98e-3 against the store's 0.91e-3). Unknowns about 20% above
  the store's (5184 against 4320 on the torus), fill 1.5 to 2 times, time
  1.5 to 2 times. Holonomy off the normal: 1e-16 clean, 0.01 to 0.06 noisy;
  about the normal 0.002 to 0.17 even clean (the `kappa_g` held at its
  start values while the target changes the intrinsic curvature).
- Frustrated twist, length memory held at 0.05: at twist 0.6 the edge store
  blows up (realized twist -286, roughness 572, mean edge 0.60x), E3
  crumples (roughness 1.74, mean edge 0.68x), E2 stays smooth (roughness
  0.30, mean edge 0.97x) and takes up the twist in the store (its twist
  +0.37 against +0.41 asked) as holonomy (off-normal 0.21). The frustration
  lives in curvature space instead of wrecking the positions.

**Decision:** E2 for the C++ contact net store: better noisy accuracy, the
only form robust under frustration, holonomy as a native frustration
measure, and an intrinsic `kappa_g` channel tied to the tangent field; at
the cost of 1.5 to 2 times the factor time. **Revised 2026-11-12:** E2 was
slower than hoped and its twist stayed latent in the store; the C++ E2 mode
is removed and the contact net moves to E3 with native rows (below).

**E3 with native rows** (`e3_native.py`, 2026-11-12). Unknowns `x` and `Y`
per vertex, no edge store. Coupling: per directed spoke the frames'
Rodrigues normal curvature `k_vj` and the tangency
`n_v . d / |d| = k_vj |d| / 2 + c0` (`c0` the frame start's offset), so the
positions' curvature follows the frames'; a unit row on `|Y|` (or the kite
rows), a weak gauge row on the turn of `T` about `N`. The vertex tensor
`S_v` in `v`'s own frame is the least-squares fit of `k = t^T S t` and
`tau = t^T S (N x t)` over the spokes (within 3% of the analytic torus on
24 x 12). Rows that compare vertices read `S` along frozen world directions
projected into each tangent plane, or carry `z` by the frames' relative
turn (`T_j` carried onto `i`'s plane by the shortest arc, `z_j e^{2 i phi}`),
so the gauge drops out. Jacobians by colored central differences (the row
pattern probed once per mesh). Native rows: shape memory (AMAP-style:
`z_v` held at its frame-start value turned by the ring's best common turn,
`H_v` held), tension (constant, frame_mean, volume controller), minimizer
(across-line Dirichlet of `k1`, `k2` on frozen principal frames), Dupin
(along-line), growth (normal curvature along / across a direction), `det`
(ring sum of face solid angles of the normals), `gradtl` (Dirichlet on `z`
with the frames' turn), twist (geodesic torsion along a winding direction,
frame start plus `tau0` per frame), kite (`|Y_v|^2 = e^{u_v}` and spin-scale
edges `|d| = l0 |Y_i| |Y_j|`), position growth (edge lengths
`l0 (1 + g (d0 . a)^2)`).

Torus `R = 1, r = 0.35`, 24 x 12, 3 frames, coupling 10, length memory
0.05 (tube 0.350, volume 2.28, roughness 0.114 at the start):

| rows | tube | volume | roughness | note |
|---|---|---|---|---|
| none | 0.350 | 2.28 | 0.114 | exact fixed point (no iterate taken) |
| growth (x 0.3 at 0.6 rad) | 0.250 | 1.07 | 0.178 | more curvature along the axis: the tube tightens |
| growth + memory 1 | 0.281 | 1.37 | 0.156 | memory resists |
| growth + det 100 | 0.227 | 0.94 | 0.127 | `det` cannot resist a uniform tube shrink: the normals, so the Gauss image, do not change |
| tension constant `H = -1.5` | 0.297 | 1.66 | 0.004 | lands on the target, very smooth |
| tension frame_mean | 0.318 | 1.91 | 0.006 | CMC fairing: smooth, loses 17% volume |
| tension volume | 0.340 | 2.18 | 0.010 | the controller recovers the volume |
| minimizer 1 | 0.309 | 1.91 | 0.093 | smooths; the torus's hoop curvature varies across its lines, so not a fixed point |
| dupin 1 | 0.394 | 2.79 | 0.120 | drifts 4% a frame from the vertex fit's discretization (0.4% at weight 0.1) |
| gradtl 1 / 0.05 | 0.462 / 0.369 | 4.03 / 2.53 | 0.352 / 0.123 | `z` varies around the tube; strong weights roughen |
| kite (alpha 0.1, cap 0.1) | 0.437 | 4.28 | 0.081 | inflates smoothly; `|Y|` 1.10 to 1.12 |
| twist `tau0 = 0.1` (absolute) | 0.493 | 4.28 | 1.58 | the start's torsion is about -1: a large frustrated jump, crumples |
| twist `tau0 = 0.1` (increment), 4 frames | 0.294 | 1.64 | 0.102 | smooth |
| twist + position growth 0.1, 4 frames | 0.398 | 3.18 | 0.102 | smooth, grows |
| twist + kite, 4 frames | 0.434 | 4.41 | 0.135 | lines start to turn at frame 3 (turn 0.18, anisotropy 0.53 -> 0.61) |

Frames take 1 to 9 s in Python (finite-difference Jacobians, 8 LM
iterates).

**C++ (2026-11-13).** Stage A: `contact_net_gn` takes per-vertex
quaternions as its second block (`shell_cr_block::frame_site`,
`curvature_config::contact_frames`, probe `cnframes`), `n = M(Y) / |Y|^2`
with analytic `dn/dY`, a unit row on `|Y|` and a gauge row. The gauge first
held the turn of `T` carried by the shortest arc from the frame-start normal;
that stalled cube -> sphere (shape error 1.0 after 355 iterates), the arc's
holonomy fighting the frames under large turns of the normal. The twist of
`conj(Y0) Y` about the body normal (its `z` component, linear in `Y`) at
weight 1e-3 lands every closed suite case (round-off on the tori and the
ellipsoid, 3e-10 on cube -> sphere in 126 iterates against 77); a strong
gauge (1) still stalls cube -> sphere. On the probe the frames track the
normals (tube 0.149 against 0.147 after three frames) at 2.3 times the frame
time.

Native rows (`frame_rows.hpp`): the rows of the prototype on `(x, Y)`, with
the vertex tensors as unknowns of their own (`vertex_shape_dof`, seeded from
the frame-start fit) tied to the frames by per-spoke fit rows. Chaining the
rows through the fitted tensors instead (the prototype's form) gave the same
iterates but 13 s a frame against 7: a row between two vertices then couples
both one-rings and the factorization dominates (3700 of 5000 samples). Probe,
three frames, `axis_ratio` off (idle: exact fixed point):

| rows | tube | roughness | note |
|---|---|---|---|
| none | 0.248 | 0.163 | fixed point |
| minimizer 0.5 | 0.251 | 0.162 | smooths slightly |
| tension 0.5 frame_mean / volume | 0.258 / 0.256 | 0.172 | the controller slows the drift |
| dupin 0.5 | 0.255 | 0.170 | the vertex fit's slow drift |
| growth 0.3 (w 0.5) | 0.192 | 0.215 | tightens the tube |
| + det 2 | 0.192 | 0.182 | `det` evens it |
| + shape memory 1 | 0.231 | 0.190 | memory resists |
| gradtl 0.1 | 0.254 | 0.170 | |
| twist 0.5, `tau0` 0.05 | 0.249 | 0.167 | |
| + stretch 0.5 | 0.318 | 0.176 | lines turn at frame 2 (anisotropy 0.70 -> 0.66) |

Next: drop the edge store from frames mode (store-free tangency as in the
prototype), the medial shear.

**Advection.** Semi-Lagrangian on the vertex `(H, z)` field: RK2 back trace
walking across faces (straightest path), barycentric sampling of the
corners' ambient tensors `[T B] S [T B]^T` (no frame angles), read back in
the destination frame. Transport `pullback` (Nabizadeh et al. 2022,
`reference/nabizadeh2022-covector-fluids.pdf`): `A <- J^T A J`,
`J = exp(-dt G)` with `G` the one-ring least-squares velocity gradient,
its normal column completed skew (`n . G e = -e . G n`), exact for rigid
motions. Finite differences of traces were tried first and blew up (1e67
in 24 steps): over a 1e-3 h offset, faceted path differences of 1e-2 give
Jacobian columns of norm 24. Schemes `first` (Stam), `bfecc` (Kim et al.
2005, `kim2005-flowfixer-bfecc.pdf`), `maccormack` with Selle's clamp on
`H` and a modulus clamp on `z` (frame invariant; `selle2008-...pdf`).

`check_curvature_advection.py` (all pass). Sphere, rigid rotation, one
full turn, Gaussian bump in `H` and `z`:

| scheme (pullback) | cube 6 (218 v, 24 steps) | cube 12 (866 v, 48 steps) | order |
|---|---|---|---|
| first | 0.582 (peak 0.47) | 0.437 (0.52) | 0.41 |
| bfecc | 0.257 (0.71) | **0.087 (0.93)** | **1.56** |
| maccormack | 0.570 (0.55) | 0.298 (0.79) | 0.94 |
| maccormack, no limiter | 0.574 | 0.306 | 0.91 |
| bfecc, parallel | 1.014 | 1.193 | - |

Parallel transport is wrong for a rotation: carrying the frame along a
latitude is not the flow's pushforward, so `z`'s phase drifts (error near
1 at any resolution). MacCormack's first order is not the limiter (free
and limited agree); on this vertex sampling BFECC's extra forward step is
what reaches second order. Torus loxodrome: peak offsets 0.28 / 0.44 /
0.87 edges, height kept 0.54 / 1.22 / 1.18 (BFECC and MacCormack
overshoot the unlimited `z` growth from the pullback). Flat shear:
pullback matches the exact `J^T A J` (0.000; undeformed 0.202). A sharp
step: MacCormack stays in [0, 1], BFECC overshoots to [-0.21, 1.20].

**Self-advection** (`burgers_curvature.py`). A perturbation `q = A sin(u)`
of `H` rides the torus ring direction at speed `c q`: Burgers on every ring,
shock at `t* = rho / (c A)`, so the inner equator (rho 1.3) creases before
the outer (2.7) from a uniform start. Scalar (parallel transport, 1024
vertices, 24 steps): every scheme follows the exact steepening until the
gradient reaches the grid (inner 1.08 / 1.80 / 4.1 against exact 1.11 /
2.00 / 10 at t 0.4 / 0.8 / 1.2), then saturates at the shock (first 7, BFECC
8, MacCormack 7.3, steady); the outer ring trails as predicted. BFECC
overshoots past the shock (q to 1.125), MacCormack keeps the range
(0.998). Tensor (pullback): `J^T A J` amplifies `H` where the flow
compresses, a faster wave compresses more, and it blows up (`|q| > 20 A`)
at a fixed 1.62 `t*` for `A` 0.25, 0.5, 1 (self-similar), generating
anisotropy `|z|` up to 3.4 A from an isotropic start. That is the curvature
instability asked for, but unbounded: a usable driver needs a cap (the
limiter on `|z|` and `H`, or a viscosity `nu L` on the store).

**As geometry** (`advect_geometry.py`). Per frame: advect the vertex field
on the live mesh, map the change by `B`, and give `ContactNetGN` increment
targets on the mesh's own read-out (`k += dk11`, `tau += dk12`, tangency
`c += |d| / 2 * d kappa_n`, `kappa_n = H + Re(z conj(w)^2)`), so the
read-out's projection loss never enters; the field is state, carried to
the next frame's vertex frames through its ambient form. Seeded from the
CR read-out: the contact read-out's hinge column spikes (150 on a torus
whose curvatures are near 1) where a bump makes a dual edge degenerate,
and the fit reproduced the spikes. Torus 48 x 24, normal bump 0.12 on the
outer equator, BFECC, 30 Gauss-Newton iterates, about 4 s per frame:

| scenario | transport | bump `u` at frame 12 (exact) | height kept | torus off |
|---|---|---|---|---|
| no flow | - | 0 (0) | 1.00 | 0 (exact) |
| rigid rotation, 0.25 / frame | pullback | 2.99 (3.0) | 0.64 | 0.009 |
| rigid rotation | frame | 2.50 at frame 10 (2.5) | 0.58 | 0.006 |
| ring speed 0.5 (shear) | pullback | wrings the torus (0.45 off), stalls | - | 0.45 |
| ring speed 0.5 | parallel | blows up from frame 1 | - | 0.53 |
| ring speed 0.5 | frame | 2.29 (2.22) | 0.35 (streaks) | 0.004 |

Transport is the deciding choice once the flow is not an isometry. The
pullback deforms the torus's own curvature (`J^T A J` of a non-isometric
flow is no surface's second form), and ambient `parallel` tilts it (each
departure tensor is projected without turning). The `frame` transport
(`(H, z e^{-2 i psi})` interpolated as numbers against a reference tangent
field, here the ring direction; `set_reference`) moves only what the flow
moves: the bump travels at the flow's rate and smears as the differential
speed shears it, and the torus holds. That reference is what the demo
already has (the shear axis, the store's principal directions).

Burgers as geometry (ring wave of 0.08 in the tube offset, `c` set so the
scalar inner shock falls near frame 6, MacCormack, frame transport): the
realized inner-equator slope doubles (0.062 -> 0.120 by frame 11, then
level), the outer barely moves (0.030 -> 0.039), `H` stays bounded (max
0.073 throughout). The crease forms, inner first, but slower than the
scalar wave: the contact net realizes each frame's increment only partly
(memory, 30 iterates), so the geometry lags the field.

Next, if this goes to the demo: `frame` transport against the shear axis
as a target driver for the store (C++: the advector, the increment targets
into `contact_net_gn`), BFECC for travelling patterns, MacCormack for
creases; a capped pullback only where the drive should stretch curvature.
