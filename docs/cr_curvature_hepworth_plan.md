# CR curvature in hepworth and kusama — implementation plan

Anisotropic curvature control on triangle meshes through the Crouzeix–Raviart
(CR) shape operator, built as composable hepworth constraints (coupled to the
curvature blocks the demo already carries) and as a standalone kusama solver
next to the Eigensatz solver (uncoupled use).

Derivations, experiments and numbers live in
`scripts/sympy/cr_poisson/MATH.md` (sections 9–12); the roadmap is
`docs/cr_poisson_shape_roadmap.md` (this is its phase C). Papers:
Stein, Wardetzky, Jacobson, Grinspun 2020 (`reference/stein2020-vector-dirichlet-energy.pdf`,
App. A.2); Wang, Liu, Tong 2012 (`reference/wang2012-discrete-fundamental-forms.pdf`);
Bouaziz et al. 2014 (projective dynamics); Overby et al. 2017 (ADMM ⊇ PD);
Eigensatz et al. 2010 (`reference/eigensatz2010-constraint-based-surface-processing.pdf`).

---

## 0. As built (part 2)

Sections 3.2a and 3.3a fixed the design; sections 4 and 5 below are the
earlier plan, and their file names were superseded as follows.

| Piece | File |
|---|---|
| Stencil walker (shell, triangle list) | `asawa/shell/stencil_walker.hpp` |
| Flat stencil table, `permuted_simplex_view<6>` gather | `kusama/stencil_table.hpp` |
| Generated edge read-out (9 readings, 9 x 18 `J`) | `scripts/sympy/generated/cr_curvature_edge_generated.hpp` |
| Builders and the uncoupled solver | `kusama/weingarten.hpp`: CR derivative, edge coupling, intrinsic frames and memory rows, spin-2 element, vertex stars (edge to vertex), `weingarten_solver` |
| Edge store | `shell_cr_block::shape` / `shape_dof`, the live edge-average reading each prepare; massless, so its rows determine it |
| Read-out and memory rows | `hepworth/block/weingarten_gn.hpp`: `weingarten_gn` (the coupling row between position and shape space: positions and store columns, rows `w k A_e`), `dirichlet_memory_gn` (rows `w sqrt(A)`) |
| Unified fairing (target rows) | `hepworth/block/edge_shape_fair.hpp` (rows `pin` / `dirichlet` / `dirichlet_trace` / `h1`, `p = A g(q)`, the same `k A_e` factor as the coupling), `edge_shape_goals.hpp` (zero, fixed, anisotropy against the frame-start store, ratio smooth, principal mean) |
| Vertex charts | `vertex_frame_cache_compute` takes `h` from the store's local polar; `frame_fair` and the principal solver read it live through `vertex_stars` |
| Demo | `monge_fairing_demo.hpp`: `curvature_config { curvature_solver kind, w_coupling, principal settings }`, `position_distortion::w_memory`, `shape_distortion` (`anisotropy` with `w_anisotropy` and `smooth`, grad, mean, ratio, `w_minimizer`), `curvature_test_target` as an edge pin |
| Tests | `test/weingarten_tests.hpp` (`tests/weingarten_tests.cpp`) |

Findings from building it:
- **The projection solver needs a line search for stiff GN rows.** With full
  steps a stiff coupling diverges; `projection_solver::gn_backtracks`
  (halve toward the previous iterate while the energy rises) fixes it.
- **Inertia is the proximal anchor.** Over 20 frames the torus tube reaches
  5.8e-13 and cube to sphere 8.6e-12 (coupling `w = 1e3` and `1e4`, store
  pin ten times that, unit vertex mass): low-frequency rounding advances a
  fraction per frame.
- **The `sobolev` step was tried and removed:** 0.2 to 0.12 in 20 frames
  against 2e-12 for Gauss–Newton, and stable only near `w = 1`.
- **Every row is `w` times its geometry; the coupling keeps shape and
  position in sync.** Shape rows act on the store `S`, position rows on `x`,
  and `weingarten_gn` couples them, like the rod rows that tie positions to
  frames. Its residual is a curvature; times the edge's dual area `A_e` it is
  an integrated curvature, the quantity the cotan `bend` row weighs, and one
  constant `k_coupling = 9` (`kusama/weingarten.hpp`) makes the two match at
  equal weights (position stiffness ratio 0.99 on the demo torus and 1.06 on
  a coarser one; both are dimensionless, so it holds across size and
  resolution). The store is a massless
  dof, the target rows use the same factor, and the least-squares solve
  combines them: `S = (w_s^2 R + sum w_k^2 g_k) / (w_s^2 + sum w_k^2)` for
  pins (exact to 1e-9 in `weingarten_store_is_the_blend_of_its_rows`), and
  positions feel a target at `w_s^2 w_k^2 / (w_s^2 + w_k^2)` (measured ratio
  1.80026 against 1.8). An earlier version scaled the rows to inertia and
  needed a store-side unit to match them; the fixed factor replaces both.
- **`h1` rows carry uniform edits.** Their pin part weighs exactly like the
  read-out; Dirichlet rows see only the variation of an edit (6 to 9 times
  less motion for a ratio edit on the torus).
- **Twist turns the curvature lines only slightly** (about 0.01 to 0.02 rad
  per frame for a 0.3 rad knob) while moving the mesh a lot: a uniformly
  turned principal field is not realizable in general (Gauss–Codazzi).
- **Edge-store charts beat the face-operator charts** on the analytic torus
  (rms `|H|` error 0.7% against 2.4%).
- **Ratio-type goals use `|k_hi| / |k_lo|`**; ordering by value flips for
  negative curvature and made the goals oscillate.
- Not ported: the ring-based Crane spoke map (the anisotropy knobs act per
  edge), and the face-site operators (`face_shape`, `face_shape_pin`,
  `face_fair`, `one_ring_operators` off) stay as they were.

### 0.1 Shape-domain controls

Each edge stores a symmetric 2x2 curvature tensor `S = (k11, k12, k22)`.
Every knob is written in one split of it:

```
S = H I + D,   H = (k0 + k1) / 2   (mean curvature)
               D = (d, tau)         (anisotropic part, |D| = (k0 - k1) / 2)
```

`k0 >= k1` are the principal curvatures (red and green in the demo). In the
frame of the principal axes `tau = 0` and `d = |D|`.

**Ratio.** `r = |k_hi| / |k_lo| >= 1`: how many times harder the surface bends
one way than the other. `q = (r - 1) / (r + 1) = min(|H|, |D|) / max(|H|, |D|)`
puts the same number in [0, 1], and `log r = 2 atanh(q)`.

| point | `r` | `q` |
|---|---|---|
| round (sphere, `k0 = k1`) | 1 | 0 |
| minimal saddle (`k0 = -k1`) | 1 | 0 |
| torus outer side, `R = 3 r_tube` | 4 | 0.6 |
| parabolic (cylinder, one direction flat) | infinity | 1 |

**Setting a ratio** keeps the larger of `|H|` and `|D|` and sets the smaller to
`q` times it, with signs and axes kept. Elliptic edges (`|D| <= |H|`) therefore
keep their mean curvature, and ratio 1 makes them round. Saddle edges keep
their anisotropic part, and ratio 1 makes them minimal. The two rules meet at
parabolic edges, so the edit is continuous; nothing is regularized and there
is no hidden scale. Round edges have no axes and ignore ratios above 1.

**Knobs** (`anisotropy_target`, fields in declaration order: `ratio`, `log`,
`relative`, `shear`, `mean`), applied
in this order against the frame-start store `h0`:

1. `shear` adds `shear * J D0` (`J` turns `(d, tau)` by 90 degrees): torsion
   on the frame-start lines. The axes turn by `atan(shear) / 2` and `|D|`
   grows by `sqrt(1 + shear^2)`.
2. `mean` sets `H = mean * H0`. Otherwise `H` is the live store's.
3. `ratio` (with `log` / `relative`) as above. Relative mode multiplies the
   frame-start ratio (a parabolic edge counts as 1000).

`shape_distortion` holds one list of these (`anisotropy`) with one weight
(`w_anisotropy`) and `smooth`: 0 pins each edge to its goal, above 0 adds the
spin-2 Dirichlet rows on the difference (the `h1` rows).

With no knob on the deviator, the deviator is the live store's, so the
defaults change nothing. Each goal depends on the live store only through
`H` (or `|D|` on saddles), so it is idempotent under the per-iterate
re-projection, and the per-frame knobs compound once per frame.

**Metrics** for gains (`curvature_metric`): `aniso = q`, `log_ratio = log r`
(capped at `log 1000`), `ratio = 1 / r`.

**Ring tilt** (`tilt_target`, position space). It turns each vertex's
frame-start one-ring by `gain * metric * |κ|^power` per frame and does not
settle. A principal axis is a line, so tilting about `a` and about `-a` go
opposite ways, and the mesh orientation (`N`) does not choose between them.
- `spin` turns about `N`, so it needs no sign.
- `twist` and `bend` turn each spoke by an angle proportional to its own
  coordinate, so flipping `a` cancels.
- `rigid` with `sign = gradient` uses the rotation vector
  `gain * l * (grad metric . (N x a)) a`, which is the same for `-a`.
- `rigid` with `sign = propagate` orients the axis field breadth-first over
  the mesh. Seams remain only around umbilics.

---

## 1. What we are building, in one paragraph

Each face reads its shape operator `S_f` from its 6-vertex stencil with the CR
derivative recomputed from the current positions (the read-out follows the
metric). A curvature row asks `S_f` to match what the curvature block (the
vertex `H` the demo already solves for) says it should be. An edge-centered ARAP
term (rotations fitted over each edge's two faces) is the memory that keeps the
metric from wandering. Both are ordinary hepworth constraints: the ARAP rows are
projection rows, the curvature rows are Gauss–Newton rows whose Jacobian is a
generated, analytic 3 x 18 block per face. Every frame the system is rebuilt
(the demo remeshes every frame anyway), factored, and iterated about 30 times;
the next frame starts from where this one stopped. Frames are the restarts that
made the Python solver converge exactly. There is no separate curvature goal:
constraints on `H` are the goal, and the positions react inside the same solve.

## 2. Findings this plan rests on

| Finding | Evidence (MATH.md) |
|---|---|
| The read-out must follow the metric. With `D` frozen at rest, non-isometric targets stall at 1.2–3.5% of size; with `D` recomputed per face stencil they converge to the exact shape. | §11 |
| Each restart (frame) cuts the error about 10x; the first solve is already within 0.2%. | §11 table |
| Edge rotations (CR-ARAP, diamond Procrustes) beat vertex rotations as the memory: 8x less bias on isometric bends, half the error on cube → sphere. | §10 |
| `S = S_g` alone does not pin the shape; the memory selects among solutions. | §9 |
| Curvature as a right-hand-side force only (`H = w L`, no `J^T J`) still converges, ~10x slower per frame, weight-sensitive. | §11b |
| The per-face read-out and its Jacobian are generated and verified (Python and C++). | §12 |

## 3. The math, precisely

### 3.1 Stencil and read-out

Face `f` with corners `P0, P1, P2` in winding order. For local edge
`i = (P_i, P_{i+1})`, `P_{3+i}` is the third vertex of the neighboring face
across that edge; `h_i = 1` if the neighbor exists, `0` on a boundary edge.
On a boundary edge the slot holds the reflected dummy
`P_{3+i} = P_i + P_{i+1} - P_{i+2}` (recomputed from the live corners each
evaluation): the masked neighbor area stays nonzero, so the masked derivative
is finite (`0 * finite`), and the dummy's Jacobian columns are exactly zero.

    s_f(P) = (k11, k12, k22) of S_f in the frame
             T = (P1 - P0) / |P1 - P0|,  N = (P1 - P0) x (P2 - P0) / |.|,  B = N x T

    S_f = K * sum_i  l_i / A * (N . g_perp,i) * v_perp,i v_perp,i^T        (K = -1)
    g_perp,i = ( D_f,i + h_i D_g,i ) P / M_e,     M_e = (A + h_i A_g) / 3

`D_f,i`, `D_g,i` are the per-face contributions to the CR perpendicular
derivative row of edge `i` (same assembly as `cr_coordinate_gradients`), so
`g_perp` uses the live lengths and areas of both faces of the edge. The edge
orientation sign cancels; the face winding is used throughout.

`cr_generated::cr_curvature_stencil(P, h, s, J)` (generated header) returns
`s` and `J = ds/dP` (3 x 18, column `3k + c` is `P[k][c]`).
Properties used as test oracles: `J` annihilates translations and
infinitesimal rotations; `J p = -s` (homogeneous of degree -1).

### 3.2 Where the curvature lives (no separate goal)

Two stores already exist as solver unknowns on `shell_cr_block`:

| Store | Layout | Already acted on by |
|---|---|---|
| `vertex_h_dof` | `(k11, k12, k22)` per vertex in its chart `(T_v, B_v)` (`vertex_frames`) | `shape_distortion` rows (anisotropy, Crane spokes, delta, grad, curvature mean / ratio Laplacians); `eigensatz_gn` targets |
| `face_dof` (`face_shape`) | `(k11, k12, k22)` per face in the face's first-edge frame | `face_shape_pin` |

**Vertex coupling.** Lift each corner's `H_v` to a 3D tangent tensor and read it
in the face frame:

    S_v = a T_v T_v^T + b (T_v B_v^T + B_v T_v^T) + c B_v B_v^T,      H_v = (a, b, c)
    P_f H = 1/3 sum_{v in f} ( T_f.S_v.T_f,  T_f.S_v.B_f,  B_f.S_v.B_f )

With `t1 = T_f.T_v, t2 = T_f.B_v, u1 = B_f.T_v, u2 = B_f.B_v`, corner `v`'s
3 x 3 block of `P_f` is

    [ t1^2      2 t1 t2            t2^2  ]
    [ t1 u1     t1 u2 + t2 u1      t2 u2 ]  / 3
    [ u1^2      2 u1 u2            u2^2  ]

`P_f` (3 x 9) is built from the frame-start charts and face frames and held for
the frame, so the coupling is linear in `H`.

**Delta form.** Averaging vertex curvature onto faces smooths it, so the plain
residual `s_f(X) - P_f H` is not zero at rest and would fair the surface every
frame. Measure changes instead:

    r_f(X, H) = s_f(X) - ( s_f(X0) + P_f (H - H0) )

`X0, H0` are the frame-start positions and `H`. At frame start `r_f = 0`; the
rows only act when something moves `H` (or the shape) away from its start.

**Face coupling.** `r_f = s_f(X) - (s_f(X0) + (h_f - h_f0))`, `dr/dh_f = -I`.
No frame transport; same frame as the read-out.

**Target mode** (tests, fixed edits, uncoupled): `r_f = s_f(X) - t_f`, with
`t_f` from a lambda (may read any block) or a frame-start cache.

### 3.2a Site decision (MATH.md section 14): the tensor lives on edges

Supersedes the vertex / face coupling of section 3.2 for the new solver.
The curvature store is per edge, `(k11, k12, k22)` in the edge's unfolded
frame `(t_e, N x t_e)` (`shell_cr_block::shape` / `shape_dof`, which
`emit_dof_blocks` must start emitting). Coupling: the mean of the edge's two
faces' live readings, each face's read-out rotated into the edge frame,

    r_e = 1/2 (s_f1^(e)(X) + s_f2^(e)(X)) - h_e,

3 rows per edge over the union of the two face stencils (8 vertices),
exactly zero at a realizable target. Smoothness rows on the store: the CR
spin-2 element, `E = sum_f sum_ab K_ab q(M_a h_a, M_b h_b)` (CR scalar
stiffness times spin-2 frame maps; Stein's vector Dirichlet construction
applied to tensors), second-order convergent on every triangulation tested.
The face site is ruled out: its dual-graph smoothness energy is singular
(circumcentric weights) or inconsistent (barycentric) on general meshes, the
face-butterfly distortion. The vertex `H` row is retired; its constraints
move to the edge store.

The generated read-out needs a 9-output variant (each face's operator in its
three edges' frames, Jacobian 9 x 18) from
`scripts/sympy/cr_poisson/cr_curvature_jacobian.py --emit`.

### 3.3a Memory decision (MATH.md section 13)

Under the working solver (live read-out, restarts), Stein's vector Dirichlet
energy on the intrinsic edge frames beats CR-ARAP: same accuracy on the tori,
20x better on an open isometric bend, and three orders of magnitude better on
cube -> sphere / ellipsoid at frame 8, converging about 10x per frame from the
first frame. It is the default memory:

    U_e = [t_e, b_e]^T [g_par, g_perp],   t = g_par / |g_par|,  b = n x t
    E = 1/2 sum_j c_j^T (L_dir + mu M2) c_j,   c_j = column j of U(D X) - U0

As hepworth rows: `L_dir` is assembled from per-face element matrices
`K_f` (6 x 6, positive semidefinite, over the face's three edges' parallel /
perpendicular components), so the energy splits by face:

    E = 1/2 sum_f sum_j c_{f,j}^T K_f c_{f,j}  (+ the mu M2 term per edge)

One Gauss–Newton constraint per face, residual `r_f = K_f^{1/2} c_f` (12
rows: 6 per column `j`), Jacobian `K_f^{1/2} dc_f/dX` over the face's three
edge diamonds, which is the same 6-vertex stencil as the curvature rows.
`K_f^{1/2}` comes from a 6 x 6 eigendecomposition at frame start. CR-ARAP
below stays as an option (constant projection rows).

### 3.3 CR-ARAP (edge ARAP)

Rotation per edge, energy over its diamond (the two incident faces):

    E = sum_e sum_{f ni e} (A_f / 3) |grad_f X - R_e grad_f X0|^2
      = sum_e sum_{f ni e} 1/3 sum_{i in f} w_f,i |d_i - R_e d0_i|^2

`d_i = x_head - x_tail` for the face's three spokes, `w_f,i = cot(angle opposite
i) / 2` (the per-face gradient energy is exactly that cotan sum). The rotation
`R_e = polar( sum_{f ni e} sum_{i in f} w_f,i d_i d0_i^T )` is the Procrustes fit
over the diamond: the center edge gets its full cotan weight, the four wing
spokes their own face's half. Each face's weights over its three edges sum to
one face's worth, so `sum A^T A = L_cot` (constant within the frame). The
gradient is exact (envelope theorem: the energy is minimized over the rotations).

As PD rows: per edge, 6 spoke rows of 3 (two faces x three spokes), row
`sqrt(w * max(w_f,i, 0) / 3) * (x_head - x_tail)`, projection
`sqrt(w * max(w_f,i, 0) / 3) * R_e d0_i`. Negative cotan weights (obtuse
angles) are clamped to zero to keep `A^T A` positive semidefinite.

Rest `d0`: the frame-start shape (default), which is the restart protocol that
converged in Python — the memory resists within a frame and moves with the
shape across frames.

### 3.4 What the solver minimizes

`projection_solver` with Gauss–Newton rows (`solver.hpp`):

    E(q) = |q - s|^2_{M/h^2} + sum_proj |A_i q - p_i(q)|^2 + sum_gn |r_j(q)|^2
    (M/h^2 + A^T A + J_ref^T J_ref) q_{k+1} = (M/h^2) s + A_proj^T p + J_ref^T J_ref q_k - J(q_k)^T r(q_k)

`J_ref` is filled at frame start and refreshed every `gn_refresh_every` iterates
(0 or >= `iterations`: frozen for the frame) or when a step on the stale factor
raises `E` (`gn_refresh_on_ascent`). The gradient is fresh every iterate.
Inertia `|q - s|^2_{M/h^2}` is the proximal anchor; each frame re-anchors at the
momentum prediction `s` (with damping this is close to the previous solution).

**Right-hand-side mode** (`mode = rhs`): the constraint reserves its rows but
emits no Jacobian triplets, so the left-hand side never sees `J`; only
`-J^T r` enters. Converges ~10x slower per frame (§11b); a knob, not the default.

### 3.5 Weights

* Curvature rows: `w` is relative to inertia, like `eigensatz_gn`. At the first
  `fill_J` of a frame, `scale = w * mean(m / h^2) / mean(diag(J^T J))` over the
  position columns, held for the frame. Per-face area weight inside:
  `r_f -> sqrt(A_f / mean A) * r_f`.
* CR-ARAP rows: raw `w` like the other hepworth position rows (`init_area`,
  `init_bending`), with the cotan and `1/3` factors inside.
* Coupling balance: `H`'s mass is ~1e-8, so `H` is held by its own rows and by
  the coupling. Strong `H` rows relative to the coupling: the shape follows your
  edits. Weak: `H` relaxes toward what the surface can do. This ratio is the
  main art knob.

---

## 4. Files

| File | Contents |
|---|---|
| `scripts/sympy/cr_poisson/cr_curvature_jacobian.py` | derivation, verification, `--emit` (exists) |
| `scripts/sympy/generated/cr_curvature_stencil_generated.hpp` | `cr_generated::cr_curvature_stencil` (exists) |
| `include/gaudi/hepworth/block/cr_curvature_stencil.hpp` | stencil topology, dummy reflection, face frame, `evaluate`, `P_f` lift, diamond tables |
| `include/gaudi/hepworth/block/cr_edge_arap.hpp` | `cr_edge_arap` projection constraint, `init_cr_edge_arap` |
| `include/gaudi/hepworth/block/cr_curvature_gn.hpp` | `cr_curvature_gn` Gauss–Newton constraint (vertex / face / target), `init_cr_curvature` |
| `include/gaudi/hepworth/block/cr_curvature_targets.hpp` | target helpers and `cr_curvature_target_compute` (`cache_compute`) |
| `include/gaudi/kusama/cr_curvature_solver.hpp` | standalone solver, API shaped like `eigensatz_solver` |
| `include/gaudi/duchamp/monge_fairing_demo.hpp` | `cr_curvature_config`, `push_cr_curvature`, test modes |
| `include/gaudi/test/cr_curvature_tests.hpp`, `tests/cr_curvature_tests.cpp` | tests, registered in `tests/CMakeLists.txt` like `eigensatz_tests` |

The generated header lives under `scripts/sympy/generated/`, which the
existing generated headers already use as an include root
(`#include "cr_curvature_stencil_generated.hpp"`).

---

## 5. Pseudocode

Pseudocode is C++-flavored; names match hepworth conventions
(`_blocks[k]->get_vec3(i, q)`, `get_offset_idx(i)`, `_id0`, `DEFINE_CREATE_FUNC`).

### 5.1 Stencils (`cr_curvature_stencil.hpp`)

```cpp
struct cr_face_stencil {
  std::array<index_t, 6> v;   // corners 0..2, opposite vertices 3..5 (-1 if boundary)
  vec3 h;                     // 1 = neighbor exists, 0 = boundary edge
  int face = -1;              // face id (for face_dof)
};

struct cr_diamond {           // one edge and its incident faces
  index_t e = -1;
  std::array<int, 2> faces{-1, -1};   // -1 if boundary
};

// Built from the current mesh each frame (topology changes with the remesher).
std::vector<cr_face_stencil> build_face_stencils(const asawa::shell::shell &M) {
  for each face f with fsize(f) == 3:
    tri = M.get_tri(f)                 // winding order
    st.v[0..2] = tri
    for i in 0..2:
      corner c = corner of edge (tri[i] -> tri[i+1]) in f
      c_opp = M.other(c)               // the twin corner across the edge
      if c_opp is a boundary / has no face: st.v[3+i] = -1; st.h[i] = 0
      else: st.v[3+i] = third vertex of face(c_opp); st.h[i] = 1
    stencils.push_back(st)
}

// Live positions for one stencil, with the reflected dummy on boundary edges.
std::array<vec3, 6> gather(const cr_face_stencil &st, const sim_block &pos, const vecX &q) {
  P[k] = pos.get_vec3(st.v[k], q) for k = 0..2
  for i in 0..2:
    P[3+i] = st.h[i] > 0 ? pos.get_vec3(st.v[3+i], q)
                         : P[i] + P[(i+1)%3] - P[(i+2)%3];
}

struct cr_face_eval { vec3 s; Eigen::Matrix<real, 3, 18> J; bool ok; };

cr_face_eval evaluate(const std::array<vec3, 6> &P, const vec3 &h) {
  if face (P0, P1, P2) is degenerate (area or an edge below eps): ok = false
  cr_generated::cr_curvature_stencil(P, h, s, J);
  ok = s.allFinite() && J.allFinite();
}

// Face frame, the same one the generated read-out uses.
void face_frame(P0, P1, P2, vec3 &T, vec3 &B, vec3 &N) {
  T = (P1 - P0).normalized(); N = (P1 - P0).cross(P2 - P0).normalized(); B = N.cross(T);
}

// Corner block of P_f (3 x 3), from the frame-start charts and face frame.
mat3 lift_block(const vec3 &Tf, const vec3 &Bf, const vertex_chart &c) {
  t1 = Tf.dot(c.T); t2 = Tf.dot(c.B); u1 = Bf.dot(c.T); u2 = Bf.dot(c.B);
  return (1/3) * [[t1*t1, 2*t1*t2,       t2*t2],
                  [t1*u1, t1*u2 + t2*u1, t2*u2],
                  [u1*u1, 2*u1*u2,       u2*u2]];
}

// Diamonds and cotan spoke weights for CR-ARAP.
std::vector<cr_diamond> build_diamonds(const asawa::shell::shell &M);
std::array<real, 3> face_spoke_weights(x0, tri);  // cot(opposite angle) / 2, spoke i = (tri[i], tri[i+1])
```

### 5.2 CR-ARAP projection constraint (`cr_edge_arap.hpp`)

```cpp
class cr_edge_arap : public block_constraint {
public:
  DEFINE_CREATE_FUNC(cr_edge_arap)
  // ids: the up-to-4 diamond vertices; spokes: (local tail, local head, weight)
  // over both faces; d0: rest spokes (frame start).
  cr_edge_arap(std::vector<index_t> ids, std::vector<spoke> spokes,
               std::vector<vec3> d0, real w, std::vector<sim_block::ptr> blocks);

  std::string name() override { return "cr_edge_arap"; }

  void fill_A(index_t &id0, std::vector<trip> &trips) override {
    _id0 = id0;
    for (k, sp) in enumerate(_spokes):
      s = sqrt(_w * max(sp.weight, 0) / 3);
      ct = _blocks[0]->get_offset_idx(_ids[sp.tail]);
      ch = _blocks[0]->get_offset_idx(_ids[sp.head]);
      for ax in 0..2:
        trips.push_back(trip(id0 + 3k + ax, ct + ax, -s));
        trips.push_back(trip(id0 + 3k + ax, ch + ax,  s));
    id0 += 3 * _spokes.size();
  }

  void project(const vecX &q, vecX &p) override {
    C = 0 (3 x 3)
    for sp in _spokes:
      d = x[head] - x[tail];  C += max(sp.weight, 0) * d * d0[k]^T;
    R = polar(C)                    // SVD, det fix (as triangle_strain / shape_match)
    for (k, sp): p.segment<3>(_id0 + 3k) = sqrt(_w * max(sp.weight, 0) / 3) * R * d0[k];
  }
};

void init_cr_edge_arap(asawa::shell::shell &M, std::vector<constraint::ptr> &out,
                       const std::vector<vec3> &x0, real w,
                       std::vector<sim_block::ptr> blocks) {
  if (!(w > 0)) return;
  for each edge e with its incident faces (1 or 2):
    collect diamond vertices -> ids (dedup), spokes of each face with
    face_spoke_weights(x0, tri), rest spokes d0 from x0
    out.push_back(cr_edge_arap::create(ids, spokes, d0, w, blocks));
}
```

### 5.3 Curvature Gauss–Newton constraint (`cr_curvature_gn.hpp`)

```cpp
enum class cr_curvature_source { vertex_h, face_shape, target };
enum class cr_curvature_mode   { gn, rhs };

struct cr_curvature_config {
  real w = 1.0;                                   // inertia-relative
  cr_curvature_source source = cr_curvature_source::vertex_h;
  cr_curvature_mode mode = cr_curvature_mode::gn;
};

// Per-face target for `source = target` (may read any captured block).
using cr_face_target_fn = std::function<bool(const vecX &q, int face, vec3 &t)>;

class cr_curvature_gn : public block_constraint {
public:
  DEFINE_CREATE_FUNC(cr_curvature_gn)
  // blocks[0]: positions; blocks[1]: vertex_h or face_dof (source-dependent).
  // Frame-start data baked by init: s0 = s_f(X0); H0 (3 corners) or h_f0;
  // Pf (3 x 9) for vertex_h; area weight a_f.
  cr_curvature_gn(cr_face_stencil st, cr_curvature_config cfg, real h,
                  baked frame-start data, std::vector<sim_block::ptr> blocks,
                  cr_face_target_fn target = {});

  bool gauss_newton() const override { return true; }
  std::string name() override { return "cr_curvature_gn"; }

  // r = sqrt(scale * a_f) * (s_f(X) - goal(H)),  goal per source (3.2).
  vec3 residual(const vecX &q, cr_face_eval &ev) const {
    P = gather(_st, *_blocks[0], q);  ev = evaluate(P, _st.h);
    if (!ev.ok) return vec3::Zero();
    switch (_cfg.source):
      vertex_h:   goal = _s0 + _Pf * (H(q) - _H0)       // H(q): 9-vector of corner H
      face_shape: goal = _s0 + (h_f(q) - _hf0)
      target:     goal = _target(q, _st.face)            // skip row if false
    return ev.s - goal;
  }

  void fill_J(const vecX &q, index_t &id0, std::vector<trip> &trips) override {
    cr_face_eval ev;  residual(q, ev);
    ensure_scale(ev)                // first fill of the frame: inertia-relative scale
    const real s = sqrt(_scale * _a);
    if (_cfg.mode == gn && ev.ok) {
      for k in 0..5 (skip k >= 3 && _st.h[k-3] == 0):
        col = _blocks[0]->get_offset_idx(_st.v[k]);
        for r in 0..2, c in 0..2: trips.push_back(trip(id0 + r, col + c, s * ev.J(r, 3k + c)));
      if source == vertex_h:
        for corner j in 0..2: colH = _blocks[1]->get_offset_idx(_st.v[j]);
          for r, c: trips.push_back(trip(id0 + r, colH + c, -s * _Pf(r, 3j + c)));
      if source == face_shape:
        colh = _blocks[1]->get_offset_idx(_st.face);
        for r: trips.push_back(trip(id0 + r, colh + r, -s));
    }
    id0 += 3;                        // same row count every call (solver contract)
  }

  void add_gradient(const vecX &q, vecX &g) override {
    cr_face_eval ev;  vec3 r = residual(q, ev);
    if (!ev.ok) return;
    const real w = _scale * _a;
    // positions: J^T r
    for k (as above): g.segment<3>(col_k) += w * ev.J.block<3,3>(0, 3k).transpose() * r;
    // curvature block: (dr/dH)^T r
    vertex_h:   for j: g.segment<3>(colH_j) += -w * _Pf.block<3,3>(0, 3j).transpose() * r;
    face_shape: g.segment<3>(colh) += -w * r;
  }

  real residual_sq(const vecX &q) override {
    cr_face_eval ev;  vec3 r = residual(q, ev);
    return ev.ok ? _scale * _a * r.squaredNorm() : 0;
  }
};

void init_cr_curvature(asawa::shell::shell &M, const shell_cr_block &cr,
                       std::vector<constraint::ptr> &out, const std::vector<vec3> &x0,
                       const cr_curvature_config &cfg, real h,
                       std::vector<sim_block::ptr> blocks, cr_face_target_fn target = {}) {
  stencils = build_face_stencils(M);
  mean_area = mean of face areas at x0
  for st in stencils:
    P0 = gather at x0; ev0 = evaluate(P0, st.h); if !ev0.ok continue;
    baked.s0 = ev0.s;  baked.a = area(P0) / mean_area;
    if source == vertex_h:
      (Tf, Bf, Nf) = face_frame(P0);
      for j in 0..2: c = cr.vertex_frames[st.v[j]]; if !c.ok -> skip face
        baked.Pf.block<3,3>(0, 3j) = lift_block(Tf, Bf, c);  baked.H0[j] = cr.vertex_h[st.v[j]];
    if source == face_shape: baked.hf0 = cr.face_shape[st.face];
    out.push_back(cr_curvature_gn::create(st, cfg, h, baked, blocks, target));
}
```

The inertia-relative scale is shared across the frame's faces (one mean over
all curvature rows), so `init_cr_curvature` hands every face a pointer to one
`cr_curvature_scale` that the first `fill_J` sets (reset when the next frame
builds new constraints).

### 5.4 Targets (`cr_curvature_targets.hpp`)

```cpp
// Frame-start targets (tests, tool strokes). Written once per frame.
class cr_curvature_target_compute : public cache_compute {
  std::vector<vec3> *targets;  mode: scale(k), umbilic(1/R), mesh(Y), lambda
  void compute(vecX &q) override {
    for each face f:
      scale:   t_f = k * s_f(X)
      umbilic: t_f = (1/R, 0, 1/R) * K-sign convention
      mesh:    t_f = s_f(Y)          // same connectivity: torus radii, cube -> sphere
      lambda:  t_f = fn(q, f)
  }
};

// World tensor field -> face-frame target (painted / analytic fields).
vec3 face_target_from_world(const mat3 &S, const vec3 &Tf, const vec3 &Bf);
// Principal values + directions at vertices -> vertex_h in its chart.
vec3 vertex_h_from_principal(real k0, real k1, const vec3 &d0, const vertex_chart &c);
```

For the demo the test targets act on `H` (a pin row on `vertex_h` toward
`vertex_h_from_principal` values), which keeps "no separate goal": the test is
just another constraint on the curvature block.

### 5.5 Demo wiring (`monge_fairing_demo.hpp`)

```cpp
struct cr_curvature_demo_config {
  bool on = false;
  real w = 1.0;                       // curvature coupling, inertia-relative
  real w_arap = 1.0;                  // CR-ARAP memory
  hepworth::block::cr_curvature_source source = hepworth::block::cr_curvature_source::vertex_h;
  hepworth::block::cr_curvature_mode mode = hepworth::block::cr_curvature_mode::gn;
};
// in shape_project:  cr_curvature_demo_config cr_curvature{};
// demo-only:         struct cr_curvature_test { enum mode { off, scale, umbilic, torus }; ... };

.with_recompute([this](solver_context &ctx) { push_cr_curvature(ctx); })

void push_cr_curvature(solver_context &ctx) {
  const auto &c = _cfg.project.cr_curvature;
  if (!c.on || !_cr || ctx.blocks.empty() || !__M->verts_are_dense_packed()) return;
  ctx.shell_mesh = __M.get();
  const auto &x0 = _shell->xs->get();
  init_cr_edge_arap(*__M, ctx.constraints, x0, c.w_arap, {ctx.blocks[0]});
  sim_block::ptr store = c.source == vertex_h ? _cr->vertex_h_dof() : _cr->face_dof();
  init_cr_curvature(*__M, *_cr, ctx.constraints, x0, {c.w, c.source, c.mode},
                    _cfg.step.h, {ctx.blocks[0], store}, test_target_or_empty());
  if (test mode != off) push the H pin toward the test values (5.4)
}
```

Step settings: `fairing_step.gn_refresh_every = 0` freezes `J` for the frame
(matches the Python restarts); keep `gn_refresh_on_ascent`. `iterations = 30`.
Turn `eigensatz.on` off when `cr_curvature.on` is on (both map `H` to positions).

### 5.6 kusama: standalone solver (`kusama/cr_curvature_solver.hpp`)

For uncoupled use (no hepworth, no curvature block): positions only, per-face
targets, CR-ARAP memory, frozen-factor Gauss–Newton with backtracking, and the
restart protocol. API shaped like `eigensatz_solver` so the demo's
`global` / `gauss_seidel` paths and probes can swap solvers.

```cpp
namespace gaudi::kusama {

struct cr_curvature_weights { real k_curv = 1; real k_arap = 1e-2; real k_close = 0; };

struct cr_curvature_report { real energy0, energy1, step; int line_halvings; bool refactored; };

class cr_curvature_solver {
public:
  cr_curvature_weights weights;
  int refresh_every = 30;   // factor reuse (frozen J) between refactors
  int backtracks = 12;

  void reset(int n, const std::vector<tri3> &tris, const std::vector<vec3> &x_rest) {
    build stencils (from tris: twin lookup by sorted edge key), diamonds,
    rest spokes and cotan weights (CR-ARAP), L_cot (kron I3), areas; drop the factor
  }
  void reset_rest(const std::vector<vec3> &x_rest);  // restart: new metric + ARAP rest

  std::vector<vec3> readout(const std::vector<vec3> &x) const;          // s_f per face
  real energy(const std::vector<vec3> &x, const std::vector<vec3> &t) const {
    E = sum_f a_f |s_f(x) - t_f|^2 + k_arap * E_crarap(x) + k_close |x - x_rest|^2
  }
  void linearize(const std::vector<vec3> &x, const std::vector<vec3> &t,
                 vecX &r, matS &J) const;                                // curvature rows (tests)

  cr_curvature_report iterate(std::vector<vec3> &x, const std::vector<vec3> &t) {
    g = sum_f a_f J_f^T r_f + k_arap (L_cot x - b(R(x))) + k_close (x - x_rest)
    if (!_factored || _reused >= refresh_every)
      factor H = sum_f a_f J_f^T J_f + k_arap L_cot + (k_close + eps) I;  _reused = 0
    d = -H^{-1} g;  backtrack a = 1, 1/2, ... until Armijo on energy
    if no step: refactor at x and retry once; else ++_reused
  }

  // Restart protocol (§11): solve to convergence, re-rest, repeat.
  void solve(std::vector<vec3> &x, const std::vector<vec3> &t,
             int restarts = 5, int max_its = 150, real tol = 1e-9) {
    for (k < restarts): reset_rest(x); iterate until step < tol or max_its
  }
};

// Eigensatz-style principal targets -> face targets, using the current
// principal directions (vertex charts) with new values.
std::vector<vec3> face_targets_from_principal(x, tris, kappa_hat, charts);
}
```

Uncoupled hepworth use of either kusama solver keeps the two existing modes:
Gauss–Newton rows (an `eigensatz_gn`-style adapter over
`cr_curvature_solver::linearize`) or Gauss–Seidel (`shape_match` pin toward one
`iterate`). The coupled path (5.3) is the recommended one in the demo.

| | coupled (H block) | uncoupled |
|---|---|---|
| CR curvature | `cr_curvature_gn` + `cr_edge_arap` (hepworth) | `kusama::cr_curvature_solver` |
| Eigensatz | `eigensatz_gn` (hepworth, exists) | `kusama::eigensatz_solver` (exists) |

---

## 6. Tests (`test/cr_curvature_tests.hpp`, `tests/cr_curvature_tests.cpp`)

| Test | Assertion |
|---|---|
| read-out vs `cr_face_operators` | random-perturbed torus: `s_f` equals `cr_face_operators` S read in the face frame (sign and frame convention), < 1e-10 |
| read-out vs stored `face_shape` / vertex `h` | after `shell_cr_block` prepare: `s_f` equals `face_shape[f]`; area-average lift of `s_f` matches `vertex_h` |
| `J` vs central differences | 200 random stencils, a third with a boundary edge: < 1e-8 relative |
| `J` invariances | translation and rotation in the null space; `J p = -s` |
| `P_f` lift | constant tensor field (`S` fixed in 3D): lifted vertex `h` reproduces the face-frame read of `S` |
| CR-ARAP idle and rigid | at rest `A q = p`; rigid motion leaves the energy 0; projection is a rotation |
| delta coupling idle | frame start: all curvature residuals 0; a solve with no `H` rows leaves `X` (and `H`) unchanged |
| target mode in `projection_solver` | torus 0.7 → 0.9 (target read from the torus mesh): shape error < 1% of size after 5 frames, decreasing every frame |
| cube → sphere | same, radial projection target: < 1% after 8 frames |
| `H` coupling | pin `vertex_h` to the torus-target values: positions follow (error < 2% after 8 frames) |
| `rhs` mode | converges on the torus target (slower), no blow-up |
| kusama solver | restart protocol reproduces Python §11: torus tube and cube → sphere to < 1e-6 of size; isometric cylinder to the CR floor |

Python oracles: `scripts/sympy/cr_poisson/` (`cr_curvature_jacobian.py`,
`cr_goal.goal_gn_loop(readout="live")`, `live_compare.py`).

---

## 7. Order of work and exit criteria

1. **Stencil + read-out** (`cr_curvature_stencil.hpp`), read-out and Jacobian
   tests. Exit: read-out matches `cr_face_operators` and the stored blocks.
2. **CR-ARAP** (`cr_edge_arap.hpp`). Exit: idle / rigid tests; with only
   CR-ARAP and inertia the demo is stable through remeshing.
3. **Curvature GN, target mode** (`cr_curvature_gn.hpp`). Exit: torus and cube
   targets in `projection_solver` meet the thresholds.
4. **Vertex / face coupling.** Exit: delta idle test; `H` pin test.
5. **Demo wiring** with test modes. Exit: torus radius change and a painted
   anisotropy edit visibly track over frames with remeshing on.
6. **kusama solver** and the uncoupled adapters. Exit: Python parity.
7. **Docs**: MATH.md pointer, roadmap phase C status, TODO item 5.

## 8. Risks and open questions

* **Vertex vs face tangent planes.** `P_f` projects the vertex tensor onto the
  face plane; on coarse, highly curved meshes the lift loses some normal
  information. The delta form makes this a second-order effect (only changes
  of `H` are transported).
* **`vertex_h` reseeding.** `seed_cr_from_dirichlet` with `reseed` rewrites
  `vertex_h` from the shape each frame; edits persist only through constraints
  re-pushed each frame (which is the design: constraints are the goal).
* **Obtuse triangles.** Clamped cotan weights in CR-ARAP; the remesher keeps
  most triangles well shaped.
* **Non-realizable `H`.** Constraints asking for curvature no surface has
  settle at a compromise set by `w`, the `H` rows and `w_arap`.
* **Scale.** The inertia-relative scale is set once per frame from the frame's
  first `J`; a large in-frame change of curvature stiffness is absorbed by the
  ascent refresh, not rescaled.
* **Cost.** ~5.7k flops per face per evaluation (value + `J`), twice per
  iterate (`add_gradient`, `residual_sq` share one evaluation if cached by
  `q` hash — optional).
* **Name.** "CR curvature" in code; the solver family nickname is open.
