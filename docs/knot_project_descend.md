# Knot project-then-descend

Design notes for knot / closed-rod fairing when tangent-point (TP) gradient descent alone is awkward: project into a **manifold family** whose intrinsic minimizer is a circle (or \(N\)-cover of a circle), then descend an energy on that projected curve — hierarchical least-squares projection instead of ADMM, with an optional chain rule through \(P\).

Related local papers / code:

- Crane et al., *Repulsive Curves* — [`crane2021-repulsive-curves.pdf`](../reference/crane2021-repulsive-curves.pdf): fractional TP + Sobolev-type metric
- Dennis & Hannay 2005, *Geometry of Călugăreanu's theorem* — [`dennis2005-geometry-calugareanu.pdf`](../reference/dennis2005-geometry-calugareanu.pdf): projection averages for Tw/Wr; **writhe framing** (\(\mathrm{Lk}=0\))
- Feng et al., *Points as Tori* — [`feng-points-as-tori.pdf`](../reference/feng-points-as-tori.pdf)
- Kidambi & Newton 2000, *Streamline topologies for integrable vortex motion on a sphere* — Physica D; DOI [10.1016/S0167-2789(99)00233-X](https://doi.org/10.1016/s0167-2789(99)00233-x) (limaçon streamline class; not yet in `reference/`)
- Braid host: [`braid_planar_demo.hpp`](../include/gaudi/duchamp/braid_planar_demo.hpp)
- Cauchy TP + line search: [`tp_energies.hpp`](../include/gaudi/duchamp/modules/tp_energies.hpp), [`tp_line_search.hpp`](../include/gaudi/duchamp/modules/tp_line_search.hpp)
- Calder MLS / LS fits (projection precedent): [`least_squares_fit.hpp`](../include/gaudi/calder/least_squares_fit.hpp)
- Circularity / asphericity energies: [`convex_medial_axis_energy.md`](convex_medial_axis_energy.md)

---

## Setup

Closed discrete rod (polyline) \(x = (x_0,\ldots,x_{n-1})\), edge lengths \(\ell_i\), unit tangents \(T_i\). Energy of interest is typically Cauchy TP

\[
K(dp,N)
  \;=\;
  \frac{|N\cdot dp|^p}{(|dp|^2+\varepsilon^2)^p},
\]

as in the braid demo (`tangent_point_type::cauchy`), optionally plus stretch / smooth / collision constraints in Hepworth.

**Observed pathology.** Minimizing TP on the bare polyline tends to produce forces concentrated in the **contact / binormal** direction of near pairs. Tangential redistribution and “become round” degrees of freedom move slowly or not at all — hence Crane’s Sobolev metric, and hence the geometric alternative below.

---

## Core idea: project, then descend

Let \(\mathcal{M}_\Theta\) be a **parametric family** of curves (or curve germs) with parameters \(\theta\in\Theta\), and let \(\Phi(\theta;\,u)\) reconstruct a point at chart coordinate \(u\).

Define a (possibly hierarchical) **projection**

\[
y \;=\; P(x),\qquad y_i \;\approx\; \Phi\bigl(\theta^\star(x);\, u_i(x)\bigr),
\]

then an energy \(\mathcal{E}(y)\) on the projected curve (TP, asphericity, holonomy, …), and descend

\[
E(x) \;=\; \mathcal{E}\bigl(P(x)\bigr).
\]

```text
x  ──►  hierarchical LS project P  ──►  y ∈ M_Θ
                                              │
                                              ▼
                                    E = TP / asphericity / …
                                              │
                                              ▼
                         g = −∇_y E   (freeze P, or chain through DP)
                                              │
                                              ▼
                                    update x, re-project
```

**Why this helps.** Choose \(\mathcal{M}_\Theta\) so that the family’s **intrinsic** minimizer is a round circle (ellipse → equal axes; limaçon → eccentricity \(0\); torus knot → \(N\)-cover of a circle as major/minor radius collapses). Shape DOFs live in \(\theta\); TP (or a cheaper proxy) on \(y\) only has to finish the job inside an already “almost circular” reduced space.

This is **not** ADMM. There is no dual ascent on constraints. Projection is ordinary weighted least squares at several bandwidths; descent is ordinary gradient (or Hepworth force) on \(\mathcal{E}\).

---

## Sobolev metric as geometry

In *Repulsive Curves*, the Sobolev / fractional metric changes the **inner product** used to convert \(\mathrm{d}E_{\mathrm{TP}}\) into a search direction: high-frequency contact spikes are damped, so the flow behaves like an \(H^{s}\) gradient rather than an \(L^2\) spike in the normal.

Geometric reading:

| Mechanism | What it changes | Roundness DOFs |
|-----------|-----------------|----------------|
| Sobolev / fractional metric | soft metric on variations of \(x\) | indirect (smoothed TP force) |
| Project to ellipse / limaçon | hard reduced space \(\mathcal{M}_\Theta\) | explicit in \(\theta\) |
| Hierarchical LS \(P\) | band-limited admissible shapes | coarse \(\theta\) locks global roundness |

Hierarchical projection is the **hard** analogue of Sobolev band-limiting: coarse levels admit only low-frequency shapes (near-circles in the family); fine levels carry residual polyline freedom. Same multigrid spirit as Calder multi-scale MLS, without dual variables.

---

## Hierarchical least-squares projection

Levels \(\ell = 0,\ldots,L\) with bandwidths

\[
\sigma_\ell \;=\; 2^\ell\,\sigma_0
\]

(\(\sigma_0\) in mean-edge or rod-radius units). Same Barnes–Hut / MLS gather pattern as Calder shell/rod fits.

### Per-level fit

At query sample \(i\) (vertex or edge midpoint), with neighbors weighted by a unit-peak window \(w_{ij}^{(\ell)}\) (e.g. \(\exp(-( |x_j-x_i|/\sigma_\ell)^2)\) or inv-dist in length scale \(\sigma_\ell\)):

\[
\theta_\ell^\star(i)
  \;=\;
  \arg\min_{\theta}
  \sum_j w_{ij}^{(\ell)}\,
  \bigl|\, x_j - \Phi(\theta;\, u_j) \,\bigr|^2.
\]

Chart coordinates \(u_j\) are arc-length, polar angle in a local frame, or spherical angle — depending on \(\Phi\).

### Reconstruction (partition of unity)

\[
y_i
  \;=\;
  \sum_{\ell=0}^{L}
  \alpha_\ell(i)\,
  \Phi\bigl(\theta_\ell^\star(i);\, u_i\bigr),
\qquad
\sum_\ell \alpha_\ell = 1,\quad \alpha_\ell \ge 0.
\]

Typical choice: fine levels dominate locally (\(\alpha_0\) large near high curvature); coarse levels get enough weight to lock **global** roundness / winding. Equivalently, treat coarse \(\Phi(\theta_L^\star)\) as a prior and fit fine residuals — “weighted LS of the higher-order idea.”

### What this is / is not

- **Is:** multigrid-style projection; each layer has its own manifold fit; local solution blends them.
- **Is not:** ADMM, Augmented Lagrangian, or a single global hard constraint \(x\in\mathcal{M}\).

---

## Chain rule through the projection

Write \(y = P(x)\) and \(E(x) = \mathcal{E}(P(x))\).

### Freeze-\(P\) step (recommended first)

Match the current TP kernel freeze in the braid loop:

1. Rebuild \(P\) from current \(x\) (hierarchical LS).
2. Hold \(P\) (and thus \(y\)) fixed for the force evaluation.
3. Compute \(g = -\nabla_y \mathcal{E}\).
4. Scatter \(g\) onto rod samples (nearest vertices / edges), or apply

\[
\nabla_x E \;\approx\; (DP)^\top g
\]

with \(DP\) from the last linearized LS if available.
5. Integrate (CCD / Hepworth velocity block), then re-project next frame.

### Differentiable LS (envelope)

If \(\theta^\star\) solves \(A(\theta)^\top A(\theta)\,\theta = A(\theta)^\top b(x)\) at a critical point,

\[
\mathrm{d}\theta
  \;=\;
  \bigl(A^\top A\bigr)^{-1}
  \bigl(\partial_x(A^\top b) - \partial_\theta(A^\top A\,\theta)\bigr)\mathrm{d}x
\]

(schematic; exact Jacobian depends on whether \(A\) depends on \(\theta\)). Then

\[
\nabla_x E
  \;=\;
  (D_x y)^\top \nabla_y \mathcal{E},
\qquad
D_x y \;=\; D_\theta\Phi\cdot D_x\theta^\star + \cdots.
\]

**Cost.** Full \(\mathrm{d}\theta\) per sample is expensive; freeze-\(P\) is enough to test whether the *family* helps. Use envelope only if freeze fights stretch/smooth constraints badly.

**Warning.** Projecting then descending TP on \(y\) without any pull-back can fight edge-stretch / smooth / min-kink already in [`braid_planar_demo`](../include/gaudi/duchamp/braid_planar_demo.hpp). Prefer scatter of \(\nabla_y\mathcal{E}\) as a soft force alongside those constraints, not a hard replace \(x\leftarrow y\).

---

## Candidate manifolds

| Family \(\Phi\) | Fit data | Natural minimizer of family | Notes |
|-----------------|----------|-----------------------------|-------|
| Ellipse (planar section) | points in local frame | circle (equal axes) | Closest to barrel \(E^{\mathrm{asp}}\) |
| **Twisted ellipse** | ellipse \((a,b)\) + twist \(\alpha\) about major spine | circle: \(a\to b\) and \(\alpha\to 0\) | Rotation DOF \(\alpha\) is the untwist direction |
| Limaçon on \(S^2\) | spherical samples | circle (\(b\to 0\)) | Kidambi–Newton; related endpoint of half-twist story |
| Cubic Hermite / Bézier | 2 ends + 2 tangents | **not** automatically a circle | Good *stencil*, bad *objective* alone |
| Sphere-constrained Bézier | same + \(\|p\|=R\) | spherical circle | “Bézier on the sphere” |
| Torus geodesic / \((p,q)\) curve | \((R,r)\), winding | circle if unknot class + \(E_\kappa\) or \(E_{r/R}\); else singular \(N\)-cover | See torus section; not geodesic length alone |

### Ellipse / ellipsoid section

In a local plane (or null plane of a coarse axis), fit ellipse parameters \(\theta = (c, R_a, R_b, \psi)\). Asphericity of the section

\[
E^{\mathrm{asp}}
  \;=\;
  \frac{(\lambda_1-\lambda_2)^2}{(\lambda_1+\lambda_2)^2}
\]

(or axis ratio \(|R_a-R_b|\)) vanishes on the circle. Same story as [`convex_medial_axis_energy.md`](convex_medial_axis_energy.md): **frame from geometry, descend circularity on the section** — do not promote \((c,R)\) to free fit parameters inside the descent energy if you want to avoid fit-then-push collapse.

### Twisted ellipse — rotation as the untwist DOF

**Yes, this is a shape (a low-dimensional family),** and it matches the “solve by calculating a rotation” instinct.

**Construction.** Take a planar ellipse with major “spine” axis \(e_a\) (semi-axis \(a\)) and minor axis \(e_b\) (semi-axis \(b\)). **Twist** means: as you travel in the ellipse parameter \(t\in[0,2\pi]\), rotate the out-of-plane lift of the minor direction about \(e_a\) by a cumulative angle \(\alpha(t)\). The simplest closed one-parameter twist is a **constant rate** with total twist \(\alpha\):

\[
\gamma_{a,b,\alpha}(t)
  \;=\;
  c
  \;+\;
  a\cos t\, e_a
  \;+\;
  b\sin t\bigl(
    \cos(\alpha\, t/2\pi)\, e_b
    \;+\;
    \sin(\alpha\, t/2\pi)\, e_a\times e_b
  \bigr).
\]

(Up to rigid frame \((c,e_a,e_b)\).) Special cases:

| \(\alpha\) | Shape |
|------------|--------|
| \(0\) | planar ellipse in \(\mathrm{span}\{e_a,e_b\}\) |
| small | mildly nonplanar “twisted ellipse” |
| \(\pi\) (half-twist) | half-turn of the minor lobe about the spine — silhouette / offset constructions often look **limaçon-like** (related to twisted / conchoidal surfaces; not identically every limaçon) |
| \(2\pi\) | full twist; returns the plane with a spin covering |

**Natural minimizer of the family.** Two coupled descents:

\[
E_{\mathrm{ell}}
  \;=\;
  \underbrace{\tfrac12\Bigl(\frac{a-b}{a+b}\Bigr)^2}_{E_{\mathrm{asp}}\text{: ellipse}\to\mathrm{circle}}
  \;+\;
  \underbrace{\tfrac12\alpha^2}_{E_\alpha\text{: untwist}},
\]

or replace \(E_\alpha\) by Cauchy TP / bending on \(\gamma_{a,b,\alpha}\). Gradient directions are clean:

- \(\partial/\partial\alpha\): **one rotation** that untwists toward the plane (unknot/uncoil in this family),
- \(\partial/\partial(a-b)\): squash to a circle once planar (or jointly).

That is exactly project-then-descend with a **rotation parameter** as the primary unknotting control — same spirit as local ARAP \(R_i\), but here \(R\) is a **global screw about the ellipse spine** instead of a per-hinge Procrustes.

**Fit.** Hierarchical LS of \((c,e_a,e_b,a,b,\alpha)\) to the rod window (nonlinear in \(\alpha\); seed \(\alpha\) from out-of-plane moment of the samples about \(e_a\)). Freeze, descend \(E_{\mathrm{ell}}\) (and/or TP), scatter to vertices; or pull \(\alpha\) with a single axis-angle about the fitted spine and Poisson-reconstruct.

**Why this is attractive vs raw TP.** The DOF that “looks like a knot” in this family is largely \(\alpha\); the minimal object is the circle; GD has an explicit rotation to turn.

**Limaçon link (honest).** Half-twist (\(\alpha=\pi\)) sits in the same zoo as twisted/conchoidal constructions whose profiles include limaçons; minimizing a limaçon \(r=a+b\cos\varphi\) by \(b\to 0\) is the planar/spherical cousin of minimizing \(\alpha\) and \(a-b\) on the twisted ellipse. Treat them as **neighboring families**, not as a theorem that \(\alpha=\pi\) *equals* a limaçon for all \((a,b)\).

### Limaçon on the sphere

Kidambi & Newton classify integrable point-vortex streamline topologies on \(S^2\); limaçon-type streamlines appear as continuous deformations of simpler closed patterns. As a **projection family**, take a spherical polar form in a chart (or Cartesian limaçon in the tangent plane of a pole), e.g.

\[
r(\varphi) \;=\; a + b\cos\varphi,
\]

with frame \(\theta = (a,b,\text{rotation on }S^2)\). Least-squares fit to the closed rod’s spherical samples; the family’s minimizer is \(b\to 0\) (circle). Local reality = limaçon (still an unknot for small \(b\)); descent drives it to the circle.

Nonlinear LS is enough for a first prototype; no need for the full vortex ODE.

### Cubic Hermite / Bézier — answers the “2 points + 2 tangents” question

A space cubic Hermite segment is uniquely determined by endpoints \(p_0,p_1\) and end tangents \(T_0,T_1\) (scale encoded in \(\|T\|\) or handle lengths). So **yes**, it is the natural 3D stencil for “two points and two tangents.”

It does **not** variationally minimize to a circle under length or bending alone. Circles appear only with extra structure, e.g.

- restrict to a plane or sphere, **and**
- equal chord / equal \(\|T\|\), planar \(T_0,T_1\), and a circularity / constant-curvature energy,

or fit a sphere-constrained Bézier and minimize curvature variance / isoperimetric deficit on that spherical curve.

**Conclusion:** Hermite/Bézier = good **local projection stencil** (hierarchical Hermite LS). Circular minimizer = separate energy on the reconstructed curve (or restrict \(\Phi\) to spherical circles).

### Torus knot / torus geodesic — is there \(E\colon \mathrm{Geo}_{\mathrm{torus}}\to\mathrm{circle}\)?

**Yes, but only for the right \((p,q)\), and the energy is on the space curve (or on \((R,r)\)), not “geodesic length on the torus.”**

Standard embedding of a \((p,q)\) curve on a torus of revolution (major \(R\), minor \(r\)):

\[
\gamma_{R,r}(t)
  \;=\;
  \bigl(
    (R + r\cos qt)\cos pt,\;
    (R + r\cos qt)\sin pt,\;
    r\sin qt
  \bigr),
  \qquad t\in[0,2\pi].
\]

(Flat-torus geodesics that close are the same winding classes; on the embedded torus of revolution, closed geodesics are a slightly thinner set, but the **family** \(\{\gamma_{R,r}\}\) is the practical projection target.)

#### Knot type and which winding can collapse

Use the usual torus-of-revolution parametrization: **\(p\)** = trips around the **major** axis (through the hole), **\(q\)** = trips around the **minor** tube.

| Collapse | What happens to \(\gamma_{R,r}\) | Unwinds to a simple circle? |
|----------|----------------------------------|------------------------------|
| \(r\to 0\) (collapse tube) | Curve → circle of radius \(R\), traversed **\(p\) times** | **Only if \(p=1\)**. Large \(q\) disappears with the tube — many minor windings are “solved” by thinning. |
| \(R\to 0\) (collapse hole) | Degenerate / singular (spindle); not a clean unknot path | Not the intended fairing |
| Stay in fixed \((p,q)\) with \(p\ge 2\) | Family cannot change covering number | **Cannot** go \(p\ge 2\to p=1\) inside this family |

So:

- **Minor winding \(q\) large, major \(p=1\)** (unknot, knot-like look): \(E_\kappa\) or \(E_{r/R}\) that drives \(r\downarrow\) is exactly the variational path \(\mathrm{Geo}_{\mathrm{torus}}\to\) **1-winding circle**. This is the attractive candidate.
- **Major winding \(p\ge 2\)**: collapsing the minor axis only yields a **\(p\)-cover** of a circle. Unwinding that cover to a simple circle is a **different problem** — you need another family (ellipse / limaçon / Hermite hierarchy, a \(p\to p-1\) surgery family, or free TP+Sobolev on the polyline). The torus geodesic family is the wrong reduced space for that goal.

**Goal clarification:** “Unwind over-2 major winding to 1” is *not* what torus projection buys you. Torus projection buys you: **absorb excess minor winding by collapsing \(r\)**, landing on a single major circle when \(p=1\).

#### Practical split

1. Fit / detect \((p,q)\). If \(p=1\), project to \(\gamma_{R,r}\) and descend \(E_{r/R}=\tfrac12(r/R)^2\) or bending \(E_\kappa=\int|\kappa|^2\mathrm{d}s\) at fixed length (+ TP barrier).
2. If \(p\ge 2\), do **not** expect torus collapse to finish the job — switch family (or allow \(p\) to change via a separate covering-reduction energy), using torus only as a local “tube” chart for the \(q\)-part.

Scale constraint: never minimize length alone in \((R,r)\) — that collapses to a point. Anchor \(R\) or total length.

---

## Energies on the projected curve

Once \(y = P(x)\) is fixed for the step:

1. **Cauchy TP** on \(y\) — same core as braid (`calc_tangent_point_radius_grad_cauchy`); gives self-avoidance / barrier in the reduced space.
2. **Asphericity / radial variance** in a section of \(y\) — cheaper roundness drive when TP is too stiff (see medial barrel note).
3. **Holonomy / twist** — see next section.

Combine (1) as barrier + (2) or (3) as fairing drive if TP alone still lacks tangential roundness content.

---

## Twist / parallel-transport as a 1D problem

For a closed framed curve, Călugăreanu / White:

\[
\mathrm{Lk} \;=\; \mathrm{Tw} + \mathrm{Wr}.
\]

- \(\mathrm{Lk}\): topological (integer) — ribbon / push-off linking; for a fixed knot+framing class this is the invariant you must **maintain**.
- \(\mathrm{Tw}\): **1D** — accumulated material rotation of the frame about the tangent along the rod,
  \[
  \mathrm{Tw}
    \;=\;
    \frac{1}{2\pi}\int_0^L \omega(s)\,\mathrm{d}s
  \]
  (twist rate \(\omega\) of the material frame relative to a reference).
- \(\mathrm{Wr}\): **3D geometry of the centerline alone** (Gauss double integral) — how much the curve coils through space.

### Converting the 3D knot to a 1D problem

Along arclength, the centerline is described by curvature/torsion \((\kappa(s),\tau(s))\), or equivalently by a **Bishop (rotation-minimizing) frame**: parallel-transport the normal plane so the connection form has **zero twist rate** by construction. Closing the loop leaves a **holonomy** — a net rotation of that frame after one circuit. That holonomy, together with material twist relative to Bishop, is the 1D bookkeeping of “how much the rod had to spin.”

Heuristic you stated: an *optimal* knot (for a fixed topological class) is a **minimum-twist** configuration — only as much parallel-transport / material rotation as the topology requires — while a repulsive energy keeps the tube from crossing.

That heuristic matches White’s identity if you phrase the variational problem as:

\[
\min\; |\mathrm{Tw}|
\quad\text{(or }\tfrac12\mathrm{Tw}^2\text{)}
\quad\text{subject to}
\quad \mathrm{Lk}=k_{\mathrm{target}},
\quad\text{thickness / no self-intersection}.
\]

Because \(\mathrm{Tw}=\mathrm{Lk}-\mathrm{Wr}\), minimizing \(|\mathrm{Tw}|\) at fixed \(\mathrm{Lk}\) pushes

\[
\mathrm{Wr} \;\to\; \mathrm{Lk}:
\]

**topology lives in the centerline’s writhe (3D shape); material twist is driven toward zero.** Soft energy on the rod:

\[
E_{\mathrm{Tw}}
  \;=\;
  \tfrac12\bigl(\mathrm{Tw}(x)-0\bigr)^2
  \;=\;
  \tfrac12\bigl(\mathrm{Lk}-\mathrm{Wr}(x)\bigr)^2
\]

(with \(\mathrm{Lk}\) held by framing / ribbon class, or replaced by a soft target \(k\)). Equivalently, maximize writhe toward the topological budget without paying material spin.

**Role of TP.** Cauchy TP (braid path) is the **thickness / barrier** constraint in that problem — not the sole fairing drive. Without it, a writhe-seeking curve can collapse through itself. Soft stretch / smooth / min-kink still regularize \(\kappa\); \(E_{\mathrm{Tw}}\) (or holonomy-to-Bishop) targets the **global 1D rotational content**.

### Soft holonomy form (same idea, Bishop language)

\[
E_{\mathrm{hol}}
  \;=\;
  \tfrac12\bigl(\mathrm{holonomy}_{\mathrm{Bishop}}(x) - 2\pi m\bigr)^2.
\]

For the round unknot, writhe vanishes and a twist-free framing is compatible with a simple circle — the joint minimum of \(E_{\mathrm{Tw}}\) + mild bending + TP. For nontrivial knot types, you cannot reach \(\mathrm{Tw}=\mathrm{Wr}=0\); the minimizer is “as little material twist as topology allows,” with writhe carrying \(\mathrm{Lk}\).

### What this is / is not

| Claim | Status |
|-------|--------|
| Net twist along the rod is a **1D** scalar (integral of \(\omega\)) | Exact |
| Fixed \(\mathrm{Lk}\), min \(\|\mathrm{Tw}\|\) ⇒ shape absorbs topology as \(\mathrm{Wr}\) | Exact (White) |
| TP / rope thickness as barrier during that flow | Sound; matches braid use of Cauchy TP |
| “Optimal knot = min parallel-transport configuration” | Good **heuristic** for fair / tight ribbons; related to (but not identical to) ideal-knot / Möbius / elastic minimizers |
| Min twist alone replaces TP fairing | No — need thickness + usually bending; twist only places topology in writhe |

### Quaternion / Bishop formulation (“reset every element”)

Yes — that picture makes sense, and it is already close to how rods store frames here (`rod::_update_frames` = Bishop / parallel transport).

**Circle → toroidal spring (the forward problem).** Put twist (or “vortex” / linking energy) into a closed elastic ring: White says excess \(\mathrm{Tw}\) can dump into \(\mathrm{Wr}\). The centerline buckles into a helical / toroidal coil (writhe) so material twist drops — classic twist-induced writhing. Torque-loop style demos in-repo (“helix / writhe dump seam”) are that direction.

**Opposite (this note’s goal).** Start from a coiled / toroidal spring and **minimize material twist** (or forbid storing it), keep thickness with TP, let geometry shed writhe toward a circle when topology allows.

**Rigid material frame vs reset frame.**

| Object | What it stores | Wrap / \(2\pi\) issues |
|--------|----------------|-------------------------|
| Material Cosserat quaternion \(q(s)\) | full SO(3) / Spin(3) along the rod | Phase can wind forever; twist is \(\omega = \mathrm{Im}(\bar q\, q')\cdot T\) |
| **Bishop / PT quaternion** | only the shortest rotation that maps \(T_i\to T_{i+1}\) | **Resets twist each step** — no accumulated material spin in the frame chain |
| Writhe | centerline Gauss integral / tangent indicatrix | Lives in **shape**, not in an ever-growing quaternion angle |

Discrete Bishop step (pure bend):

\[
q_{i\to i+1}
  \;=\;
  \mathrm{quat}\bigl(\mathrm{shortest\ rotation}\; T_i \mapsto T_{i+1}\bigr),
\qquad
u_{i+1}
  \;=\;
  q_{i\to i+1}\, u_i
\]

(no extra \(\exp(\tfrac12\tau\, T)\) factor). That is “**not** maintaining a rigid material frame; reset the rotation at every element.” Writhe never has to “wrap the quaternion infinitely” — multiple helical turns show up as the **path of \(T\) on \(S^2\)** (and in \(\mathrm{Wr}\)), while each hinge only stores a small PT update. Closing the loop: compare \(u_{\mathrm{end}}\) to \(u_{\mathrm{start}}\); the mismatch is **holonomy** (one \(S^1\) angle), not an unbounded phase along the chain.

Material twist, if needed, is a **separate 1D scalar** \(\theta_i\) (or hinge twist strain) on top of Bishop — exactly \(\mathrm{Tw}=\frac1{2\pi}\sum \Delta\theta\), which you can soft-minimize with \(E_{\mathrm{Tw}}\) while TP blocks intersection.

**Toroidal quaternion configuration (practical recipe).**

1. Centerline \(x_i\) on (or near) \(\gamma_{R,r}\).
2. Rebuild \(u_i\) with Bishop PT from tangents (as `_update_frames`).
3. Optional: set hinge twist rest to \(0\) (min-twist ribbon); topology sits in \(\mathrm{Wr}(x)\).
4. Forward buckling: increase prescribed \(\mathrm{Tw}\) or applied \(M_t\) → writhe dump into a spring.
5. Reverse: \(E_{\mathrm{Tw}}\to 0\) + bending + Cauchy TP → shed coil toward a circle when \(p=1\).

### Local rotations → global vertices (without manifold projection)

Two pipelines for “apply a local rotation, update the curve”:

```text
A) Projection family (this note’s default)
   x → Bishop PT frames (ephemeral) → fit Φ(θ) → E on y=P(x)
   → ∇_y E → scatter to x   [no stored Cosserat chain]

B) Rotation / strain reconstruct (ARAP-rod / Cosserat style)
   local R_i (or hinge quaternion) → prescribed edge ê_i = R_i (x_{i+1}-x_i)
   → recover x by least squares / Poisson → optional Bishop rebuild
```

**Pipeline A (your sketch).** Correct: do **not** need stored Cosserat quats; Bishop from \(0..L\) is rebuilt from \(x\); “local rotate” lives in the **manifold parameters** \(\theta\) or in forces on \(y=P(x)\); vertices update by scattering \(\nabla_y E\) (or a soft pull \(x\to y\)), not by integrating a new frame chain.

**Pipeline B — recover vertices from local rotations.**  
After each element gets a local rotation \(R_i\in\mathrm{SO}(3)\) (torque, \(E_{\mathrm{Tw}}\), bend, …), set target edge

\[
\hat e_i \;=\; R_i\,(x_{i+1}-x_i)
\]

(or \(\hat e_i=\ell_i R_i T_i\)). New positions solve

\[
\min_{x'}
\sum_i \bigl| x'_{i+1} - x'_i - \hat e_i \bigr|^2
\quad\text{(+ pin / kill rigid motion)}.
\]

Normal equations on the cycle are a **weighted discrete Laplacian**:

\[
L x' \;=\; b(\hat e),
\qquad
(L x')_i \;=\; 2x'_i - x'_{i-1} - x'_{i+1}
\]

(uniform weights). So yes: **Laplacian / Poisson reconstruct** — the 1D-cycle ARAP global step. Closed-rod compatibility: \(\sum_i\hat e_i=0\) (else best \(L^2\) fit leaves a closure residual).

Naive alternative: integrate \(x_{i+1}=x_i+\hat e_i\) from a seed, then rigid-align / project to close — only exact when \(\hat e\) already closes.

**Frames in B:** store Cosserat quats and update with \(R_i\), *or* discard and rebuild Bishop from \(x'\) after Poisson (`rod::_update_frames`). Writhe in \(x'\); material twist in quats or a hinge scalar.

| | A — project-then-descend | B — local \(R_i\) + Laplacian |
|--|--------------------------|------------------------------|
| Primary DOF | \(x\); \(\theta\) of \(\Phi\) ephemeral | \(x\) and/or hinge \(R_i\) |
| Cosserat quats | not stored; Bishop from \(x\) | optional; or rebuild after Poisson |
| Global update | scatter \(\nabla_y E\) | \(L x'=b(\hat e)\) |
| Natural for | ellipse / limaçon / torus family | elastic twist→writhe dump |

Hepworth rod quaternion + bend/twist ≈ **B**. Manifold \(P\) in this note ≈ **A**. Hybrid: \(R_i\) from \(E_{\mathrm{Tw}}\), Poisson to \(x'\), optional \(P(x')\) for roundness.

### What this is called

There isn’t one exotic name for “toroidal quaternion turbulence,” but the pieces already have standard names:

| Informal picture | Standard name |
|------------------|---------------|
| Closed curve + continuous frame / ribbon | **Framed knot** (or **ribbon**) |
| Integer topological invariant of curve + push-off | **Linking number** \(\mathrm{Lk}\in\mathbb{Z}\) |
| Real “coiling of the centerline through space” | **Writhe** \(\mathrm{Wr}\in\mathbb{R}\) |
| Real accumulated material spin along the rod | **Twist** \(\mathrm{Tw}\in\mathbb{R}\) |
| White / Călugăreanu | \(\mathrm{Lk}=\mathrm{Tw}+\mathrm{Wr}\) |
| Twist-free frame rebuilt by shortest \(T_i\to T_{i+1}\) | **Bishop frame** / **parallel transport (PT) frame** / **rotation-minimizing frame** |
| Net mismatch of PT frame after one loop | **Holonomy** (a **geometric phase** / **anholonomy**) — a real angle, not forced to \(2\pi\mathbb{Z}\) |
| Full rigid-section Cosserat orientation along \(s\) | **Cosserat / Kirchhoff directors** (quaternions common) |
| \(2\pi\) vs \(4\pi\) / sign flip in quaternion lift | **Spin cover** \(\mathrm{Spin}(3)\to\mathrm{SO}(3)\) (belt trick) |

**Why non-integer multiples of \(\pi\) show up.**  
\(\mathrm{Lk}\) is integer; \(\mathrm{Tw}\) and \(\mathrm{Wr}\) are **real**. A closed curve can store a non-integer amount of geometric coiling as writhe while the ribbon’s linking stays integral — that is exactly the “internal” budget that is *not* a Cosserat twist angle wrapping forever, and *not* “PT has zero content.” PT kills *local* twist rate; **closing** still leaves a real holonomy, and the centerline still has real writhe. Quaternion \(2\pi\sim -1\) adds another half-turn bookkeeping layer; it is not the same as writhe, but it is why frame phase feels “toroidal / π-odd.”

Colloquial “internal turbulence” ≈ **writhe + frame holonomy** (geometric phase in the framed curve). The theorem that ties them is **Călugăreanu–White–Fuller**.

Dennis & Hannay (2005) sharpen the geometry: \(\mathrm{Wr}\) is the direction-average of signed axis self-crossings; \(2\,\mathrm{Tw}\) is the direction-average of signed *local* (ribbon edge-on) crossings. They also define the **writhe framing** — a canonical ribbon with \(\mathrm{Tw}=-\mathrm{Wr}\) so \(\mathrm{Lk}=0\). That is the mathematical name for the “twist cancels writhe / no net link” ribbon we want as a min-twist reference (alongside Bishop PT for the discrete frame chain).

---

## Spine + variable circle (canal) and “dual” manifold

**Your picture.** A thick closed knot ≈ **convolution / sweep of a circle of radius \(r(s)\) along a closed spine \(\gamma(s)\)** — a **canal surface** / pipe surface (torus = circular spine + constant \(r\)).

\[
S(\gamma,r)
  \;=\;
  \bigl\{\,
    \gamma(s) + r(s)\bigl(\cos\theta\, N(s) + \sin\theta\, B(s)\bigr)
  \,\bigr\}.
\]

So the “manifold” attached to the knot is really a **pair** \((\gamma, r)\), or the solid tube they bound — not only the 1D curve.

### What “dual” can mean here (no classical dual knot)

Knot theory does not give a unique “dual knot” for a single knot the way planar graphs have duals. Useful dualities for minimization:

| Dual partner | Object | Minimize what? |
|--------------|--------|----------------|
| Spine \(\gamma\) | 1D core | length, bending, \(E_{\mathrm{Tw}}\), TP on \(\gamma\) |
| Radius \(r(s)\) | thickness field | \(\int r\), variance of \(r\), or fix \(r\ge r_{\min}\) |
| Tube surface \(S\) | 2D canal | Willmore / elastic surface energy; mesh ARAP on \(S\) |
| Seifert surface \(\Sigma\) | spanning 2-manifold with \(\partial\Sigma=\gamma\) | area(\(\Sigma\)) (Seifert energy) — topology carrier |
| Tube solid | 3-volume | volume / ropelength duals |

**Geodesic on a manifold.** Two precise readings:

1. \(\gamma\) is a **geodesic of a torus / canal metric** (already covered: torus \((p,q)\) family).
2. \(\gamma\) is a critical point of an energy on the space of closed curves (length with thickness = **ideal knot** / ropelength) — “geodesic” in a Finsler-type configuration metric, not a surface geodesic.

### Minimization in \((\gamma,r)\) space

**Ropelength / ideal knot (standard):**

\[
\min_\gamma\; \mathrm{Length}(\gamma)
\quad\text{subject to}\quad
\mathrm{reach}(\gamma) \ge r_{\min}
\]

(thickness = inscribed-ball / self-distance constraint). Equivalent spirit to **Cauchy TP as soft thickness** + length or bending on the spine. The “manifold” (tube) is the constraint set, not a separate free surface.

**Variable-radius canal energy (your convolution view):**

\[
E(\gamma,r)
  \;=\;
  \underbrace{\int |\kappa_\gamma|^2\,\mathrm{d}s}_{E_\kappa(\gamma)}
  \;+\;
  \underbrace{\tfrac12\int (r-r_\star)^2\,\mathrm{d}s}_{E_r}
  \;+\;
  E_{\mathrm{TP}}(\gamma;\,r)
  \;+\;
  \text{optional } E_{\mathrm{Tw}}(\gamma).
\]

- Collapse “knotty thickness variation” by driving \(r\to r_\star\) (or \(r\) small on minor coils when \(p=1\)).
- Fair the spine with \(E_\kappa\) / ARAP-rod Laplacian as in pipeline B.
- TP uses local \(r(s)\) as the forbidden approach distance (variable-radius barrier).

**ARAP in this space — where it lives**

| ARAP on… | Role |
|----------|------|
| Spine polyline (1D cycle) | local \(R_i\) + Laplacian → fair \(\gamma\) (previous section) |
| Canal surface mesh \(S\) | shell ARAP / Willmore-like fairing of the tube skin |
| Not natural | ARAP on “the knot” as an abstract group element |

Minimizing **only** the tube surface (ARAP on \(S\)) without controlling \(\gamma\) can inflate/slide the skin while the knot type of the core stays; usually you minimize **spine + thickness**, and treat \(S\) as derived.

**Seifert dual (optional heavier machinery):** minimize area of a spanning surface \(\Sigma\) with \(\partial\Sigma=\gamma\). Drives \(\gamma\) toward unknots that bound small disks (round unknot). Real knot types have a positive area lower bound (genus); you cannot area-minimize to a flat disk without changing knot type.

### Practical takeaway

- Yes: model the thick knot as **spine \(\gamma\) + radius \(r(s)\)** (canal / convolution).
- “Minimize the manifold” ≈ minimize **energies on \((\gamma,r)\)** (and maybe Seifert area), with TP as thickness — not a mysterious dual knot.
- **ARAP** applies cleanly to \(\gamma\) (1D) or to the tube mesh \(S\) (2D); the joint control variables are still \((\gamma,r)\).
- Links: torus section (\(r\) collapse for \(p=1\)); medial/barrel notes (circle along an axis); Crane TP as soft reach constraint.

---

## Descent protocol (like braid TP freeze)

1. Hierarchical LS → \(\theta_\ell^\star\), reconstruct \(y = P(x)\).
2. Freeze \(P\).
3. Evaluate \(\mathcal{E}(y)\) (Cauchy TP and/or \(E^{\mathrm{asp}}\) and/or \(E_{\mathrm{hol}}\)).
4. Scatter \(-\nabla_y\mathcal{E}\) onto the rod; run CCD / velocity block as in [`braid_planar_demo`](../include/gaudi/duchamp/braid_planar_demo.hpp).
5. Remesh / dynamic step; repeat.

Do **not** differentiate through coarse SVD/LS on the first pass.

---

## Recommended first experiments (later; not this doc)

1. **Ellipse hierarchy only** — planar ellipse LS on the sphere-rod polyline; force \(-\nabla E^{\mathrm{asp}}\) in the section; no TP.
2. **Ellipse + Cauchy TP on \(y\)** — freeze \(P\); compare \(E_{\mathrm{TP}}\) and roundness to the current braid Cauchy + CCD path.
3. **Limaçon on \(S^2\)** — same hierarchy with \(r=a+b\cos\varphi\); watch \(b\to 0\).
4. **Holonomy + TP barrier** — \(E_{\mathrm{hol}}\) with existing Cauchy TP as collision/self-avoidance only.

---

## Contrast summary

| Approach | Mechanism | Failure mode |
|----------|-----------|--------------|
| Bare Cauchy TP on polyline | \(L^2\)-ish contact forces | Normal-dominated; slow roundness |
| Sobolev / fractional TP (Crane) | change metric on variations | Needs good \(s\); still no explicit circle family |
| Fit cylinder/circle then push residual | model residual | Fit noise ↔ force fight |
| Hierarchical project \(P\) + descend on \(y\) | reduced space whose min is a circle | Bad family / over-coarse \(\sigma\) kills detail |
| Holonomy + TP barrier | global twist + local avoidance | Framing / \(k\) choice; needs stable frames |

---

## Open implementation notes

- Start with freeze-\(P\) and ellipse or limaçon; Hermite only as a stencil, not as the circularity objective.
- Keep Hepworth stretch/smooth; add projection forces as soft extras.
- Kidambi–Newton PDF is not yet in [`reference/`](../reference/); cite by title + DOI until fetched into MANIFEST.
- BH: same integrators as TP / MLS; projection fits are ordinary weighted reductions at each \(\sigma_\ell\).
