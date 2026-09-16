# Convex-weighted medial axis from normals

Design notes for recovering a local **center** \(c\), **axis** \(a\), and optionally **radius** \(\rho\) from surface samples \((p_j, N_j)\), grounded in an energy / supporting-halfspace argument. Intended as an alternative (or prior) to post-hoc cylinder MLS + \(\tfrac12(\rho-r)^2\) barreling.

Related code / ideas: pierce cylinder medial, convexity gates in Calder MLS, solid-angle / winding via hierarchical sums, trunk spoke geometry.

---

## Setup

At a query sample \(i\) (face or vertex):

- position \(p_i\), unit normal \(N_i\) (inward for the formulas below; flip signs if the mesh stores outward normals),
- neighbors \(j\) with \((p_j, N_j)\),
- displacement \(d_{ij} = p_j - p_i\).

The classical medial / MAT picture: centers of maximal inscribed balls lie where **spokes** (lines \(p - t N\)) meet. On a tube, normals are orthogonal to the ruling direction \(a\).

---

## Convexity weights

The medial axis of a region is determined by **supporting half-spaces** of a locally convex body. Integrate only neighbors that participate in that support (same spirit as cylinder MLS convexity gates).

Example weight (product of gates × length scale):

\[
w_{ij}
  \;\propto\;
  \max\bigl(0,\; -\,d_{ij}\cdot N_i\bigr)
  \;\cdot\;
  \max\bigl(0,\; \widehat{d_{ij}}\cdot N_j\bigr)
  \;\cdot\;
  G(d_{ij};\,\sigma).
\]

- First factor: \(j\) lies on the inward side of the tangent plane at \(i\).
- Second: \(N_j\) faces toward the segment (local convex agreement).
- \(G\): Gaussian / MLS bandwidth (edge length or radius units depending on phase).

Green / divergence intuition: only the **convex envelope** of samples contributes to the inscribed-ball geometry; occluded or reflex neighbors are suppressed.

---

## Phase A — axis from normals

### Energy

On a cylinder (or tube), \(N_j \perp a\). Penalize alignment of normals with a candidate axis:

\[
E(a)
  \;=\;
  \sum_j w_{ij}\,(a\cdot N_j)^2
  \;=\;
  a^{\mathsf T}
  \underbrace{\Bigl(\sum_j w_{ij}\, N_j N_j^{\mathsf T}\Bigr)}_{M_i}
  a.
\]

### Stationarity

Minimize \(E(a)\) subject to \(|a|=1\):

\[
a^\star
  \;=\;
  \arg\min_{|a|=1} E(a)
  \;=\;
  \text{smallest eigenvector of } M_i.
\]

Equivalently: SVD / eigendecomposition of \(M_i\); axis ↔ singular vector for the smallest singular value.

Optional coupling terms of the form \(N_i N_j^{\mathsf T}+N_j N_i^{\mathsf T}\) can be added; the minimal energy above is already well-posed.

### Geometric reading

\(M_i\) is a **covariance of convex-weighted normals**. Large eigenvalues span the spoke / radial plane; the near-null direction is the ruling (medial tangent).

---

## Phase B — center (and radius) from half-spaces

Each sample defines a supporting plane. An inscribed ball of center \(c\) and radius \(\rho\) satisfies

\[
(p_j - c)\cdot N_j \;\ge\; \rho
\qquad\text{for all }j\text{ in the convex neighborhood}
\]

(with \(N\) inward). The medial ball is the LP

\[
\max_{c,\,\rho}\;\rho
\quad\text{subject to}\quad
(p_j - c)\cdot N_j \;\ge\; \rho.
\]

At optimality, active constraints are equalities: spokes \(c = p_j - \rho N_j\).

### Two-plane bisector

For two planes \((p_i,N_i)\), \((p_j,N_j)\), points equidistant to both lie on the **angle-bisector plane** (the “half-plane between the two half-planes”). Local MLS / averaging of bisectors, or a small weighted max-\(\rho\) solve, yields \(c\) (and \(\rho\) as the dual / slack).

### Spoke form

Once \(\rho\) is known (or taken from pierce half-travel),

\[
c_i \;=\; p_i - \rho_i\, N_i.
\]

Neighbor agreement \(\|c_i - c_j\|^2\) is a discrete **spoke / cut-locus** energy (medial coherence without an explicit volume grid).

---

## Single-phase vs two-phase

| Mode | Procedure |
|------|-----------|
| **Single-phase** | Build \(w_{ij}\), form \(M_i\), take \(a^\star\), then one small inscribed-ball / bisector fit for \((c,\rho)\). |
| **Two-phase** | Fair / average \(a\) over a neighborhood first (e.g. SVD of summed \(M\)), then fit \((c,\rho)\mid a\). |

Two-phase helps when per-face axes flicker; single-phase is simpler for debugging.

---

## Optional — winding / solid angle as an axis score

Full winding about an interior point is \(\sim\pm 1\) and does **not** encode an axis. With a **guess** \((c,a)\):

1. Work in the plane \(\Pi\) through \(c\) with normal \(a\).
2. Project a ring (or local surface patch) into \(\Pi\).
3. Evaluate planar winding \(W(a;c)\) of that curve about the projected center.

For a clean tube section, \(|W|\approx 1\). A tilted axis pinches the projection and destabilizes \(W\).

**Axis refinement:**

\[
a^{\star\star}
  \;=\;
  \arg\max_{|a|=1}
  \;
  \mathrm{coherence}\bigl(W(a;c)\bigr)
\]

e.g. maximize proximity of \(|W|\) to \(1\), or minimize variance of local solid angle around the ring.

This is harmonic in \(\Pi\) (argument / solid angle ↔ conjugate of a 2D harmonic potential). Hierarchical double-layer / solid-angle sums (Barnes–Hut) evaluate the score; they do not “cancel” the enclosed winding so much as **measure section quality** for a candidate ruling.

Projected “winding along \(x,y,z\)” is less natural than this **scalar functional on \(a\)**.

---

## Relation to distance, divergence, and cut locus

Let \(d(x)=\mathrm{dist}(x,\partial\Omega)\) inside the solid. Then \(|\nabla d|=1\) a.e., integral curves of \(\nabla d\) are spokes, and

\[
\Delta d \;=\; \mathrm{div}(\nabla d)
\]

is smooth on parallel surfaces and **blows up / concentrates** at the cut locus (medial axis). The convex-weighted normal energy is a **surface-only surrogate** for “normals should meet on a coherent medial structure,” without solving eikonal in the volume.

True eikonal is nonlinear and not a Barnes–Hut sum. Harmonic / Yukawa / solid-angle fields *are* BH-native and can supply soft ridges related to (but not identical to) the MAT.

---

## Soft medial-field energies — colocate centers, align axes

Local SVD / inscribed-ball gives a raw field \((c_i, a_i)\) per sample. Fair that field with **soft** energies (forces / Hepworth residuals), not hard equality. Same freeze pattern: rebuild \((c,a)\) from geometry (or treat them as auxiliary variables), then take gradients.

Notation: undirected neighbors \(i\sim j\) with Gaussian × convex weights \(w_{ij}\). Axes are **unsigned** (\(a\sim -a\)).

### Center colocation (Dirichlet on \(c\))

\[
E^{\mathrm{cen}}
  \;=\;
  \tfrac12\sum_{i\sim j} w_{ij}\,|c_i - c_j|^2.
\]

**Soft push:** \(\partial E/\partial c_i = \sum_j w_{ij}(c_i-c_j)\). Pull neighboring centers together.

**Onto the surface** (spoke map \(c = p - \rho N\), \(\rho\) frozen for the step):

\[
\nabla_{p_i} E^{\mathrm{cen}}
  \;=\;
  \sum_j w_{ij}\,(c_i - c_j)
  \quad\text{(same vector; spoke translation moves }c\text{ with }p\text{)}.
\]

Optional: only the part that preserves the ball to first order, or project the force onto \(N_i\) if the surface may only move along spokes.

### Axis alignment — projector Dirichlet (ARAP for line fields)

Do **not** use \(|a_i-a_j|^2\) (sign flips). Use projectors \(P_i = a_i a_i^{\mathsf T}\in\mathrm{Sym}\), the natural \(\mathbb{RP}^{2}\) / line-field metric:

\[
E^{\mathrm{ax}}
  \;=\;
  \tfrac12\sum_{i\sim j} w_{ij}\,
  \|a_i a_i^{\mathsf T} - a_j a_j^{\mathsf T}\|_F^2
  \;=\;
  \sum_{i\sim j} w_{ij}\,\bigl(1 - (a_i\cdot a_j)^2\bigr).
\]

Vanishes iff neighboring axes are parallel (either orientation). This is the discrete Dirichlet energy on unsigned directions — the same object as ARAP when the “local frame” is only required to match along one axis (Chao / spokes-and-rims style, reduced to a line field; cf. `chao2010-geometric-elastic-deformations.pdf`, cubic stylization’s orientation locals).

**Gradient on the sphere** \(|a_i|=1\) (treat neighbors frozen, or simultaneous):

\[
\nabla_{a_i}^{\mathrm{tan}} E^{\mathrm{ax}}
  \;=\;
  (I - a_i a_i^{\mathsf T})
  \Bigl(
    -2\sum_j w_{ij}\,(a_i\cdot a_j)\,a_j
  \Bigr).
\]

Descent rotates \(a_i\) toward the weighted neighbor mean in the unsigned sense (equivalent to the leading eigenvector of \(\sum_j w_{ij}\,a_j a_j^{\mathsf T}\)).

**ARAP / Procrustes local–global variant** (often stabler than raw gradient on \(a\)):

1. **Local:** for each \(i\), set
   \[
   a_i^\star
     \;=\;
     \text{largest eigenvector of }
     \sum_j w_{ij}\, a_j a_j^{\mathsf T}
   \]
   (or blend with the data term from Phase A: \(\sum w N N^{\mathsf T}\) nullspace vs neighbor projectors).
2. **Global:** if \(a\) is tied to the mesh through normals / spokes, one Poisson / ARAP step on positions; if \(a\) is an auxiliary field, just replace \(a_i\leftarrow a_i^\star\) and re-orthonormalize.

That local step **is** ARAP for a single axis: Procrustes matching of projectors, same spirit as Hepworth `arap` / normal-align locals.

### Spine coupling — centers lie along the axis

Colocation alone collapses centers to a point; alignment alone ignores where the axis sits. Couple them so \(\Delta c\) is parallel to \(a\):

\[
E^{\mathrm{spine}}
  \;=\;
  \tfrac12\sum_{i\sim j} w_{ij}\,
  \bigl|(I - \bar a_{ij}\bar a_{ij}^{\mathsf T})(c_i - c_j)\bigr|^2,
\qquad
\bar a_{ij}
  \;=\;
  \frac{a_i\pm a_j}{|a_i\pm a_j|}
\]

(pick the sign that maximizes \(a_i\cdot(\pm a_j)\)). Vanishes iff neighboring centers differ only along the shared ruling — a soft “medial space curve” with tangent \(a\).

**Joint soft medial energy:**

\[
E^{\mathrm{med}}
  \;=\;
  \alpha\, E^{\mathrm{spine}}
  \;+\;
  \beta\, E^{\mathrm{ax}}
  \;+\;
  \gamma\, E^{\mathrm{cen}}_{\perp}
\]

where \(E^{\mathrm{cen}}_{\perp}\) can be omitted if \(E^{\mathrm{spine}}\) already handles transverse colocation; keep a tiny isotropic \(E^{\mathrm{cen}}\) only if you want mild longitudinal smoothing of arc-length spacing.

Gradients:

- \(\partial E^{\mathrm{spine}}/\partial c_i\): transverse spring — pull \(c_i\) toward the axis through \(c_j\) (and vice versa).
- \(\partial E^{\mathrm{spine}}/\partial a\): weak; usually freeze \(a\) from \(E^{\mathrm{ax}}\) / SVD, then push only on \(c\) (or on \(p\) via spokes).

### Data attachment (keep recovery honest)

Soft fairing without a data term drifts. Attach to the Phase A/B fit \((c_i^0, a_i^0)\):

\[
E^{\mathrm{data}}
  \;=\;
  \tfrac12\sum_i \mu_i\,|c_i - c_i^0|^2
  \;+\;
  \tfrac12\sum_i \nu_i\,\bigl(1 - (a_i\cdot a_i^0)^2\bigr).
\]

Full field step: minimize \(E^{\mathrm{med}}+E^{\mathrm{data}}\) over \((c,a)\) (or local–global ARAP on \(a\), Jacobi/Gauss–Seidel on \(c\)), then map residuals to surface forces through spokes / normals.

### What to implement first

| Term | Role | Gradient / update |
|------|------|-------------------|
| \(E^{\mathrm{ax}}\) | Align axes | Projector Dirichlet or local eig of \(\sum w aa^{\mathsf T}\) |
| \(E^{\mathrm{spine}}\) | Colocate in the cross-section | Transverse springs on \(c\); push \(p\) along spokes |
| \(E^{\mathrm{data}}\) | Stay near SVD / ball fit | Soft anchors |
| \(E^{\mathrm{var}}\) | Round the ring | After \((c,a)\) faired |

Axis alignment really is the ARAP/line-field piece; center colocation is ordinary soft Dirichlet, refined by spine so you do not squash the trunk lengthwise.

---

## Barreling / circularity — differentiable energies

With \((c,a)\) held **fixed** for one force step (same local-global pattern as TP: rebuild kernels / fields from current geometry, then take \(\nabla\) w.r.t. positions only), circularity is a smooth energy on the **projected ring** in \(\Pi\perp a\).

Define the planar foot and radius of sample \(j\) about query \(i\)'s medial frame:

\[
q_j \;=\; (I - a a^{\mathsf T})(p_j - c),
\qquad
\rho_j \;=\; |q_j|,
\qquad
\hat\rho_j \;=\; q_j / \rho_j.
\]

Use the same Gaussian × convex weights \(w_{ij}\) that isolated the nearby axis. Restrict the sum to the “ball” in center / surface space so the energy is **local to this trunk segment**.

### Recommended first try: radial variance (TP-style local energy)

\[
\bar\rho_i
  \;=\;
  \frac{\sum_j w_{ij}\,\rho_j}{\sum_j w_{ij}},
\qquad
E^{\mathrm{var}}_i
  \;=\;
  \tfrac12 \sum_j w_{ij}\,(\rho_j - \bar\rho_i)^2.
\]

**No fitted cylinder \(r\)** — the target radius is the weighted mean of the local section. Vanishes iff all weighted samples share one radius (circle about \(c\) in \(\Pi\)).

**Gradient** (treat \(w_{ij}\), \(c\), \(a\) as frozen; same freeze as TP kernels):

\[
\nabla_{p_k} E^{\mathrm{var}}_i
  \;=\;
  w_{ik}\,(\rho_k - \bar\rho_i)\,\hat\rho_k
  \;-\;
  \Bigl(\sum_j w_{ij}(\rho_j-\bar\rho_i)\Bigr)
  \frac{w_{ik}}{\sum_\ell w_{i\ell}}\,\hat\rho_k.
\]

The second term is the mean-shift correction from \(\partial\bar\rho/\partial p_k\). A common TP-like approximation drops it (treat \(\bar\rho\) as constant for the force), leaving the spoke force

\[
f_k \;\approx\; -\,w_{ik}\,(\rho_k - \bar\rho_i)\,\hat\rho_k
\]

(sign: descent). Optional surface clamp: \(f \leftarrow (f\cdot N)N\).

Per-vertex energy: evaluate \(E^{\mathrm{var}}\) in the convex Gaussian neighborhood of that vertex (or its incident faces), accumulate forces on the same samples — analogous to a local TP stencil, not a global double-integral.

### One traversal for \((c,a,r)\) and \(\nabla E\)?

**Dependency order** (current geometry, no lag):

1. Weights \(w_{ij}\) from \((p,N)\).
2. \(a\) from \(M=\sum w\,NN^{\mathsf T}\) (3×3 eig **after** the sum).
3. \(c\) (and optional MLS \(r\)) need \(a\) (cylinder rows / projection to \(\Pi\)), or a separate inscribed-ball solve.
4. \(\rho_j=|(I-aa^{\mathsf T})(p_j-c)|\) and \(\bar\rho\) need \((c,a)\).
5. \(\nabla_p E^{\mathrm{var}}\) needs \(\bar\rho\) (and frozen \(c,a\)).

So a **streaming** pass cannot emit the force until \((c,a)\) exist. You cannot get a consistent current-frame \((c,a,r)\) and \(\nabla E\) from one pure online reduction without either storing the neighbor list or lagging the frame.

**Practical “one traversal” (TP / ARAP freeze):** use **lagged** medial frame \((c^{-},a^{-})\) from the previous step.

Single gather over \(j\):

\[
\begin{aligned}
M &\mathrel{+}= w\,N_j N_j^{\mathsf T}, \\
\text{(optional cylinder moments for next }c,r\text{)}, \\
\rho_j &= |(I-a^{-}(a^{-})^{\mathsf T})(p_j-c^{-})|, \\
W &\mathrel{+}= w,\quad S \mathrel{+}= w\,\rho_j, \\
f_i &\mathrel{+}= \text{spoke residual vs running or deferred }\bar\rho.
\end{aligned}
\]

After the gather: \(\bar\rho=S/W\), finish \(f_i\propto -(\rho_i-\bar\rho)\hat\rho_i\), then \(a\leftarrow\mathrm{eigvec}_{\min}(M)\), refresh \(c\) (spoke or moments). Next step’s energy uses that frame. Same pattern as TP kernels and ARAP local rotations: **field from \(x^{n}\), force on \(x^{n}\), field update for \(n+1\)**.

**One gather, two touches of a small list** (no lag, still one tree query):

1. BH / knn once → buffer \(\{(p_j,N_j,w_j)\}\).
2. Reduce \(M\) → \(a\); reduce moments → \(c,r\) (or \(\bar\rho\) after \(c,a\)).
3. Replay buffer → \(E\), \(\nabla_p E\).

Cost is one spatial traversal; the replay is \(O(k)\) in the stencil size.

**Joint energy, alternating (one energy, still two phases):**

\[
E(p;\,c,a,r)
  \;=\;
  \underbrace{\tfrac12\sum_j w_j\bigl(|(I-aa^{\mathsf T})(p_j-c)| - r\bigr)^2}_{E^{\mathrm{cyl}}\;\sim\;E^{\mathrm{var}}\text{ if }r=\bar\rho}
  \;+\;
  \underbrace{\tfrac12\,a^{\mathsf T} M(p)\,a}_{E^{\mathrm{ax\text{-}data}}}
  \;+\;\text{optional spine / data}.
\]

- Fix \(p\): one stencil reduction → \(M\), then \(a\); then \(c,r=\arg\min E^{\mathrm{cyl}}\) (or \(r=\bar\rho\)).
- Fix \((c,a,r)\): same stencil (replay or lag) → \(\nabla_p E^{\mathrm{cyl}}\).

You do **not** need to differentiate through the 3×3 eig for a usable barrel force. Differentiating through \(a(p)\) and \(c(p)\) is possible but noisy and unnecessary if the soft spine / data terms already fair the field.

**Bottom line:** one spatial traversal yes (gather once); current-frame \((c,a,r)\) and \(\nabla E\) need a replay or a lagged frame — same as TP. Prefer lag for growth timesteps; prefer gather+replay for debugging a single static evaluation.

### Pairwise form (even closer to TP)

\[
E^{\mathrm{pair}}_i
  \;=\;
  \tfrac14\sum_{j,k} w_{ij}\,w_{ik}\,(\rho_j - \rho_k)^2
  \;=\;
  \Bigl(\sum_j w_{ij}\Bigr)\, E^{\mathrm{var}}_i.
\]

Same critical points as variance; gradient is pairwise “radius Dirichlet.” Implement whichever is cheaper; they are equivalent up to the total-weight factor.

### No-mean alternative: planar asphericity (SVD-native)

Second-moment tensor of the projected cloud:

\[
Q_i
  \;=\;
  \sum_j w_{ij}\, q_j q_j^{\mathsf T}
  \qquad\text{(lives in }\Pi\text{; }Qa=0\text{)}.
\]

For a circle about \(c\), \(Q_i\) is **isotropic in \(\Pi\)**. Let \(\lambda_1\ge\lambda_2\) be the nonzero eigenvalues of \(Q_i\) (or eigenvalues of the \(2\times2\) reduction in an orthonormal frame of \(\Pi\)). Asphericity energy:

\[
E^{\mathrm{asp}}_i
  \;=\;
  \frac{(\lambda_1 - \lambda_2)^2}{(\lambda_1 + \lambda_2)^2}
  \;=\;
  1 - \frac{4\det_{\Pi} Q_i}{(\mathrm{tr}\,Q_i)^2}.
\]

Vanishes iff the weighted gyration ellipse is a circle. Differentiable in \(q_j\) (hence \(p_j\)) wherever \(\mathrm{tr}\,Q>0\). Gradient via \(\partial Q/\partial p_k = w_{ik}(\hat\rho_k q_k^{\mathsf T} + q_k \hat\rho_k^{\mathsf T})\) projected to \(\Pi\), then chain rule through \(\det\) and \(\mathrm{tr}\).

This never introduces an explicit \(r\); it only kills elliptical / oval sections. Combine with variance if you also care about absolute radius uniformity vs a mean (asphericity alone allows a thin elliptical shell that is still “round” only in second moment — for dense ring samples they agree in practice).

### Isoperimetric deficit (true circle, harder stencil)

Ordered projected polygon (needs a ring, not a raw MLS cloud):

\[
A = \tfrac12\oint(x\,dy-y\,dx),
\qquad
P = \oint ds,
\qquad
E^{\mathrm{iso}} = \frac{P^2}{4\pi A} - 1 \;\ge\; 0.
\]

Vanishes iff the section is a circle. Differentiable w.r.t. polygon vertices via Green, but requires a coherent loop (winding slice or mesh belt). Prefer after variance/asphericity are working.

### What not to use as the primary descent energy

\[
\tfrac12(\rho_i - r_{\mathrm{MLS}})^2
\]

with \(r\) from a separate cylinder fit: the target and the geometry fight each other when the fit is noisy. Variance / pairwise / asphericity keep the target **inside** the same weighted neighborhood that defined \((c,a)\).

---

## Descent protocol (like TP)

1. Rebuild convex-Gaussian neighborhood; form \(M\), take \(a^\star\); fit / average \(c\) (inscribed ball or spoke) → data \((c^0,a^0)\).
2. Soft-fair the medial field: local–global \(E^{\mathrm{ax}}\) (ARAP/projector) + \(E^{\mathrm{spine}}\) (+ \(E^{\mathrm{data}}\)).
3. Freeze faired \((c,a,w)\).
4. Evaluate \(E^{\mathrm{var}}\) or \(E^{\mathrm{asp}}\) and accumulate \(\nabla_p E\) (spoke directions in \(\Pi\)).
5. Integrate / project (normal clamp, CFL, …).
6. Remesh / topology only after the field is stable.

Do **not** differentiate through the axis SVD on the first pass — same freeze as TP kernels and ARAP local frames.

---

## Pipeline summary

```text
samples (p, N)
    │
    ▼
convex × Gaussian weights w_ij  (isolate nearby axis)
    │
    ├──────────────────────────────┐
    ▼                              ▼
M = Σ w N Nᵀ                 inscribed-ball / spokes
a⁰ = argmin aᵀ M a                (c⁰, ρ̄)
    │                              │
    └──────────┬───────────────────┘
               ▼
        soft fair: E^ax (projector / ARAP) + E^spine (+ E^data)
               ▼
        freeze faired (c, a, w)
               ▼
        E^var or E^asp  →  ∇_p  (spoke forces)
               ▼
        optional later: E^iso on ordered ring / winding slice
               ▼
        soft-slice TP in Π ⊥ a (Gaussian axial window) → ∇_p
```

Implemented as `calder::tangent_point_gradient_axis_slice` and growth `_slice_tp_weight`.

---

## Contrast with fitted cylinder residual

| Approach | Mechanism | Failure mode |
|----------|-----------|--------------|
| MLS \(r\) then \(\tfrac12(\rho-r)^2\) | Fit then push | Fit noise ↔ force fight |
| Convex normal PCA + inscribed \(c\) | Energy on half-spaces | Needs good \(w\) |
| \(E^{\mathrm{ax}}\) + \(E^{\mathrm{spine}}\) | Soft line-field ARAP + transverse center springs | Over-smooth branches if \(w\) too wide |
| \(E^{\mathrm{var}}\) / \(E^{\mathrm{pair}}\) | Local radius Dirichlet about frozen \(c,a\) | Needs stable axis; empty nbhd |
| \(E^{\mathrm{asp}}\) | Isotropic gyration in \(\Pi\) | Thin / sparse samples |
| \(E^{\mathrm{iso}}\) | True circle (Green) | Needs ordered ring |

---

## Open implementation notes

- Sign of \(N\): pick inward consistently with pierce / solid-angle.
- Bandwidth \(\sigma\): edge length for surface MLS; radius units if operating in center space — same ball that isolates the axis should bound circularity **and** medial fairing stencils.
- BH: solid-angle / double-layer for optional \(W(a;c)\); \(Q\), \(\bar\rho\), and projector sums are ordinary weighted reductions.
- Freeze \((c,a,w)\) during \(\nabla E^{\mathrm{var}}\); fair \((c,a)\) softly before circularity forces.
- Axis: prefer \(1-(a_i\cdot a_j)^2\) / local eig of \(\sum w aa^{\mathsf T}\), not \(|a_i-a_j|^2\).
- Start with \(E^{\mathrm{ax}}+E^{\mathrm{spine}}+E^{\mathrm{data}}\), then \(E^{\mathrm{var}}\); add asphericity if ovals persist.
