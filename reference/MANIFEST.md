# Reference paper index

Local PDFs live in this directory. Prefer this manifest for tag/keyword search before opening a PDF.
Page text is cached by [`pdf_text.py`](pdf_text.py) (`pages`, `search`) under `reference/.cache/`.
`python3 reference/pdf_text.py chunk reference/<file>.pdf` writes one markdown file per table-of-contents section, plus `INDEX.md`, in that cache. Open the index, then one section.
Provenance URLs/DOIs: [SOURCES.txt](SOURCES.txt). PDFs may move to external storage later; keep this index in-repo.

## How to add a paper

1. Drop or fetch the PDF into `reference/` (stable filename: `authorYEAR-short-slug.pdf` when renaming).
2. Append an entry below (title, file, url, tags, summary, maps_to).
3. Add URL/DOI line in `SOURCES.txt`.

---

## Barill et al. 2018 — Fast Winding Numbers for Soups and Clouds
- **file:** `fast-winding-numbers-for-soups-and-clouds-siggraph-2018-barill-et-al.pdf`
- **url:** project / SIGGRAPH 2018 (local copy)
- **tags:** `[winding, inside-outside, point-cloud, triangle-soup, hierarchy, asawa, geometry-processing]`
- **summary:** Hierarchical fast evaluation of generalized winding numbers on triangle soups and point clouds for robust inside/outside queries and reconstruction pipelines.
- **maps_to:** Solid queries, containment, mesh cleanup / segmentation prep.

## Jacobson, Kavan, Sorkine-Hornung 2013 — Robust Inside-Outside Segmentation using Generalized Winding Numbers
- **file:** `robust-inside-outside-segmentation-using-generalized-winding-numbers-siggraph-2013-jacobson-et-al.pdf`
- **url:** SIGGRAPH 2013 (local copy)
- **tags:** `[winding, inside-outside, tetrahedral-meshing, segmentation, asawa, geometry-processing]`
- **summary:** Generalized winding number on imperfect meshes (holes, intersections, non-manifold) to segment inside/outside and drive tetrahedralization.
- **maps_to:** Robust volume/inside tests before remesh or collision volumes.

## Kim & Delaney 2013 — Subspace Fluid Re-Simulation
- **file:** `KimDelaney2013.pdf`
- **url:** (local copy; SIGGRAPH / TOG era)
- **tags:** `[fluids, subspace, model-reduction, re-simulation, duchamp]`
- **summary:** Subspace integration for fast re-simulation of novel fluid dynamics with large speedups vs full solve (optional velocity reconstruction).
- **maps_to:** Reduced-order / subspace ideas for time-dependent fields (not shell PD).

## Cui, Langlois, Sen, Kim — Spiral-Spectral Fluid Simulation
- **file:** `SpiralSpectralFluids.pdf`
- **url:** (local copy)
- **tags:** `[fluids, spectral, spiral, radial-domain, sphere, torus, duchamp]`
- **summary:** Fast spectral fluid method on radial domains (discs, spheres, cylinders, ellipses, spheroids, tori) for surfaces and volumes.
- **maps_to:** Spectral / radial domain simulation patterns.

## Liu, Bargteil, O'Brien, Kavan 2013 — Fast Simulation of Mass-Spring Systems
- **file:** `Liu-FSM-2013-11.pdf`
- **url:** (local copy; SIGGRAPH Asia / TOG)
- **tags:** `[mass-spring, implicit-euler, optimization, cloth, hepworth, projection-dynamics]`
- **summary:** Reformulates implicit Euler mass-spring as an optimization; few iterations give real-time cloth comparable to slower Newton.
- **maps_to:** Historical lineage for projective / constraint-based cloth; related spirit to Hepworth PD loops.

## Liu & Jacobson 2019 — Cubic Stylization
- **file:** `liu2019-cubic-stylization.pdf` (arxiv); hi-res: `liu2019-cubic-stylization-hires.pdf`
- **url:** https://arxiv.org/pdf/1910.02926.pdf
- **project:** https://www.dgp.toronto.edu/projects/cubic-stylization/
- **hires:** https://www.dgp.toronto.edu/projects/cubic-stylization/cubicStyle_high.pdf
- **tags:** `[stylization, ARAP, L1-normals, ADMM, local-global, axis-alignment, hepworth, shell, spokes-and-rims]`
- **summary:** Cubify a mesh via ARAP (Chao spokes-and-rims) plus ℓ¹ on rotated vertex normals; local ADMM soft-threshold on `z = R n̂`, global ARAP Poisson; artistic λ / orientation controls.
- **maps_to:** Hepworth `arap` / `arap_normal_align` + `init_cubification`; demo `cubic_stylization`. Hard `N_goal` tip-in = folding-ish mode; paper cubeness = soft ‖z‖₁ ADMM.
- **cites:** Stein et al. 2018a (developability); Chao et al. 2010 (spokes-and-rims).

## Chao, Pinkall, Sanan, Schröder 2010 — A Simple Geometric Model for Elastic Deformations
- **file:** `chao2010-geometric-elastic-deformations.pdf`
- **url:** https://multires.caltech.edu/pubs/GeomElastic.pdf
- **doi:** https://doi.org/10.1145/1833349.1778775
- **tags:** `[ARAP, spokes-and-rims, elasticity, hepworth, shell, geometry-processing]`
- **summary:** Geometric elasticity as distance of the deformation gradient from SO(3); clarifies spoke vs spokes-and-rims ARAP stencils on triangle meshes.
- **maps_to:** Local Procrustes edge sets in `arap` / `arap_normal_align` / `init_arap` (rim edges + cotan at center).

## Crane, Pinkall, Schröder / collaborators — Repulsive Curves
- **file:** `crane2021-repulsive-curves.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/RepulsiveCurves.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/
- **doi:** https://doi.org/10.1145/3439429
- **tags:** `[repulsive, curves, tangent-point, Sobolev, rods, calder, TP, geometry-processing]`
- **summary:** Fractional tangent-point energies for curve fairing / self-avoidance. The Sobolev–Slobodeckij metric is \(A = B + B_0\) (§4.2, discrete §5.2): \(B\) compares tangential derivatives \(\mathcal{D}u\), and \(B_0\) is the low-order term weighted by the tangent-point kernel. The gradient solves \(A g = \mathrm{d}E\).
- **maps_to:** Calder TP rod energies / gradients; Sobolev Jacobi metric in `tp_sobolev_jacobi.hpp` / [`docs/sobolev_B_tangent_filter.md`](../docs/sobolev_B_tangent_filter.md); braid / rod demos that need self-repulsion.
- **cites / related:** Supplemental evaluation `crane2021-repulsive-curves-supplemental.pdf` (which inner product actually untangles). Dennis & Hannay 2005 for ribbon twist / writhe, which this energy does not include.

## Yu, Schumacher, Crane 2021 — Repulsive Curves, supplemental (evaluation)
- **file:** `crane2021-repulsive-curves-supplemental.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/RepulsiveCurves_Supplemental.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveCurves/
- **doi:** https://doi.org/10.1145/3439429
- **tags:** `[repulsive, curves, Sobolev, preconditioner, untangle, knots, optimization, rods, TP]`
- **summary:** TOG supplemental to *Repulsive Curves* (6 pages). Not the derivation of \(B\) and \(B_0\) — that is the main paper, §4.2 and Appendix A. This note compares 19 optimizers on KnotPlot embeddings (128 knot types + 100 trefoils, \(\mathcal{E}^2_{4.5}\)). Inner products tried as drop-in replacements: \(L^2\) (mass), \(H^1\) (stiffness), \(H^2 = LM^{-1}L\), fractional \(H^s\), and a convexified edge-pair Hessian. Hard constraints (fixed barycenter, fixed edge length) vs soft Hencky strain. \(H^s\) projected gradient gets stuck on 6.5% of plotted hard-constraint examples; \(H^1\) fails about 64% (stuck about 47%). Iteration count for \(H^s\) stays flat under refinement; \(H^1\) does not, because the inner-product order no longer matches the energy.
- **maps_to:** Evidence for keeping the fractional metric rather than a plain Laplacian preconditioner on the braid solve. The untangling motion in the main paper’s Figure 7 is a near-tangential slide of the centerline, produced by \(A^{-1}\mathrm{d}E\), not a material-frame twist.

## Crane et al. — Repulsive Surfaces
- **file:** `crane2021-repulsive-surfaces.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveSurfaces/RepulsiveSurfaces.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveSurfaces/
- **tags:** `[repulsive, surfaces, tangent-point, shell, calder, TP, geometry-processing]`
- **summary:** Surface analogue of repulsive / tangent-point energies for self-avoiding surface flows.
- **maps_to:** Shell TP modules; soft contact / untangling-adjacent ideas.

## Gillespie et al. — Harnack Tracing
- **file:** `gillespie-harnack-tracing.pdf`
- **url:** https://markjgillespie.com/Research/harnack-tracing/HarnackTracing.pdf
- **tags:** `[tracing, implicit, contours, asawa, geometry-processing]`
- **summary:** Robust tracing of algebraic / implicit curves and surfaces via Harnack-type guarantees.
- **maps_to:** Contour extraction, level-set / implicit debugging, medial or ridge tracing pipelines.

## Gillespie et al. — Integer Coordinates
- **file:** `gillespie-integer-coordinates.pdf`
- **url:** https://markjgillespie.com/Research/integer-coordinates/integer_coordinates.pdf
- **tags:** `[coordinates, integers, parameterization, asawa, geometry-processing]`
- **summary:** Integer-valued coordinate constructions on meshes (quantized / discrete coordinate frames).
- **maps_to:** Discrete frames, polycube-adjacent orientation, integer constraints in Hepworth-style locals.

## Feng et al. — Points as Tori
- **file:** `feng-points-as-tori.pdf`
- **url:** https://nzfeng.github.io/research/PointsAsTori/PointsAsTori.pdf
- **tags:** `[tori, points, medial, cyclide, albers, calder, geometry-processing]`
- **summary:** Represent / fit toroidal structure to point geometry; links point clouds to torus / cyclide models.
- **maps_to:** Darboux / cyclide / medial fits; albers torus-constrained shape operators.

## Minarčík — Surface Untangling
- **file:** `minarcik-surface-untangling.pdf`
- **url:** https://minarcik.com/SurfaceUntangling.pdf
- **tags:** `[untangle, intersection, shell, remesh, asawa, geometry-processing]`
- **summary:** Untangle self-intersecting or folded triangle meshes toward valid embeddings.
- **maps_to:** Mesh repair before remesh / TP; shell collision and inversion recovery.

## Shen, O'Brien, Shewchuk 2004 — Interpolating and Approximating Implicit Surfaces from Polygon Soup
- **file:** `Shen-IAI-2004-08.pdf`
- **url:** SIGGRAPH 2004 (local copy)
- **tags:** `[implicit-surface, polygon-soup, MLS, smoothing, envelopes, asawa, calder]`
- **summary:** Build interpolating or approximating implicits from polygon soup; smooth small features, simplify topology, simulation envelopes.
- **maps_to:** Implicit / envelope construction from messy surface data.

## Stein, Grinspun, Crane 2018a — Developability of Triangle Meshes
- **file:** `stein2018-developability-triangle-meshes.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/DiscreteDevelopable/paper.pdf
- **arxiv:** https://arxiv.org/pdf/1804.07232.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/DiscreteDevelopable/
- **doi:** https://doi.org/10.1145/3197517.3201303
- **tags:** `[developable, flatten, rulings, seams, vertex-star, shell, asawa, hepworth, paper-folding]`
- **summary:** Discrete developability for triangle meshes (flattenability + straight rulings). Variational flow concentrates Gaussian curvature onto seam curves; energy supported on the vertex star; no explicit patch clustering. Cited by Cubic Stylization as a related geometric flow.
- **maps_to:** Paper/cardboard nets and crease/seam formation; shell energies that drive ruled, flattenable patches. Complements stretch (`triangle_strain` / `edge_strain`) + crease hinges rather than cubic axis-alignment.

## Stein, Wardetzky, Jacobson, Grinspun 2020 — A Simple Discretization of the Vector Dirichlet Energy
- **file:** `stein2020-vector-dirichlet-energy.pdf`
- **url:** https://ddg.math.uni-goettingen.de/pub/vector-dirichlet.pdf
- **doi:** https://doi.org/10.1111/cgf.14070
- **tags:** `[vector-dirichlet, covariant-derivative, Crouzeix-Raviart, direction-fields, asawa, geometry-processing]`
- **summary:** Concise triangle-mesh discretization of the covariant vector Dirichlet energy using Crouzeix–Raviart elements (SGP 2020 / CGF).
- **maps_to:** Direction / vector-field smoothness energies on shells.

## Turk & O'Brien 1999 — Variational Implicit Surfaces
- **file:** `turk1999-variational-implicit-surfaces.pdf`
- **url:** https://repository.gatech.edu/server/api/core/bitstreams/ebf3251d-156b-4cbc-8b49-0317de363ed7/content
- **note:** Georgia Tech tech report GIT-GVU-99-15; see also Turk & O'Brien, *Modelling with implicit surfaces that interpolate* (TOG 2002).
- **tags:** `[implicit-surface, RBF, variational, scattered-data, biharmonic, thin-plate, level-set, sculpting, calder, albers, hepworth]`
- **summary:** Build a smooth 3D implicit whose zero set *interpolates* prescribed surface points, with sign fixed by interior (+) or exterior (−) constraints, by solving a variational scattered-data problem in biharmonic RBFs φ(r)=r³ plus a degree-1 polynomial. One dense symmetric system in the constraint count (not a grid), so no discretization artifacts and no resolution parameter. Surface *normals* are prescribed directly via small offset constraints, which blobbies cannot do. Includes polygon-mesh → implicit conversion and CSG / ray-trace / particle-sculpting uses.
- **maps_to:** RBF field family in `rbf_displacement_basis.hpp` / `local_slab_field_fit.hpp` (`local_rbf_mls_field`) — the paper is the interpolating, normal-constrained variant of what we fit approximately; `implicit_field.hpp` / `level_set_field_bundle.hpp` level-set evaluation; albers normal-constrained fits (quadric / Darboux) as the low-order algebraic counterpart; Harnack tracing for rendering the zero set without a voxel grid.
- **cites / related:** Shen, O'Brien, Shewchuk 2004 (same lineage, *approximating* implicits from polygon soup); Barill et al. 2018 (winding number as the sign source instead of ± constraints).

## Dennis & Hannay 2005 — Geometry of Călugăreanu's theorem
- **file:** `dennis2005-geometry-calugareanu.pdf`
- **url:** https://arxiv.org/pdf/math-ph/0503012
- **doi:** https://doi.org/10.1098/rspa.2005.1527
- **tags:** `[framed-knot, ribbon, linking, twist, writhe, Calugareanu, White, Fuller, geometric-phase, Bishop, rods, duchamp]`
- **summary:** Geometric reading of Călugăreanu–White–Fuller \(\mathrm{Lk}=\mathrm{Tw}+\mathrm{Wr}\) for closed ribbons. Writhe = average signed self-crossings of the axis over projection directions; **twice the twist** = average signed *local* (edge-on) crossings of the two ribbon edges. Constructs the **writhe framing**: a natural ribbon on any closed space curve with \(\mathrm{Tw}=-\mathrm{Wr}\) hence \(\mathrm{Lk}=0\) always — twist exactly cancels writhe.
- **maps_to:** [`docs/knot_project_descend.md`](../docs/knot_project_descend.md) min-twist / Bishop / holonomy story; rod Bishop frames in `asawa/rod/rod.hpp` (`_update_frames`); braid / torque-loop twist–writhe dump; Cauchy TP as thickness while descending \(E_{\mathrm{Tw}}\).
- **cites / related:** Crane repulsive curves (TP barrier); White 1969 / Fuller DNA writhe.

## Peng, Deng, Zhang, Geng, Qin, Liu 2018 — Anderson Acceleration for Geometry Optimization and Physics Simulation
- **file:** `peng2018-anderson-acceleration.pdf`
- **url:** https://arxiv.org/pdf/1805.05715
- **doi:** https://doi.org/10.1145/3197517.3201290
- **tags:** `[anderson, acceleration, local-global, projective-dynamics, ARAP, fixed-point, solver, hepworth]`
- **summary:** Treats a local-global solver (Shape-Up, ARAP, projective dynamics) as a fixed-point iteration `Q_LG = G(Q)`: local projections `P(Q)`, then the global solve. Anderson acceleration takes the current and last `m` iterates and picks the affine combination whose residual `F = G(Q) - Q` is smallest: `theta = argmin |F_k - sum theta_j dF_{k-j}|`, `Q_{k+1} = G_k - sum theta_j dG_{k-j}` (an `m x m` least-squares solve). Safeguard (Algorithm 1): compute the local step at the accelerated iterate; if its energy is not below the previous iterate's, fall back to the plain local-global iterate. The projections from that check are reused as the next global step's input, so the check costs about one distance evaluation per iterate. The default is `m = 5`; past `m = 6` gains are small. They relate it to L-BFGS on projective dynamics (Liu et al. 2017) and Chebyshev acceleration (Wang 2015), with similar per-iteration decrease to L-BFGS.
- **maps_to:** `projection_solver::step` in `hepworth/block/solver.hpp` (projections then the cached LDLT solve are `G`; `gn_energy` is the target energy). Gauss–Newton rows are part of `G`, so a `J` refresh changes the map and should reset the history. As a composable iterate corrector beside the planned modal (low-mode Newton) correction.
- **cites / related:** Anderson 1965; Walker and Ni 2011; Bouaziz et al. 2014 (projective dynamics); Liu, Bouaziz, Kavan 2017 (quasi-Newton projective dynamics); Wang 2015 (Chebyshev semi-iterative PD/PBD).

## Crane, Pinkall, Schröder 2013 — Robust Fairing via Conformal Curvature Flow
- **file:** `crane2013-conformal-willmore-flow.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/ConformalWillmoreFlow/paper.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/ConformalWillmoreFlow/
- **doi:** https://doi.org/10.1145/2461912.2461986
- **tags:** `[willmore, conformal, fairing, curvature-flow, CMC, shell, hepworth]`
- **summary:** Willmore flow on a triangle mesh, integrated in curvature space so a single explicit step can be huge and triangle aspect ratios stay put. Mean curvature is \(H = \tfrac12\langle N, \Delta f\rangle\). The flow direction is the \(L^2\) gradient of Willmore energy on a curvature variable \(\rho\), then projected orthogonal to the linear constraints that keep the result integrable: constant total curvature, exactness against a harmonic-vector-field basis (non-simply-connected surfaces), and a Möbius/inversion constraint that limits area distortion. Positions are recovered by a Dirac eigenvalue problem for new tangents, then a Poisson solve \(\Delta\tilde f = \nabla\cdot\tilde T\). Same loop as their length-preserving curve flow (curvature → tangents → positions). Also builds CMC surfaces.
- **maps_to:** Shell fairing beside `edge_willmore` / cotangent Laplacian in `shell_constraints.hpp`. The dihedral projection zeros hinge angle; this flow changes \(H\) and rebuilds the immersion conformally. Heller & Pedit (arXiv:1705.03217) use this flow for their constrained-Willmore experiments.
- **cites / related:** Bobenko & Schröder 2005 and Wardetzky et al. 2007 (discrete Willmore); Desbrun et al. 1999 implicit fairing; mean-curvature flow as the membrane analogue. Soliman et al. 2021 replace the Möbius penalty with an exact length-cross-ratio constraint.

## Soliman, Chern, Diamanti, Knöppel, Pinkall, Schröder 2021 — Constrained Willmore Surfaces
- **file:** `soliman2021-constrained-willmore-surfaces.pdf`
- **url:** https://www.yousufsoliman.com/projects/download/ConstrainedWillmoreSurfaces.pdf
- **doi:** https://doi.org/10.1145/3450626.3459759
- **tags:** `[willmore, conformal, fairing, CETM, shell, hepworth]`
- **summary:** Minimizes discrete Willmore \(W(f)=\tfrac14\langle Lf,\star Lf\rangle\) on a triangle mesh subject to an exact conformal-class constraint. Two metrics are discretely conformally equivalent iff their length cross ratios agree (Springborn, Pinkall, Schröder 2008): \(\log c = C\log\ell\), and conformal classes are \(\log\ell\) modulo vertex scale factors. The Lagrange multipliers are discrete quadratic differentials, read as edge tension that balances \(\operatorname{grad} W\). The saddle point is integrated by competitive gradient descent in a Sobolev metric, not by projecting an \(L^2\) step back onto the constraint. Positions may be vertex coordinates or differential coordinates, so an abstract metric can be realized with no initial immersion. Point, area, and volume constraints sit in the same Lagrangian. Without the conformal constraint every torus collapses to the Clifford torus and four point constraints on a sphere yield the round sphere.
- **maps_to:** Shell fairing beside Crane, Pinkall, Schröder 2013 and `edge_willmore` in `shell_constraints.hpp`. Crane’s Möbius penalty limits area distortion and can drift in conformal class; this paper holds the class exactly via length cross ratios.
- **cites / related:** Crane et al. 2013 (conformal Willmore flow; they note its conformal-class drift); Springborn, Pinkall, Schröder 2008 (CETM); Canham 1970 / Helfrich 1973 (area–volume vesicles).

## Knöppel, Crane, Pinkall, Schröder 2015 — Stripe Patterns on Surfaces
- **file:** `knoeppel2015-stripe-patterns.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/paper.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/StripePatterns/
- **doi:** https://doi.org/10.1145/2767000
- **tags:** `[stripes, direction-field, eigen, laplacian, shell, hepworth]`
- **summary:** Evenly spaced stripes on a triangle mesh, with singularities inserted so a prescribed orientation and spacing can be met. The pattern is the global minimizer of a convex quadratic energy. That minimizer is the smallest generalized eigenpair \(Ax=\lambda Bx\): \(A\) has the sparsity of the cotangent Laplacian, \(B\) is the mass matrix. They factor \(A\) once and run inverse power iteration (Algorithm 6), a random start times a few backsolves. Texture coordinates are then read off the eigenvector in each face.
- **maps_to:** Eigenpatterns beside `duchamp/arnoldi.hpp` and the projection solve in `hepworth/block/solver.hpp`. The monge demo's edge block (`shell_cr_block`) is the same kind of tangent unknown. Inverse power matches this paper; Arnoldi is the Krylov form when several modes are wanted. Neither replaces the projection backsolve.
- **cites / related:** Mullen et al. 2008 (spectral coordinates); Crane et al. 2013 conformal Willmore (Dirac eigenproblem for tangents, a different operator).

## Knöppel, Crane, Pinkall, Schröder 2013 — Globally Optimal Direction Fields
- **file:** `knoeppel2013-globally-optimal-direction-fields.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/
- **doi:** https://doi.org/10.1145/2461912.2462005
- **tags:** `[direction-field, connection-laplacian, principal-curvature, eigen, holonomy, shell]`
- **summary:** Smooth n-direction fields stored as the n-th power of one complex number per vertex, so line fields (n = 2) and crosses (n = 4) need no sign or period jumps. A quadratic connection Dirichlet energy (with a holomorphic / anti-holomorphic blend `s`); its smallest generalized eigenvector `A u = lambda M u` is the globally smoothest field, and an alignment term with a guidance field turns it into one shifted linear solve `(A - lambda_t M) u = M q` (Eq. 16), `lambda_t` trading alignment (toward minus infinity) against smoothness (toward `lambda_1`). Principal-curvature guidance is the Hopf differential with magnitude `(k1 - k2)^2`, so umbilics defer to smoothness automatically. Transport angles from angle-rescaled vertex charts, holonomy per triangle.
- **maps_to:** `kusama/direction_field.hpp` (aligned solve on the vertex charts, `lambda_t = -1 / l^2`, plain cotan connection Laplacian through `kusama/twisted_laplacian.hpp`), demo `position_distortion::align_length` for the area_ratio shear direction. Stripe patterns (Knöppel 2015) share the twisted-Laplacian assembly.
- **cites / related:** Knöppel et al. 2015 (stripes, which take a 2013 field as input); Crane, Desbrun, Schröder 2010 (trivial connections); Stein et al. 2020 (vector Dirichlet on CR edges).

## Knöppel, Crane, Pinkall, Schröder 2013 — Globally Optimal Direction Fields
- **file:** `knoeppel2013-globally-optimal-direction-fields.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf
- **doi:** https://doi.org/10.1145/2461912.2462005
- **tags:** `[direction-field, parallel-transport, connection-laplacian, eigen, curvature, shell]`
- **summary:** Smooth n-direction fields as one complex value per vertex raised to the n-th power, so the energy is quadratic. The smoothest field is the smallest eigenvector of a connection Laplacian against the mass matrix (inverse power on one Cholesky factor); a family `s` of energies trades straighter lines (Dirichlet) for fewer singularities (holomorphic). Alignment with a guidance field (principal curvature: the Hopf differential, weight `(k1 - k2)^2`, so umbilics defer to smoothness) is one shifted linear solve `(A - lambda_t M) u = M q`, `lambda_t` below the smallest eigenvalue trading alignment for smoothness.
- **maps_to:** `kusama/direction_field.hpp` (aligned solve, `lambda_t = -1 / l^2`, on the vertex charts, plain cotan form via `kusama/twisted_laplacian.hpp`), demo `position_distortion::align_length` for the area_ratio shear direction. Same assembly as the stripe solve (Knöppel 2015), with the transport angle in place of the phase advance.
- **cites / related:** Knöppel et al. 2015 (stripes take this field as input); Crane, Desbrun, Schröder 2010 (trivial connections).

## Sharp, Soliman, Crane 2019 — The Vector Heat Method
- **file:** `sharp2019-vector-heat-method.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/paper.pdf
- **project:** https://www.cs.cmu.edu/~kmcrane/Projects/VectorHeatMethod/
- **doi:** https://doi.org/10.1145/3243651
- **tags:** `[parallel-transport, vector-dirichlet, direction-field, heat, laplacian, shell]`
- **summary:** Parallel transport of tangent vectors (and other vector-valued data) by short-time heat flow with the connection Laplacian. One implicit step \((M - tL^\nabla)X = X_0\) diffuses vectors along the surface; the magnitudes are restored from a scalar heat solve with the same \(t\), so the directions come out as parallel transport along shortest geodesics without tracing any. Gives the logarithmic map, Karcher means, and vector field extension from sparse samples, all by a few prefactored sparse solves.
- **maps_to:** `kusama/material_direction.hpp::smooth`: one implicit step of the Stein et al. vector Dirichlet energy (`kusama/vector_dirichlet.hpp`, same CR edge unknowns as `shell_cr_block`) guided by last step's carried field, the per-step smoother for the `area_ratio_target::material` directions. Their stable alternative to shooting geodesics.
- **cites / related:** Crane, Weischedel, Wardetzky 2013 (scalar heat method); Knöppel et al. 2013 (globally optimal direction fields); Stein et al. 2020 (vector Dirichlet).

## Crane, Desbrun, Schröder 2010 — Trivial Connections on Discrete Surfaces
- **file:** `crane2010-trivial-connections.pdf`
- **url:** https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/paper.pdf
- **doi:** https://doi.org/10.1111/j.1467-8659.2010.01761.x
- **tags:** `[direction-field, parallel-transport, holonomy, singularities, shell]`
- **summary:** Smoothest direction fields with prescribed singularities. Choose the singular vertices and their indices (summing to the Euler characteristic, Poincaré–Hopf); one sparse linear solve then gives the connection closest to Levi-Civita whose holonomy cancels the Gaussian curvature around every cycle except at the singularities, plus prescribed holonomy around the nontrivial generators. Transport by that connection is path independent, so a vector at one vertex extends to a global field.
- **maps_to:** The closure step after the material directions (TODO 17): a direction field that comes back to itself around every loop (closed winding families), on any genus, once the singularities are fixed. On a torus with no singularities it reduces to prescribing the two generator holonomies.
- **cites / related:** Knöppel et al. 2013 (globally optimal direction fields); Knöppel et al. 2015 (stripe patterns place singularities automatically).

## Pinkall 1985 — Hopf tori in S³
- **file:** `pinkall1985-hopf-tori.pdf`
- **url:** https://www3.math.tu-berlin.de/geometrie/Lehre/SS12/MathVisII/Pinkall%20-%201985%20-%20Hopf%20tori%20in%20S3.pdf
- **doi:** https://doi.org/10.1007/BF01389060
- **tags:** `[willmore, hopf, elastic, shell, hepworth]`
- **summary:** A Hopf torus is the preimage of a closed curve on \(S^2\) under the Hopf map \(S^3\to S^2\). If the curve has length \(L\) and encloses oriented area \(A\), the torus is isometric to \(\mathbb{R}^2/\Gamma\) with lattice generated by \((2\pi,0)\) and \((A/2, L/2)\). That area is the holonomy of the Hopf connection; it is how the fiber closes, not a local tilt of the curvature circle. The mean curvature of the torus equals the curvature \(\kappa\) of the spherical curve, and the Willmore energy reduces to \(\int(1+\kappa^2)\,ds\). Critical points are the Hopf tori over elastic curves on \(S^2\). A great circle gives the Clifford torus. Stereographic projection sends these to Willmore tori in \(\mathbb{R}^3\).
- **maps_to:** `hopf_sweep` in `include/gaudi/hepworth/block/hopf_sweep.hpp`. The surface in this paper is a circle bundle over a curve on \(S^2\). The local section modulus \(\kappa=\mathrm{cn}(s,m)\) is not derived here; Pinkall cites Langer and Singer for the elastic curves.
- **cites / related:** Langer and Singer, *The total squared curvature of closed curves*, J. Differential Geom. 20 (1984), DOI `10.4310/jdg/1214438990` — elastic curves on \(S^2\), not fetched (no open PDF). Ferus and Pedit, Math. Z. 204 (1990), DOI `10.1007/BF02570873` — \(S^1\)-equivariant Willmore tori, not fetched (no open PDF). The arXiv preprint `dg-ga/9507001` is their later paper *Curved flats in symmetric spaces*, not this article.

## Eigensatz 2010 — Constraint-Based Surface Processing for Geometric Modeling and Architecture
- **file:** `eigensatz2010-constraint-based-surface-processing.pdf`
- **url:** https://diglib.eg.org/bitstream/handle/10.2312/8240/eigensatz.pdf?sequence=1
- **doi:** https://doi.org/10.2312/diss.20108240
- **tags:** `[curvature, metric, fairing, least-squares, Levenberg-Marquardt, shell, hepworth, paneling]`
- **summary:** ETH dissertation (Diss. ETH No. 19140). Geometric constraints on position, metric, and curvature, solved by Levenberg–Marquardt. Chapter 2 is the surface toolkit: first fundamental form, curvature, fairness, Gauss–Newton, and LM. Chapter 3 is constraint-based deformation, including curvature-domain shape processing (the Eurographics 2008 paper) and separate positional, metric, and curvature controls (Eurographics 2009). Chapter 4 is architectural paneling. Appendix A differentiates area, length, angle, and curvature on a triangle mesh.
- **maps_to:** `kusama/eigensatz_curvature.hpp` (eq. 3.13 tensor on the barycentric region, A.4 eigenvalue derivative, A.6 umbilic fallback) and `kusama/eigensatz_solver.hpp` (eq. 3.12 with eq. 3.11 conformal, areal, closeness; §2.4 Levenberg–Marquardt). Wrapped as the goal of `hepworth/block/shape_match.hpp`. Section text is in `reference/.cache/eigensatz2010-constraint-based-surface-processing/INDEX.md`.
- **cites / related:** Eigensatz, Sumner, Pauly 2008, *Curvature-Domain Shape Processing*; Eigensatz and Pauly 2009, *Positional, Metric, and Curvature Control*; Eigensatz et al. 2010, *Paneling Architectural Freeform Surfaces*.

## Bouaziz et al. 2014 — Projective Dynamics: Fusing Constraint Projections for Fast Simulation
- **file:** not fetched
- **doi:** https://doi.org/10.1145/2601097.2601116
- **tags:** `[hepworth, projective-dynamics, local-global, shape-matching]`
- **summary:** Constraints as \(\tfrac{w}{2}\|A q - p(q)\|^2\) with \(p\) a projection onto the constraint set. A local step projects every constraint, a global step solves a constant prefactored system with inertia. Generalizes shape matching and position-based dynamics.
- **maps_to:** `hepworth/block/solver.hpp` (`projection_solver`); every `block_constraint`. `shape_match` is a constraint whose projection is a full solve.

## Müller et al. 2005 — Meshless Deformations Based on Shape Matching
- **file:** not fetched
- **doi:** https://doi.org/10.1145/1073204.1073216
- **tags:** `[shape-matching, polar-decomposition, hepworth]`
- **summary:** Goal positions are the rest shape moved by the best-fit rigid transform (polar factor of \(A_{pq}\)); particles move a fraction \(\alpha\) toward them. Extensions: best-fit linear map blended with \(R\) by \(\beta\), volume-normalized; quadratic deformations; overlapping clusters; plasticity.
- **maps_to:** `hepworth/block/shape_match.hpp` (`shape_align::rigid` / `linear`, `rest_shape_goal`).

## Overby, Brown, Li, Narain 2017 — ADMM ⊇ Projective Dynamics
- **file:** not fetched
- **doi:** https://doi.org/10.1109/TVCG.2017.2730875
- **tags:** `[projective-dynamics, ADMM, hepworth]`
- **summary:** Projective dynamics is ADMM without the dual update. Adding a scaled dual per constraint converges to the exact joint minimizer and admits arbitrary nonlinear local solves.
- **maps_to:** Follow-up for `shape_match`: a dual on the goal rows if the inertia / goal balance point is wrong.

## Wang, Liu, Tong 2012 — Linear Surface Reconstruction from Discrete Fundamental Forms on Triangle Meshes
- **file:** `wang2012-discrete-fundamental-forms.pdf`
- **url:** https://www.cse.msu.edu/~ytong/DiscreteFundamentalForms.pdf
- **doi:** https://doi.org/10.1111/j.1467-8659.2012.03153.x
- **tags:** `[fairing, conformal, ARAP, curvature, hepworth]`
- **summary:** Edge lengths as the discrete first fundamental form and dihedral angles as the second; reconstructs vertex positions (unique up to rigid motion) by one sparse linear solve, from a convex quadratic energy on per-face frames lifted to the 6D rigid-motion space. Gives a local integrability condition (discrete Gauss / Mainardi–Codazzi). Demonstrates rotation-invariant deformation through point and orientation constraints.
- **maps_to:** `scripts/sympy/cr_poisson/face_loops.py` (`run_rotations`: CR shape-operator targets to dihedrals, rotation sync, Poisson); a route to curvature-domain (conformal) processing in hepworth.

## Dynamic mesh: correctness and quality (not fetched)

Background for `docs/dynamic_mesh_correctness_plan.md` (TODO item 9). Titles and summaries only; fetch when implementing.

- **Lawson 1977 — Software for C1 surface interpolation.** Introduces the edge-flip algorithm; flips by the opposite-angle criterion reach the planar Delaunay triangulation. Tags `[delaunay, remesh, asawa]`. Maps to `dynamic::flip_edges`.
- **Rosenberg, Stenger 1975 — A lower bound on the angles of triangles constructed by bisecting the longest side.** Repeated longest-side bisection keeps every angle at least half the smallest initial angle. Tags `[remesh, delaunay]`. Maps to the split order in `dynamic::subdivide_edges`.
- **Rivara 1984 — Algorithms for refining triangular grids suitable for adaptive and multigrid techniques; Rivara 1997 — New longest-edge algorithms for the refinement and/or improvement of unstructured triangulations.** Longest-edge refinement; the longest-edge propagation path (Lepp) ends at a terminal edge (longest in both incident triangles), and Lepp-Delaunay refinement splits terminal edges then restores Delaunay. Tags `[remesh, delaunay, parallel]`. Maps to terminal-edge (local-maximum) split selection.
- **Dey, Edelsbrunner, Guha, Nekhayev 1999 — Topology preserving edge contraction.** The link condition `Lk(a) ∩ Lk(b) = Lk(ab)` for edge contraction to preserve the topology of a triangulated manifold. Tags `[remesh, asawa, verification]`. Maps to `collapse_edge` guards.
- **Bobenko, Springborn 2007 — A discrete Laplace–Beltrami operator for simplicial surfaces; Fisher, Springborn, Schröder, Bobenko 2007 — An algorithm for the construction of intrinsic Delaunay triangulations with applications to digital geometry processing.** Intrinsic Delaunay triangulations of a surface mesh, reached by intrinsic flips. Tags `[delaunay, intrinsic, remesh]`. Maps to the flip-convergence question (H3).
- **Sharp, Soliman, Crane 2019 — Navigating intrinsic triangulations.** Signpost data structure for intrinsic triangulations, intrinsic Delaunay flips and refinement on surfaces. Tags `[delaunay, intrinsic, remesh]`. Maps to an intrinsic-flip reference in the quality experiment.
- **Botsch, Kobbelt 2004 — A remeshing approach to multiresolution modeling.** Isotropic remeshing by split long, collapse short, flip, tangential smoothing. Tags `[remesh, asawa]`. Maps to `dynamic::step`.
- **Gonthier 2008 — Formal proof: the four-color theorem.** Coq proof built on hypermaps (the combinatorial-map structure of the shell). Tags `[verification, lean, asawa]`. Maps to the Lean spec (track d).

## Angelidis, Cani, Wyvill, King 2006 — Swirling-Sweepers: Constant-Volume Modeling
- **file:** `angelidis2006-swirling-sweepers.pdf`
- **url:** https://dgp.toronto.edu/~silex/Publications/ACWK04/ACWK04b.pdf (HAL: https://inria.hal.science/inria-00402557, behind a bot check)
- **doi:** https://doi.org/10.1016/j.gmod.2005.11.001 (Graphical Models 68(4), 2006; Pacific Graphics 2004)
- **tags:** `[volume-preserving, space-deformation, swirl, sweepers, sculpting, conservative-distortion]`
- **summary:** A swirl is a point tool `c` with a rotation `M` of angle `theta` about axis `v`; a point moves by `M^phi(p)` with `phi = mu(|p - c|)` a smooth decreasing influence of radius `lambda` and `M^alpha = exp(alpha log M)`. Twisting space about an axis by a falloff angle has Jacobian determinant 1 (proof in Appendix A), so a swirl preserves volume. A swirling-sweeper realizes a drag `t`: n swirls (8 in practice) on a circle of radius `r` perpendicular to `t`, tangent axes, radii `2r`, angles chosen so the circle center maps exactly to `h + t` (Appendix B), blended by `M + N = exp(log M + log N)`. Large drags are split into `s` steps (proportional to speed, inversely to tool size) so the map stays foldover free and volume preserving; volume errors of 0.01 to 3.6% come from the finite steps.
- **maps_to:** TODO item 14 (sculpting demo). For TODO 13(c): a swirl's Jacobian away from its center is a rotation times a tangential shear of strength `r theta'(r)` (area preserving in the plane), the kind of per-ring target a conservative one-ring constraint holds.
- **cites / related:** Angelidis, Wyvill, Cani 2004 (Sweepers); Alexa 2002 (linear combination of transformations); von Funck et al. 2006.

## von Funck, Theisel, Seidel 2006 — Vector Field Based Shape Deformations
- **file:** `vonfunck2006-vector-field-shape-deformations.pdf`
- **url:** https://graphics.stanford.edu/courses/cs468-07-winter/Papers/fts-vfbsd-06.pdf
- **doi:** https://doi.org/10.1145/1141911.1142002 (SIGGRAPH 2006)
- **tags:** `[divergence-free, volume-preserving, deformation, stream-function, conservative-distortion]`
- **summary:** Deformation as path-line integration of a C1 time-dependent divergence-free field, so it preserves volume and cannot self-intersect (path lines do not cross in space-time). 2D: `v = (-p_y, p_x)` from one scalar `p` (a stream function); 3D: `v = grad p x grad q`. Locality: inner region (`r(x) < r_i`) holds a chosen simple field `grad e x grad f` (constant: translation; linear: rotation; quadratic), the outer region is zero, and a Bernstein blend `b(r)` on the scalar fields keeps `v` divergence free and C1 in between. Modeling metaphors drag, twist and bend through these choices.
- **maps_to:** TODO 13(c): the stream-function family of area-preserving one-ring maps (`psi = uv` squeeze, `r^2 / 2` rotation, `f(r)` twirl, `f(u)` swept shear) is the 2D case of this construction, applied to rest spokes in the principal chart.
- **cites / related:** Angelidis et al. 2004/2006 (swirls); Theisel et al. 2005 (path lines in space-time); Davis 1967 (div-free fields from potentials).

## Combaz, Neyret 2002 — Painting Folds using Expansion Textures
- **file:** `combaz2002-expansion-textures.pdf` (to fetch; the host times out from the agent)
- **url:** http://evasion.imag.fr/Publications/2002/CN02/pg.pdf (mirror: http://evasion.inrialpes.fr/Publications/2002/CN02/pg.pdf)
- **doi:** https://doi.org/10.1109/pccga.2002.1167853 (Pacific Graphics 2002)
- **tags:** `[growth, expansion, rest-length, folds, buckling, remesh, morphogenesis]`
- **summary:** A painted expansion tensor field lengthens each edge's rest length (integrated along the edge); a surface-equilibrium solver (Green–Lagrange strain per triangle, plus curvature smoothing) then finds the folded shape. Fold wavelength and regularity are set by an anisotropic smoothing kernel. Small expansions are one step; large ones alternate partial expansion and equilibrium, and very large growth needs plasticity: relax the rest lengths toward the current ones (fade the memory of the initial state) so stress does not accumulate. The mesh is kept by Delaunay-style flips, splits on length or normal angle, and collapses at smaller limits for hysteresis.
- **maps_to:** TODO 13(a) (curvature-driven growth: rest growth, then ARAP to the grown rest, with rest relaxation as plasticity) and TODO 9(g) (split / collapse hysteresis).
- **cites / related:** Rodriguez et al. 1994 (growth as rest-state change); Baraff and Witkin 1998; Desbrun et al. 1999.

## Rodriguez, Hoger, McCulloch 1994 — Stress-dependent finite growth in soft elastic tissues
- **file:** `rodriguez1994-finite-growth.pdf`
- **url:** https://biomechanics.stanford.edu/me338/me338_project04.pdf
- **doi:** https://doi.org/10.1016/0021-9290(94)90021-3 (J. Biomechanics 27(4), 455–467)
- **tags:** `[growth, morphogenesis, rest-state, residual-stress, elasticity]`
- **summary:** Volumetric growth changes the local zero-stress reference state: a growth deformation gradient `F_g` (taken as a pure stretch `U_g`) maps the old stress-free state to a new, possibly incompatible one; an elastic `F_e` restores continuity, so the observed growth is `F = F_e F_g` and residual stress comes from `F_e` alone. The growth rate tensor `D_g` (symmetric part of the growth velocity gradient, `tr D_g` the dilation rate) can follow a constitutive law in the stress.
- **maps_to:** TODO 13(a): per-ring growth is `F_g` on the rest spokes; the spoke ARAP is `F_e`; incompatible growth (curvature-dependent `F_g`) is what buckles and branches the surface.
- **cites / related:** Skalak 1981 (incompatible growth); Fung 1990 (residual strain, opening angles).

## Merks 2003 — Branching Growth in Stony Corals: a Modelling Approach (thesis)
- **file:** `merks2003-branching-coral-thesis.pdf` (University of Amsterdam thesis; personal-use download)
- **url:** https://pure.uva.nl/ws/files/3391622/26478_Merks.pdf
- **related paper:** Merks, Hoekstra, Kaandorp, Sloot 2003, *Models of coral growth: spontaneous branching, compactification and the Laplacian growth assumption*, J. Theor. Biol. 224(2), 153–166, https://doi.org/10.1016/S0022-5193(03)00140-1
- **tags:** `[growth, branching, accretive-growth, curvature-rule, Laplacian-growth, morphogenesis]`
- **summary:** Accretive growth: layers deposited normal to the surface, thickness set by growth functions of local shape and environment. Earlier models (Kaandorp 1991–1995) needed a curvature rule (growth highest where curvature is high, suppressed on flat or hollow sections) and explicit branch initiation; the thesis derives branching instead from resource transport (diffusion-limited aggregation and Laplacian growth, both unstable by the Mullins–Sekerka mechanism, where protruding tips grow faster and split).
- **maps_to:** TODO 13(a): the curvature rule is the direct recipe (normal growth or rest growth scaled by curvature); the Laplacian-growth chapters explain when tips split rather than smooth out.
- **cites / related:** Kaandorp and Kubler 2001; Witten and Sander 1981 (DLA); Mullins and Sekerka 1964.

---

## Tag cheat sheet

| Tag | Use when discussing… |
|-----|----------------------|
| `hepworth` | Projection dynamics / shell-rod constraints |
| `asawa` | Mesh topology, winding, surface ops |
| `calder` | Integrators / MLS / field fits |
| `albers` | Normal-constrained LS shape fits (quadric, cyclide) |
| `duchamp` | Demos / fluid / experiment hosts |
| `ARAP` / `ADMM` / `L1-normals` / `spokes-and-rims` | Cubic stylization & Chao ARAP |
| `repulsive` / `TP` | Tangent-point / self-avoiding curves & surfaces |
| `framed-knot` / `writhe` / `Calugareanu` | Ribbon Lk=Tw+Wr, writhe framing, rod twist |
| `untangle` | Self-intersection / fold repair |
| `developable` / `paper-folding` | Flattenable / ruled shells, nets, seams |
| `willmore` / `conformal` / `fairing` | Conformal curvature flow, Willmore fairing |
| `anderson` / `acceleration` / `fixed-point` | Accelerating local-global / projective-dynamics solvers |
| `growth` / `branching` / `morphogenesis` / `rest-state` | Growth by changing the rest metric, curvature-driven branching |
| `volume-preserving` / `divergence-free` / `swirl` / `conservative-distortion` | Area- or volume-preserving maps, stream functions, swirling-sweepers |
| `hopf` / `elastic` | Hopf tori, elastic curves on \(S^2\) |
| `winding` | Inside-outside / solid angle |
| `fluids` / `spectral` / `subspace` | Fluid papers |
| `mass-spring` | Implicit cloth / optimization time integration |
| `implicit-surface` | Polygon-soup implicits, variational / RBF implicits
| `RBF` / `variational` / `scattered-data` | Interpolating implicit fields, RBF mode bases |
| `vector-dirichlet` | Smooth vector fields on meshes |
| `tracing` / `coordinates` | Harnack tracing / integer coords |
| `delaunay` / `intrinsic` / `remesh` | Flips, Delaunay quality, split / collapse order |
| `verification` / `lean` | Topology invariants, link condition, formal specs of mesh operations |
