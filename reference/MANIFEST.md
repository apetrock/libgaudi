# Reference paper index

Local PDFs live in this directory. Prefer this manifest for tag/keyword search before opening a PDF.
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
- **tags:** `[repulsive, curves, tangent-point, rods, calder, TP, geometry-processing]`
- **summary:** Fractional tangent-point energies for curve fairing / self-avoidance without collisions as hard constraints.
- **maps_to:** Calder TP rod energies / gradients; braid / rod demos that need self-repulsion.

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
| `winding` | Inside-outside / solid angle |
| `fluids` / `spectral` / `subspace` | Fluid papers |
| `mass-spring` | Implicit cloth / optimization time integration |
| `implicit-surface` | Polygon-soup implicits, variational / RBF implicits
| `RBF` / `variational` / `scattered-data` | Interpolating implicit fields, RBF mode bases |
| `vector-dirichlet` | Smooth vector fields on meshes |
| `tracing` / `coordinates` | Harnack tracing / integer coords |
