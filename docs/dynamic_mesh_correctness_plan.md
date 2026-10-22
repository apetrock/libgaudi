# Dynamic mesh: correctness, quality, parallelism

Scope: the corner shell (`include/gaudi/asawa/shell/shell.hpp`), its topology
operations (`operations.hpp`), and the remeshing pass (`dynamic.hpp`,
`dynamic::step`). Four tracks, in order of cost: (a) guard tightness and
invariants by exhaustive small-mesh checking, (b) Delaunay quality versus
edge order, (c) parallel splits on the CPU and then the GPU, (d) a Lean spec
used as a proof object and a differential-testing oracle. TODO item 9.

---

## 0. The shell as a combinatorial map

Corners are darts. Slot `c` is live when `next(c) >= 0`.

| Structure | In the code | As a map |
|---|---|---|
| edge | corner pair `{2k, 2k + 1}`, `other(c) = 2 (c / 2) + (c + 1) % 2` | fixed-point-free involution `alpha` |
| face | cycle of `next` | cycles of `phi = next` |
| vertex | cycle of `vnext(c) = other(prev(c))` (`for_each_vertex`) | cycles of `alpha . phi^{-1}` |
| labels | `corners_vert`, `corners_face` | constant on the matching cycles |
| heads | `vert_begin`, `face_begin` | a dart in the matching cycle (cache) |

### Invariants (I1–I9)

- I1 `next` and `prev` are inverse on live corners.
- I2 Live is closed under `other`: both corners of a pair are live or both dead.
- I3 A dead corner has `next = prev = vert = face = -1`.
- I4 `face(c)` is constant on every `next` cycle; distinct cycles have distinct face ids.
- I5 `vert(c)` is constant on every `vnext` cycle; distinct cycles have distinct vertex ids (a vertex owns exactly one ring).
- I6 A live face `f` has `fbegin(f)` in its cycle; a dead face has `fbegin(f) = -1`. Same for vertices with `vbegin`.
- I7 Every face cycle has length at least 3 after `delete_degenerates` (and exactly 3 after `triangulate` / on the remesher's input).
- I8 `vert(c) != vert(other(c))` (no loop edges) after `delete_degenerates`.
- I9 Euler bookkeeping: `chi = V - E + F` and the component count change exactly as the operation table says.

`shell::audit_links` checks part of I1–I2 today; this work adds
`asawa::shell::check_invariants(const shell &)` returning the first violated
invariant and the offending ids.

### Operations and their expected effect

| Operation (`operations.hpp`) | Guards in the code | dV | dE | dF | d chi |
|---|---|---|---|---|---|
| `split_edge` | none | +1 | +1 | 0 | 0 |
| `split_face(c0, c1)` | caller | 0 | +1 | +1 | 0 |
| `merge_face` (general) | `assert(c1 == other(c0))` | 0 | -1 | -1 | 0 |
| `merge_face` -> `remove_dangling` | same face, `c0 == next(c1)` | -1 | -1 | 0 | 0 |
| `merge_face` -> `remove_cruft` | same face, digon | derive | -1 | -1 or -2 | derive |
| `merge_edge(c)` (contraction step) | none | -1 | -1 | 0 | 0 |
| `subdivide_edge` | both faces triangles | +1 | +3 | +2 | 0 |
| `collapse_edge` | `count_cycle <= 1`, `vsize >= 3` at ends, `>= 4` at opposite vertices | -1 | -3 | -2 | 0 iff the link condition holds |
| `flip_edge` | both faces triangles | 0 | 0 | 0 | 0 |
| `merge_edge(cA0, cB0)` (bridge) | valences `>= 4`, `adjacent` -> weld | derive per branch | | | handle or component join |
| `dissolve_vertex(v)` (`element_qc.hpp`) | fan of triangles, ring vertices distinct with valence `>= 4`, no diagonal already an edge, best triangulation above a quality floor | -1 | -3 | -2 | 0 (reuses the freed pairs and faces) |
| `retriangulate_face(f)` (`element_qc.hpp`) | consistent cycle, no diagonal already an edge | 0 | +(n-3) | +(n-3) | 0 |
| `pack` | none | relabeling | | | 0 (isomorphism; datums permuted consistently) |

"derive" entries are part of the work: compute them from the code, then check.

---

## (a) Guard tightness and invariants (C++ first, then Python model)

**Goal.** Every operation preserves I1–I9 on every small closed
triangulation, and `collapse_edge`'s guards are compared with the exact
topological condition.

**Link condition** (Dey, Edelsbrunner, Guha, Nekhayev 1999): contracting
edge `ab` of a closed triangulated 2-manifold preserves its topology iff
`Lk(a) ∩ Lk(b) = Lk(ab)`. For a triangle mesh: `a` and `b` have exactly two
common neighbors (the two opposite vertices), and the mesh is not a
tetrahedron.

**Meshes.** Exhaustive or near-exhaustive families, all closed:
- spheres: the tetrahedron, octahedron, icosahedron, and every
  triangulation reachable by random edge flips and splits up to 12 vertices
  (or `plantri` output if installed);
- tori: 3 x 3 to 5 x 5 grids, plus random flips;
- genus 2: two tori bridged by `merge_edge`;
- degenerate-adjacent: vertices of valence 3, 4; edges whose opposite
  vertices have valence 3 (the `vsize < 4` guard's case).

**Checks.**
1. `check_invariants` after every operation on every edge (or edge pair for
   the bridge) of every mesh, in a new `include/gaudi/test/shell_topology_tests.hpp`.
2. For `collapse_edge`, classify every edge: guard accepts / rejects versus
   link condition true / false. Report counts per family. Any
   "accepted, link condition false" case is a bug: save the mesh (OBJ) as a
   regression fixture. "Rejected, link condition true" is the conservatism
   that leaves short edges and slivers behind.
3. Euler bookkeeping per operation against the table.

**Exit.** A table of guard outcomes per family; zero invariant violations
(or each one filed with a fixture); a proposed exact guard if the current
one is too strict or too loose.

---

## (b) Delaunay quality versus edge order

**Question.** Does longest-first splitting keep an already well-triangulated
(Delaunay) mesh Delaunay, and does the global length sort matter?

**Known.**
- Planar Lawson flips reach the Delaunay triangulation from any
  triangulation (Lawson 1977). Inserting a point into a Delaunay
  triangulation and flipping the suspect edges around it restores Delaunay
  (incremental insertion). A split without flips does not guarantee it.
- On surfaces, intrinsic flips reach the intrinsic Delaunay triangulation
  (Bobenko & Springborn 2007; Fisher et al. 2007); Sharp, Soliman & Crane
  2019 do intrinsic Delaunay refinement.
- Longest-edge bisection bounds angles below (Rosenberg & Stenger 1975);
  longest-edge refinement and Lepp-Delaunay improvement (Rivara 1984, 1997)
  work on terminal edges, edges that are the longest in both incident
  triangles: local length maxima, not a global sort. That is the same rule
  as a parallel priority-based selection (track c).

**What the remesher does** (`dynamic::step`): split long edges (sorted
longest first, greedy with `one_ring_dedup`, threshold `_Cs`), break
cycles, collapse short edges (sorted shortest first), merge, then one
`flip_edges` sweep in random order. The flip test compares the sums of
opposite angles for the two diagonals (Lawson's criterion on the surface's
own angles) plus a dihedral penalty `10 (theta_flip - theta_same)^2`.
Because a flip changes the surface, the planar convergence argument does not
apply.

**Experiment** (Python, `scripts/mesh_quality/`):
- meshes: planar Delaunay of random points; sphere and torus at the
  remesher's target edge length; a remeshed demo frame dumped from duchamp;
- selection: (1) current sorted greedy + `one_ring_dedup`; (2) terminal-edge
  selection (local maxima, hash tie-break), one round and rounds to maximal;
- flips: none; one random sweep (current); sweeps to convergence (count
  sweeps, detect cycling); intrinsic flips (reference);
- metrics: fraction of locally Delaunay edges (opposite angles sum to at
  most pi), minimum-angle histogram, valence histogram, edges split per
  pass, passes to reach the target length, flips per sweep.

**Hypotheses.**
- H1: midpoint split of a terminal edge plus local flips around the new
  vertex keeps a Delaunay mesh Delaunay (planar: a theorem; surface: measure).
- H2: terminal-edge selection gives quality at least equal to the sorted
  greedy order, with more splits per pass.
- H3: the current flip criterion terminates when swept to convergence (or
  find the cycling configuration the dihedral penalty creates).

**Proof targets after the data.** A potential that every accepted flip
strictly decreases (termination of H3); a surface version of the
Rosenberg–Stenger / Rivara angle bound for well-sampled meshes.

**Exit.** A written result for H1–H3 with plots; a recommendation on the
split order and the flip schedule.

---

## (c) Parallel splits, then GPU

**Profile first.** Time `dynamic::step` (per phase) against the PD solve per
frame on the monge fairing and knotted surface demos. The GPU track is
worth it only if remeshing is a large share or the mesh must live on the GPU
for rendering anyway.

**What is already parallel-shaped** (`op_edges` in `dynamic.hpp`): all new
corners, vertices, faces and datum slots are allocated up front; each
subdivide gets deterministic slots (`vs + i`, `cs + 6 i`, `fs + 2 i`); the
OpenMP loop is commented out.

**What blocks it.**
1. Selection: `gather_edges` sorts every edge and marks neighborhoods
   greedily (serial).
2. Writes outside the op's two faces: `split_face` calls `vupdate` /
   `fupdate` (whole-ring label rewrites, idempotent values) and touches
   `vbegin`.
3. Datum `subdivide` per op: audit that each datum writes only its own
   preallocated slots.

**CPU design.**
1. Candidates: edges over threshold with two triangular faces.
2. Terminal-edge selection: an edge is selected if it has the largest
   (length, hash) among candidates sharing a face with it; repeat for a
   maximal set. Deterministic, no sort.
3. Split in parallel (OpenMP) into the preallocated slots with ring and head
   updates disabled inside the op.
4. Datums in parallel.
5. One fix-up pass: rebuild `corners_vert` / `corners_face` labels and
   `vbegin` / `fbegin` (parallel per corner; `update_head` is the serial
   version).
Conflict zone: the two incident faces, instead of both endpoints' one-rings.
Correctness argument: (d) proves splits on disjoint face pairs commute and
deferring the label and head updates is exact.

**GPU design** (WebGPU compute, as vermeer already uses WebGPU): the shell is
structure-of-arrays integers (`corners_next`, `corners_prev`,
`corners_vert`, `corners_face`, `vert_begin`, `face_begin`), one buffer each.
Kernels: mark, select (1–3 rounds), scan (prefix sum for slot offsets),
split, datum interpolation, fix-up, pack (stream compaction). Flips are the
same shape (two faces). Collapses need the link-condition check (2-ring
reads) and a larger conflict zone. Bridging merges (BVH proximity, non-local)
stay on the CPU.

**Exit.** CPU: split pass parallel with identical topology to the serial
pass on the (a) meshes (up to relabeling), speedup measured. GPU: split /
flip / pack kernels matching the CPU on the same meshes.

---

## (d) Lean spec

**Goal.** An executable, proved specification of the shell operations, used
to (1) prove invariant preservation and the Euler table, (2) state the
collapse and bridge theorems with the link condition, (3) export test
vectors that the C++ (and a future Rust) shell must reproduce. Precedent:
the four-colour theorem in Coq is built on hypermaps, the same structure
(Gonthier 2008).

**Model** (`lean/Asawa/`, Lean 4 + Mathlib):
- `Shell`: a finite type of dart slots, a `Finset` of live darts, `alpha`
  (involution), `phi : Equiv.Perm` (next), vertex and face labels as
  functions, heads as caches.
- `Invariants`: I1–I9 as `Prop`s.
- `Ops`: each operation as `Shell -> Option Shell` (none when the guard
  rejects), written at the permutation level: most are a transposition
  composed with `phi` (merging two cycles or splitting one), which is the
  common rule behind `split_face`, `merge_face`, its `remove_dangling` /
  `remove_cruft` branches, vertex split and merge, and the bridge's
  `swap_rows`.
- `Theorems`: preservation for split, flip, subdivide (first); the Euler
  table; commutativity of operations with disjoint face supports and
  exactness of deferred label / head updates (for track c); collapse under
  the link condition; bridge Euler change.
- `Enumerate`: bounded checking with `decide` over small cases before
  proofs.
- `Export`: JSON test vectors (mesh, operation, expected mesh up to
  relabeling) replayed by a C++ test against `asawa::shell`.

**Milestones.**
1. Model, invariants, `split_edge` / `split_face` / `flip_edge` with
   preservation proofs and the exporter; C++ replay test green.
2. `merge_face` branches and `subdivide_edge`; Euler table.
3. Commutativity and deferred updates (unblocks c).
4. `collapse_edge` with the link condition; the bridge.

---

## References

Lawson 1977 (flip algorithm); Rosenberg & Stenger 1975 (longest-side
bisection angle bound); Rivara 1984, 1997 (longest-edge refinement, Lepp,
terminal edges); Dey, Edelsbrunner, Guha, Nekhayev 1999 (link condition);
Bobenko & Springborn 2007, Fisher, Springborn, Schröder, Bobenko 2007
(intrinsic Delaunay); Sharp, Soliman, Crane 2019 (intrinsic triangulations);
Botsch & Kobbelt 2004 (isotropic remeshing); Gonthier 2008 (four-colour
theorem, hypermaps). Parallel mesh adaptation with independent sets and
deferred updates: PRAgMaTIc (Rokos, Gorman; to verify). Entries in
`reference/MANIFEST.md`.
