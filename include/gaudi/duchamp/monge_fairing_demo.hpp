#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gaudi/asawa/primitive_objects.hpp"
#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/dynamic.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/duchamp/fields.hpp"
#include "gaudi/geometry_logger.hpp"
#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/face_element.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/hepworth/block/shape_det_conserve.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/kusama/material_direction.hpp"
#include "gaudi/kusama/stripes.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

namespace gaudi {
namespace duchamp {

struct fairing_step {
  real h = 0.1;
  real damping = 0.1;
  int iterations = 45;
  /// Gauss–Newton rows keep their factor this many iterates (fresh gradient
  /// every iterate); 1 refactors every iterate, 0 only on ascent.
  int gn_refresh_every = 15;
  /// Step halvings when an iterate raises the energy (0: full GN steps).
  /// Stiff curvature rows diverge without it.
  int gn_backtracks = 0;
  /// On an iterate that raises the energy, refactor at the previous iterate
  /// and keep going (true), or stop the step there with an imperfect solve
  /// and no refactor (false). Either way the step returns its best iterate.
  bool gn_refresh_on_ascent = true;
  /// Relative energy rise an iterate may take before it counts as an
  /// ascent. A frozen factor converges to its model's minimum in a few
  /// iterates and then wobbles a few percent; this keeps iterating through
  /// that (the step still returns its best iterate). 0 stops on any rise.
  real gn_ascent_tol = 0.05;
  /// Log which constraints minimize an energy at the start of each step
  /// (`hepworth::block::audit_variational`); slow, for diagnosis.
  bool audit_variational = false;
  /// The step: `ldlt` (the factored solve) or `modal` (the same with the
  /// low-mode Newton / buckling correction, configured by `modal`).
  hepworth::block::step_preconditioner preconditioner = hepworth::block::step_preconditioner::ldlt;
  /// The buckling preconditioner, fields in order: `modes` (smallest
  /// generalized eigenvalues corrected per step), `ncv`, `lambda_max` (Newton
  /// step below it), `lambda_neg` (negative curvature at or below it),
  /// `tau` (negative-mode step: lowers the model by `tau` times the energy),
  /// `fd_eps`, `max_iterations`.
  hepworth::block::modal_config modal{};
};

/// Dirichlet bake onto edge tangents (golden field). `lambda` is the guidance
/// weight. `reseed` starts square one each frame: replace carried tangents,
/// rewrite vertex frames and shape operators (including pin rest). Off fills
/// only unmarked zeros and leaves an existing rest. `frame_angle` turns every
/// stored frame by this many radians once on write. Positive turns the
/// tangent toward N × T. `w_pin` holds the operator at rest. 0 leaves the pin off.
struct fairing_dirichlet {
  real lambda = 50.0;
  bool reseed = true;
  real frame_angle = 0.0;
  real w_pin = 0.0;
};

/// Vertex frames from the golden edge-tangent weingarten (`vertex_frames`),
/// drawn after `relax` and before remesh. Red is k0 (larger eigenvalue),
/// green is k1 (smaller), blue the normal. `edge_tangents` draws the stored
/// tangents after the split, cyan, through the midpoint.
struct fairing_draw {
  bool frames = true;
  /// Scale red/green by raw `|κ|`: half-length is `nscale * |κ|` for that
  /// principal value (k0 and k1 each get their own length). Off keeps fixed
  /// `nscale` on both. Always reads the solved shape-block operator at vertices.
  bool curvature = false;
  bool edge_tangents = false;
  /// Axis half-length in mean edge lengths (`nscale = length * mean_edge`).
  /// With `curvature`, red uses `nscale * |k0|` and green `nscale * |k1|`.
  real length = 1.5;
  /// Orange: the carried `material` stretch axes, after
  /// `relax` and before remesh.
  bool material = false;
};

/// What couples the edge curvature store to positions.
///  off        — nothing.
///  weingarten — the store against the mean of each edge's two faces' live CR
///               read-outs (`weingarten_gn`), with the intrinsic vector
///               Dirichlet memory (`position_distortion.w_memory`).
enum class curvature_solver { off, weingarten };

/// `w_coupling` is the coupling row between position space and shape space
/// (`weingarten_gn`), like the rod rows that keep positions and frames in
/// sync: `shape_distortion` rows act on the edge store, `position_distortion`
/// rows on the positions, and the coupling keeps the two in sync. Its weight
/// is a row weight like every other (`A = w * factor`, the factor the edge's
/// dual area, so its residual is an integrated curvature, the quantity
/// `bend` weighs). A shape row at the same weight lands the store halfway;
/// keep it at or above the position rows that hold curvature (`w_bend`) for
/// shape edits to show in the mesh. Without it (`off`) the store reads its
/// frame-start value at weight 1.
struct curvature_config {
  curvature_solver kind = curvature_solver::weingarten;
  real w_coupling = 1.0;
  /// Damping on the curvature store: every iterate it is also pulled toward
  /// its previous iterate at this weight (the `w_coupling` scale). Like mass,
  /// it slows the store's steps, but it vanishes once the iterates settle, so
  /// the frame's result is unchanged. Helps the solve on large meshes. 0 is
  /// off.
  real w_store_damping = 0.0;
};

/// Shape → position projectors: Monge, Albers, spine, frame fair, and the
/// curvature solver. These land or couple stored curvature into immersion.
/// `one_ring_operators` builds spine, Monge, and stretch on the vertex
/// one-ring. Off keeps those on the face. Albers stays on the one-ring
/// either way.
struct shape_project {
  /// One shape per line. Comment a line to drop it. `quadric` matches the
  /// fitted jet, so it wins every spoke it is allowed into.
  /// `cylinder_minor` projects onto the tight principal circle.
  /// `cylinder_major` projects onto the loose one. The other axis is the ruling.
  std::vector<hepworth::block::monge_shape> shapes{
      //hepworth::block::monge_shape::cylinder_minor,
      //hepworth::block::monge_shape::cylinder_major,
      //hepworth::block::monge_shape::sphere,
      //hepworth::block::monge_shape::osculating_torus,
      //hepworth::block::monge_shape::hyperboloid,
      //hepworth::block::monge_shape::quadric,
  };
  /// Closest-point landings on the one ring. Comment a line to drop it.
  /// `w_fair_albers` scales the row. 0 leaves it off.
  std::vector<hepworth::block::albers_shape> albers{
      //hepworth::block::albers_shape::quadric,
      //hepworth::block::albers_shape::sphere,
      //hepworth::block::albers_shape::darboux_cyclide,
  };
  real w_fair_albers = 0.0e-1;
  /// Canal spines. Comment a line to drop it. The smaller ring residual wins.
  /// `w_frame_fair` scales the row.
  std::vector<hepworth::block::spine_kind> spine{
      hepworth::block::spine_kind::cylinder,
      //hepworth::block::spine_kind::sphere,
      hepworth::block::spine_kind::torus,
      //hepworth::block::spine_kind::ellipsoid,
  };
  bool one_ring_operators = true;
  real w_frame_fair = 0.0;
  /// Canal-radial residual bite in [0, 1]. `1` lands fully each project;
  /// `0.25` takes a quarter of `r_tgt − |p − C|` along the canal radial.
  real spine_step = 0.0;
  /// Edge curvature store → positions.
  curvature_config curvature{};
};

/// Curvature edits on the edge curvature store (one 2x2 curvature tensor per
/// edge). The mesh follows through `shape_project.curvature.w_coupling`: at
/// equal weights the store lands halfway between the mesh and the edit.
struct shape_distortion {
  /// One entry per edit, applied every frame against the start of the frame.
  /// Fields, in the order they must be written (leave out what you do not
  /// use); ratio means |larger curvature| / |smaller curvature|:
  ///   .ratio        ratio to reach; 1 round (saddles: minimal). NaN off.
  ///   .log          ratio is given as log(ratio); 0 round.
  ///   .relative     ratio multiplies the frame-start ratio (log: adds).
  ///   .mean         mean curvature times this, per frame. NaN off.
  /// Examples: {.ratio = 1.0} round; {.ratio = 0.9, .relative = true} 10%
  /// less anisotropic per frame. The lines' shear is its own row (`shear`).
  std::vector<hepworth::block::anisotropy_target> anisotropy{
      hepworth::block::anisotropy_target{
          //.ratio = 1.6180339887498949, // golden ratio
          //.ratio = 0.5, // golden ratio
          .ratio = 1.0,
          .log = false,
      },
  };
  /// Weight of the `anisotropy` rows. 0 is off.
  real w_anisotropy = 0.0;
  /// Shears the curvature pattern of each vertex ring per frame, the
  /// shape-space match of `position_distortion::shear`: the ring's
  /// frame-start tensors pushed through the area-preserving shear by
  /// `1 + shear` along `.angle` (radians from the tight tangent; NaN is the
  /// 45-degree bisector) and `1 - shear` across, then turned by the ring's
  /// best common turn (ARAP in shape space); mean curvature free.
  /// `.winding = {.ring, .tube, .handedness}` sets the angle at startup like
  /// stretch. Comment a line to drop it. `w_shear` scales every entry; 0
  /// leaves the row off.
  std::vector<hepworth::block::curvature_shear_target> shear{
      hepworth::block::curvature_shear_target{.shear = 0.5},
  };
  real w_shear = 0.0;
  /// Smoothing of the edits: 0 pulls each edge to its goal on its own; above
  /// 0 also makes an edit that changes from edge to edge cost more, so it
  /// lands smoothly (the old Crane spoke map).
  real smooth = 0.0;

  /// Willmore flow on the store with a feature-size equalizer. Fields, in
  /// the order they must be written:
  ///   .bands  {radius, gain} pairs in ascending radius (mean edge lengths).
  ///           A band covers sizes from its radius to the next band's; the
  ///           last runs to the whole surface. gain is the fraction of the
  ///           flow per frame in that band; 0 leaves the band free (the row
  ///           ignores it, so position rows can bend and buckle it).
  ///   .mean   0 rounds only (anisotropic part); 1 also the mean curvature.
  /// One band, or equal gains everywhere, is plain Willmore (no filtering).
  /// Several entries can be active at once; each is its own row.
  std::vector<hepworth::block::willmore_target> willmore{
      // Plain Willmore, 30% per frame, every feature size:
      //hepworth::block::willmore_target{.bands = {{0.0, 0.3}}},
      // Plain Willmore, 10% per frame, every feature size:
      //   hepworth::block::willmore_target{.bands = {{0.0, 0.1}}},
      // Gentle plain Willmore, 2% per frame:
      //   hepworth::block::willmore_target{.bands = {{0.0, 0.02}}},
      // Full Willmore, mean curvature too (rounds and shrinks):
      //   hepworth::block::willmore_target{.bands = {{0.0, 0.1}}, .mean = 1.0},
      // Mids only: 2 to 6 edge lengths flow, finer and coarser free:
      //   hepworth::block::willmore_target{
      //       .bands = {{0.0, 0.0}, {2.0, 0.1}, {6.0, 0.0}}},
      // Crane's high-pass: small features flow, the large shape is free:
      //   hepworth::block::willmore_target{.bands = {{0.0, 0.1}, {6.0, 0.0}}},
      // Fillet creases only, everything else free:
      //   hepworth::block::willmore_target{.bands = {{0.0, 0.3}, {1.5, 0.0}}},
      // Crease -> cylinder -> round, mids free to buckle:
      //   hepworth::block::willmore_target{
      //       .bands = {{0.0, 0.2}, {2.0, 0.0}, {8.0, 0.1}}, .mean = 0.5},
      // Balloon only: the large scale rounds and evens out, the rest free:
      //   hepworth::block::willmore_target{
      //       .bands = {{0.0, 0.0}, {8.0, 0.1}}, .mean = 1.0},
      // Strong crease plumping, slow balloon, wide free middle:
       hepworth::block::willmore_target{
      .bands = {{0.0, 1.0}, {3.0, 0.0}, {6.0, 1.0}}, .mean = 0.5},
  };
  /// Weight of the `willmore` rows. 0 is off.
  real w_willmore = 0.1;

  /// Shape ARAP: holds each edge's curvature at its frame-start value (how
  /// anisotropic it is, and the pattern of its lines) but lets each patch
  /// turn its curvature lines together. Edits pushed harder than this have
  /// to give somewhere. 0 is off.
  real w_shape_arap = 0.0;
  /// Shape ARAP also holds the mean curvature; off holds only the
  /// anisotropic part.
  bool shape_arap_mean = false;

  /// Spin-2 Dirichlet rows toward zero: fairs the store. Do not run at full
  /// weight with `w_curvature_ratio`.
  real w_grad = 0.0;

  /// Dirichlet rows on the trace `k11 + k22` toward zero (fairs the mean
  /// curvature). 0 leaves it off.
  real w_curvature_mean = 0.0;
  /// Evens out the anisotropy: pins each edge toward the mean anisotropy
  /// `q = (ratio - 1) / (ratio + 1)` of itself and its face neighbors,
  /// keeping its mean curvature (saddles: the anisotropic part) and axes.
  /// 0 leaves it off.
  real w_curvature_ratio = 0.0;
  /// The frame minimizer: pin toward each edge's principal curvatures set to
  /// the neighbors' normal curvature along the same directions. 0 leaves it
  /// off.
  real w_minimizer = 0.0;
  /// Dupin rows: each principal curvature held constant along its own
  /// direction, so curvature lines become circles and the surface is locally
  /// a Dupin cyclide (torus, cylinder, cone and their Moebius images). A
  /// fixed quadratic per frame (frames frozen at the frame start). 0 leaves
  /// it off.
  real w_dupin = 0.0;
  /// Skew rows: the Dupin rows holding each principal curvature's gradient
  /// along its own direction (how the curvature changes from one side of
  /// the ring to the other, per axis) at `skew_scale` times its frame-start
  /// value instead of zero. 1 maintains the skew against the other edits,
  /// above 1 grows it each frame, below 1 relaxes it toward Dupin. Shape
  /// space only; the mesh follows through `w_coupling`. 0 leaves it off.
  real w_skew = 0.0;
  real skew_scale = 0.0;
  /// Directional curvature growth, the shape-space match of position
  /// `stretch`: each frame the normal curvature along `.angle` (radians from
  /// the tight direction toward the other) times `1 + x`, across it times
  /// `1 + y`, against the frame-start store; negative bends less that way.
  /// `.power` scales each fraction by `(|k| / mean |k|)^power` of that axis's
  /// frame-start curvature, so positive grows the more curved edges faster.
  /// `.winding` sets the angle at startup like stretch. Comment a line to
  /// drop it. `w_growth` scales every entry; 0 leaves the row off.
  std::vector<hepworth::block::curvature_growth_target> growth{
      hepworth::block::curvature_growth_target{
        .angle = 1.0, .x = 1.0, .y = 0.0,
        .winding = {.ring = 1, .tube = 3, .handedness = 1},
        .power = 0.0
      },
  };
  real w_growth = 0.0;
  /// Holds each vertex ring's total Gauss-image area `sum K dA` of the store
  /// at its frame-start value and lets it redistribute between the ring's
  /// edges and between their mean curvature and anisotropy
  /// (`ring_det_conserve`), the shape-space match of `area`. `K` is
  /// intrinsic, so this leaves bending free and resists only what would
  /// change the metric. 0 leaves it off.
  real w_det = 0.0;
  /// Holds each vertex's Gaussian curvature `det h_v` of the store at a rest
  /// `K0` baked on its first frame and carried through remeshing
  /// (`vertex_gauss_conserve`). Mean curvature and anisotropy are free, so
  /// more curvature has to come out as anisotropy, not scale (pointwise `K`
  /// changes as `1 / s^2`). 0 leaves it off.
  real w_gauss = 0.5;
  /// Rest relaxation of `w_gauss`, per frame `K0 += gauss_relax (K - K0)`:
  /// 0 keeps the first frame's `K0`, 1 resets it every frame (like `w_det`),
  /// small values let a pattern that persists (a buckle) creep into the rest.
  real gauss_relax = 0.5;
};

/// Position-space distortion: stretch, area ratio, area, bend, smooth.
struct position_distortion {
  /// Position spokes scaled on the live minor and major tangents. `angle` is
  /// signed radians from the minor tangent toward the major. 0 leaves `x` on
  /// the minor and `y` on the major. Those are extra fractions: 0 keeps that
  /// axis, 0.1 makes it 1.1 times as long, and a negative value shrinks it.
  /// `power` multiplies each fraction by `|κ|^power` of that axis. 0 leaves
  /// the fractions. Negative pushes a flat direction more.
  /// (0, 0) adds no row. Comment a line to drop it. `w_stretch` scales every
  /// entry. 0 leaves the row off. `shape_project.one_ring_operators` builds
  /// this on the vertex ring.
  /// `.winding = {.ring, .tube, .handedness}` (last field) replaces `angle` at
  /// startup with the angle of a closed line that winds `ring` times around
  /// the ring and `tube` times around the tube of the fitted torus; {1, 1} is
  /// the Hopf fiber (Villarceau circle). Both 0 keeps `angle`.
  /// `.direction` (`hepworth::block::direction_source`) is where the axis
  /// comes from each step: `angle` off the live principal frame, `material`
  /// (set from the angle once, then carried with the surface;
  /// `material_smooth`, `material_follow`), or `stripes` (the angle's
  /// direction closed by a warm-started stripe pattern; `stripes`).
  std::vector<hepworth::block::stretch_target> stretch{
      hepworth::block::stretch_target{
          .angle = 1.1,
          .x = 0.2,
          .y = -0.04,
          .power = 0.5,
          .winding = {.ring = 1, .tube = 7, .handedness = 1},
      },
      // Along the Hopf fiber, computed from the torus at startup:
      //   hepworth::block::stretch_target{
      //       .angle = 0.0, .x = 0.2, .y = -0.04, .power = 0.45,
      //       .winding = {.ring = 1, .tube = 1, .handedness = 1}},
      // Two turns around the tube per trip around the ring:
      //   hepworth::block::stretch_target{
      //       .angle = 0.0, .x = 0.2, .y = -0.04, .power = 0.45,
      //       .winding = {.ring = 1, .tube = 2, .handedness = 1}},
  };
  real w_stretch = 0.0;
  /// Ring tilt: each vertex's one-ring, as it was at the start of the frame,
  /// is turned by `gain * metric * |κ|^power` radians per frame and held
  /// there, so it never settles (the row also holds the ring's pose against
  /// other drift). Fields (`hepworth::block::tilt_target`):
  ///  axis   — k0 (red) or k1 (green) principal tangent.
  ///  mode   — rigid: tilt the whole ring about the axis. spin: turn it
  ///           about the normal (counter-clockwise from outside). twist:
  ///           warp it into a saddle along the axis. bend: curl it across
  ///           the axis.
  ///  sign   — rigid only. gradient: tilt toward rising `metric`, scaled by
  ///           how fast it rises across the axis (no seams). propagate: one
  ///           axis direction for the whole mesh, seams only around round
  ///           points.
  ///  gain   — radians per frame. 0 leaves the entry off.
  ///  metric — what multiplies `gain`: aniso `(ratio - 1) / (ratio + 1)`,
  ///           log_ratio, or 1 / ratio.
  ///  power  — exponent on `|κ|` along the axis; 0 leaves it out.
  /// `w_tilt` scales every entry; 0 leaves the row off. One-ring operators
  /// only.
  std::vector<hepworth::block::tilt_target> tilt{
      hepworth::block::tilt_target{
          .axis = hepworth::block::principal_axis::k0,
          .mode = hepworth::block::tilt_mode::rigid,
          .gain = 0.5,
          .metric = hepworth::block::curvature_metric::aniso,
          .power = 1.0,
      },
  };
  real w_tilt = 0.0;
  /// Locks each one ring's covariant anisotropy: spoke rms along the tight
  /// principal tangent over that along the other held at `ratio` (1 is
  /// round). Area, skew and rotation are free, so it composes with `area`
  /// and `shear`. Comment a line to drop it. `w_axis_ratio` scales every
  /// entry; 0 leaves the row off.
  std::vector<hepworth::block::axis_ratio_target> axis_ratio{
      hepworth::block::axis_ratio_target{.ratio = 0.5},
  };
  real w_axis_ratio = 0.2;
  /// Shears each one ring every frame: stretch by `1 + shear` along `.angle`
  /// (radians from the tight tangent, as stretch; NaN is the 45-degree
  /// bisector) and `1 - shear` across, against the frame-start ring. Only
  /// that ratio is held; area, skew and rotation are free. At 45 degrees it
  /// is independent of `axis_ratio`; nearer the principal tangents the two
  /// overlap. `.winding = {.ring, .tube, .handedness}` sets the angle at
  /// startup like stretch (on the demo torus the bisector is about {1, 3}).
  /// Comment a line to drop it. `w_shear` scales every entry; 0 leaves the
  /// row off.
  std::vector<hepworth::block::shear_target> shear{
      hepworth::block::shear_target{
          .shear = 0.125,
          .winding = {.ring = 1, .tube = 3, .handedness = 1},
      },
      // Along the Hopf fiber, computed from the torus at startup:
      //   hepworth::block::shear_target{
      //       .shear = 0.125, .winding = {.ring = 1, .tube = 1, .handedness = 1}},
  };
  real w_shear = 1.0;
  /// Shape-to-position shear on each triangle: the face's curvature tensor
  /// (its edges' frame-start store tensors) is sheared as `shape.shear`
  /// shears it (by `1 + shear` along `.angle` from its tight axis and
  /// `1 - shear` across, `.winding` setting the angle at startup), and the
  /// triangle is held at the in-plane stretch that carries its CR tensor
  /// there with the edge normals held (`F = S*^-1 S`, area kept), any
  /// rotation. Faces near parabolic are skipped. Comment a line to drop it.
  /// `w_face_shear` scales every entry; 0 leaves the row off.
  std::vector<hepworth::block::curvature_shear_target> face_shear{
      hepworth::block::curvature_shear_target{.shear = 0.1},
  };
  real w_face_shear = 0.0;
  /// Smoothing length, in mean edge lengths, of the vector heat step applied
  /// to the `material` stretch axes at each step start. 0 only transports
  /// them. The steps compound (about `sqrt(steps)` times this length), so
  /// keep it small.
  real material_smooth = 0.0;
  /// Fraction each step by which the `material` stretch axes lean back
  /// toward the direction `angle` reads from the step-start charts, so the
  /// drive can follow the buckle. 0 is pure material; 1 follows the
  /// principal frame as `angle` does.
  real material_follow = 0.0;
  /// Stripe period, in mean edge lengths, for `stripes` stretch axes
  /// (`kusama::stripes`): the direction `angle` reads from the step-start
  /// charts becomes the level lines of the stripe pattern running across it,
  /// so the axis follows closed curves. The pattern is carried and
  /// warm-started, so it keeps its winding from step to step. Below about 8
  /// the mesh cannot resolve it.
  real stripes = 1.0;
  /// Cotangent one-ring. Holds the mean-curvature vector baked on first use,
  /// so it resists curvature edits. 0 leaves it off.
  real w_bend = 1.0e-2;
  /// Per-triangle rest area: weight `.w` (0 leaves it off) and `.band`, the
  /// range of its step-start area each face may take freely (`{1, 1}`
  /// conserves it; `{0.5, 1.5}` lets it drift 50% per step before pushing
  /// back).
  hepworth::block::area_config area{.w = 0.5, .band = {0.9, 1.1}};
  /// Cotangent Laplacian on positions, target Δx = 0. 0 leaves it off.
  real w_smooth = 0.0e-1;
  /// Per-triangle conformal (closest similarity) against the start-of-step
  /// shape. Scale is free, shear pays. 0 leaves it off.
  real w_conformal = 0.0;
  /// Intrinsic vector Dirichlet memory of the metric (`dirichlet_memory_gn`,
  /// rest at the frame start; factor the mean edge length). Weingarten solver
  /// only.
  real w_memory = 0.5;
  /// Holds the area-weighted centroid where the scene started (translation
  /// only; shape and rotation are free). 0 is off.
  real w_center = 0.5;
  /// After each step, move the mesh back so its centroid is exactly where it
  /// started, and remove the mean velocity. Exact and free, outside the solve.
  bool recenter = false;
};

/// Demo-only check of the coupled loop: pins the edge store to a known
/// shape. `scale` asks for `scale` times the frame-start curvature; `torus`
/// for the analytic torus (`torus_R`, `torus_tube`) about +Z. `w` is
/// a store pin like the `shape_distortion` rows.
/// Off adds no pin.
struct curvature_test_target {
  enum class mode { off, scale, torus };
  mode kind = mode::off;
  real scale = 1.0;
  real torus_R = 1.0;
  real torus_tube = 0.4;
  real w = 1.0;
};

struct monge_fairing_demo_config {
  /// Empty builds a torus. `.obj` and ASCII `.ply` load from this path.
  std::string mesh_path;
  /// The built torus's vertex counts around the ring and around the tube.
  int torus_ring = 48;
  int torus_tube = 24;
  /// Normal offset applied once at startup, in mean edge lengths. 0 is off.
  /// Each live vertex draws a uniform sample in [-1, 1] along its normal.
  real perturb = 0.0;
  fairing_step step{};
  fairing_dirichlet dirichlet{};
  fairing_draw draw{};
  shape_project project{};
  shape_distortion shape{};
  position_distortion position{};
  curvature_test_target curvature_test{};
  /// The remesher's element QC (`asawa::shell::qc_config`, fields in order):
  /// `datum`, `q_warn` (repair triangles below this quality), `q_min` (hard
  /// failure), `val_lo` / `val_hi` (interior valence band), `rounds`,
  /// `enabled`. Example: `.qc = {.q_warn = 0.15, .val_hi = 7}`.
  asawa::shell::qc_config qc{};
  /// Write bad triangles to `dump/bad_triangles.log` and an OBJ of their
  /// neighborhood: after each solve those below the solver's degenerate
  /// quality (the worst few when the solve blew up or did not converge), and
  /// after each remesh those the QC left below `qc.q_warn`.
  bool dump_bad_triangles = true;
};

inline const monge_fairing_demo_config k_monge_fairing_demo_default{};

/// Shell positions and edge tangents. The frame is derived from that field.
/// No rod, so nothing welds the surface.
class monge_fairing_demo {
public:
  using ptr = std::shared_ptr<monge_fairing_demo>;

  static ptr create(const monge_fairing_demo_config &cfg =
                        k_monge_fairing_demo_default) {
    return std::make_shared<monge_fairing_demo>(cfg);
  }

  explicit monge_fairing_demo(const monge_fairing_demo_config &cfg = {})
      : _cfg(cfg) {
    reset_scene();
  }

  void reset_scene() {
    __M = make_mesh(_cfg.mesh_path, _cfg.torus_ring, _cfg.torus_tube);
    asawa::shell::triangulate(*__M);
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);
    perturb_vertices_noise(*__M, x, _cfg.perturb);

    _center0 = area_centroid(x);
    set_stretch_angles(x);

    asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    const size_t n_stretch = _cfg.position.stretch.size();
    _material_slot.assign(n_stretch, -1);
    _gauss_slot = -1;
    _stripe_slot.assign(n_stretch, -1);
    _stripe_dirs.assign(n_stretch, {});
    for (size_t t = 0; t < n_stretch; ++t) {
      const hepworth::block::direction_source d = _cfg.position.stretch[t].direction;
      if (d == hepworth::block::direction_source::material)
        _material_slot[t] = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
      if (d == hepworth::block::direction_source::stripes)
        _stripe_slot[t] = asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    }
    _material_seeded.assign(_material_slot.size(), 0);
    _material_x0.clear();
    auto xs = std::make_shared<shell_vert_positions>(__M, 0);
    auto vs = std::make_shared<shell_vert_velocities>(__M, 1);
    _shell =
        std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);
    _cr = std::make_shared<hepworth::block::shell_cr_block>(__M);
    _cr->face_site = !_cfg.project.one_ring_operators || _cfg.dirichlet.w_pin > 0.0;
    real C = 1.5;
    const real l0 = asawa::shell::avg_length(*__M, xs->get());
    __surf = asawa::shell::dynamic::create(__M, C * l0, real(1.25) * C * l0,
                                           C * real(0.5) * l0);
    __surf->set_qc(_cfg.qc);
    __surf->set_qc_handler([this](asawa::shell::shell &, const asawa::shell::qc_report &rep) {
      if (!rep.bad_faces.empty())
        dump_bad_triangles("qc_left", false, 8, _cfg.qc.q_warn);
    });

    _solver_config =
        hepworth::block::block_solver_builder<
            hepworth::block::shell_position_block,
            hepworth::block::shell_cr_block>::create()
            .with_blocks(_shell, _cr)
            .with_presolve([this](hepworth::block::solver_context &) {
              if (_cr)
                hepworth::block::seed_cr_from_dirichlet(
                    *__M, *_cr, _cfg.dirichlet.reseed, _cfg.dirichlet.lambda,
                    _cfg.dirichlet.frame_angle);
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.area.w > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_area(*__M, ctx.constraints, _shell->xs->get(),
                                         _cfg.position.area,
                                         hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_center > 0.0) || !__M->verts_are_dense_packed())
                return;
              const std::vector<vec3> &x = _shell->xs->get();
              std::vector<index_t> ids(x.size());
              for (size_t i = 0; i < x.size(); ++i)
                ids[i] = static_cast<index_t>(i);
              ctx.constraints.push_back(hepworth::block::centroid_pin::create(
                  ids, vertex_areas(x), _center0, _cfg.position.w_center,
                  hepworth::block::select_blocks<0>(ctx)));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_conformal > 0.0) ||
                  !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_triangle_conformal(
                  *__M, ctx.constraints, _shell->xs->get(),
                  _cfg.position.w_conformal,
                  hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (_cfg.position.w_smooth <= 0.0)
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_laplacian(
                  *__M, ctx.constraints, _shell->xs->get(),
                  hepworth::block::laplacian_mode::null,
                  hepworth::block::laplacian_stencil::cotan, _cfg.position.w_smooth,
                  hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              ctx.shell_mesh = __M.get();
              if (!__M->verts_are_dense_packed())
                return;
              if (_cr && !ctx.blocks.empty()) {
                hepworth::block::init_shell_caches(
                    *__M, *_cr, ctx.blocks[0], ctx.computes);
                hepworth::block::init_face_shape_pin(
                    *__M, *_cr, ctx.blocks[0], ctx.constraints,
                    _cfg.dirichlet.w_pin);
                if (_cfg.project.one_ring_operators) {
                  std::vector<char> live(static_cast<size_t>(__M->vert_count()), 1);
                  hepworth::block::init_directed_stretch(
                      *__M, ctx.constraints, _shell->xs->get(), live,
                      {ctx.blocks[0]}, _cfg.position.w_stretch, stretch_now(),
                      &_cr->vertex_frames,
                      _cfg.position.w_stretch > 0.0
                          ? stretch_dirs_now()
                          : std::vector<const std::vector<vec3> *>{});
                  hepworth::block::init_directed_tilt(
                      *__M, ctx.constraints, _shell->xs->get(), live,
                      {ctx.blocks[0]}, _cfg.position.w_tilt, _cfg.position.tilt,
                      &_cr->vertex_frames);
                  hepworth::block::init_axis_ratio(
                      *__M, ctx.constraints, _shell->xs->get(), live, {ctx.blocks[0]},
                      _cfg.position.w_axis_ratio, _cfg.position.axis_ratio,
                      &_cr->vertex_frames);
                  hepworth::block::init_shear(
                      *__M, ctx.constraints, _shell->xs->get(), live, {ctx.blocks[0]},
                      _cfg.position.w_shear, shear_now(), &_cr->vertex_frames);
                } else {
                  hepworth::block::init_face_stretch(
                      *__M, *_cr, ctx.blocks[0], ctx.constraints, _cfg.position.w_stretch,
                      stretch_now());
                }
              }
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_bend > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_bending(
                  *__M, ctx.constraints, _shell->xs->get(), _cfg.position.w_bend,
                  hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              const bool face = _cfg.project.w_frame_fair > 0.0;
              const bool albers = _cfg.project.w_fair_albers > 0.0 && !_cfg.project.albers.empty();
              if (!face && !albers)
                return;
              ctx.shell_mesh = __M.get();
              if (face && _cr && !ctx.blocks.empty()) {
                std::vector<hepworth::block::objective> diamond;
                diamond.reserve(_cfg.project.shapes.size() + _cfg.project.spine.size());
                for (const hepworth::block::monge_shape shape : _cfg.project.shapes)
                  diamond.push_back(hepworth::block::monge_fairing{shape});
                for (const hepworth::block::spine_kind shape : _cfg.project.spine) {
                  if (shape == hepworth::block::spine_kind::none)
                    continue;
                  diamond.push_back(hepworth::block::spine_fairing{shape});
                }
                if (_cfg.project.one_ring_operators && !diamond.empty()) {
                  std::vector<char> live(static_cast<size_t>(__M->vert_count()), 1);
                  hepworth::block::init_frame_fair(
                      *__M, ctx.constraints, _shell->xs->get(), live,
                      _cfg.project.w_frame_fair, {ctx.blocks[0]}, diamond, -1,
                      &_cr->vertex_frames, _cr->shape_dof(), &_cr->stars,
                      _cfg.project.spine_step);
                } else {
                  hepworth::block::init_face_fair(
                      *__M, *_cr, ctx.blocks[0], ctx.constraints,
                      _cfg.project.w_frame_fair, diamond);
                }
              }
              if (albers && !ctx.blocks.empty()) {
                std::vector<hepworth::block::objective> spokes;
                spokes.reserve(_cfg.project.albers.size());
                for (const hepworth::block::albers_shape shape : _cfg.project.albers)
                  spokes.push_back(hepworth::block::albers_fairing{shape});
                std::vector<char> live(static_cast<size_t>(__M->vert_count()), 1);
                hepworth::block::init_albers_fair(
                    *__M, ctx.constraints, _shell->xs->get(), live,
                    _cfg.project.w_fair_albers, {ctx.blocks[0]}, spokes);
              }
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              push_curvature(ctx);
            })
            .dt(_cfg.step.h)
            .damping(_cfg.step.damping)
            .iterations(_cfg.step.iterations)
            .build();
  }

  /// Couples the edge store to positions with the configured solver, then
  /// adds the target rows on the store at the read-out row's unit.
  void push_curvature(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    const curvature_config &c = _cfg.project.curvature;
    if (!_cr || !_shell || ctx.blocks.empty() || !_cr->shape_dof() ||
        !__M->verts_are_dense_packed() || _cr->table.n_edges <= 0)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    if (static_cast<int>(x.size()) != _cr->table.n_verts)
      return;
    ctx.shell_mesh = __M.get();
    hepworth::sim_block::ptr pos = ctx.blocks[0];
    if (c.kind == curvature_solver::weingarten && c.w_coupling > 0.0) {
      ctx.constraints.push_back(hb::weingarten_gn::create(
          &_cr->table, x, c.w_coupling,
          std::vector<hepworth::sim_block::ptr>{pos, _cr->shape_dof()}));
      if (_cfg.position.w_memory > 0.0)
        ctx.constraints.push_back(hb::dirichlet_memory_gn::create(
            &_cr->table, x, _cfg.position.w_memory, std::vector<hepworth::sim_block::ptr>{pos}));
      push_edge_fairing(ctx);
      return;
    }
    // No live read-out row: the store reads its frame-start value at weight 1.
    ctx.constraints.push_back(hb::edge_shape_fair::create(
        &_cr->table, x, hb::edge_fair_rows::pin, real(1),
        std::vector<hepworth::sim_block::ptr>{_cr->shape_dof()}, hb::goal_fixed(_cr->shape)));
    push_edge_fairing(ctx);
  }

  /// `shape_distortion`, the frame minimizer and the test target as
  /// `edge_shape_fair` rows on the edge store.
  void push_edge_fairing(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    const std::vector<vec3> &x = _shell->xs->get();
    hepworth::sim_block::ptr store = _cr->shape_dof();
    const std::vector<hepworth::sim_block::ptr> held{store};
    auto fair = [&](hb::edge_fair_rows rows, real w, hb::edge_goal_fn goal,
                    real lambda = real(1), const std::vector<vec3> *frames = nullptr,
                    bool tensor_metric = false) {
      if (w > 0.0 && std::isfinite(w))
        ctx.constraints.push_back(hb::edge_shape_fair::create(&_cr->table, x, rows, w, held,
                                                              std::move(goal), lambda, frames,
                                                              tensor_metric));
    };
    const shape_distortion &s = _cfg.shape;
    // The store is the live reading at prepare: the rest the knobs measure against.
    const std::vector<vec3> &h0 = _cr->shape;
    const bool smooth = s.smooth > 0.0;
    for (const hb::anisotropy_target &spec : s.anisotropy)
      fair(smooth ? hb::edge_fair_rows::h1 : hb::edge_fair_rows::pin, s.w_anisotropy,
           hb::goal_anisotropy(store, spec, h0), smooth ? s.smooth : real(1), nullptr, true);
    fair(hb::edge_fair_rows::pin, s.w_shape_arap,
         hb::goal_shape_arap(store, &_cr->table, h0, s.shape_arap_mean));
    fair(hb::edge_fair_rows::pin, _cfg.project.curvature.w_store_damping, hb::goal_live(store));
    if (s.w_willmore > 0.0)
      for (const hb::willmore_target &spec : s.willmore)
        if (hb::edge_goal_fn goal = hb::goal_willmore_eq(store, _cr->table, x, h0, spec))
          fair(hb::edge_fair_rows::pin, s.w_willmore, std::move(goal));
    fair(hb::edge_fair_rows::dirichlet, s.w_grad, hb::goal_zero());
    fair(hb::edge_fair_rows::dupin, s.w_dupin, hb::goal_zero(), real(1), &h0);
    if (s.w_skew > 0.0 && std::isfinite(s.w_skew))
      ctx.constraints.push_back(hb::edge_shape_fair::create(
          &_cr->table, x, hb::edge_fair_rows::dupin, s.w_skew, held, hb::goal_zero(), real(1), &h0,
          false, s.skew_scale));
    fair(hb::edge_fair_rows::dirichlet_trace, s.w_curvature_mean, hb::goal_zero());
    std::shared_ptr<const hb::edge_frames> frames;
    const bool shearing = s.w_shear > 0.0 && !s.shear.empty();
    if (s.w_curvature_ratio > 0.0 || s.w_minimizer > 0.0 || shearing ||
        _cfg.curvature_test.kind == curvature_test_target::mode::torus)
      frames = hb::build_edge_frames(_cr->table, x);
    const position_distortion &pd = _cfg.position;
    if (pd.w_face_shear > 0.0 && !pd.face_shear.empty()) {
      if (!frames)
        frames = hb::build_edge_frames(_cr->table, x);
      std::vector<hb::curvature_shear_target> fs = pd.face_shear;
      for (hb::curvature_shear_target &spec : fs)
        spec.angle = winding_angle(spec.winding, spec.angle);
      hb::init_face_shear(_cr->table, *frames, h0, ctx.constraints, x, {ctx.blocks[0]},
                          pd.w_face_shear, fs);
    }
    for (hb::curvature_growth_target spec : s.growth) {
      spec.angle = winding_angle(spec.winding, spec.angle);
      fair(hb::edge_fair_rows::pin, s.w_growth, hb::goal_curvature_growth(spec, h0), real(1),
           nullptr, true);
    }
    hb::init_ring_det(_cr->table, x, h0, s.w_det, store, ctx.constraints);
    if (s.w_gauss > 0.0 && std::isfinite(s.w_gauss) &&
        x.size() == static_cast<size_t>(__M->vert_count())) {
      const hb::vertex_star_charts charts = hb::build_vertex_star_charts(_cr->table, x);
      std::vector<real> &K0 = gauss_rest_mut();
      const std::vector<real> K = hb::vertex_gauss(charts, h0);
      const real a = std::clamp(s.gauss_relax, real(0), real(1));
      for (size_t v = 0; v < K0.size() && v < K.size(); ++v) {
        if (!std::isfinite(K[v]))
          continue;
        K0[v] = std::isfinite(K0[v]) ? K0[v] + a * (K[v] - K0[v]) : K[v];
      }
      hb::init_vertex_gauss(_cr->table, x, charts, K0, s.w_gauss, store, ctx.constraints);
    }
    if (shearing)
      for (hb::curvature_shear_target spec : s.shear) {
        spec.angle = winding_angle(spec.winding, spec.angle);
        fair(hb::edge_fair_rows::pin, s.w_shear,
             hb::goal_curvature_shear(store, &_cr->table, *frames, _cr->vertex_frames, spec, h0));
      }
    if (s.w_curvature_ratio > 0.0)
      fair(hb::edge_fair_rows::pin, s.w_curvature_ratio, hb::goal_ratio_smooth(store, frames));
    fair(hb::edge_fair_rows::pin, s.w_minimizer, hb::goal_principal_mean(store, frames));
    const curvature_test_target &test = _cfg.curvature_test;
    if (test.kind != curvature_test_target::mode::off)
      fair(hb::edge_fair_rows::pin, test.w, hb::goal_fixed(test_targets(x, frames.get())));
  }

  /// `curvature_test` targets per edge, in the store's frames and sign.
  std::vector<vec3> test_targets(const std::vector<vec3> &x,
                                 const hepworth::block::edge_frames *frames) const {
    const curvature_test_target &test = _cfg.curvature_test;
    std::vector<vec3> h = _cr->shape;
    if (test.kind == curvature_test_target::mode::scale) {
      for (vec3 &v : h)
        v *= test.scale;
      return h;
    }
    if (!frames)
      return h;
    const kusama::stencil_table &t = _cr->table;
    std::vector<vec3> mid(static_cast<size_t>(t.n_edges), vec3::Zero());
    for (int f = 0; f < t.n_faces; ++f)
      for (int i = 0; i < 3; ++i)
        mid[static_cast<size_t>(t.face_edges[f][i])] =
            real(0.5) * (x[static_cast<size_t>(t.ids[f][i])] +
                         x[static_cast<size_t>(t.ids[f][(i + 1) % 3])]);
    real agree = 0;
    for (int e = 0; e < t.n_edges; ++e) {
      const vec3 &p = mid[static_cast<size_t>(e)];
      const real rho = std::hypot(p[0], p[1]);
      const vec3 hoop = rho > real(0) ? vec3(-p[1] / rho, p[0] / rho, 0) : vec3::UnitX();
      const vec3 c(p[0] / std::max(rho, real(1e-12)) * test.torus_R,
                   p[1] / std::max(rho, real(1e-12)) * test.torus_R, 0);
      const vec3 n = (p - c).normalized();
      const vec3 tube = n.cross(hoop);
      const real cphi = n.dot(vec3(c.normalized()));
      const real k_tube = real(1) / test.torus_tube;
      const real k_hoop = cphi / (test.torus_R + test.torus_tube * cphi);
      const mat3 S = k_tube * tube * tube.transpose() + k_hoop * hoop * hoop.transpose();
      const vec3 &te = frames->t[static_cast<size_t>(e)], &be = frames->b[static_cast<size_t>(e)];
      const vec3 he(te.dot(S * te), te.dot(S * be), be.dot(S * be));
      agree += (he[0] + he[2]) * (h[static_cast<size_t>(e)][0] + h[static_cast<size_t>(e)][2]);
      h[static_cast<size_t>(e)] = he;
    }
    if (agree < real(0))
      for (vec3 &v : h)
        v = -v;
    return h;
  }

  /// One solver pass. Does not remesh, so a meridional seam stays a seam.
  void relax() {
    _solver_config.dt = _cfg.step.h;
    _solver_config.damping = _cfg.step.damping;
    _solver_config.iterations = _cfg.step.iterations;
    _solver.gn_refresh_every = _cfg.step.gn_refresh_every;
    _solver.gn_backtracks = _cfg.step.gn_backtracks;
    _solver.gn_ascent_tol = _cfg.step.gn_ascent_tol;
    _solver.run_variational_audit = _cfg.step.audit_variational;
    _solver.preconditioner = _cfg.step.preconditioner;
    _solver.modal = _cfg.step.modal;
    _solver.gn_on_ascent = _cfg.step.gn_refresh_on_ascent
                               ? hepworth::block::projection_solver::ascent_policy::refresh
                               : hepworth::block::projection_solver::ascent_policy::stop;
    if (_cr)
      _cr->frame_angle = _cfg.dirichlet.frame_angle;
    hepworth::block::run_solver_step(_solver_config, _solver);
    transport_material_dirs();
  }

  void step(int frame, bool remesh = true) {
    _frame = frame;
    relax();
    recenter_mesh();
    transport_tangents();
    if (_cr)
      _cr->commit_edge_tangents();
    // Frames before remesh: projection field from vertex_weingarten.
    draw_solved_frames();
    measure_material_deflection();
    draw_material_dirs();
    // A solve that blew up, went non-finite or stopped short of converging
    // dumps its worst triangles even when none is below the threshold.
    const auto &r = _solver.last_step;
    dump_bad_triangles("after_solve", !std::isfinite(r.E_last) ||
                                          r.E_last > real(1e3) * std::max(r.E_first, real(1e-30)) ||
                                          r.dq_last > real(1e-3) * r.q_norm);
    if (!remesh)
      return;
    __surf->step(true);
    if (!_cfg.qc.enabled)
      dump_bad_triangles("after_remesh");
    draw_remeshed_edge_tangents();
  }

  /// Faces below `q_below` (quality `4 sqrt3 A / sum l^2`; default the
  /// solver's degenerate threshold), or the `n_worst` worst when `force`: one
  /// block per call
  /// appended to `dump/bad_triangles.log` (frame, stage, face, vertex ids,
  /// positions, edge lengths, area, quality), and
  /// `dump/bad_triangles_f<frame>_<stage>.obj` with those faces (group
  /// `bad`) and every face sharing a vertex with one (group `ring`). Off with
  /// `monge_fairing_demo_config::dump_bad_triangles`. Returns the number
  /// below the threshold.
  int dump_bad_triangles(const char *stage, bool force = false, int n_worst = 8,
                         real q_below = kusama::weingarten::k_degenerate_quality) {
    if (!_cfg.dump_bad_triangles || !__M)
      return 0;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    struct tri {
      int f;
      std::array<int, 3> v;
      real q;
    };
    std::vector<tri> all;
    for (auto fi : __M->get_face_range()) {
      const asawa::shell::FaceId fid = asawa::shell::face_id(fi);
      if (__M->fsize(fid) != 3)
        continue;
      const auto t = __M->get_tri(fid);
      const std::array<int, 3> v{static_cast<int>(t[0]), static_cast<int>(t[1]),
                                 static_cast<int>(t[2])};
      if (v[0] < 0 || v[1] < 0 || v[2] < 0 || static_cast<size_t>(v[0]) >= x.size() ||
          static_cast<size_t>(v[1]) >= x.size() || static_cast<size_t>(v[2]) >= x.size())
        continue;
      all.push_back({static_cast<int>(fi), v,
                     kusama::weingarten::face_quality(x[static_cast<size_t>(v[0])],
                                                      x[static_cast<size_t>(v[1])],
                                                      x[static_cast<size_t>(v[2])])});
    }
    std::sort(all.begin(), all.end(), [](const tri &a, const tri &b) { return a.q < b.q; });
    int n_bad = 0;
    while (n_bad < static_cast<int>(all.size()) &&
           !(all[static_cast<size_t>(n_bad)].q >= q_below))
      ++n_bad;
    const int n_out =
        force ? std::max(n_bad, std::min(n_worst, static_cast<int>(all.size()))) : n_bad;
    if (n_out == 0)
      return 0;
    std::error_code ec;
    std::filesystem::create_directories("dump", ec);
    std::ofstream log("dump/bad_triangles.log", std::ios::app);
    log << "frame " << _frame << " " << stage << ": " << n_bad << " of " << all.size()
        << " faces below quality " << q_below
        << (force ? " (solve blew up or did not converge: worst faces)" : "") << "\n";
    std::vector<char> hot(x.size(), 0);
    for (int i = 0; i < n_out; ++i) {
      const tri &t = all[static_cast<size_t>(i)];
      const vec3 &a = x[static_cast<size_t>(t.v[0])], &b = x[static_cast<size_t>(t.v[1])],
                 &c = x[static_cast<size_t>(t.v[2])];
      log << "  face " << t.f << " quality " << t.q << " area "
          << real(0.5) * (b - a).cross(c - a).norm() << " verts " << t.v[0] << " " << t.v[1]
          << " " << t.v[2] << " lengths " << (b - a).norm() << " " << (c - b).norm() << " "
          << (a - c).norm() << " x (" << a.transpose() << ") (" << b.transpose() << ") ("
          << c.transpose() << ")\n";
      for (int k : t.v)
        hot[static_cast<size_t>(k)] = 1;
    }
    const std::string obj =
        "dump/bad_triangles_f" + std::to_string(_frame) + "_" + std::string(stage) + ".obj";
    std::ofstream out(obj);
    std::vector<int> local(x.size(), -1);
    int nv = 0;
    auto vid = [&](int v) {
      if (local[static_cast<size_t>(v)] < 0) {
        local[static_cast<size_t>(v)] = ++nv;
        out << "v " << x[static_cast<size_t>(v)].transpose() << "\n";
      }
      return local[static_cast<size_t>(v)];
    };
    std::vector<char> is_bad(all.size(), 0);
    for (int i = 0; i < n_out; ++i)
      is_bad[static_cast<size_t>(i)] = 1;
    std::ostringstream bad_faces, ring_faces;
    for (size_t i = 0; i < all.size(); ++i) {
      const tri &t = all[i];
      if (!is_bad[i] && !hot[static_cast<size_t>(t.v[0])] && !hot[static_cast<size_t>(t.v[1])] &&
          !hot[static_cast<size_t>(t.v[2])])
        continue;
      std::ostringstream &o = is_bad[i] ? bad_faces : ring_faces;
      o << "f " << vid(t.v[0]) << " " << vid(t.v[1]) << " " << vid(t.v[2]) << "\n";
    }
    out << "g bad\n" << bad_faces.str() << "g ring\n" << ring_faces.str();
    log << "  wrote " << obj << "\n";
    std::cout << "[monge] frame " << _frame << " " << stage << ": " << n_bad
              << " triangle(s) below " << q_below << ", wrote dump/bad_triangles.log and " << obj
              << std::endl;
    return n_bad;
  }

  /// Same charts projection uses: `vertex_frames` with `ok`. No invent path.
  void draw_solved_frames() {
    if (!_cfg.draw.frames || !_cr || !__M || !__M->verts_are_dense_packed())
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    const int nV = __M->vert_count();
    if (static_cast<int>(x.size()) < nV)
      return;
    const std::vector<hepworth::block::vertex_chart> &frames = _cr->vertex_frames;
    if (frames.empty())
      return;
    const real l0 = asawa::shell::avg_length(*__M, x);
    const real nscale = std::max(l0 * _cfg.draw.length, real(1e-3));
    const vec4 red(1, 0, 0, 1);
    const vec4 green(0, 1, 0, 1);
    const vec4 blue(0, 0, 1, 1);
    for (int vi = 0; vi < nV; ++vi) {
      if (static_cast<size_t>(vi) >= frames.size())
        continue;
      const hepworth::block::vertex_chart &c = frames[static_cast<size_t>(vi)];
      hepworth::block::shape_w W;
      vec3 T, B, N;
      if (!hepworth::block::vertex_chart_at(&frames, vi, W, T, B, N))
        continue;
      const vec3 &xi = x[static_cast<size_t>(vi)];
      if (!xi.allFinite())
        continue;
      N = N.normalized();
      geometry_logger::line(xi - nscale * N, xi + nscale * N, blue);
      hepworth::block::shape_w axes = hepworth::block::shape_with_axes(W);
      if (!axes.ok) {
        geometry_logger::line(xi - nscale * T, xi + nscale * T, red);
        geometry_logger::line(xi - nscale * B, xi + nscale * B, green);
        continue;
      }
      auto draw_principal = [&](real ex, real ey, real k, const vec4 &col) {
        vec3 e = ex * T + ey * B;
        e -= N * e.dot(N);
        if (!(e.squaredNorm() > real(1e-16)) || !e.allFinite())
          return;
        e.normalize();
        real half = nscale;
        if (_cfg.draw.curvature) {
          const real ak = std::abs(k);
          if (!std::isfinite(ak))
            return;
          half = nscale * ak;
        }
        if (!(half > real(0)))
          return;
        geometry_logger::line(xi - half * e, xi + half * e, col);
      };
      draw_principal(axes.e1x, axes.e1y, axes.k1, red);
      draw_principal(axes.e2x, axes.e2y, axes.k2, green);
    }
  }

  /// World tangents on the edge datum after the split. Cyan, through the
  /// midpoint. Empty slots are the ones the next seed still has to fill.
  void draw_remeshed_edge_tangents() {
    if (!_cfg.draw.edge_tangents || !_cr || !_cr->tangents || !__M)
      return;
    const std::vector<vec3> &Tx = _cr->tangents->data();
    const std::vector<vec3> &px = asawa::const_get_vec_data(*__M, 0);
    const real nscale = std::max(
        asawa::shell::avg_length(*__M, px) * _cfg.draw.length * real(0.5),
        real(1e-3));
    const vec4 cyan(0, 0.85, 1, 1);
    for (const index_t c : __M->get_edge_range()) {
      const int slot = static_cast<int>(c) / 2;
      if (slot < 0 || static_cast<size_t>(slot) >= Tx.size())
        continue;
      vec3 dir = Tx[static_cast<size_t>(slot)];
      if (!dir.allFinite() || !(dir.squaredNorm() > real(1e-16)))
        continue;
      dir.normalize();
      const int v0 = static_cast<int>(__M->vert(asawa::shell::corner_id(c)));
      const int v1 =
          static_cast<int>(__M->vert(__M->other(asawa::shell::corner_id(c))));
      if (v0 < 0 || v1 < 0 || v0 >= static_cast<int>(px.size()) ||
          v1 >= static_cast<int>(px.size()))
        continue;
      const vec3 mid = real(0.5) * (px[static_cast<size_t>(v0)] +
                                    px[static_cast<size_t>(v1)]);
      geometry_logger::line(mid - nscale * dir, mid + nscale * dir, cyan);
    }
  }

  /// Edge-tangent polar on the solved positions, then rebuild golden
  /// `vertex_frames` from that field (same weingarten as projection).
  void transport_tangents() {
    if (!_shell || !_shell->dof_vec3())
      return;
    vecX q;
    _shell->dof_vec3()->map_to_x(q);
    for (const hepworth::cache_compute::ptr &compute : _solver._computes) {
      const auto tang =
          std::dynamic_pointer_cast<hepworth::block::tangent_cache_compute>(
              compute);
      if (tang)
        tang->compute(q);
    }
    if (!_cr || !__M || _cr->topology.nE <= 0 ||
        static_cast<int>(_cr->u.size()) < 2 * _cr->topology.nE)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    const int n = __M->vert_count();
    std::vector<vec3> H(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> wrote(static_cast<size_t>(n), 0);
    const std::vector<vec3> *gpar =
        _cr->g_par.empty() ? nullptr : &_cr->g_par;
    const std::vector<vec3> *gperp =
        _cr->g_perp.empty() ? nullptr : &_cr->g_perp;
    hepworth::block::vertex_weingarten(*__M, x, _cr->u, H, tang, axes, wrote,
                                       &_cr->vertex_frames, gpar, gperp);
  }

  /// |ωz| against the tilt |(ωx, ωy)| along one tube of the built-in torus,
  /// after relax(). ωz is the spin about the frame normal.
  void report_tube_spin(std::ostream &out) {
    const real R = 1;
    const real half = real(0.75) * (real(2) * M_PI / real(48));
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    struct seam_vert {
      int i = -1;
      real phi = 0;
    };
    std::vector<seam_vert> seam;
    const int n = __M->vert_count();
    for (int i = 0; i < n; ++i) {
      if (i >= static_cast<int>(x.size()))
        continue;
      const vec3 &p = x[static_cast<size_t>(i)];
      const real rho = std::hypot(p.x(), p.y());
      if (!(rho > real(1e-8)))
        continue;
      const real theta = std::atan2(p.y(), p.x());
      if (std::abs(theta) > half)
        continue;
      seam.push_back(
          seam_vert{i, std::atan2(p.z(), rho - R)});
    }
    std::sort(seam.begin(), seam.end(),
              [](const seam_vert &a, const seam_vert &b) { return a.phi < b.phi; });
    out << "tube seam vertices=" << seam.size()
        << "  w_z is spin about N, tilt is hypot(w_x, w_y)\n";
    int edges = 0;
    int flips = 0;
    real spin_abs = 0;
    real tilt_abs = 0;
    real prev_z = 0;
    bool have_prev = false;
    const int m = static_cast<int>(seam.size());
    for (int k = 0; k < m; ++k) {
      const int i = seam[static_cast<size_t>(k)].i;
      const int j = seam[static_cast<size_t>((k + 1) % m)].i;
      const auto ring = __M->get_one_ring(asawa::shell::vert_id(i));
      bool linked = false;
      for (asawa::shell::VertId vj : ring) {
        if (static_cast<int>(vj) == j) {
          linked = true;
          break;
        }
      }
      if (!linked)
        continue;
      vec3 Ti, Bi, Ni, Tj, Bj, Nj;
      if (!vertex_dirichlet_frame(i, x, Ti, Bi, Ni, nullptr) ||
          !vertex_dirichlet_frame(j, x, Tj, Bj, Nj, nullptr))
        continue;
      const quat qi =
          hepworth::block::frame_quat_from_axes(Ti, Bi, Ni, quat::Identity());
      const quat qj =
          hepworth::block::frame_quat_from_axes(Tj, Bj, Nj, qi);
      const hepworth::block::frame_body_omega w =
          hepworth::block::frame_body_omega_of(qi, qj);
      if (!w.ok)
        continue;
      const real tilt = std::hypot(w.x, w.y);
      out << "  phi " << seam[static_cast<size_t>(k)].phi << " -> "
          << seam[static_cast<size_t>((k + 1) % m)].phi << "  w_z " << w.z
          << "  tilt " << tilt << "\n";
      spin_abs += std::abs(w.z);
      tilt_abs += tilt;
      if (have_prev && w.z * prev_z < real(0))
        ++flips;
      prev_z = w.z;
      have_prev = true;
      ++edges;
    }
    out << "edges=" << edges;
    if (edges > 0) {
      out << "  mean |w_z|=" << spin_abs / real(edges)
          << "  mean tilt=" << tilt_abs / real(edges)
          << "  sign flips=" << flips;
    }
    out << "\n";
  }

  const hepworth::block::shell_cr_block::ptr &cr() const { return _cr; }
  /// Convergence of the last solver step (iterates, energies, last step).
  const hepworth::block::projection_solver::step_report &last_solve() const {
    return _solver.last_step;
  }
  /// Mean angle (radians, sign free) between the first carried stretch axis
  /// and the one its target's `angle` reads from the solved charts,
  /// after the last solve: how far the principal frame turned away from the
  /// material. NaN without a carried direction.
  real material_deflection() const { return _material_deflection; }

  /// The `w_gauss` rest `K0` per vertex (NaN until baked).
  const std::vector<real> &gauss_rest() { return gauss_rest_mut(); }
  /// Rebake `K0` from the live surface on the next frame.
  void reset_gauss() {
    for (real &k : gauss_rest_mut())
      k = std::numeric_limits<real>::quiet_NaN();
  }
  /// The last `stripes` solve: energy per unit mass (`lambda`),
  /// faces its read-out skipped, and whether it started cold.
  struct stripe_report {
    real lambda = std::numeric_limits<real>::quiet_NaN();
    int skipped = 0;
    bool cold = false;
  };
  const stripe_report &last_stripes() const { return _stripe_report; }

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;

private:
  /// One normal wave. Normals are taken before any vertex moves.
  static void perturb_vertices(asawa::shell::shell &M, std::vector<vec3> &x,
                               real amp_edges) {
    if (!(std::abs(amp_edges) > real(0)))
      return;
    const real l0 = asawa::shell::avg_length(M, x);
    const real amp = amp_edges * l0;
    const int n = static_cast<int>(x.size());
    std::vector<vec3> off(static_cast<size_t>(std::max(n, 0)), vec3::Zero());
    for (int i = 0; i < n; ++i) {
      if (static_cast<int>(M.vbegin(asawa::shell::vert_id(i))) < 0)
        continue;
      vec3 nrm = asawa::shell::vert_normal(M, asawa::shell::vert_id(i), x);
      if (!nrm.allFinite() || !(nrm.squaredNorm() > real(1e-16)))
        continue;
      nrm.normalize();
      const vec3 &p = x[static_cast<size_t>(i)];
      const real wave = std::sin(real(3) * p.x()) * std::cos(real(2) * p.y()) *
                        std::sin(real(5) * p.z());
      off[static_cast<size_t>(i)] = amp * nrm;
    }
    for (int i = 0; i < n; ++i)
      x[static_cast<size_t>(i)] += off[static_cast<size_t>(i)];
  }

  /// Independent normal offset. Normals are taken before any vertex moves.
  static void perturb_vertices_noise(asawa::shell::shell &M, std::vector<vec3> &x,
                                    real amp_edges) {
    if (!(std::abs(amp_edges) > real(0)))
      return;
    const real l0 = asawa::shell::avg_length(M, x);
    const real amp = amp_edges * l0;
    const int n = static_cast<int>(x.size());
    std::vector<vec3> off(static_cast<size_t>(std::max(n, 0)), vec3::Zero());
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<real> pick(real(-1), real(1));
    for (int i = 0; i < n; ++i) {
      if (static_cast<int>(M.vbegin(asawa::shell::vert_id(i))) < 0)
        continue;
      vec3 nrm = asawa::shell::vert_normal(M, asawa::shell::vert_id(i), x);
      if (!nrm.allFinite() || !(nrm.squaredNorm() > real(1e-16)))
        continue;
      nrm.normalize();
      off[static_cast<size_t>(i)] = amp * pick(rng) * nrm;
    }
    for (int i = 0; i < n; ++i)
      x[static_cast<size_t>(i)] += off[static_cast<size_t>(i)];
  }

  static bool ends_with(const std::string &s, const std::string &suf) {
    return s.size() >= suf.size() &&
           s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
  }

  static asawa::shell::shell::ptr mesh_from_tables(
      std::vector<vec3> vertices, std::vector<std::vector<int>> faces) {
    std::vector<index_t> corners_next, corners_vert, corners_face;
    asawa::assemble_table(vertices, faces, corners_next, corners_vert,
                          corners_face);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(
        corners_next, corners_vert, corners_face);
    auto vdata = asawa::datum_t<vec3>::create(asawa::prim_type::VERTEX, vertices);
    M->insert_datum(vdata);
    return M;
  }

  static asawa::shell::shell::ptr load_ascii_ply(const std::string &path) {
    std::ifstream in(path);
    if (!in)
      throw std::runtime_error("PLY asset '" + path + "' did not open");
    std::string line;
    std::getline(in, line);
    int nvert = 0;
    int nface = 0;
    int vert_props = 0;
    bool in_vert = false;
    while (std::getline(in, line)) {
      if (line == "end_header")
        break;
      if (line.rfind("element vertex", 0) == 0) {
        nvert = std::stoi(line.substr(std::string("element vertex").size()));
        in_vert = true;
      } else if (line.rfind("element face", 0) == 0) {
        nface = std::stoi(line.substr(std::string("element face").size()));
        in_vert = false;
      } else if (in_vert && line.rfind("property", 0) == 0) {
        ++vert_props;
      }
    }
    if (nvert <= 0 || nface <= 0 || vert_props < 3)
      throw std::runtime_error("PLY asset '" + path + "' has no triangles");
    std::vector<vec3> vertices;
    vertices.reserve(static_cast<size_t>(nvert));
    for (int i = 0; i < nvert; ++i) {
      double x = 0, y = 0, z = 0;
      in >> x >> y >> z;
      for (int p = 3; p < vert_props; ++p) {
        double skip = 0;
        in >> skip;
      }
      vertices.emplace_back(x, y, z);
    }
    std::vector<std::vector<int>> faces;
    faces.reserve(static_cast<size_t>(nface));
    for (int i = 0; i < nface; ++i) {
      int n = 0;
      in >> n;
      std::vector<int> f(static_cast<size_t>(std::max(n, 0)));
      for (int k = 0; k < n; ++k)
        in >> f[static_cast<size_t>(k)];
      if (n >= 3)
        faces.push_back(std::move(f));
    }
    if (!in)
      throw std::runtime_error("PLY asset '" + path + "' ended early");
    return mesh_from_tables(std::move(vertices), std::move(faces));
  }

  static asawa::shell::shell::ptr make_mesh(const std::string &path, int ring = 48,
                                            int tube = 24) {
    if (path.empty()) {
      const real R = 1.0;
      const real r = 0.33;
      auto M = asawa::shell::load_torus(R, r, std::max(ring, 3), std::max(tube, 3));
      return M;
    }
    if (ends_with(path, ".ply") || ends_with(path, ".PLY"))
      return load_ascii_ply(path);
    return asawa::shell::load_obj(path);
  }

  /// Fan normal and the average incident edge vector. `lock` keeps the red
  /// axis from flipping between draws.
  bool vertex_dirichlet_frame(int i, const std::vector<vec3> &x, vec3 &T,
                              vec3 &B, vec3 &N, const vec3 *lock) const {
    if (!_cr || !__M)
      return false;
    const int nE = _cr->count();
    const std::vector<real> &coeff = _cr->u;
    if (nE <= 0 || static_cast<int>(coeff.size()) < 2 * nE)
      return false;
    const asawa::shell::VertId vi = asawa::shell::vert_id(i);
    if (static_cast<int>(__M->vbegin(vi)) < 0)
      return false;
    N = asawa::shell::vert_normal(*__M, vi, x);
    if (!N.allFinite() || !(N.squaredNorm() > real(1e-12)))
      return false;
    N.normalize();
    const std::vector<int> &slot = _cr->topology.slot_map;
    T = vec3::Zero();
    __M->const_for_each_vertex(
        vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &S) {
          const int s = static_cast<int>(ci) / 2;
          if (s < 0 || s >= static_cast<int>(slot.size()))
            return;
          const int e = slot[static_cast<size_t>(s)];
          if (e < 0 || e >= nE)
            return;
          vec3 W = kusama::edge_guidance_vector_3d(
              S, ci, x, N, coeff[static_cast<size_t>(e)],
              coeff[static_cast<size_t>(e + nE)]);
          W -= N * W.dot(N);
          T += W;
        });
    if (!T.allFinite() || !(T.squaredNorm() > real(1e-16)))
      return false;
    if (lock && lock->allFinite() && T.dot(*lock) < real(0))
      T = -T;
    T.normalize();
    B = N.cross(T);
    if (!(B.squaredNorm() > real(1e-24)))
      return false;
    B.normalize();
    T = B.cross(N);
    return T.allFinite() && B.allFinite();
  }

  /// Per-vertex dual areas (a third of each incident face) and their
  /// area-weighted centroid.
  std::vector<real> vertex_areas(const std::vector<vec3> &x) const {
    std::vector<real> a(x.size(), real(0));
    for (auto fi : __M->get_face_range()) {
      const asawa::shell::FaceId fid = asawa::shell::face_id(fi);
      if (__M->fsize(fid) != 3)
        continue;
      const auto t = __M->get_tri(fid);
      const real A = asawa::shell::face_area(*__M, fid, x) / real(3);
      for (int k = 0; k < 3; ++k) {
        const int vi = static_cast<int>(t[k]);
        if (vi >= 0 && static_cast<size_t>(vi) < a.size())
          a[static_cast<size_t>(vi)] += A;
      }
    }
    return a;
  }

  vec3 area_centroid(const std::vector<vec3> &x) const {
    const std::vector<real> a = vertex_areas(x);
    vec3 c = vec3::Zero();
    real A = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      c += a[i] * x[i];
      A += a[i];
    }
    return A > real(0) ? vec3(c / A) : c;
  }

  /// `position_distortion::recenter`: translate back to `_center0` and drop
  /// the area-weighted mean velocity.
  void recenter_mesh() {
    if (!_cfg.position.recenter || !__M)
      return;
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    std::vector<vec3> &v = asawa::get_vec_data(*__M, 1);
    const std::vector<real> a = vertex_areas(x);
    vec3 c = vec3::Zero(), vm = vec3::Zero();
    real A = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      c += a[i] * x[i];
      if (i < v.size())
        vm += a[i] * v[i];
      A += a[i];
    }
    if (!(A > real(0)))
      return;
    const vec3 d = c / A - _center0;
    vm /= A;
    for (vec3 &p : x)
      p -= d;
    for (vec3 &u : v)
      u -= vm;
  }

  /// Torus of revolution through `x`: axis the least-spread direction about
  /// the centroid, `R` and `r` from the distances to it. False unless
  /// `0 < r < R`.
  static bool fit_torus(const std::vector<vec3> &x, real &R, real &r) {
    if (x.size() < 8)
      return false;
    vec3 c = vec3::Zero();
    for (const vec3 &p : x)
      c += p;
    c /= real(x.size());
    mat3 C = mat3::Zero();
    for (const vec3 &p : x)
      C += (p - c) * (p - c).transpose();
    Eigen::SelfAdjointEigenSolver<mat3> es(C);
    const vec3 axis = es.eigenvectors().col(0).normalized();
    real lo = std::numeric_limits<real>::max(), hi = 0;
    for (const vec3 &p : x) {
      const vec3 d = p - c;
      const real rho = (d - axis * d.dot(axis)).norm();
      lo = std::min(lo, rho);
      hi = std::max(hi, rho);
    }
    R = real(0.5) * (hi + lo);
    r = real(0.5) * (hi - lo);
    return r > real(0) && r < R;
  }

  /// Fits the startup torus once (`_torus_R`, `_torus_r`) for the stretch
  /// windings and logs the angles it gives.
  void set_stretch_angles(const std::vector<vec3> &x) {
    _torus_fit = fit_torus(x, _torus_R, _torus_r);
    auto log = [&](const char *what, const hepworth::block::stretch_winding &w, real angle) {
      if (w.ring != 0 || w.tube != 0)
        std::cout << "[monge] " << what << " winding (" << w.ring << ", " << w.tube << "): "
                  << (_torus_fit ? "torus R " : "no torus fit, R ") << _torus_R << " r "
                  << _torus_r << ", angle " << angle << std::endl;
    };
    for (const hepworth::block::stretch_target &s : stretch_now())
      log("stretch", s.winding, s.angle);
    for (const hepworth::block::shear_target &s : shear_now())
      log("shear", s.winding, s.angle);
  }

  /// The angle from the tube circle (the tight tangent) of a line winding
  /// `w` on the startup torus, `handedness * atan2(ring * sqrt(R^2 - r^2),
  /// tube * r)`; `fallback` without a fit or windings.
  real winding_angle(const hepworth::block::stretch_winding &w, real fallback) const {
    if (!_torus_fit)
      return fallback;
    return hepworth::block::loxodrome_angle(w, _torus_R, _torus_r, fallback);
  }

  /// The config's stretch list with each `winding` entry's angle set from the
  /// startup torus.
  std::vector<hepworth::block::stretch_target> stretch_now() const {
    std::vector<hepworth::block::stretch_target> out = _cfg.position.stretch;
    for (hepworth::block::stretch_target &s : out)
      s.angle = winding_angle(s.winding, s.angle);
    return out;
  }

  /// The config's shear list with each `winding` entry's angle set from the
  /// startup torus.
  std::vector<hepworth::block::shear_target> shear_now() const {
    std::vector<hepworth::block::shear_target> out = _cfg.position.shear;
    for (hepworth::block::shear_target &s : out)
      s.angle = winding_angle(s.winding, s.angle);
    return out;
  }

  /// Per stretch target its axis field for the step (`direction`), null for
  /// `angle`. A `material` target carries its field: seeded from the
  /// step-start charts at the target's angle and oriented, later leaned
  /// toward those charts by `material_follow` and smoothed, with the
  /// step-start positions kept for `transport_material_dirs`. A `stripes`
  /// target takes the oriented chart directions and closes them with its
  /// stripe pattern.
  std::vector<const std::vector<vec3> *> stretch_dirs_now() {
    std::vector<const std::vector<vec3> *> out(_material_slot.size(), nullptr);
    if (!_cr || !_shell || !__M->verts_are_dense_packed())
      return out;
    const std::vector<vec3> &x = _shell->xs->get();
    const std::vector<hepworth::block::stretch_target> targets = stretch_now();
    const real t = std::pow(_cfg.position.material_smooth * asawa::shell::avg_length(*__M, x), 2);
    for (size_t k = 0; k < _material_slot.size() && k < targets.size(); ++k) {
      if (_stripe_slot[k] >= 0 && _cfg.position.stripes > 0.0) {
        std::vector<vec3> &u = _stripe_dirs[k];
        u = hepworth::block::chart_directions(_cr->vertex_frames, targets[k].angle);
        u.resize(x.size(), vec3::Zero());
        kusama::material_direction::orient(*__M, u);
        close_with_stripes(u, _stripe_slot[k], x);
        out[k] = &u;
        continue;
      }
      if (_material_slot[k] < 0)
        continue;
      std::vector<vec3> &u = asawa::get_vec_data(*__M, _material_slot[k]);
      if (u.size() != x.size())
        continue;
      if (!_material_seeded[k]) {
        u = hepworth::block::chart_directions(_cr->vertex_frames, targets[k].angle);
        u.resize(x.size(), vec3::Zero());
        kusama::material_direction::orient(*__M, u);
        _material_seeded[k] = 1;
      } else if (_cfg.position.material_follow > 0.0) {
        kusama::material_direction::lean_toward(
            u, hepworth::block::chart_directions(_cr->vertex_frames, targets[k].angle),
            _cfg.position.material_follow);
      }
      kusama::material_direction::smooth(*__M, x, u, t);
      out[k] = &u;
      _material_x0 = x;
    }
    return out;
  }

  /// Replaces `u` by the level lines of the stripe pattern across it
  /// (`stripes`), warm-started from the pattern in vertex datum `slot`
  /// (`(Re, Im, 0)`; all zero is a cold start). Reversing `u` conjugates the
  /// pattern, so the warm start is the carried pattern or its conjugate,
  /// whichever has the lower energy. Vertices the pattern cannot read (zeros,
  /// field flips) keep `u`; a failed solve keeps all of it and drops the
  /// pattern.
  void close_with_stripes(std::vector<vec3> &u, index_t slot, const std::vector<vec3> &x) {
    if (slot < 0)
      return;
    const size_t n = x.size();
    std::vector<vec3> across(n, vec3::Zero());
    for (size_t i = 0; i < n && i < u.size(); ++i) {
      const vec3 nrm = kusama::material_direction::unit_vert_normal(*__M, static_cast<int>(i), x);
      const vec3 c = nrm.cross(u[i]);
      if (c.squaredNorm() > real(1e-24))
        across[i] = c.normalized();
    }
    const real nu = real(2) * real(M_PI) /
                    (_cfg.position.stripes * asawa::shell::avg_length(*__M, x));
    const kusama::stripes::system S = kusama::stripes::build(*__M, x, across, nu);
    std::vector<vec3> &P = asawa::get_vec_data(*__M, slot);
    P.resize(n, vec3::Zero());
    vecX psi(2 * static_cast<int>(n));
    bool cold = true;
    for (size_t i = 0; i < n; ++i) {
      psi[2 * i] = P[i].x();
      psi[2 * i + 1] = P[i].y();
      cold &= !(P[i].squaredNorm() > real(0));
    }
    if (cold) {
      psi.resize(0);
    } else {
      vecX conj = psi;
      for (size_t i = 0; i < n; ++i)
        conj[2 * i + 1] = -conj[2 * i + 1];
      if (kusama::stripes::rayleigh(S, conj) < kusama::stripes::rayleigh(S, psi))
        psi = conj;
    }
    _stripe_report = {};
    _stripe_report.cold = cold;
    _stripe_report.lambda = kusama::stripes::solve(S, psi, cold ? 40 : 2);
    if (!std::isfinite(_stripe_report.lambda)) {
      std::fill(P.begin(), P.end(), vec3::Zero());
      return;
    }
    for (size_t i = 0; i < n; ++i)
      P[i] = vec3(psi[2 * i], psi[2 * i + 1], 0);
    const std::vector<vec3> d =
        kusama::stripes::level_directions(*__M, x, S, psi, &_stripe_report.skipped);
    for (size_t i = 0; i < n && i < u.size(); ++i)
      if (d[i].squaredNorm() > real(0))
        u[i] = d[i].dot(u[i]) < real(0) ? vec3(-d[i]) : d[i];
  }

  /// Turns the carried directions with the solve (`stretch_dirs_now` kept
  /// the step-start positions; same topology until the remesh).
  void transport_material_dirs() {
    if (_material_x0.empty() || !__M)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    if (x.size() == _material_x0.size())
      for (const index_t slot : _material_slot)
        if (slot >= 0)
          kusama::material_direction::transport(*__M, _material_x0, x,
                                                asawa::get_vec_data(*__M, slot));
    _material_x0.clear();
  }

  /// `material_deflection` after the solve, before remesh.
  void measure_material_deflection() {
    _material_deflection = std::numeric_limits<real>::quiet_NaN();
    if (!_cr || !__M)
      return;
    const std::vector<hepworth::block::stretch_target> targets = stretch_now();
    for (size_t k = 0; k < _material_slot.size() && k < targets.size(); ++k) {
      if (_material_slot[k] < 0)
        continue;
      const std::vector<vec3> &u = asawa::get_vec_data(*__M, _material_slot[k]);
      const std::vector<vec3> c =
          hepworth::block::chart_directions(_cr->vertex_frames, targets[k].angle);
      real sum = 0;
      int n = 0;
      for (size_t i = 0; i < u.size() && i < c.size(); ++i) {
        if (!(u[i].squaredNorm() > real(0)) || !(c[i].squaredNorm() > real(0)))
          continue;
        sum += std::acos(std::clamp(std::abs(u[i].dot(c[i])), real(0), real(1)));
        ++n;
      }
      if (n > 0)
        _material_deflection = sum / real(n);
      return;
    }
  }

  /// Orange segments along each carried stretch axis
  /// (`fairing_draw::material`), through the vertex.
  void draw_material_dirs() {
    if (!_cfg.draw.material || !__M)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    const real half = std::max(asawa::shell::avg_length(*__M, x) * _cfg.draw.length, real(1e-3));
    const vec4 orange(1, 0.55, 0, 1);
    for (const index_t slot : _material_slot) {
      if (slot < 0)
        continue;
      const std::vector<vec3> &u = asawa::get_vec_data(*__M, slot);
      for (size_t i = 0; i < x.size() && i < u.size(); ++i)
        if (u[i].allFinite() && u[i].squaredNorm() > real(0))
          geometry_logger::line(x[i] - half * u[i], x[i] + half * u[i], orange);
    }
  }

  monge_fairing_demo_config _cfg;
  /// Vertex datum per stretch target with `material` (-1 without), its
  /// seeded flag, and the step-start positions while a solve is open.
  std::vector<real> &gauss_rest_mut() {
    if (_gauss_slot < 0)
      _gauss_slot = asawa::init_vert_datum<real>(*__M, std::numeric_limits<real>::quiet_NaN());
    return asawa::get_real_data(*__M, _gauss_slot);
  }
  index_t _gauss_slot = -1;
  std::vector<index_t> _material_slot;
  /// Vertex datum per `stripes` stretch target holding its stripe pattern
  /// (-1 without), and its closed axes for the step.
  std::vector<index_t> _stripe_slot;
  std::vector<std::vector<vec3>> _stripe_dirs;
  stripe_report _stripe_report;
  std::vector<char> _material_seeded;
  std::vector<vec3> _material_x0;
  real _material_deflection = std::numeric_limits<real>::quiet_NaN();
  bool _torus_fit = false;
  real _torus_R = 0, _torus_r = 0;
  vec3 _center0 = vec3::Zero();
  int _frame = 0;
  hepworth::block::shell_position_block::ptr _shell;
  hepworth::block::shell_cr_block::ptr _cr;
  hepworth::block::block_solver_config<hepworth::block::shell_position_block,
                                       hepworth::block::shell_cr_block>
      _solver_config;
  hepworth::block::projection_solver _solver;
};

} // namespace duchamp
} // namespace gaudi
