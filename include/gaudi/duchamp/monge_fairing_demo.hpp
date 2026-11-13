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
#include "gaudi/hepworth/block/amap_gn.hpp"
#include "gaudi/hepworth/block/contact_net_gn.hpp"
#include "gaudi/hepworth/block/dupin_circle.hpp"
#include "gaudi/hepworth/block/edge_shape_fair.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"
#include "gaudi/hepworth/block/face_element.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/kite_growth.hpp"
#include "gaudi/hepworth/block/osculating_torus_ring.hpp"
#include "gaudi/hepworth/block/shape_det_conserve.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/hepworth/blocks/shell_cr_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/blocks/shell_shape_block.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"
#include "gaudi/kusama/material_direction.hpp"
#include "gaudi/kusama/stripes.hpp"

namespace gaudi {
namespace duchamp {

struct fairing_step {
  real h = 0.1;
  real damping = 0.1;
  int iterations = 10;
  /// Gauss–Newton rows keep their factor this many iterates (fresh gradient
  /// every iterate); 1 refactors every iterate, 0 only on ascent.
  int gn_refresh_every = 15;
  /// Step halvings when an iterate raises the energy (0: full GN steps).
  /// Stiff curvature rows diverge without it.
  int gn_backtracks = 1;
  /// On an iterate that raises the energy, refactor at the previous iterate
  /// and keep going (true), or stop the step there with an imperfect solve
  /// and no refactor (false). Either way the step returns its best iterate.
  bool gn_refresh_on_ascent = true;
  /// With `gn_refresh_on_ascent`: a step from a fresh factor that still goes
  /// uphill ends the frame at its best iterate (the Gauss-Newton model is
  /// wrong there; more steps compound it).
  bool gn_stop_on_fresh_ascent = true;
  /// With `gn_refresh_on_ascent`: refactors on ascent per frame; past it an
  /// ascent ends the frame at its best iterate (bounds the frame at about
  /// `iterations / gn_refresh_every + 1 + this` factorizations). Negative: no
  /// cap.
  int gn_max_ascent_refreshes = 1;
  /// Relative energy rise an iterate may take before it counts as an
  /// ascent. A frozen factor converges to its model's minimum in a few
  /// iterates and then wobbles a few percent; this keeps iterating through
  /// that (the step still returns its best iterate). 0 stops on any rise.
  real gn_ascent_tol = 0.1;
  /// Log which constraints minimize an energy at the start of each step
  /// (`hepworth::block::audit_variational`); slow, for diagnosis.
  bool audit_variational = false;
  /// The step: `ldlt` (the factored solve) or `modal` (the same with the
  /// low-mode Newton / buckling correction, configured by `modal`).
  hepworth::block::step_preconditioner preconditioner = hepworth::block::step_preconditioner::ldlt;
  /// The iterate corrector after each global step (`projection_solver::accel`):
  /// `none` (default), `anderson`, `chebyshev`, `anderson_chebyshev`
  /// (Chebyshev only while Anderson's history is short). A corrected iterate
  /// is taken only if it does not raise the energy. Measured on the probe it
  /// did not change the frames and doubled the time of the stateful couplings
  /// (MATH.md 19.12), so it is off. `anderson_m` the window, `cheb_delay` /
  /// `cheb_window` when and over how many iterates Chebyshev estimates its
  /// rate.
  hepworth::block::step_accel accel = hepworth::block::step_accel::anderson;
  int anderson_m = 5;
  int cheb_delay = 10;
  int cheb_window = 5;
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
  real lambda = 200.0;
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
///  off         — nothing; the store holds its frame-start value.
///  contact_net — the principal contact element net (MATH.md 19.10 to 19.12):
///                positions and normals solved together as Gauss-Newton rows
///                in the projection solver, one sphere per edge through the
///                Rodrigues relation with twist, the store's `(k11, k12)` per
///                edge as its targets.
enum class curvature_solver { off, contact_net };

/// `w_coupling` is the contact net's row weight between position space and
/// shape space: `shape_distortion` rows act on the edge store,
/// `position_distortion` rows on the positions, and the coupling keeps the
/// two in sync. A shape row at the same weight lands the store halfway; keep
/// it at or above the position rows that hold curvature (`w_bend`) for shape
/// edits to show in the mesh. Without it (`off`) the store reads its
/// frame-start value at weight 1.
struct curvature_config {
  /// `curvature_solver::off` or `curvature_solver::contact_net`.
  curvature_solver kind = curvature_solver::contact_net;
  real w_coupling = 0.5;
  /// Tangency target: `vertex` (default) from the store's vertex tensor;
  /// `edge` the edge's `k` plus the frame start's one-sided offset.
  hepworth::block::contact_net_gn::tangency contact_tangency = hepworth::block::contact_net_gn::tangency::vertex;
  /// Normals' mass per unit vertex dual area: a pull toward the frame-start
  /// normals that damps steps on a stale factor (0: none).
  real contact_normal_mass = 0.0;
  /// One quaternion frame per vertex in place of the normals (MATH.md 19.21,
  /// E3): the normal is `Y e3 conj(Y) / |Y|^2`. Lands where the normals do
  /// on the contact-net suite. Off by default.
  bool contact_frames = false;
};

/// Shape → position: the osculating-torus ring row and the curvature solver.
/// `one_ring_operators` builds stretch, tilt and axis ratio on the vertex
/// one-ring. Off keeps stretch on the face.
struct shape_project {
  bool one_ring_operators = true;
  /// Osculating-torus ring row. The store sets both curvatures and the
  /// principal angle at each foot. The 1-ring and 2-ring fit only the pose,
  /// the foot's normal offset and the normal's tilt. Every ring point lands
  /// on that posed torus. 0 is off.
  real w_osculating_torus = 0.0;
  /// Row scale (`two_ring`), pose iterations, and pose prior. `.volume`
  /// `corotational` keeps each one-ring fan prism's rest volume (to its
  /// tube-axis points) while landing; `off` lands only.
  hepworth::block::osculating_torus_ring_config osculating_torus{
      .two_ring = true,
      .volume = {.mode = hepworth::block::ring_volume_mode::corotational,
                 .corotational = {.w_proj = 1.0, .w_rest = 0, .w_stay = 0, .w_volume = 5.0e-1}},
  };
  /// Edge curvature store → positions.
  curvature_config curvature{};
};

/// Curvature edits on the edge curvature store (one 2x2 curvature tensor per
/// edge). The mesh follows through `shape_project.curvature.w_coupling`: at
/// equal weights the store lands halfway between the mesh and the edit.
struct shape_distortion {
  /// Shape ARAP: holds each edge's curvature at its frame-start value (how
  /// anisotropic it is, and the pattern of its lines) but lets each patch
  /// turn its curvature lines together. Edits pushed harder than this have
  /// to give somewhere. 0 is off.
  real w_shape_arap = 0.0;
  /// Shape ARAP also holds the mean curvature; off holds only the
  /// anisotropic part.
  bool shape_arap_mean = true;

  /// The traceless spin-2 Dirichlet rows toward zero: smooths the principal
  /// direction field weighted by anisotropy, so directions turn freely near
  /// umbilics; mean curvature free. The Dupin rows are its directional
  /// restriction. 0 leaves it off.
  real w_grad_traceless = 0.0;

  /// The frame minimizer: pin toward each edge's principal curvatures set to
  /// the neighbors' normal curvature along the same directions. 0 leaves it
  /// off.
  real w_minimizer = 0.5;
  /// Surface tension in curvature space (`goal_mean_tension`): each edge's
  /// mean curvature pulled to `tension.h0` (`constant`, Young-Laplace:
  /// pressure over twice the tension; outward normals make a sphere's `H`
  /// negative) or to the store's mean (`mean`, constant-mean-curvature
  /// fairing), trace-free part kept. 0 leaves it off.
  real w_tension = 0.15;
  hepworth::block::tension_target tension{};
  /// Dupin rows: each principal curvature held constant along its own
  /// direction, so curvature lines become circles and the surface is locally
  /// a Dupin cyclide (torus, cylinder, cone and their Moebius images). A
  /// fixed quadratic per frame (frames frozen at the frame start). 0 leaves
  /// it off.
  real w_dupin = 0.0;
  /// The Dupin rows' stencil: 0 reads each principal curvature's derivative
  /// on one face; above 0 traces each curvature line `dupin_reach` mean edge
  /// lengths both ways and holds the residual of a degree-`dupin_degree` fit
  /// of its normal curvature along the line (0: constant, the line a circle;
  /// 1: linear change allowed). Longer reaches settle lower-frequency
  /// departures from Dupin.
  real dupin_reach = 0.0;
  int dupin_degree = 0;
  /// Directional curvature growth, the shape-space match of position
  /// `stretch`: each frame the normal curvature along `.angle` (radians from
  /// the tight direction toward the other) times `1 + x`, across it times
  /// `1 + y`, against the frame-start store; negative bends less that way.
  /// `.power` scales each fraction by `(|k| / mean |k|)^power` of that axis's
  /// frame-start curvature, so positive grows the more curved edges faster.
  /// `.winding` sets the angle at startup like stretch. Comment a line to
  /// drop it. `w_growth` scales every entry; 0 leaves the row off.
  std::vector<hepworth::block::curvature_growth_target> growth{
      hepworth::block::curvature_growth_target{.angle = 1.0,
                                               .x = 1.0,
                                               .y = 0.0,
                                               .winding = {.ring = 1, .tube = 5, .handedness = 1},
                                               .power = 0.5},
  };
  real w_growth = 0.0;
  /// Holds each vertex ring's total Gauss-image area `sum K dA` of the store
  /// at its frame-start value and lets it redistribute between the ring's
  /// edges and between their mean curvature and anisotropy
  /// (`ring_det_conserve`), the shape-space match of `area`. `K` is
  /// intrinsic, so this leaves bending free and resists only what would
  /// change the metric. 0 leaves it off.
  real w_det = 0.0;
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
          .winding = {.ring = 1, .tube = 5, .handedness = 1},
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
          .gain = 0.1,
          .metric = hepworth::block::curvature_metric::aniso,
          .power = 1.0,
      },
  };
  real w_tilt = 0.0;
  /// Locks each one ring's covariant anisotropy: spoke rms along the tight
  /// principal tangent over that along the other held at `ratio` (1 is
  /// round). Area, skew and rotation are free, so it composes with `area`.
  /// Comment a line to drop it. `w_axis_ratio` scales every entry; 0 leaves
  /// the row off.
  std::vector<hepworth::block::axis_ratio_target> axis_ratio{
      hepworth::block::axis_ratio_target{.ratio = 0.5},
  };
  real w_axis_ratio = 0.55;
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
  real w_bend = 1.0e-1;
  /// Per-triangle rest area: weight `.w` (0 leaves it off) and `.band`, the
  /// range of its step-start area each face may take freely (`{1, 1}`
  /// conserves it; `{0.5, 1.5}` lets it drift 50% per step before pushing
  /// back).
  hepworth::block::area_config area{.w = 1.0e-1, .band = {0.9, 1.1}};
  /// Cotangent Laplacian on positions, target Δx = 0. 0 leaves it off.
  real w_smooth = 0.0e-1;
  /// Per-triangle conformal (closest similarity) against the start-of-step
  /// shape. Scale is free, shear pays. 0 leaves it off.
  real w_conformal = 0.0;
  /// As-Moebius-as-possible rows (Vaxman, Mueller, Weber 2015; `amap_gn`):
  /// per-vertex quaternions `Y_i` and, per edge, `x_k - x_i` against
  /// `conj(Y_i) q_ik Y_k` (`q_ik` the frame-start edge). Vertex-based and
  /// Moebius-exact, where `w_conformal` is per face and first order.
  /// `amap_inversion` damps the inversions (scale variation; the paper's
  /// `alpha_INV`, relative). 0 leaves it off.
  real w_amap = 0.0;
  real amap_inversion = 0.1;
  /// Position-space Dupin rows (`dupin_circle.hpp`): four consecutive points
  /// along a curvature line fitted to a planar circle and projected onto it.
  /// One entry per row set; `.kind` picks where on the circle the points go:
  ///   closest     the plain Dupin row (closest point on the live circle)
  ///   slide       advance `.delta` mean edge lengths of arc per frame
  ///   even        ends kept, inner points at thirds of the arc
  ///   graded      ends kept, arcs in ratio 1 : `.g` : `.g`^2
  ///   hyperbolic  Moebius regrading toward / away from the group's middle
  ///               (`.lambda`), cross ratio kept
  ///   parabolic   Moebius slide (`.beta`), cross ratio kept
  ///   loxodromic  the line turned about its middle vertex by `.psi` on its
  ///               curvature sphere (and scaled by `.lambda`)
  /// `.family` 0 both, 1 tight, 2 loose. Direction-dependent kinds (slide,
  /// graded, hyperbolic, parabolic) take their sign from `dupin_circle_ref`;
  /// without it they
  /// add no rows. `w_dupin_circle` scales every entry; 0 leaves them off.
  /// One preset per mode; keep the one you want (several can be on at once).
  /// Fields must be written in this order: `.kind`, `.family`, `.w`, `.delta`,
  /// `.g`, `.lambda`, `.beta`, `.psi`, `.max_angle`, `.min_aniso`.
  std::vector<hepworth::block::dupin_circle_spec> dupin_circle{
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::closest},
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::slide, .delta = 0.2},
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::even},
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::graded, .g = 1.5},
      hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::hyperbolic, .lambda = 2.2},
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::parabolic, .beta = 0.2},
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::loxodromic, .psi = 0.3},
      // Tight family only (1) or loose only (2), e.g.:
      //hepworth::block::dupin_circle_spec{.kind = hepworth::block::dupin_circle_spec::mode::loxodromic, .family = 2, .psi = 0.1},
  };
  real w_dupin_circle = 0.2;
  vec3 dupin_circle_ref = vec3::Zero();
  /// Surface tension in position space: per-triangle area rows with a zero
  /// target (`area_mode::zero`, the area gradient is the mean-curvature
  /// normal). 0 leaves it off.
  real w_tension = 0.0;
  /// Moebius kites as an edge growth (`kite_growth.hpp`): a vertex log-scale
  /// `u` per frame, then per edge a length row at `e^{(u_i + u_j) / 2}` times
  /// its frame-start length (triangles change shape where their vertices'
  /// `u` differ). `.source`:
  ///   mean_curvature  u = alpha (-H - level): each element's Moebius map on
  ///                   its mean-curvature sphere; convex grows, concave
  ///                   shrinks, saddles (H ~ 0) stay
  ///   anisotropy      u = alpha (|z| - level)
  ///   tangent         the Dirichlet edge tangents turned by `.angle` (0
  ///                   along them, trapezoids; pi / 4, kites), times
  ///                   alpha |z|, through one Poisson solve (always mean zero;
  ///                   on a torus they circulate and give almost no u)
  /// `.balanced` removes the mean of a scalar source (redistributes scale, no
  /// net growth). `.max_u` caps the per-frame `|u|`. `w_kite` scales the
  /// rows; 0 leaves them off.
  /// Fields must be written in this order: `.alpha`, `.angle`, `.max_u`,
  /// `.source`, `.level`, `.balanced`.
  hepworth::block::kite_target kite{
      .alpha = 0.1,              // log-scale per unit curvature per frame (tangent: per unit |z|)
      .angle = real(M_PI) / 4.0, // tangent source only: 0 along them (trapezoids), pi / 4 kites
      .max_u = 0.1,              // cap on the per-frame |u|
      .source = hepworth::block::kite_source::anisotropy, // mean_curvature, anisotropy, or tangent
      .level = 0.0,              // curvature at which nothing changes (scalar sources)
      .balanced = false,         // remove the mean: no net growth (scalar sources)
  };
  real w_kite = 0.0;
  /// Intrinsic vector Dirichlet memory of the metric (`dirichlet_memory_gn`,
  /// rest at the frame start; factor the mean edge length). Contact net only.
  real w_memory = 0.5;
  /// Holds the area-weighted centroid where the scene started (translation
  /// only; shape and rotation are free). 0 is off.
  real w_center = 0.5;
  /// After each step, move the mesh back so its centroid is exactly where it
  /// started, and remove the mean velocity. Exact and free, outside the solve.
  bool recenter = false;
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

  static ptr create(const monge_fairing_demo_config &cfg = k_monge_fairing_demo_default) {
    return std::make_shared<monge_fairing_demo>(cfg);
  }

  explicit monge_fairing_demo(const monge_fairing_demo_config &cfg = {}) : _cfg(cfg) { reset_scene(); }

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
    // Persistent one-ring medial depth for the torus row's volume term,
    // carried through the remesher (NaN until first seen).
    _torus_depth_slot = asawa::init_vert_datum<real>(*__M, std::numeric_limits<real>::quiet_NaN());
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
    _shell = std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);
    _cr = std::make_shared<hepworth::block::shell_cr_block>(__M);
    _cr->face_site = !_cfg.project.one_ring_operators || _cfg.dirichlet.w_pin > 0.0;
    const bool net = _cfg.project.curvature.kind == curvature_solver::contact_net;
    _cr->normal_site = net && !_cfg.project.curvature.contact_frames;
    _cr->frame_site = net && _cfg.project.curvature.contact_frames;
    _cr->normal_mass_scale = _cfg.project.curvature.contact_normal_mass;
    _cr->amap_site = _cfg.position.w_amap > 0.0;
    real C = 3.0;
    const real l0 = asawa::shell::avg_length(*__M, xs->get());
    __surf = asawa::shell::dynamic::create(__M, C * l0, real(1.5) * C * l0, C * real(1.0) * l0);
    __surf->set_qc(_cfg.qc);
    __surf->set_qc_handler([this](asawa::shell::shell &, const asawa::shell::qc_report &rep) {
      if (!rep.bad_faces.empty())
        dump_bad_triangles("qc_left", false, 8, _cfg.qc.q_warn);
    });

    _solver_config =
        hepworth::block::block_solver_builder<hepworth::block::shell_position_block,
                                              hepworth::block::shell_cr_block>::create()
            .with_blocks(_shell, _cr)
            .with_presolve([this](hepworth::block::solver_context &) {
              if (_cr)
                hepworth::block::seed_cr_from_dirichlet(*__M, *_cr, _cfg.dirichlet.reseed,
                                                        _cfg.dirichlet.lambda, _cfg.dirichlet.frame_angle);
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.area.w > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_area(*__M, ctx.constraints, _shell->xs->get(), _cfg.position.area,
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
              if (!(_cfg.position.w_conformal > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_triangle_conformal(*__M, ctx.constraints, _shell->xs->get(),
                                                       _cfg.position.w_conformal,
                                                       hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_amap > 0.0) || !_cr || !_cr->amap_dof() || ctx.blocks.empty() ||
                  !__M->verts_are_dense_packed())
                return;
              const std::vector<vec3> &x = _shell->xs->get();
              if (_cr->table.n_verts != static_cast<int>(x.size()))
                return;
              ctx.constraints.push_back(hepworth::block::amap_gn::create(
                  &_cr->table, x, _cfg.position.w_amap, _cfg.position.amap_inversion,
                  std::vector<hepworth::sim_block::ptr>{ctx.blocks[0], _cr->amap_dof()}));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_dupin_circle > 0.0) || !_cr || ctx.blocks.empty() ||
                  !__M->verts_are_dense_packed())
                return;
              const std::vector<vec3> &x = _shell->xs->get();
              if (_cr->table.n_verts != static_cast<int>(x.size()) || _cr->vertex_frames.size() != x.size())
                return;
              std::vector<vec3> ref;
              if (_cfg.position.dupin_circle_ref.squaredNorm() > 0.0) {
                ref.assign(x.size(), vec3::Zero());
                for (size_t i = 0; i < x.size(); ++i) {
                  const vec3 n = _cr->vertex_frames[i].N.normalized();
                  ref[i] = _cfg.position.dupin_circle_ref - n * n.dot(_cfg.position.dupin_circle_ref);
                  if (ref[i].squaredNorm() > 0.0)
                    ref[i].normalize();
                }
              }
              hepworth::block::init_dupin_circles(_cr->table, x, _cr->vertex_frames, _cfg.position.dupin_circle,
                                                  _cfg.position.w_dupin_circle, ref.empty() ? nullptr : &ref,
                                                  ctx.constraints, hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_kite > 0.0) || !_cr || ctx.blocks.empty() || !__M->verts_are_dense_packed())
                return;
              const std::vector<vec3> &x = _shell->xs->get();
              if (_cr->table.n_verts != static_cast<int>(x.size()) ||
                  _cr->shape.size() != static_cast<size_t>(_cr->table.n_edges))
                return;
              const std::vector<real> u = hepworth::block::init_kite_growth(
                  _cr->table, x, _cr->edge_tangent_dirs(x), _cr->carried_edge_tangents(), _cr->shape,
                  _cfg.position.kite, _cfg.position.w_kite, ctx.constraints,
                  hepworth::block::select_blocks<0>(ctx));
              if (!u.empty()) {
                const auto mm = std::minmax_element(u.begin(), u.end());
                _kite_u_range = *mm.second - *mm.first;
              }
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_tension > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_area(*__M, ctx.constraints, _shell->xs->get(), _cfg.position.w_tension,
                                         hepworth::block::select_blocks<0>(ctx), hepworth::block::area_mode::zero);
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (_cfg.position.w_smooth <= 0.0)
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_laplacian(*__M, ctx.constraints, _shell->xs->get(),
                                              hepworth::block::laplacian_mode::null,
                                              hepworth::block::laplacian_stencil::cotan,
                                              _cfg.position.w_smooth, hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              ctx.shell_mesh = __M.get();
              if (!__M->verts_are_dense_packed())
                return;
              if (_cr && !ctx.blocks.empty()) {
                hepworth::block::init_shell_caches(*__M, *_cr, ctx.blocks[0], ctx.computes);
                hepworth::block::init_face_shape_pin(*__M, *_cr, ctx.blocks[0], ctx.constraints,
                                                     _cfg.dirichlet.w_pin);
                if (_cfg.project.one_ring_operators) {
                  std::vector<char> live(static_cast<size_t>(__M->vert_count()), 1);
                  hepworth::block::init_directed_stretch(
                      *__M, ctx.constraints, _shell->xs->get(), live, {ctx.blocks[0]},
                      _cfg.position.w_stretch, stretch_now(), &_cr->vertex_frames,
                      _cfg.position.w_stretch > 0.0 ? stretch_dirs_now()
                                                    : std::vector<const std::vector<vec3> *>{});
                  hepworth::block::init_directed_tilt(*__M, ctx.constraints, _shell->xs->get(), live,
                                                      {ctx.blocks[0]}, _cfg.position.w_tilt,
                                                      _cfg.position.tilt, &_cr->vertex_frames);
                  hepworth::block::init_axis_ratio(*__M, ctx.constraints, _shell->xs->get(), live,
                                                   {ctx.blocks[0]}, _cfg.position.w_axis_ratio,
                                                   _cfg.position.axis_ratio, &_cr->vertex_frames);
                } else {
                  hepworth::block::init_face_stretch(*__M, *_cr, ctx.blocks[0], ctx.constraints,
                                                     _cfg.position.w_stretch, stretch_now());
                }
              }
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.position.w_bend > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_bending(*__M, ctx.constraints, _shell->xs->get(), _cfg.position.w_bend,
                                            hepworth::block::select_blocks<0>(ctx));
            })
            .with_recompute([this](hepworth::block::solver_context &ctx) { push_curvature(ctx); })
            .with_recompute([this](hepworth::block::solver_context &ctx) {
              if (!(_cfg.project.w_osculating_torus > 0.0) || !_cr || ctx.blocks.empty() ||
                  !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hepworth::block::init_osculating_torus_ring(
                  *__M, ctx.constraints, _cfg.project.w_osculating_torus, {ctx.blocks[0], _cr->shape_dof()},
                  &_cr->vertex_frames, &_cr->stars, _cfg.project.osculating_torus,
                  &asawa::get_real_data(*__M, _torus_depth_slot));
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
    if (!_cr || !_shell || ctx.blocks.empty() || !_cr->shape_dof() || !__M->verts_are_dense_packed() ||
        _cr->table.n_edges <= 0)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    if (static_cast<int>(x.size()) != _cr->table.n_verts)
      return;
    ctx.shell_mesh = __M.get();
    hepworth::sim_block::ptr pos = ctx.blocks[0];
    hepworth::sim_block::ptr normals = c.contact_frames ? hepworth::sim_block::ptr(_cr->frame_dof())
                                                        : hepworth::sim_block::ptr(_cr->normal_dof());
    if (c.kind == curvature_solver::contact_net && c.w_coupling > 0.0 && normals && _cr->normal.size() == x.size()) {
      // The net's rows couple positions, normals (or frames) and the store in
      // the one solve; the store answers to them and its shape rows.
      std::vector<hepworth::sim_block::ptr> net_blocks{pos, normals, _cr->shape_dof()};
      auto net = hb::contact_net_gn::create(&_cr->table, x, _cr->normal, _cr->shape, c.w_coupling, net_blocks);
      net->tangency_target = c.contact_tangency;
      ctx.constraints.push_back(net);
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

  /// `shape_distortion` as `edge_shape_fair` rows on the edge store.
  void push_edge_fairing(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    const std::vector<vec3> &x = _shell->xs->get();
    hepworth::sim_block::ptr store = _cr->shape_dof();
    const std::vector<hepworth::sim_block::ptr> held{store};
    auto fair = [&](hb::edge_fair_rows rows, real w, hb::edge_goal_fn goal, real lambda = real(1),
                    const std::vector<vec3> *frames = nullptr, bool tensor_metric = false) {
      if (w > 0.0 && std::isfinite(w))
        ctx.constraints.push_back(hb::edge_shape_fair::create(&_cr->table, x, rows, w, held, std::move(goal),
                                                              lambda, frames, tensor_metric));
    };
    const shape_distortion &s = _cfg.shape;
    // The store is the live reading at prepare: the rest the knobs measure against.
    const std::vector<vec3> &h0 = _cr->shape;
    fair(hb::edge_fair_rows::pin, s.w_shape_arap,
         hb::goal_shape_arap(store, &_cr->table, h0, s.shape_arap_mean));
    if (s.dupin_reach > 0.0) {
      if (s.w_dupin > 0.0 && std::isfinite(s.w_dupin))
        ctx.constraints.push_back(hb::edge_shape_fair::create(
            &_cr->table, x, hb::edge_fair_rows::dupin_line, s.w_dupin, held, hb::goal_zero(), real(1), &h0, false,
            real(0), hb::dupin_line_spec{s.dupin_reach, s.dupin_degree}));
    } else {
      fair(hb::edge_fair_rows::dupin, s.w_dupin, hb::goal_zero(), real(1), &h0);
    }
    fair(hb::edge_fair_rows::dirichlet_traceless, s.w_grad_traceless, hb::goal_zero());
    for (hb::curvature_growth_target spec : s.growth) {
      spec.angle = winding_angle(spec.winding, spec.angle);
      fair(hb::edge_fair_rows::pin, s.w_growth, hb::goal_curvature_growth(spec, h0), real(1), nullptr, true);
    }
    hb::init_ring_det(_cr->table, x, h0, s.w_det, store, ctx.constraints);
    if (s.w_minimizer > 0.0)
      fair(hb::edge_fair_rows::pin, s.w_minimizer,
           hb::goal_principal_mean(store, hb::build_edge_frames(_cr->table, x)));
    fair(hb::edge_fair_rows::pin, s.w_tension, hb::goal_mean_tension(store, s.tension));
  }

  /// One solver pass. Does not remesh, so a meridional seam stays a seam.
  void relax() {
    _solver_config.dt = _cfg.step.h;
    _solver_config.damping = _cfg.step.damping;
    _solver_config.iterations = _cfg.step.iterations;
    _solver.gn_refresh_every = _cfg.step.gn_refresh_every;
    _solver.gn_backtracks = _cfg.step.gn_backtracks;
    _solver.gn_ascent_tol = _cfg.step.gn_ascent_tol;
    _solver.gn_stop_on_fresh_ascent = _cfg.step.gn_stop_on_fresh_ascent;
    _solver.gn_max_ascent_refreshes = _cfg.step.gn_max_ascent_refreshes;
    _solver.run_variational_audit = _cfg.step.audit_variational;
    _solver.preconditioner = _cfg.step.preconditioner;
    _solver.modal = _cfg.step.modal;
    _solver.accel = _cfg.step.accel;
    _solver.anderson_m = _cfg.step.anderson_m;
    _solver.cheb_delay = _cfg.step.cheb_delay;
    _solver.cheb_window = _cfg.step.cheb_window;
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
    //step N times to relax topology
    // I suspect that bad sheet breakage
    for (int i = 0; i < 2; ++i) {
      //ideally step should report whether or not it did anywork...
      __surf->step(true);
    }
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
      const std::array<int, 3> v{static_cast<int>(t[0]), static_cast<int>(t[1]), static_cast<int>(t[2])};
      if (v[0] < 0 || v[1] < 0 || v[2] < 0 || static_cast<size_t>(v[0]) >= x.size() ||
          static_cast<size_t>(v[1]) >= x.size() || static_cast<size_t>(v[2]) >= x.size())
        continue;
      all.push_back(
          {static_cast<int>(fi), v,
           kusama::weingarten::face_quality(x[static_cast<size_t>(v[0])], x[static_cast<size_t>(v[1])],
                                            x[static_cast<size_t>(v[2])])});
    }
    std::sort(all.begin(), all.end(), [](const tri &a, const tri &b) { return a.q < b.q; });
    int n_bad = 0;
    while (n_bad < static_cast<int>(all.size()) && !(all[static_cast<size_t>(n_bad)].q >= q_below))
      ++n_bad;
    const int n_out = force ? std::max(n_bad, std::min(n_worst, static_cast<int>(all.size()))) : n_bad;
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
      log << "  face " << t.f << " quality " << t.q << " area " << real(0.5) * (b - a).cross(c - a).norm()
          << " verts " << t.v[0] << " " << t.v[1] << " " << t.v[2] << " lengths " << (b - a).norm() << " "
          << (c - b).norm() << " " << (a - c).norm() << " x (" << a.transpose() << ") (" << b.transpose()
          << ") (" << c.transpose() << ")\n";
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
    std::cout << "[monge] frame " << _frame << " " << stage << ": " << n_bad << " triangle(s) below "
              << q_below << ", wrote dump/bad_triangles.log and " << obj << std::endl;
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
    const real nscale =
        std::max(asawa::shell::avg_length(*__M, px) * _cfg.draw.length * real(0.5), real(1e-3));
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
      const int v1 = static_cast<int>(__M->vert(__M->other(asawa::shell::corner_id(c))));
      if (v0 < 0 || v1 < 0 || v0 >= static_cast<int>(px.size()) || v1 >= static_cast<int>(px.size()))
        continue;
      const vec3 mid = real(0.5) * (px[static_cast<size_t>(v0)] + px[static_cast<size_t>(v1)]);
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
      const auto tang = std::dynamic_pointer_cast<hepworth::block::tangent_cache_compute>(compute);
      if (tang)
        tang->compute(q);
    }
    if (!_cr || !__M || _cr->topology.nE <= 0 || static_cast<int>(_cr->u.size()) < 2 * _cr->topology.nE)
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    const int n = __M->vert_count();
    std::vector<vec3> H(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> tang(static_cast<size_t>(n), vec3::Zero());
    std::vector<vec3> axes(static_cast<size_t>(n), vec3::Zero());
    std::vector<char> wrote(static_cast<size_t>(n), 0);
    const std::vector<vec3> *gpar = _cr->g_par.empty() ? nullptr : &_cr->g_par;
    const std::vector<vec3> *gperp = _cr->g_perp.empty() ? nullptr : &_cr->g_perp;
    hepworth::block::vertex_weingarten(*__M, x, _cr->u, H, tang, axes, wrote, &_cr->vertex_frames, gpar,
                                       gperp);
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
      seam.push_back(seam_vert{i, std::atan2(p.z(), rho - R)});
    }
    std::sort(seam.begin(), seam.end(), [](const seam_vert &a, const seam_vert &b) { return a.phi < b.phi; });
    out << "tube seam vertices=" << seam.size() << "  w_z is spin about N, tilt is hypot(w_x, w_y)\n";
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
      const quat qi = hepworth::block::frame_quat_from_axes(Ti, Bi, Ni, quat::Identity());
      const quat qj = hepworth::block::frame_quat_from_axes(Tj, Bj, Nj, qi);
      const hepworth::block::frame_body_omega w = hepworth::block::frame_body_omega_of(qi, qj);
      if (!w.ok)
        continue;
      const real tilt = std::hypot(w.x, w.y);
      out << "  phi " << seam[static_cast<size_t>(k)].phi << " -> "
          << seam[static_cast<size_t>((k + 1) % m)].phi << "  w_z " << w.z << "  tilt " << tilt << "\n";
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
      out << "  mean |w_z|=" << spin_abs / real(edges) << "  mean tilt=" << tilt_abs / real(edges)
          << "  sign flips=" << flips;
    }
    out << "\n";
  }

  const hepworth::block::shell_cr_block::ptr &cr() const { return _cr; }
  /// Convergence of the last solver step (iterates, energies, last step).
  const hepworth::block::projection_solver::step_report &last_solve() const { return _solver.last_step; }
  /// Mean angle (radians, sign free) between the first carried stretch axis
  /// and the one its target's `angle` reads from the solved charts,
  /// after the last solve: how far the principal frame turned away from the
  /// material. NaN without a carried direction.
  real material_deflection() const { return _material_deflection; }

  /// Mean distance of the curvature-line groups' points to their fitted
  /// circles, over the mean edge length (`dupin_circle_residual`): 0 when every
  /// curvature line is a circle (a Dupin cyclide). NaN without charts.
  real kite_u_range() const { return _kite_u_range; }

  real dupin_circle_residual() const {
    if (!_cr || !_shell)
      return std::numeric_limits<real>::quiet_NaN();
    const std::vector<vec3> &x = _shell->xs->get();
    if (_cr->table.n_verts != static_cast<int>(x.size()) || _cr->vertex_frames.size() != x.size())
      return std::numeric_limits<real>::quiet_NaN();
    return hepworth::block::dupin_circle_residual(_cr->table, x, _cr->vertex_frames);
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
  static void perturb_vertices(asawa::shell::shell &M, std::vector<vec3> &x, real amp_edges) {
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
      const real wave = std::sin(real(3) * p.x()) * std::cos(real(2) * p.y()) * std::sin(real(5) * p.z());
      off[static_cast<size_t>(i)] = amp * nrm;
    }
    for (int i = 0; i < n; ++i)
      x[static_cast<size_t>(i)] += off[static_cast<size_t>(i)];
  }

  /// Independent normal offset. Normals are taken before any vertex moves.
  static void perturb_vertices_noise(asawa::shell::shell &M, std::vector<vec3> &x, real amp_edges) {
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
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
  }

  static asawa::shell::shell::ptr mesh_from_tables(std::vector<vec3> vertices,
                                                   std::vector<std::vector<int>> faces) {
    std::vector<index_t> corners_next, corners_vert, corners_face;
    asawa::assemble_table(vertices, faces, corners_next, corners_vert, corners_face);
    asawa::shell::shell::ptr M = asawa::shell::shell::create(corners_next, corners_vert, corners_face);
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

  static asawa::shell::shell::ptr make_mesh(const std::string &path, int ring = 48, int tube = 24) {
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
  bool vertex_dirichlet_frame(int i, const std::vector<vec3> &x, vec3 &T, vec3 &B, vec3 &N,
                              const vec3 *lock) const {
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
    __M->const_for_each_vertex(vi, [&](asawa::shell::CornerId ci, const asawa::shell::shell &S) {
      const int s = static_cast<int>(ci) / 2;
      if (s < 0 || s >= static_cast<int>(slot.size()))
        return;
      const int e = slot[static_cast<size_t>(s)];
      if (e < 0 || e >= nE)
        return;
      vec3 W = kusama::edge_guidance_vector_3d(S, ci, x, N, coeff[static_cast<size_t>(e)],
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
        std::cout << "[monge] " << what << " winding (" << w.ring << ", " << w.tube
                  << "): " << (_torus_fit ? "torus R " : "no torus fit, R ") << _torus_R << " r " << _torus_r
                  << ", angle " << angle << std::endl;
    };
    for (const hepworth::block::stretch_target &s : stretch_now())
      log("stretch", s.winding, s.angle);
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
    const real nu = real(2) * real(M_PI) / (_cfg.position.stripes * asawa::shell::avg_length(*__M, x));
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
    const std::vector<vec3> d = kusama::stripes::level_directions(*__M, x, S, psi, &_stripe_report.skipped);
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
          kusama::material_direction::transport(*__M, _material_x0, x, asawa::get_vec_data(*__M, slot));
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
      const std::vector<vec3> c = hepworth::block::chart_directions(_cr->vertex_frames, targets[k].angle);
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
  /// Vertex datum holding the torus row's persistent one-ring medial depth.
  index_t _torus_depth_slot = -1;
  /// Vertex datum per stretch target with `material` (-1 without), its
  /// seeded flag, and the step-start positions while a solve is open.
  std::vector<index_t> _material_slot;
  /// Vertex datum per `stripes` stretch target holding its stripe pattern
  /// (-1 without), and its closed axes for the step.
  std::vector<index_t> _stripe_slot;
  std::vector<std::vector<vec3>> _stripe_dirs;
  stripe_report _stripe_report;
  std::vector<char> _material_seeded;
  std::vector<vec3> _material_x0;
  real _material_deflection = std::numeric_limits<real>::quiet_NaN();
  /// Range of the last frame's kite log-scale `u` (0 when off).
  real _kite_u_range = 0;
  bool _torus_fit = false;
  real _torus_R = 0, _torus_r = 0;
  vec3 _center0 = vec3::Zero();
  int _frame = 0;
  hepworth::block::shell_position_block::ptr _shell;
  hepworth::block::shell_cr_block::ptr _cr;
  hepworth::block::block_solver_config<hepworth::block::shell_position_block, hepworth::block::shell_cr_block>
      _solver_config;
  hepworth::block::projection_solver _solver;
};

} // namespace duchamp
} // namespace gaudi
