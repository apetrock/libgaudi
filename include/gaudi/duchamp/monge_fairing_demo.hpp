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
#include "gaudi/hepworth/block/chiral_drive.hpp"
#include "gaudi/hepworth/block/contact_net_gn.hpp"
#include "gaudi/hepworth/block/edge_shape_goals.hpp"
#include "gaudi/hepworth/block/flat_growth.hpp"
#include "gaudi/hepworth/block/frame_rows.hpp"
#include "gaudi/hepworth/block/generic_constraints.hpp"
#include "gaudi/hepworth/block/mobius_growth.hpp"
#include "gaudi/hepworth/block/osculating_torus_ring.hpp"
#include "gaudi/hepworth/block/shape_objectives.hpp"
#include "gaudi/hepworth/block/shape_projection.hpp"
#include "gaudi/hepworth/block/shell_constraints_init.hpp"
#include "gaudi/hepworth/block/solver_composition.hpp"
#include "gaudi/hepworth/block/weingarten_gn.hpp"
#include "gaudi/hepworth/blocks/shell_e3_block.hpp"
#include "gaudi/hepworth/blocks/shell_position_block.hpp"
#include "gaudi/hepworth/nodes/solver_builder.hpp"

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
  hepworth::block::modal_config modal{
      .modes = 4,             // smallest generalized eigenvalues corrected per step
      .ncv = 0,               // Lanczos subspace (0: max(2 modes + 1, 20))
      .lambda_max = 0.9,      // Newton step below this
      .lambda_neg = 1e-3,     // negative / vanishing curvature at or below this
      .tau = 0.05,            // negative-mode trust radius (relative P-norm)
      .fd_eps = 1e-6,         // central-difference step for H v
      .max_iterations = 200,  // Lanczos restarts
  };
};

/// The frame-start charts the position rows read (`shell_e3_block`), drawn
/// after `relax` and before remesh: red the tight principal direction, green
/// the loose one, blue the normal.
struct fairing_draw {
  bool frames = true;
  /// Scale red/green by raw `|κ|` (half-length `nscale * |κ|`); off keeps
  /// fixed `nscale` on both.
  bool curvature = false;
  /// Axis half-length in mean edge lengths (`nscale = length * mean_edge`).
  real length = 1.5;
};

/// The curvature solver.
///  off         — no curvature rows.
///  contact_net — the principal contact element net on vertex quaternion
///                frames (MATH.md 19.10 to 19.12, 19.21 E3): positions and one
///                quaternion `Y` per vertex solved together as Gauss-Newton
///                rows in the projection solver; the normal is
///                `Y e3 conj(Y) / |Y|^2` and the curvature is the frames' own
///                (no edge store). Each directed spoke's tangency reads the
///                frames' normal curvature; `shape_distortion` acts through
///                the native frame rows (`frame_rows`) on per-vertex tensors.
enum class curvature_solver { off, contact_net };

/// `w_coupling` is the contact net's row weight between positions and frames:
/// `shape_distortion` rows act on the frames' vertex tensors,
/// `position_distortion` rows on the positions, and the coupling keeps the
/// two in sync. Keep it at or above the position rows that hold curvature
/// (`w_bend`) for shape edits to show in the mesh.
struct curvature_config {
  /// `curvature_solver::off` or `curvature_solver::contact_net`.
  curvature_solver kind = curvature_solver::contact_net;
  real w_coupling = 0.5;
  /// Frames' mass per unit vertex dual area: a pull toward the frame-start
  /// frames that damps steps on a stale factor (0: none).
  real contact_normal_mass = 0.0;
};

/// Shape → position: the osculating-torus ring row and the curvature solver.
struct shape_project {
  /// Osculating-torus ring row. The frame-start charts set both curvatures
  /// and the principal angle at each foot. The 1-ring and 2-ring fit only the pose,
  /// the foot's normal offset and the normal's tilt. Every ring point lands
  /// on that posed torus. 0 is off.
  real w_osculating_torus = 0.0;
  /// Row scale (`two_ring`), pose iterations, and pose prior. `.volume`
  /// `corotational` keeps each one-ring fan prism's rest volume (to its
  /// tube-axis points) while landing; `off` lands only.
  hepworth::block::osculating_torus_ring_config osculating_torus{
      .two_ring = 1.0,          // 2-ring row scale against the 1-ring
      .iterations = 4,          // pose Gauss-Newton iterations per projection
      .prior = 1e-3,            // pull of the pose toward the chart pose
      .relative = false,        // target the frame-start fit plus the change since
      .store_gain = 1.0,        // scale on that change under `relative`
      .closest_point = false,   // land on the closest torus point (10x cost)
      .volume =
          {
              .mode = hepworth::block::ring_volume_mode::corotational,   // off, corotational
              .apex_cap = 10.0,                                          // farthest tube-axis point, in ring reaches
              .foot_weight = 1.0,                                        // the foot's absolute row
              .corotational =
                  {
                      .w_proj = 1.0,            // each point onto its landing
                      .w_rest = 0.0,            // to its corotated rest position
                      .w_stay = 0.0,            // to where it is now
                      .w_volume = 5.0e-1,       // each fan prism to its rest volume
                      .w_spine = 1.0,           // tube-axis points in the Procrustes fit
                      .iterations = 3,          // volume Gauss-Newton steps
                      .foot_stiffness = -1.0,   // negative: the ring's size
                      .persistent = false,      // rest volumes from the persistent depth
                  },
          },
  };
  /// Frames ↔ positions.
  curvature_config curvature{
      .kind = curvature_solver::contact_net,   // off, contact_net
      .w_coupling = 0.5,                        // positions <-> frames
      .contact_normal_mass = 0.0,               // frames' mass per unit dual area
  };
};

/// Curvature edits as native rows on the contact net's frames (`frame_rows`):
/// each vertex carries a curvature tensor in its own frame, tied to the frames
/// by fit rows. The mesh follows through `shape_project.curvature.w_coupling`.
struct shape_distortion {
  /// Shape ARAP (AMAP-style memory): holds each vertex's curvature at its
  /// frame-start value (how anisotropic it is, and the pattern of its lines)
  /// but lets each ring turn its curvature lines together. Edits pushed harder than this have
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

  /// The frame minimizer: each principal curvature differenced across its
  /// own lines between neighbours (across-line Dirichlet, frame-start
  /// principal frames). 0 leaves it off.
  real w_minimizer = 0.0;
  /// Surface tension in curvature space: each vertex's mean curvature pulled
  /// to `tension.h0` (`constant`, Young-Laplace: pressure over twice the
  /// tension; outward normals make a sphere's `H` negative), to the
  /// frame-start mean (`mean`, constant-mean-curvature fairing), or to that
  /// mean corrected toward the first frame's volume (`volume`,
  /// `.volume_gain`); trace-free part kept. 0 leaves it off.
  real w_tension = 0.0;
  hepworth::block::tension_target tension{
      .kind = hepworth::block::tension_target::mode::mean,   // constant, mean, volume
      .h0 = 0.0,                                             // `constant` target
      .volume_gain = 1.0,                                    // `volume`: flattening per relative deficit
  };
  /// Dupin rows: each principal curvature held constant along its own
  /// direction, so curvature lines become circles and the surface is locally
  /// a Dupin cyclide (torus, cylinder, cone and their Moebius images); per
  /// spoke, along the frame-start principal frames. 0 leaves it off.
  real w_dupin = 0.0;
  /// The Dupin rows' stencil: 0 differences each principal curvature across
  /// one spoke; above 0 traces each curvature line `dupin_reach` mean edge
  /// lengths both ways from every vertex (one mean edge length per step,
  /// about that many hops) and holds the residual of a degree-`dupin_degree`
  /// fit of its normal curvature along the line (0: constant, the line a
  /// circle; 1: linear change allowed). Traced at the frame start and frozen
  /// for the frame.
  real dupin_reach = 0.0;
  int dupin_degree = 0;
  /// Directional curvature growth, the shape-space match of position
  /// `stretch`: each frame the normal curvature along `.angle` (radians from
  /// the tight direction toward the other) times `1 + x`, across it times
  /// `1 + y`, against the frame start; negative bends less that way; the
  /// torsion along it held. `.winding` sets the angle at startup like
  /// stretch. Only the first entry is used (`.power` is not yet read).
  /// `w_growth` scales it; 0 leaves the row off.
  std::vector<hepworth::block::curvature_growth_target> growth{
      hepworth::block::curvature_growth_target{
          .angle = 1.0,                                                // radians from the tight direction
          .x = 1.0,                                                    // along: normal curvature times 1 + x
          .y = 0.0,                                                    // across: times 1 + y
          .power = 0.5,                                                // not yet read
          .winding = {.ring = 1, .tube = 5, .handedness = 1},    // {ring, tube} replaces angle on the torus
      },
  };
  real w_growth = 0.0;
  /// Holds each vertex ring's Gauss-image area (the solid angle of its
  /// faces' normals) at its frame-start value, the shape-space match of
  /// `area`. It does not resist a uniform tube shrink (the normals do not
  /// change). 0 leaves it off.
  real w_det = 0.0;
  /// The geodesic torsion along `twist_angle` (radians from the tight direction; winding
  /// like growth) raised by `tau0` per frame, a uniform twist. Pair it with a
  /// position growth (stretch, kite). 0 leaves it off.
  real w_twist = 0.0;
  real tau0 = 0.05;
  real twist_angle = 0.6;
  hepworth::block::stretch_winding twist_winding{.ring = 0, .tube = 0, .handedness = 1};
  /// Chiral curvature drive (Armon's chiral seed pods, MATH.md 19.22): each
  /// vertex's target tensor is its reference plus a saddle `k0 w` at
  /// `+-.opening` (45 degrees) about an axis, pinned through the frames' spoke
  /// readings. Held against a metric (`w_chiral_rest`) narrow regions twist
  /// (a helicoid, `H ~ 0`) and wide ones coil (a helix, `K ~ 0`); which one
  /// is set by `width sqrt(k0 / thickness)`, the solver chooses. `k0`'s sign
  /// is the handedness.
  ///   .mode     rate: reference the frame-start tensor (chirality added
  ///             every frame); held: the first frame's tensor, carried
  ///   .axis     as the shear's: principal (tight direction), frame (carried
  ///             tangent), automatic (principal where q >= .axis_q)
  ///   .angle    radians from the axis; `.winding` sets it on the torus
  ///   .measure, .power   where it acts, as the shear's (power 0: uniform)
  /// `w_chiral` scales the spoke rows; 0 leaves it off. Fields in order:
  /// `.k0`, `.angle`, `.opening`, `.winding`, `.measure`, `.power`, `.axis`,
  /// `.axis_q`, `.mode`.
  hepworth::block::chiral_target chiral{
      .k0 = 4.0,                                             // saddle per unit length; sign: handedness
      .angle = 0.0,                                          // radians from the axis
      .opening = 0.25 * M_PI,                                // saddle half-angle about the axis
      .winding = {.ring = 1, .tube = 5, .handedness = 1},    // {ring, tube} replaces angle on the torus
      .measure = hepworth::block::growth_measure::off_umbilic,   // where it acts
      .power = 0.0,                                          // 0: uniform
      .axis = hepworth::block::shear_axis::principal,        // principal, frame, automatic
      .axis_q = 0.2,                                         // automatic: below this q the frame tangent
      .mode = hepworth::block::chiral_mode::held,            // rate, held
      .keep_eigenvalues = false,                             // only turn the principal axes (no tighter tube)
      .k_max = std::numeric_limits<real>::infinity(),        // rate: cap on the change from the first frame
  };
  real w_chiral = 0.0;
  /// How the chiral spoke rows measure a miss: `curvature` (scale-free) or
  /// `chord` (an unmet target is partly met by shortening spokes, so the
  /// drive also shrinks the mesh).
  hepworth::block::spoke_target_form chiral_form = hepworth::block::spoke_target_form::curvature;
  /// The metric the chiral drive pushes against: every edge held at its
  /// length on the first frame (rest positions carried through remeshing).
  /// Without it the saddle (which asks for a tighter tube everywhere) is
  /// delivered by shrinking: the torus thins until it vanishes. Only active
  /// with the chiral drive on (it would freeze every metric drive, the shear
  /// included). 0 leaves it off.
  real w_chiral_rest = 1.0;
};

/// Position-space rows: the classic bend and area, the ring axis ratio and
/// shear (axes from the frame charts), the two growths, the metric memory and
/// the centroid pin.
struct position_distortion {
  /// Locks each one ring's covariant anisotropy: spoke rms along the tight
  /// principal tangent over that along the other held at `ratio` (1 is
  /// round). Area, skew and rotation are free, so it composes with `area`.
  /// Comment a line to drop it. `w_axis_ratio` scales every entry; 0 leaves
  /// the row off.
  std::vector<hepworth::block::axis_ratio_target> axis_ratio{
      hepworth::block::axis_ratio_target{
          .ratio = 0.5,                       // tight over loose spoke rms (1 round)
          .area = {.lo = 1.0, .hi = 1.0},       // ring area band against the frame start
      },
  };
  real w_axis_ratio = 0.45;
  /// Shears each one ring every frame: stretch by `1 + shear` along `.angle`
  /// (radians from the axis; NaN is the 45-degree bisector) and `1 - shear`
  /// across, against the frame-start ring. Only that ratio is held; area,
  /// skew and rotation are free. `.winding = {.ring, .tube, .handedness}`
  /// sets the angle at startup from the fitted torus (on the demo torus the
  /// bisector is about {1, 3}).
  /// Where: each ring's shear is scaled in log space, ratio `r^w`, with
  /// `w = m^power` (the measures in [0, 1]) or `(m / mean m)^power`, `m` the
  /// `.measure` of the frame chart's tensor (`power` 0: uniform everywhere):
  ///   off_umbilic     1 - (1 - q) min(1, |H| / mean |H|), q = |z| / (|H| + |z|):
  ///                   0 at curved round points (left to tension), 1 on flats
  ///                   and anisotropic regions
  ///   flatness        1 - |k_lo| / |k_hi|: 1 where one curvature vanishes
  ///   mean_curvature  |H| (negative power favours flatter regions)
  ///   anisotropy      |z|
  ///   k_lo            |k_lo|
  /// The axis `.angle` turns from (`.axis`):
  ///   principal       the tight principal direction (noise at round points)
  ///   frame           the carried frame tangent (coherent everywhere, but
  ///                   not tied to the curvature lines)
  ///   automatic       principal where `q >= .axis_q`, frame below it
  /// Comment a line to drop it. `w_shear` scales every entry; 0 leaves the
  /// row off. Fields in order: `.shear`, `.angle`, `.winding`, `.area`,
  /// `.measure`, `.power`, `.axis`, `.axis_q`.
  std::vector<hepworth::block::shear_target> shear{
      hepworth::block::shear_target{
          .shear = 0.125,                                        // (1 + s) / (1 - s) per frame along the axis
          .angle = std::numeric_limits<real>::quiet_NaN(),       // radians from the axis; NaN the bisector
          .winding = {.ring = 1, .tube = 5, .handedness = 1},   // sets the angle at startup on the torus
          .area = {.lo = 1.0, .hi = 1.0},                        // ring area band against the frame start
          .measure = hepworth::block::growth_measure::flatness,   // where the shear acts (table above)
          .power = 0.0,                                          // 0: uniform; higher: sharper
          .axis = hepworth::block::shear_axis::principal,       // principal, frame, automatic
          .axis_q = 0.2,                                         // automatic: below this q the frame tangent
      },
  };
  real w_shear = 1.0;
  /// Cotangent one-ring. Holds the mean-curvature vector baked on first use,
  /// so it resists curvature edits. 0 leaves it off.
  real w_bend = 1.0e-1;
  /// Per-triangle rest area: weight `.w` (0 leaves it off) and `.band`, the
  /// range of its step-start area each face may take freely (`{1, 1}`
  /// conserves it; `{0.5, 1.5}` lets it drift 50% per step before pushing
  /// back).
  hepworth::block::area_config area{
      .w = 1.0e-1,                                // 0 leaves it off
      .band = {.lo = 0.9, .hi = 1.1},          // free range against the step start
  };
  /// Moebius growth (`mobius_growth.hpp`): a vertex log-scale `u` per frame,
  /// then per edge a length row at `e^{(u_i + u_j) / 2}` times its
  /// frame-start length. `.source` (the frame charts' tensors):
  ///   mean_curvature  u = alpha (-H - level): convex grows, concave shrinks
  ///   anisotropy      u = alpha (|z| - level)
  ///   flat            u = alpha (1 - |k_lo| / |k_hi| - level): grows around
  ///                   parabolic lines, not at round points
  /// `.balanced` removes the mean; `.max_u` caps the per-frame `|u|`.
  /// Fields in order: `.alpha`, `.max_u`, `.source`, `.level`, `.balanced`.
  hepworth::block::mobius_target mobius{
      .alpha = 0.1,
      .max_u = 0.1,
      .source = hepworth::block::mobius_source::anisotropy,
      .level = 0.0,
      .balanced = false,
  };
  real w_mobius = 0.0;
  /// Directional growth with a curvature policy (`flat_growth.hpp`): per
  /// edge, goal length `l0 (1 + g w ((1 - iso) (dhat0 . a)^2 + iso))`.
  /// Where: `w = flatness^power`, or `(m / mean m)^power` for the other
  /// measures, `m` the `.measure` of the frame charts' tensor (`power` 0
  /// grows everywhere alike):
  ///   flatness        1 - |k_lo| / |k_hi|: 1 where one curvature vanishes
  ///                   (power 8 on the torus: top 1, outer 0.09, inner 0.004)
  ///   mean_curvature  |H| (negative power favours flatter regions)
  ///   anisotropy      |z|
  ///   k_lo            |k_lo| (negative power favours a flat loose direction)
  /// Along what: `a` at `.angle` radians from the loose principal direction
  /// toward the tight one (0 the loose direction, pi / 2 the tight one) or
  /// from `.winding` (a closed line on the fitted torus, as the shear's).
  /// `.iso` mixes in direction-free growth (iso 1 with power 0: a plain
  /// uniform scale growth). With flatness and the loose direction the top
  /// and bottom of the torus tube lengthen along the hoop and ruffle.
  /// Fields in order: `.g`, `.measure`, `.power`, `.angle`, `.winding`,
  /// `.iso`.
  hepworth::block::flat_growth_target flat_growth{
      .g = 0.05,                                              // length added per frame at full weight
      .measure = hepworth::block::growth_measure::mean_curvature,   // flatness, mean_curvature, anisotropy, k_lo
      .power = 0.0,                                           // 0: uniform; high: only the extremes grow
      .angle = 0.0,                                           // 0 loose direction, pi / 2 tight direction
      .winding = {.ring = 0, .tube = 0, .handedness = 1},     // {ring, tube} replaces angle on the torus
      .iso = 0.0,                                             // 0 directional, 1 direction-free
  };
  real w_flat_growth = 0.0;
  /// Intrinsic vector Dirichlet memory of the metric (`dirichlet_memory_gn`,
  /// rest at the frame start; factor the mean edge length). Required by the contact net (tangency
  /// alone leaves slides and scale free).
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
  fairing_draw draw{};
  shape_project project{};
  shape_distortion shape{};
  position_distortion position{};
  /// The remesher's element QC (`asawa::shell::qc_config`, fields in order):
  /// `datum`, `q_warn` (repair triangles below this quality), `q_min` (hard
  /// failure), `val_lo` / `val_hi` (interior valence band), `rounds`,
  /// `enabled`. Example: `.qc = {.q_warn = 0.15, .val_hi = 7}`.
  asawa::shell::qc_config qc{
      .datum = 0,         // vertex vec3 datum the metric reads
      .q_warn = 0.1,      // repair faces below this quality
      .q_min = 1e-4,      // hard failures below this
      .val_lo = 4,        // interior valence band
      .val_hi = 8,
      .rounds = 3,        // repair rounds per remesh step
      .enabled = true,
  };
  /// Write bad triangles to `dump/bad_triangles.log` and an OBJ of their
  /// neighborhood: after each solve those below the solver's degenerate
  /// quality (the worst few when the solve blew up or did not converge), and
  /// after each remesh those the QC left below `qc.q_warn`.
  bool dump_bad_triangles = true;
};

inline const monge_fairing_demo_config k_monge_fairing_demo_default{};

/// The fairing demo on the E3 contact net: positions and one quaternion frame
/// per vertex, native curvature rows on the frames, and a few position rows.
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
    _volume_ref = std::numeric_limits<real>::quiet_NaN();
    _torus_fit = fit_torus(x, _torus_R, _torus_r);
    asawa::init_vert_datum<vec3>(*__M, vec3::Zero());
    // Persistent one-ring medial depth for the torus row's volume term,
    // carried through the remesher (NaN until first seen).
    _torus_depth_slot = asawa::init_vert_datum<real>(*__M, std::numeric_limits<real>::quiet_NaN());
    // The chiral drive's held reference tensor and rest positions (NaN until
    // first seen).
    _chiral_h_slot = asawa::init_vert_datum<vec3>(*__M, vec3::Constant(std::numeric_limits<real>::quiet_NaN()));
    _chiral_x_slot = asawa::init_vert_datum<vec3>(*__M, vec3::Constant(std::numeric_limits<real>::quiet_NaN()));
    auto xs = std::make_shared<shell_vert_positions>(__M, 0);
    auto vs = std::make_shared<shell_vert_velocities>(__M, 1);
    _shell = std::make_shared<hepworth::block::shell_position_block>(__M, xs, vs);
    _e3 = std::make_shared<hepworth::block::shell_e3_block>(__M);
    _e3->frame_mass_scale = _cfg.project.curvature.contact_normal_mass;
    real C = 3.0;
    const real l0 = asawa::shell::avg_length(*__M, xs->get());
    __surf = asawa::shell::dynamic::create(__M, C * l0, real(1.5) * C * l0, C * real(1.0) * l0);
    __surf->set_qc(_cfg.qc);
    __surf->set_qc_handler([this](asawa::shell::shell &, const asawa::shell::qc_report &rep) {
      if (!rep.bad_faces.empty())
        dump_bad_triangles("qc_left", false, 8, _cfg.qc.q_warn);
    });

    namespace hb = hepworth::block;
    _solver_config =
        hb::block_solver_builder<hb::shell_position_block, hb::shell_e3_block>::create()
            .with_blocks(_shell, _e3)
            .with_recompute([this](hb::solver_context &ctx) {
              if (!(_cfg.position.area.w > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hb::init_area(*__M, ctx.constraints, _shell->xs->get(), _cfg.position.area, hb::select_blocks<0>(ctx));
            })
            .with_recompute([this](hb::solver_context &ctx) {
              if (!(_cfg.position.w_center > 0.0) || !__M->verts_are_dense_packed())
                return;
              const std::vector<vec3> &x = _shell->xs->get();
              std::vector<index_t> ids(x.size());
              for (size_t i = 0; i < x.size(); ++i)
                ids[i] = static_cast<index_t>(i);
              ctx.constraints.push_back(hb::centroid_pin::create(ids, vertex_areas(x), _center0,
                                                                 _cfg.position.w_center, hb::select_blocks<0>(ctx)));
            })
            .with_recompute([this](hb::solver_context &ctx) {
              if (!(_cfg.position.w_bend > 0.0) || !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hb::init_bending(*__M, ctx.constraints, _shell->xs->get(), _cfg.position.w_bend,
                               hb::select_blocks<0>(ctx));
            })
            .with_recompute([this](hb::solver_context &ctx) { push_position_rows(ctx); })
            .with_recompute([this](hb::solver_context &ctx) { push_curvature(ctx); })
            .with_recompute([this](hb::solver_context &ctx) {
              if (!(_cfg.project.w_osculating_torus > 0.0) || !_e3 || ctx.blocks.empty() ||
                  !__M->verts_are_dense_packed())
                return;
              ctx.shell_mesh = __M.get();
              hb::init_osculating_torus_ring(*__M, ctx.constraints, _cfg.project.w_osculating_torus, {ctx.blocks[0]},
                                             &_e3->vertex_frames, nullptr, _cfg.project.osculating_torus,
                                             &asawa::get_real_data(*__M, _torus_depth_slot));
            })
            .dt(_cfg.step.h)
            .damping(_cfg.step.damping)
            .iterations(_cfg.step.iterations)
            .build();
  }

  /// One solver pass. Does not remesh.
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
    hepworth::block::run_solver_step(_solver_config, _solver);
    if (_e3)
      _e3->commit_gauge();
  }

  void step(int frame, bool remesh = true) {
    _frame = frame;
    relax();
    recenter_mesh();
    draw_solved_frames();
    // A solve that blew up, went non-finite or stopped short of converging
    // dumps its worst triangles even when none is below the threshold.
    const auto &r = _solver.last_step;
    dump_bad_triangles("after_solve", !std::isfinite(r.E_last) ||
                                          r.E_last > real(1e3) * std::max(r.E_first, real(1e-30)) ||
                                          r.dq_last > real(1e-3) * r.q_norm);
    if (!remesh)
      return;
    for (int i = 0; i < 2; ++i)
      __surf->step(true);
    if (!_cfg.qc.enabled)
      dump_bad_triangles("after_remesh");
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

  /// The frame-start charts (red tight, green loose, blue normal).
  void draw_solved_frames() {
    if (!_cfg.draw.frames || !_e3 || !__M || !__M->verts_are_dense_packed())
      return;
    const std::vector<vec3> &x = asawa::const_get_vec_data(*__M, 0);
    const int nV = __M->vert_count();
    if (static_cast<int>(x.size()) < nV)
      return;
    const std::vector<hepworth::block::vertex_chart> &frames = _e3->vertex_frames;
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

  const hepworth::block::shell_e3_block::ptr &e3() const { return _e3; }
  /// Convergence of the last solver step (iterates, energies, last step).
  const hepworth::block::projection_solver::step_report &last_solve() const { return _solver.last_step; }
  /// Range of the last frame's Moebius log-scale `u` (0 when off).
  real mobius_u_range() const { return _mobius_u_range; }

  asawa::shell::shell::ptr __M;
  asawa::shell::dynamic::ptr __surf;

private:
  /// The ring axis ratio and shear, and the two growths, on the frame charts.
  void push_position_rows(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    if (!_e3 || ctx.blocks.empty() || !__M->verts_are_dense_packed())
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    if (_e3->vertex_frames.size() != x.size() || _e3->table.n_verts != static_cast<int>(x.size()))
      return;
    ctx.shell_mesh = __M.get();
    const position_distortion &p = _cfg.position;
    const std::vector<char> live(x.size(), 1);
    hb::init_axis_ratio(*__M, ctx.constraints, x, live, {ctx.blocks[0]}, p.w_axis_ratio, p.axis_ratio,
                        &_e3->vertex_frames);
    hb::init_shear(*__M, ctx.constraints, x, live, {ctx.blocks[0]}, p.w_shear, shear_now(), &_e3->vertex_frames);
    _mobius_u_range = 0;
    if (p.w_mobius > 0.0) {
      const std::vector<real> u = hb::init_mobius_growth(_e3->table, x, _e3->vertex_frames, p.mobius, p.w_mobius,
                                                         ctx.constraints, hb::select_blocks<0>(ctx));
      if (!u.empty()) {
        const auto mm = std::minmax_element(u.begin(), u.end());
        _mobius_u_range = *mm.second - *mm.first;
      }
    }
    hb::flat_growth_target fg = p.flat_growth;
    // A winding gives the angle from the tight direction; the growth axis is
    // measured from the loose one.
    if (_torus_fit && (fg.winding.ring != 0 || fg.winding.tube != 0))
      fg.angle = real(0.5 * M_PI) - winding_angle(fg.winding, real(0));
    hb::init_flat_growth(_e3->table, x, _e3->vertex_frames, fg, p.w_flat_growth, ctx.constraints,
                         hb::select_blocks<0>(ctx));
    push_chiral_rest(ctx, x);
  }

  /// The chiral drive's metric: each edge at its first-frame length, read
  /// from rest positions carried as a vertex datum.
  void push_chiral_rest(hepworth::block::solver_context &ctx, const std::vector<vec3> &x) {
    namespace hb = hepworth::block;
    const real w = _cfg.shape.w_chiral_rest;
    if (!(w > 0.0) || !(_cfg.shape.w_chiral > 0.0) || _cfg.shape.chiral.k0 == 0.0 || _chiral_x_slot < 0)
      return;
    std::vector<vec3> &x0 = asawa::get_vec_data(*__M, _chiral_x_slot);
    if (x0.size() != x.size())
      return;
    for (size_t i = 0; i < x.size(); ++i)
      if (!x0[i].allFinite())
        x0[i] = x[i];
    for (int e = 0; e < _e3->table.n_edges; ++e) {
      const std::array<int, 2> v = hb::growth_detail::edge_ends(_e3->table, e);
      if (v[0] < 0)
        continue;
      const real l0 = (x0[static_cast<size_t>(v[1])] - x0[static_cast<size_t>(v[0])]).norm();
      if (!(l0 > real(0)))
        continue;
      auto row = hb::edge_strain::create(std::vector<index_t>{v[0], v[1]}, x, w, hb::select_blocks<0>(ctx));
      row->set_goal_length(l0);
      ctx.constraints.push_back(row);
    }
  }

  /// The contact net on positions and frames, the metric memory and the
  /// native frame rows.
  void push_curvature(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    const curvature_config &c = _cfg.project.curvature;
    if (c.kind != curvature_solver::contact_net || !(c.w_coupling > 0.0) || !_e3 || !_shell ||
        ctx.blocks.empty() || !_e3->frame_dof() || !__M->verts_are_dense_packed() || _e3->table.n_edges <= 0)
      return;
    const std::vector<vec3> &x = _shell->xs->get();
    if (static_cast<int>(x.size()) != _e3->table.n_verts || _e3->normal.size() != x.size())
      return;
    ctx.shell_mesh = __M.get();
    hepworth::sim_block::ptr pos = ctx.blocks[0];
    ctx.constraints.push_back(hb::contact_net_gn::create(&_e3->table, x, _e3->normal, std::vector<vec3>{},
                                                         c.w_coupling,
                                                         std::vector<hepworth::sim_block::ptr>{pos, _e3->frame_dof()}));
    if (_cfg.position.w_memory > 0.0)
      ctx.constraints.push_back(hb::dirichlet_memory_gn::create(&_e3->table, x, _cfg.position.w_memory,
                                                                std::vector<hepworth::sim_block::ptr>{pos}));
    push_frame_rows(ctx);
  }

  /// `shape_distortion` as native rows on positions and frames (`frame_rows`).
  void push_frame_rows(hepworth::block::solver_context &ctx) {
    namespace hb = hepworth::block;
    const std::vector<vec3> &x = _shell->xs->get();
    if (!_e3->frame_dof() || _e3->frame.size() != x.size())
      return;
    const shape_distortion &s = _cfg.shape;
    hb::frame_rows_config fc;
    fc.w_memory = s.w_shape_arap;
    fc.memory_mean = s.shape_arap_mean;
    fc.w_tension = s.w_tension;
    fc.tension = s.tension.kind == hb::tension_target::mode::constant ? hb::frame_tension::constant
                 : s.tension.kind == hb::tension_target::mode::volume ? hb::frame_tension::volume
                                                                      : hb::frame_tension::frame_mean;
    fc.h0 = s.tension.h0;
    fc.volume_gain = s.tension.volume_gain;
    if (!std::isfinite(_volume_ref))
      _volume_ref = enclosed_volume(x);
    fc.volume_ref = _volume_ref;
    fc.w_minimizer = s.w_minimizer;
    fc.w_dupin = s.w_dupin;
    fc.dupin_reach = s.dupin_reach;
    fc.dupin_degree = s.dupin_degree;
    fc.w_det = s.w_det;
    fc.w_gradtl = s.w_grad_traceless;
    fc.w_twist = s.w_twist;
    fc.tau0 = s.tau0;
    std::vector<real> growth_angle, twist_angle;
    if (!s.growth.empty()) {
      const hb::curvature_growth_target &g = s.growth.front();
      fc.w_growth = s.w_growth;
      fc.growth_x = g.x;
      fc.growth_y = g.y;
      growth_angle.assign(x.size(), winding_angle(g.winding, g.angle));
    }
    twist_angle.assign(x.size(), winding_angle(s.twist_winding, s.twist_angle));
    const bool chiral = s.w_chiral > 0 && s.chiral.k0 != 0 && _e3->vertex_frames.size() == x.size();
    if (chiral) {
      fc.w_target = s.w_chiral;
      fc.target_form = s.chiral_form;
    }
    const bool any = fc.w_memory > 0 || fc.w_tension > 0 || fc.w_minimizer > 0 || fc.w_dupin > 0 ||
                     fc.w_growth > 0 || fc.w_det > 0 || fc.w_gradtl > 0 || fc.w_twist > 0 || fc.w_target > 0;
    if (!any)
      return;
    auto rows = hb::frame_rows::create(&_e3->table, x, _e3->frame, fc,
                                       std::vector<hepworth::sim_block::ptr>{ctx.blocks[0], _e3->frame_dof()},
                                       growth_angle, twist_angle);
    if (chiral) {
      hb::chiral_target ct = s.chiral;
      ct.angle = winding_angle(ct.winding, ct.angle);
      const std::vector<vec3> *h_first = nullptr;
      if (_chiral_h_slot >= 0) {
        std::vector<vec3> &h0 = asawa::get_vec_data(*__M, _chiral_h_slot);
        if (h0.size() == x.size()) {
          for (size_t i = 0; i < x.size(); ++i)
            if (!h0[i].allFinite() && _e3->vertex_frames[i].ok)
              h0[i] = _e3->vertex_frames[i].h;
          h_first = &h0;
        }
      }
      const std::vector<vec3> *href = ct.mode == hb::chiral_mode::held ? h_first : nullptr;
      std::vector<real> kt, taut;
      hb::chiral_spoke_targets(kusama::contact_net::build_topology(_e3->table), x, _e3->vertex_frames, href, ct,
                               kt, taut, h_first);
      rows->set_spoke_targets(std::move(kt), std::move(taut));
    }
    ctx.constraints.push_back(rows);
  }

  real enclosed_volume(const std::vector<vec3> &x) const {
    real v = 0;
    for (const auto &id : _e3->table.ids)
      v += x[static_cast<size_t>(id[0])].dot(x[static_cast<size_t>(id[1])].cross(x[static_cast<size_t>(id[2])])) /
           real(6);
    return v;
  }

  /// The config's shear list with each `winding` entry's angle set from the
  /// startup torus.
  std::vector<hepworth::block::shear_target> shear_now() const {
    std::vector<hepworth::block::shear_target> out = _cfg.position.shear;
    for (hepworth::block::shear_target &s : out)
      s.angle = winding_angle(s.winding, s.angle);
    return out;
  }

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

  /// The angle from the tube circle (the tight tangent) of a line winding
  /// `w` on the startup torus, `handedness * atan2(ring * sqrt(R^2 - r^2),
  /// tube * r)`; `fallback` without a fit or windings.
  real winding_angle(const hepworth::block::stretch_winding &w, real fallback) const {
    if (!_torus_fit)
      return fallback;
    return hepworth::block::loxodrome_angle(w, _torus_R, _torus_r, fallback);
  }



  monge_fairing_demo_config _cfg;
  /// First frame's enclosed volume (the `volume` tension's reference).
  real _volume_ref = std::numeric_limits<real>::quiet_NaN();
  /// Vertex datum holding the torus row's persistent one-ring medial depth.
  index_t _torus_depth_slot = -1;
  /// Vertex datums: the chiral drive's held reference tensor (chart frame)
  /// and rest positions.
  index_t _chiral_h_slot = -1;
  index_t _chiral_x_slot = -1;
  /// Range of the last frame's Moebius log-scale `u` (0 when off).
  real _mobius_u_range = 0;
  bool _torus_fit = false;
  real _torus_R = 0, _torus_r = 0;
  vec3 _center0 = vec3::Zero();
  int _frame = 0;
  hepworth::block::shell_position_block::ptr _shell;
  hepworth::block::shell_e3_block::ptr _e3;
  hepworth::block::block_solver_config<hepworth::block::shell_position_block, hepworth::block::shell_e3_block>
      _solver_config;
  hepworth::block::projection_solver _solver;
};

} // namespace duchamp
} // namespace gaudi
