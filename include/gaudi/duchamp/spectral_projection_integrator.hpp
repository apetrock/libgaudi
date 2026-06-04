#ifndef GAUDI_DUCHAMP_SPECTRAL_PROJECTION_INTEGRATOR_HPP
#define GAUDI_DUCHAMP_SPECTRAL_PROJECTION_INTEGRATOR_HPP

/// Subspace tracking for **vertex-sized** symmetric sparse operators: inject how `A` is built
/// from `(shell, positions)`, optionally swap the partial eigensolver, then rebuild, partial
/// spectrum, project a vertex field, optional greedy/Hungarian mode alignment, truncation,
/// and optional normal motion driven by a per-vertex scalar.
///
/// **Timestep loop (cotan / custom vertex Laplacian path):**
/// 1. `remember_embedding_before_deformation()` — snapshot `evecs` / coeffs.
/// 2. Deform `x` (e.g. `step_vertices_along_normals_scaled_by_field` using `vertex_field()`).
/// 3. `rebuild_operator()` — your builder must yield `|V|×|V|` symmetric sparse `A`.
/// 4. `compute_modes(warm_start, false)` — optional Lanczos warm-start from the previous field.
/// 5. `project_embedding_after_new_spectrum()` — project the remembered field and align modes.
///
/// **Vertex Laplacian vs vector Dirichlet (Stein edge stiffness):** this integrator operates in
/// **dense-packed vertex index space** (`A.rows() == vert_count`). The Stein operator from
/// `build_vector_dirichlet_energy` is `(2|E|)×(2|E|)` on edge DOFs — use
/// `kusama::vector_dirichlet_partial_spectrum` for that spectrum; do not feed raw Stein `L`
/// here unless you first reduce/prolong to `|V|` (not provided).
///
/// **GL / manual validation:** `projects/spectral_modes_test` (keys **B**, **N**, **W**).

#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/shell/datum_x.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/laplace_spectrum.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/kusama/laplacian_anisotropic.hpp"
#include "gaudi/common.h"
#include "gaudi/duchamp/laplace_modes_band.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

namespace gaudi {
namespace duchamp {

using kusama::index_t;
using kusama::real;

/// Same mid-band semantics as `spectral_modes_demo`.
struct spectral_projection_config {
  int requested_modes = 32;
  double epsilon_shift = 1e-6;
  laplace_modes_band band = laplace_modes_band::low_frequency;
  /// Shift–invert: if `< 0`, auto σ from `mid_slider_t` and `‖A‖_∞` (see `resolve_sigma_mid_band`).
  double mid_band_shift = -1.0;
  double mid_slider_t = -2.0;
  index_t ncv = 0;
  /// If true, `project_and_align_to_previous` uses a global optimum assignment (Hungarian);
  /// default greedy is faster for large `k`.
  bool use_hungarian_alignment = false;
  /// If true, after projecting onto the new (truncated) basis we renormalize `_coeffs`
  /// so the embedded field has unit L2 norm. Without this, each projection is a
  /// contraction (`|coeffs|_2 ≤ 1` since truncating a complete basis leaks energy),
  /// so the displacement field `f = evecs * coeffs` decays geometrically across
  /// frames — you see a big first step, then tiny steps. Useful for visualization
  /// / mode-tracking; disable for a physically decaying integrator.
  bool renormalize_projected_coeffs = false;
  /// If true (and `band == shift_invert_middle`), the next shift-invert solve uses
  /// `σ = λ̃_tracked`, the first-order Rayleigh–Schrödinger prediction for the
  /// tracked eigenvalue on the new operator. Keeps the tracked mode from jumping
  /// into a different part of the spectrum when many near-degenerate modes exist.
  bool retarget_sigma_to_tracked_mode = false;
  /// If true, `advance_projection_step` halves `dt` until the predicted perturbation
  /// is small compared to the local spectral gap (`r_tracked < cfl_alpha`), so the
  /// step doesn't cross modes. Up to `max_dt_backtracks` halvings, then proceeds anyway.
  bool adaptive_dt_cfl = false;
  double cfl_alpha = 0.5;
  int max_dt_backtracks = 5;
};

/// First-order eigen-perturbation prediction for `A_new = A_old + ΔA`, using
/// `φ` and `λ` of `A_old`. All arrays are indexed by column in the tracked basis.
struct spectrum_prediction {
  Eigen::VectorXd lambda_tilde; ///< λ̃_j = λ_j + <φ_j, ΔA φ_j>
  Eigen::VectorXd gap;          ///< gap_j = min_{i≠j} |λ_i − λ_j| (original spectrum)
  Eigen::VectorXd robustness;   ///< r_j = max_{i≠j} |<φ_i, ΔA φ_j>| / |λ_i − λ_j|
  int tracked_index = -1;       ///< argmax_j |coeffs_j|
  double tracked_lambda_tilde = 0.0;
  double tracked_gap = 0.0;
  double tracked_robustness = 0.0;
  bool valid = false;
};

namespace detail {

/// Minimum-cost perfect matching on an `n×n` matrix (rows ↔ columns). Returns `perm[r]` = column
/// assigned to row `r`. O(n³); costs must be finite.
inline std::vector<int> hungarian_min_cost_square(const std::vector<std::vector<double>> &a) {
  const int n = static_cast<int>(a.size());
  if (n == 0)
    return {};
  const double INF = 1e200;
  std::vector<double> u(static_cast<size_t>(n + 1));
  std::vector<double> v(static_cast<size_t>(n + 1));
  std::vector<int> p(static_cast<size_t>(n + 1));
  std::vector<int> way(static_cast<size_t>(n + 1));
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(static_cast<size_t>(n + 1), INF);
    std::vector<char> used(static_cast<size_t>(n + 1), false);
    do {
      used[static_cast<size_t>(j0)] = true;
      const int i0 = p[static_cast<size_t>(j0)];
      double delta = INF;
      int j1 = 0;
      for (int j = 1; j <= n; ++j) {
        if (!used[static_cast<size_t>(j)]) {
          const double cur =
              a[static_cast<size_t>(i0 - 1)][static_cast<size_t>(j - 1)] -
              u[static_cast<size_t>(i0)] - v[static_cast<size_t>(j)];
          if (cur < minv[static_cast<size_t>(j)]) {
            minv[static_cast<size_t>(j)] = cur;
            way[static_cast<size_t>(j)] = j0;
          }
          if (minv[static_cast<size_t>(j)] < delta) {
            delta = minv[static_cast<size_t>(j)];
            j1 = j;
          }
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[static_cast<size_t>(j)]) {
          u[static_cast<size_t>(p[static_cast<size_t>(j)])] += delta;
          v[static_cast<size_t>(j)] -= delta;
        } else
          minv[static_cast<size_t>(j)] -= delta;
      }
      j0 = j1;
    } while (p[static_cast<size_t>(j0)] != 0);
    do {
      const int j1 = way[static_cast<size_t>(j0)];
      p[static_cast<size_t>(j0)] = p[static_cast<size_t>(j1)];
      j0 = j1;
    } while (j0 != 0);
  }
  std::vector<int> ans(static_cast<size_t>(n), -1);
  for (int j = 1; j <= n; ++j) {
    if (p[static_cast<size_t>(j)] != 0)
      ans[static_cast<size_t>(p[static_cast<size_t>(j)] - 1)] = j - 1;
  }
  return ans;
}

} // namespace detail

/// Injected builder: must return a **square** `|V|×|V|` sparse matrix for current `x` (datum 0).
using vertex_operator_builder_fn = std::function<Eigen::SparseMatrix<double>(
    asawa::shell::shell::ptr M, std::vector<vec3> &x, const spectral_projection_config &cfg)>;

/// Partial symmetric eigensolve on `A` (same contract as `laplace_eigs_*`). Sets `sigma_used` when
/// the mid band uses shift–invert; otherwise may leave it unchanged.
using partial_vertex_eigensolve_fn = std::function<bool(
    const Eigen::SparseMatrix<double> &A, const spectral_projection_config &cfg,
    double &sigma_used, Eigen::VectorXd &evals, Eigen::MatrixXd &evecs,
    const Eigen::VectorXd *warm_start)>;

/// Default cotan stiffness: `A = -sym(L) + εI` (same as `spectral_modes_demo`).
inline Eigen::SparseMatrix<double>
make_default_cotan_vertex_operator(asawa::shell::shell::ptr M, std::vector<vec3> &x,
                                   const spectral_projection_config &cfg) {
  if (!M)
    return {};
  const int nv = static_cast<int>(x.size());
  if (nv <= 0)
    return {};
  kusama::laplacian lap(M, x);
  Eigen::SparseMatrix<double> Ls = kusama::symmetrize_sparse(lap.stiffness());
  Eigen::SparseMatrix<double> Lpos = Ls;
  Lpos *= -1.0;
  Lpos.makeCompressed();
  return kusama::regularize_stiffness(Lpos, static_cast<real>(cfg.epsilon_shift));
}

/// Curvature-guided anisotropic cotan stiffness: `A = -sym(L_aniso) + εI` (same post-process
/// as \ref make_default_cotan_vertex_operator).
inline Eigen::SparseMatrix<double>
make_anisotropic_cotan_vertex_operator(asawa::shell::shell::ptr M, std::vector<vec3> &x,
                                       const spectral_projection_config &cfg,
                                       real vd_lambda, real sigma_u, real sigma_v,
                                       asawa::shell::face_curvature_stencil stencil,
                                       kusama::anisotropic_laplacian_kind kind) {
  if (!M)
    return {};
  const int nv = static_cast<int>(x.size());
  if (nv <= 0)
    return {};
  Eigen::SparseMatrix<real> C = kusama::build_curvature_aligned_laplacian(
      *M, x, vd_lambda, sigma_u, sigma_v, stencil, kind);
  Eigen::SparseMatrix<double> Ls = kusama::symmetrize_sparse(C);
  Eigen::SparseMatrix<double> Lpos = Ls;
  Lpos *= -1.0;
  Lpos.makeCompressed();
  return kusama::regularize_stiffness(Lpos, static_cast<real>(cfg.epsilon_shift));
}

/// Default: Spectra via `laplace_eigs_*` / `spectrum.hpp` (symmetric `A`).
inline bool default_vertex_eigensolve_for_projection(const Eigen::SparseMatrix<double> &A,
                                                      const spectral_projection_config &cfg,
                                                      double &sigma_used,
                                                      Eigen::VectorXd &evals,
                                                      Eigen::MatrixXd &evecs,
                                                      const Eigen::VectorXd *warm_start) {
  const index_t k_req = std::max(1, cfg.requested_modes);
  sigma_used = 0.0;
  if (cfg.band == laplace_modes_band::low_frequency) {
    return kusama::laplace_eigs_low_frequency(A, k_req, cfg.ncv, evals, evecs, warm_start);
  }
  if (cfg.band == laplace_modes_band::largest_magnitude) {
    return kusama::laplace_eigs_largest_magnitude(A, k_req, cfg.ncv, evals, evecs, warm_start);
  }
  double sigma = cfg.mid_band_shift;
  if (sigma < 0.0) {
    const double bound = kusama::sparse_sym_max_row_sum_abs(A);
    if (cfg.mid_slider_t >= 0.0 && cfg.mid_slider_t <= 1.0) {
      const double t = std::clamp(cfg.mid_slider_t, 0.0, 1.0);
      sigma = t * std::max(bound, 1e-12);
    } else {
      // Deprecated legacy path: no slider, fall back to 0.35 of bound.
      sigma = 0.35 * std::max(bound, 1e-12);
    }
  }
  sigma_used = sigma;
  return kusama::laplace_eigs_shift_invert_nearest(
      A, static_cast<real>(sigma), k_req, cfg.ncv, evals, evecs, warm_start);
}

/// Greedy overlap assignment: tracked column `j` best matches `prev(:,j)` ↔ raw column `raw_col_for_tracked[j]`.
struct mode_alignment {
  std::vector<index_t> raw_col_for_tracked;
  std::vector<real> sign;
};

inline mode_alignment compute_mode_alignment_greedy(const Eigen::MatrixXd &prev,
                                                    const Eigen::MatrixXd &raw_new) {
  const index_t k = static_cast<index_t>(prev.cols());
  mode_alignment a;
  a.raw_col_for_tracked.resize(static_cast<size_t>(k));
  a.sign.resize(static_cast<size_t>(k));
  std::vector<char> used(static_cast<size_t>(k), 0);
  for (index_t j = 0; j < k; ++j) {
    index_t best_i = -1;
    real best_dot = real(-1);
    for (index_t i = 0; i < k; ++i) {
      if (used[static_cast<size_t>(i)])
        continue;
      const real d = std::abs(prev.col(j).dot(raw_new.col(i)));
      if (d > best_dot) {
        best_dot = d;
        best_i = i;
      }
    }
    if (best_i < 0)
      best_i = j;
    used[static_cast<size_t>(best_i)] = 1;
    a.raw_col_for_tracked[static_cast<size_t>(j)] = best_i;
    const real dot = prev.col(j).dot(raw_new.col(best_i));
    a.sign[static_cast<size_t>(j)] = (dot >= real(0)) ? real(1) : real(-1);
  }
  return a;
}

/// Global optimum for total overlap: minimize `\sum_j ( -|prev(:,j)·raw(:,i_j)| )` over permutations `i_j`.
inline mode_alignment compute_mode_alignment_hungarian(const Eigen::MatrixXd &prev,
                                                       const Eigen::MatrixXd &raw_new) {
  const index_t k = static_cast<index_t>(prev.cols());
  mode_alignment a;
  a.raw_col_for_tracked.assign(static_cast<size_t>(k), 0);
  a.sign.assign(static_cast<size_t>(k), real(1));
  if (k == 0)
    return a;
  std::vector<std::vector<double>> cost(static_cast<size_t>(k),
                                         std::vector<double>(static_cast<size_t>(k)));
  for (index_t j = 0; j < k; ++j)
    for (index_t i = 0; i < k; ++i)
      cost[static_cast<size_t>(j)][static_cast<size_t>(i)] =
          -std::abs(prev.col(j).dot(raw_new.col(i)));
  const std::vector<int> perm = detail::hungarian_min_cost_square(cost);
  for (index_t j = 0; j < k; ++j) {
    const int pi = perm[static_cast<size_t>(j)];
    const index_t i = (pi >= 0) ? static_cast<index_t>(pi) : j;
    a.raw_col_for_tracked[static_cast<size_t>(j)] = i;
    const real dot = prev.col(j).dot(raw_new.col(i));
    a.sign[static_cast<size_t>(j)] = (dot >= real(0)) ? real(1) : real(-1);
  }
  return a;
}

/// Reorder/sign `evecs` columns and `evals`, map coefficient vector from **raw** Spectra order to aligned order.
inline void apply_mode_alignment(const mode_alignment &al, Eigen::MatrixXd &evecs,
                                 Eigen::VectorXd &evals, Eigen::VectorXd &coeffs_raw_to_aligned) {
  const index_t k = static_cast<index_t>(evecs.cols());
  const Eigen::MatrixXd tmp = evecs;
  const Eigen::VectorXd etmp = evals;
  const Eigen::VectorXd ctmp = coeffs_raw_to_aligned;
  for (index_t j = 0; j < k; ++j) {
    const index_t i = al.raw_col_for_tracked[static_cast<size_t>(j)];
    const real s = al.sign[static_cast<size_t>(j)];
    evecs.col(j) = s * tmp.col(i);
    evals[j] = etmp(i);
    coeffs_raw_to_aligned[j] = s * ctmp(i);
  }
}

/// Keep the `keep_count` coefficients with largest magnitude; zero the rest.
inline void truncate_coeffs_by_magnitude(Eigen::VectorXd &coeffs, int keep_count) {
  const int k = static_cast<int>(coeffs.size());
  if (keep_count <= 0 || k <= 0)
    return;
  keep_count = std::min(keep_count, k);
  std::vector<int> ix(static_cast<size_t>(k));
  for (int i = 0; i < k; ++i)
    ix[static_cast<size_t>(i)] = i;
  std::sort(ix.begin(), ix.end(), [&](int a, int b) {
    return std::abs(coeffs[a]) > std::abs(coeffs[b]);
  });
  std::vector<char> keep(static_cast<size_t>(k), 0);
  for (int t = 0; t < keep_count; ++t)
    keep[static_cast<size_t>(ix[static_cast<size_t>(t)])] = 1;
  for (int i = 0; i < k; ++i) {
    if (!keep[static_cast<size_t>(i)])
      coeffs[i] = real(0);
  }
}

/// Integrate vertex positions along unit normals scaled by `amplitude * scalar_per_vertex[i]`.
inline void step_vertices_along_normals_scaled_by_field(
    asawa::shell::shell &M, std::vector<vec3> &x,
    const std::vector<real> &scalar_per_vertex, real dt, real amplitude = real(1)) {
  const int nv = static_cast<int>(x.size());
  if (nv <= 0 || static_cast<int>(scalar_per_vertex.size()) != nv)
    return;
  for (int vi = 0; vi < nv; ++vi) {
    const vec3 n =
        asawa::shell::vert_normal(M, asawa::shell::vert_id(vi), x).normalized();
    x[static_cast<size_t>(vi)] +=
        static_cast<real>(dt) * amplitude * static_cast<real>(scalar_per_vertex[static_cast<size_t>(vi)]) * n;
  }
}

/// Vertex Laplacian spectrum pipeline with injectable operator and optional custom eigensolver.
class spectral_projection_integrator {
public:
  explicit spectral_projection_integrator(spectral_projection_config cfg = {},
                                        vertex_operator_builder_fn vertex_builder = {},
                                        partial_vertex_eigensolve_fn eigensolve = {})
      : _cfg(std::move(cfg)), _vertex_builder(std::move(vertex_builder)),
        _eigensolve(std::move(eigensolve)) {}

  const spectral_projection_config &config() const { return _cfg; }
  spectral_projection_config &config() { return _cfg; }

  void set_vertex_operator_builder(vertex_operator_builder_fn f) { _vertex_builder = std::move(f); }

  void set_partial_eigensolve(partial_vertex_eigensolve_fn f) { _eigensolve = std::move(f); }

  void set_mesh(asawa::shell::shell::ptr M) { _M = std::move(M); }

  asawa::shell::shell::ptr mesh() const { return _M; }

  /// Build `_A` from current vertex positions (datum 0). Returns false if mesh unset or
  /// `A` is not `|V|×|V|`. Ensures the shell is dense-packed first — `laplacian(M, x)`
  /// and the resulting stiffness matrix assume `vert_id == index_in_x`, so a non-packed
  /// shell (e.g. one that has seen edge collapses without a follow-up pack) would
  /// silently produce a misaligned operator and garbage modes.
  bool rebuild_operator() {
    if (!_M)
      return false;
    if (!_M->verts_are_dense_packed()) {
      std::cerr << "[spectral_projection_integrator] shell not dense-packed; "
                   "calling pack(M) before rebuild (syncs VERTEX datum).\n";
      asawa::shell::pack(*_M);
      if (!_M->verts_are_dense_packed()) {
        std::cerr << "[spectral_projection_integrator] ERROR: pack(M) did not "
                     "yield dense-packed shell — aborting rebuild_operator.\n";
        return false;
      }
    }
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    const int nv = static_cast<int>(x.size());
    if (nv <= 0)
      return false;
    if (_vertex_builder)
      _A = _vertex_builder(_M, x, _cfg);
    else
      _A = make_default_cotan_vertex_operator(_M, x, _cfg);
    if (static_cast<int>(_A.rows()) != nv || static_cast<int>(_A.cols()) != nv ||
        _A.rows() != _A.cols()) {
      std::cerr << "[spectral_projection_integrator] rebuild_operator: expected " << nv << "×" << nv
                << " sparse A, got " << _A.rows() << "×" << _A.cols() << "\n";
      return false;
    }
    return true;
  }

  /// Partial eigen-decomposition for the current `_A` and config band.
  /// When `reset_coefficients_to_default` is true (first solve), picks a default one-hot
  /// mode like `spectral_modes_demo`. Set false when you will call
  /// `project_embedding_after_new_spectrum()` immediately after.
  bool compute_modes(const Eigen::VectorXd *warm_start = nullptr,
                     bool reset_coefficients_to_default = true) {
    if (_A.rows() < 2)
      return false;
    const index_t k_req = std::max(1, _cfg.requested_modes);
    _spectrum_ok = false;
    if (_eigensolve)
      _spectrum_ok = _eigensolve(_A, _cfg, _sigma_used, _evals, _evecs, warm_start);
    else
      _spectrum_ok =
          default_vertex_eigensolve_for_projection(_A, _cfg, _sigma_used, _evals, _evecs, warm_start);
    if (!_spectrum_ok) {
      _coeffs.resize(0);
      return false;
    }
    const int k = static_cast<int>(_evals.size());
    _coeffs.resize(k);
    if (reset_coefficients_to_default) {
      _coeffs.setZero();
      if (k > 1 && _cfg.band == laplace_modes_band::low_frequency)
        _coeffs[1] = 1.0; // skip constant mode — match spectral_modes_demo
      else if (k > 0)
        _coeffs[0] = 1.0;
    } else
      _coeffs.setZero();
    return true;
  }

  double sigma_used() const { return _sigma_used; }

  const Eigen::SparseMatrix<double> &operator_matrix() const { return _A; }
  const Eigen::VectorXd &evals() const { return _evals; }
  const Eigen::MatrixXd &evecs() const { return _evecs; }
  Eigen::VectorXd &coeffs() { return _coeffs; }
  const Eigen::VectorXd &coeffs() const { return _coeffs; }

  int mode_count() const { return static_cast<int>(_evals.size()); }
  bool spectrum_ok() const { return _spectrum_ok; }

  void set_coeffs_one_hot(int mode_index) {
    const int k = mode_count();
    if (k <= 0)
      return;
    _coeffs.setZero();
    const int j = std::clamp(mode_index, 0, k - 1);
    _coeffs[j] = 1.0;
  }

  Eigen::VectorXd vertex_field() const { return _evecs * _coeffs; }

  void project_field_onto_current_basis(const Eigen::VectorXd &f_vertex) {
    if (!_spectrum_ok || _evecs.size() == 0)
      return;
    _coeffs = _evecs.transpose() * f_vertex;
  }

  /// After a new `_evecs` solve, set `coeffs` to the projection of the previous embedded field
  /// `prev_evecs * prev_coeffs`, then align columns to `prev_evecs` and reorder `coeffs` / `evals`.
  void project_and_align_to_previous(const Eigen::MatrixXd &prev_evecs,
                                     const Eigen::VectorXd &prev_coeffs) {
    if (!_spectrum_ok || prev_evecs.cols() != _evecs.cols() ||
        prev_evecs.rows() != _evecs.rows() ||
        prev_coeffs.size() != prev_evecs.cols())
      return;
    const Eigen::VectorXd f = prev_evecs * prev_coeffs;
    Eigen::VectorXd c = _evecs.transpose() * f;
    const mode_alignment al =
        _cfg.use_hungarian_alignment
            ? compute_mode_alignment_hungarian(prev_evecs, _evecs)
            : compute_mode_alignment_greedy(prev_evecs, _evecs);
    apply_mode_alignment(al, _evecs, _evals, c);
    _coeffs = c;
  }

  /// Snapshot current `_A` / `_evecs` / `_evals` / `_coeffs` (embedding in the old basis).
  /// Call **before** `rebuild_operator` when the mesh geometry changes. `_prev_A` and
  /// `_prev_evals` feed the Rayleigh–Schrödinger prediction machinery.
  void remember_embedding_before_deformation() {
    _prev_evecs = _evecs;
    _prev_evals = _evals;
    _prev_A = _A;
    _coeffs_snapshot = _coeffs;
    _have_embedding_snapshot = _spectrum_ok && _prev_evecs.size() > 0 &&
                             _coeffs_snapshot.size() == _prev_evecs.cols();
  }

  /// First-order Rayleigh–Schrödinger prediction of the spectrum of the current `_A`
  /// using the pre-deformation basis (`_prev_A`, `_prev_evecs`, `_prev_evals`).
  /// Call **after** `rebuild_operator` and **before** `compute_modes` to decide:
  ///   (a) whether `dt` was too big (high `robustness` ⇒ mode will cross), and
  ///   (b) where to retarget σ for shift-invert (`λ̃_tracked`).
  /// Cost: two sparse·dense matmuls (`A · Φ`) plus a `k×k` dense product. Much
  /// cheaper than the Lanczos/ARPACK solve.
  spectrum_prediction predict_after_rebuild() const {
    spectrum_prediction p;
    if (!_have_embedding_snapshot)
      return p;
    const int k = static_cast<int>(_prev_evecs.cols());
    const int n = static_cast<int>(_A.rows());
    if (k <= 0 || n <= 0)
      return p;
    if (_prev_A.rows() != n || _prev_A.cols() != n)
      return p;
    if (_prev_evals.size() != k || _prev_evecs.rows() != n ||
        _coeffs_snapshot.size() != k)
      return p;
    // M(i,j) = <φ_i, ΔA φ_j> computed without forming ΔA as a sparse matrix:
    //   (A_new Φ) − (A_old Φ) is a dense n×k, then Φ^T times that.
    const Eigen::MatrixXd dB =
        (_A * _prev_evecs) - (_prev_A * _prev_evecs);
    const Eigen::MatrixXd M = _prev_evecs.transpose() * dB;
    p.lambda_tilde.resize(k);
    p.gap.resize(k);
    p.robustness.resize(k);
    for (int j = 0; j < k; ++j) {
      p.lambda_tilde[j] = _prev_evals[j] + M(j, j);
      double gmin = std::numeric_limits<double>::infinity();
      double rmax = 0.0;
      for (int i = 0; i < k; ++i) {
        if (i == j)
          continue;
        const double d = std::abs(_prev_evals[i] - _prev_evals[j]);
        if (d < gmin)
          gmin = d;
        const double r = std::abs(M(i, j)) / std::max(d, 1e-30);
        if (r > rmax)
          rmax = r;
      }
      p.gap[j] = std::isfinite(gmin) ? gmin : 0.0;
      p.robustness[j] = rmax;
    }
    int tracked = 0;
    double best = -1.0;
    for (int j = 0; j < k; ++j) {
      const double m = std::abs(_coeffs_snapshot[j]);
      if (m > best) {
        best = m;
        tracked = j;
      }
    }
    p.tracked_index = tracked;
    p.tracked_lambda_tilde = p.lambda_tilde[tracked];
    p.tracked_gap = p.gap[tracked];
    p.tracked_robustness = p.robustness[tracked];
    p.valid = true;
    return p;
  }

  const spectrum_prediction &last_prediction() const { return _last_prediction; }

  /// After `rebuild_operator` + `compute_modes`, project the remembered field and align modes.
  void project_embedding_after_new_spectrum() {
    if (!_have_embedding_snapshot)
      return;
    project_and_align_to_previous(_prev_evecs, _coeffs_snapshot);
    if (_cfg.renormalize_projected_coeffs) {
      const double n = _coeffs.norm();
      if (n > 1e-12)
        _coeffs /= n;
    }
    _have_embedding_snapshot = false;
  }

  const Eigen::MatrixXd &previous_evecs() const { return _prev_evecs; }
  const Eigen::VectorXd &previous_evals() const { return _prev_evals; }
  const Eigen::SparseMatrix<double> &previous_operator_matrix() const { return _prev_A; }

  /// One projected step:
  ///   1. snapshot (A, φ, λ, coeffs)
  ///   2. displace along normals by `|f_i| * dt * amplitude`
  ///   3. rebuild operator
  ///   4. **predict** the new spectrum via Rayleigh–Schrödinger; if `adaptive_dt_cfl`
  ///      and `r_tracked > cfl_alpha`, revert positions, halve `dt`, retry
  ///   5. **retarget** σ to `λ̃_tracked` if `retarget_sigma_to_tracked_mode` (shift-invert)
  ///   6. solve, align to previous basis, optionally renormalize coeffs
  bool advance_projection_step(real dt, real amplitude,
                               bool use_vertex_field_as_warm_start = true) {
    if (!_M || !_spectrum_ok)
      return false;
    remember_embedding_before_deformation();
    std::vector<vec3> &x = asawa::get_vec_data(*_M, 0);
    const int nv = static_cast<int>(x.size());
    if (nv <= 0)
      return false;

    // Displacement direction is fixed once per step (from the pre-step basis);
    // we only scale it by the (possibly backtracked) dt when retrying.
    // Pass the SIGNED vertex field so negative regions of the mode pull inward
    // (along −n) and positive regions push outward (along +n). Using |f| would
    // collapse both lobes of the eigenmode into a one-sided outward flow.
    const Eigen::VectorXd f = vertex_field();
    std::vector<real> scalars(static_cast<size_t>(nv));
    for (int i = 0; i < nv; ++i)
      scalars[static_cast<size_t>(i)] = static_cast<real>(f[i]);

    // Snapshot positions so we can revert between backtracks.
    const std::vector<vec3> x_saved = x;

    real trial_dt = dt;
    const int max_bt = std::max(0, _cfg.max_dt_backtracks);
    int backtracks = 0;
    spectrum_prediction pred;
    while (true) {
      x = x_saved;
      step_vertices_along_normals_scaled_by_field(*_M, x, scalars, trial_dt,
                                                  amplitude);
      if (!rebuild_operator())
        return false;
      pred = predict_after_rebuild();
      // Accept unconditionally if CFL is off, prediction invalid, or the mode is
      // isolated in its band (gap==0 means all other modes equal — nothing to
      // cross into meaningfully, so backing off further won't help).
      const bool cfl_on = _cfg.adaptive_dt_cfl && pred.valid &&
                          pred.tracked_gap > 0.0;
      if (!cfl_on || pred.tracked_robustness < _cfg.cfl_alpha)
        break;
      if (backtracks >= max_bt) {
        std::cerr << "[spectral_projection] CFL: backtrack limit ("
                  << max_bt << ") reached with r_tracked="
                  << pred.tracked_robustness
                  << " gap_tracked=" << pred.tracked_gap
                  << "; proceeding with dt=" << trial_dt << "\n";
        break;
      }
      trial_dt *= real(0.5);
      ++backtracks;
    }
    if (backtracks > 0) {
      std::cerr << "[spectral_projection] CFL: backtracked " << backtracks
                << " time(s); dt " << dt << " -> " << trial_dt
                << " (r_tracked=" << pred.tracked_robustness
                << ", gap=" << pred.tracked_gap << ")\n";
    }
    _last_prediction = pred;

    // σ retargeting: force shift-invert to center on the predicted tracked
    // eigenvalue. Only meaningful for `shift_invert_middle`. We mutate the cfg
    // copy for one solve and restore.
    const bool retarget =
        _cfg.retarget_sigma_to_tracked_mode && pred.valid &&
        pred.tracked_index >= 0 &&
        _cfg.band == laplace_modes_band::shift_invert_middle;
    const double saved_mid_shift = _cfg.mid_band_shift;
    if (retarget)
      _cfg.mid_band_shift = pred.tracked_lambda_tilde;

    Eigen::VectorXd warm;
    const Eigen::VectorXd *ws = nullptr;
    if (use_vertex_field_as_warm_start && _prev_evecs.size() > 0 &&
        _coeffs_snapshot.size() == _prev_evecs.cols()) {
      warm = _prev_evecs * _coeffs_snapshot;
      ws = &warm;
    }
    const bool ok = compute_modes(ws, false);

    if (retarget)
      _cfg.mid_band_shift = saved_mid_shift;

    if (!ok)
      return false;
    project_embedding_after_new_spectrum();
    return true;
  }

private:
  spectral_projection_config _cfg;
  asawa::shell::shell::ptr _M;
  vertex_operator_builder_fn _vertex_builder;
  partial_vertex_eigensolve_fn _eigensolve;
  Eigen::SparseMatrix<double> _A;
  Eigen::VectorXd _evals;
  Eigen::MatrixXd _evecs;
  Eigen::VectorXd _coeffs;
  Eigen::MatrixXd _prev_evecs;
  Eigen::VectorXd _prev_evals;
  Eigen::SparseMatrix<double> _prev_A;
  Eigen::VectorXd _coeffs_snapshot;
  spectrum_prediction _last_prediction;
  bool _have_embedding_snapshot = false;
  bool _spectrum_ok = false;
  double _sigma_used = 0.0;
};

} // namespace duchamp
} // namespace gaudi

#endif
