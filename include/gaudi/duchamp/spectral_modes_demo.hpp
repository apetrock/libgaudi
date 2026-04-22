#ifndef __GAUDI_SPECTRAL_MODES_DEMO__
#define __GAUDI_SPECTRAL_MODES_DEMO__

#include "gaudi/common.h"

#include "gaudi/asawa/shell/asset_loader.hpp"
#include "gaudi/asawa/shell/operations.hpp"
#include "gaudi/asawa/shell/shell.hpp"

#include "gaudi/bontecou/laplacian.hpp"
#include "gaudi/bontecou/laplace_spectrum.hpp"
#include "gaudi/duchamp/laplace_modes_band.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace duchamp {

using namespace asawa;

/// Map a scalar vertex field `phi` to a per-vertex `vec4` color by mixing `col_a` / `col_b`.
/// Uses the same min-max / z-score / 2-98% percentile auto-scaling as
/// `spectral_modes_demo::get_mesh_colors` so sparse / localized / high-frequency modes stay
/// visible. Output size matches `phi.size()`; empty input returns empty output.
inline std::vector<vec4> scalar_field_to_mesh_colors(const Eigen::VectorXd &phi,
                                                      const vec4 &col_a,
                                                      const vec4 &col_b) {
  const int n = static_cast<int>(phi.size());
  std::vector<vec4> colors;
  if (n <= 0)
    return colors;
  colors.resize(static_cast<size_t>(n), col_a);

  double mn = phi[0], mx = phi[0];
  for (int i = 1; i < n; ++i) {
    const double v = phi[i];
    mn = std::min(mn, v);
    mx = std::max(mx, v);
  }
  const double span = mx - mn;

  double mean = 0.0;
  for (int i = 0; i < n; ++i)
    mean += phi[i];
  mean /= static_cast<double>(std::max(1, n));

  double var = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = phi[i] - mean;
    var += d * d;
  }
  var /= static_cast<double>(std::max(1, n - 1));
  const double sigma = std::sqrt(std::max(0.0, var));

  enum class scale_kind { min_max, z_score, percentile };
  scale_kind kind = scale_kind::min_max;
  if (sigma > 1e-14) {
    if (span < 1e-8 || span < 4.0 * sigma)
      kind = scale_kind::z_score;
    else if (span > 8.0 * sigma)
      kind = scale_kind::percentile;
  }

  double lo = mn, hi = mx;
  if (kind == scale_kind::percentile && n >= 5) {
    std::vector<double> w;
    w.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
      w.push_back(phi[i]);
    std::sort(w.begin(), w.end());
    const size_t nu = static_cast<size_t>(n);
    const size_t il =
        std::min(nu - 1, static_cast<size_t>(std::floor(0.02 * static_cast<double>(nu - 1))));
    const size_t ih =
        std::min(nu - 1, static_cast<size_t>(std::ceil(0.98 * static_cast<double>(nu - 1))));
    if (ih > il && w[ih] > w[il]) {
      lo = w[il];
      hi = w[ih];
    } else {
      kind = scale_kind::min_max;
    }
  }

  const double pspan = hi - lo;

  for (int k = 0; k < n; ++k) {
    const double v = phi[k];
    double t = 0.5;
    if (kind == scale_kind::z_score && sigma > 1e-14)
      t = (v - mean) / (3.0 * sigma) + 0.5;
    else if (kind == scale_kind::percentile && pspan > 1e-20)
      t = (v - lo) / pspan;
    else if (span > 1e-20)
      t = (v - mn) / span;
    else if (sigma > 1e-14)
      t = (v - mean) / (3.0 * sigma) + 0.5;
    t = std::clamp(t, 0.0, 1.0);
    colors[static_cast<size_t>(k)] = (1.0 - t) * col_a + t * col_b;
  }
  return colors;
}

/// Loads a mesh, builds cotan stiffness, symmetrizes for Spectra (lower-tri
/// matvec), computes Laplacian eigenmodes, colors vertices like growth_study.
class spectral_modes_demo {
public:
  typedef std::shared_ptr<spectral_modes_demo> ptr;

  /// @param mid_slider_t  Used only for shift–invert auto \f$\sigma\f$ when
  /// `mid_band_shift < 0`: \f$-1\f$ = upper target only (legacy `--mid`); \f$[0,1]\f$ =
  /// \f$\sigma = (1-t)\sigma_{\mathrm{lo}} + t\,\sigma_{\mathrm{hi}}\f$ with
  /// \f$\sigma_{\mathrm{hi}}=0.35\|A\|_\infty\f$, \f$\sigma_{\mathrm{lo}}\ll\sigma_{\mathrm{hi}}\f$.
  /// Ignored when band is not `shift_invert_middle`.
  static ptr create(int num_modes = 32, double epsilon_shift = 1e-6,
                    laplace_modes_band band = laplace_modes_band::low_frequency,
                    double mid_band_shift = -1.0, double mid_slider_t = -2.0) {
    return std::make_shared<spectral_modes_demo>(num_modes, epsilon_shift, band,
                                                 mid_band_shift, mid_slider_t);
  }

  spectral_modes_demo(int num_modes = 32, double epsilon_shift = 1e-6,
                      laplace_modes_band band = laplace_modes_band::low_frequency,
                      double mid_band_shift = -1.0, double mid_slider_t = -2.0)
      : _epsilon(epsilon_shift), _requested_modes(num_modes), _band(band),
        _mid_band_shift(mid_band_shift), _mid_slider_t(mid_slider_t) {
    __M = asawa::shell::load_bunny();
    asawa::shell::triangulate(*__M);
    if (!__M->verts_are_dense_packed()) {
      std::cerr << "[spectral_modes_demo] Shell not dense-packed; calling pack(M) "
                   "before Laplacian / eigen analysis (permutes VERTEX datum in sync).\n";
      asawa::shell::pack(*__M);
    }
    if (!__M->verts_are_dense_packed()) {
      std::cerr << "[spectral_modes_demo] ERROR: pack(M) did not yield a dense-packed "
                   "shell — verts_are_dense_packed() still false (pack or predicate bug).\n";
      throw std::runtime_error(
          "spectral_modes_demo: pack failed to dense-pack mesh");
    }
    std::vector<vec3> &x = asawa::get_vec_data(*__M, 0);
    asawa::center(x);

    const int nv = static_cast<int>(x.size());

    bontecou::laplacian lap(__M, x);
    Eigen::SparseMatrix<double> Ls = bontecou::symmetrize_sparse(lap.stiffness());
    Eigen::SparseMatrix<double> Lpos = Ls;
    Lpos *= -1.0;
    Lpos.makeCompressed();
    _L_eig = bontecou::regularize_stiffness(Lpos, _epsilon);

    if (static_cast<int>(_L_eig.rows()) != nv) {
      std::cerr << "[spectral_modes_demo] WARNING: operator size " << _L_eig.rows()
                << " != vertex data size " << nv << "\n";
    }

    if (_band == laplace_modes_band::low_frequency)
      _spectrum_ok = bontecou::laplace_eigs_low_frequency(_L_eig, _requested_modes, 0,
                                                         _evals, _evecs);
    else if (_band == laplace_modes_band::largest_magnitude)
      _spectrum_ok = bontecou::laplace_eigs_largest_magnitude(
          _L_eig, _requested_modes, 0, _evals, _evecs);
    else {
      double sigma = _mid_band_shift;
      if (sigma < 0.0) {
        const double bound = bontecou::sparse_sym_max_row_sum_abs(_L_eig);
        if (_mid_slider_t >= 0.0 && _mid_slider_t <= 1.0) {
          const double t = std::clamp(_mid_slider_t, 0.0, 1.0);
          sigma = t * std::max(bound, 1e-12);
          std::cerr << "[spectral_modes_demo] mid-band auto shift sigma=" << sigma
                    << " (t=" << t << " * |A|_inf bound " << bound << ")\n";
        } else {
          // Deprecated legacy path: no slider, fall back to 0.35 of bound.
          sigma = 0.35 * std::max(bound, 1e-12);
          std::cerr << "[spectral_modes_demo] mid-band auto shift sigma=" << sigma
                    << " (legacy 0.35 * |A|_inf bound " << bound << ")\n";
        }
      } else {
        std::cerr << "[spectral_modes_demo] mid-band shift sigma=" << sigma
                  << " (explicit; mid slider ignored)\n";
      }
      _sigma_used = sigma;
      _spectrum_ok = bontecou::laplace_eigs_shift_invert_nearest(
          _L_eig, sigma, _requested_modes, 0, _evals, _evecs);
    }

    if (!_spectrum_ok) {
      std::cerr << "[spectral_modes_demo] Spectra eigen solve failed\n";
    } else {
      _num_modes = static_cast<int>(_evals.size());
      const char *band_s = (_band == laplace_modes_band::low_frequency)
                               ? "smallest algebraic lambda on (-L_sym+eps*I)"
                               : (_band == laplace_modes_band::largest_magnitude)
                                     ? "largest |lambda| on (-L_sym+eps*I)"
                                     : "shift-invert nearest to sigma on (-L_sym+eps*I)";
      std::cerr << "[spectral_modes_demo] computed " << _num_modes
                << " eigenmodes (" << band_s << "), matrix nnz=" << _L_eig.nonZeros()
                << "\n";
      if (static_cast<int>(_evecs.rows()) != nv) {
        std::cerr << "[spectral_modes_demo] WARNING: eigenvector length "
                  << _evecs.rows() << " != " << nv << " (colors misaligned)\n";
      }
      const int ncheck = std::min(3, _num_modes);
      for (int j = 0; j < ncheck; ++j)
        bontecou::log_laplace_eigen_stats(_L_eig, j, _evals[j], _evecs.col(j));
      if (_num_modes > 1 && _band == laplace_modes_band::low_frequency)
        _mode_idx = 1;
    }
  }

  int mode_count() const { return _spectrum_ok ? _num_modes : 1; }

  void set_mode_index(int i) {
    if (!_spectrum_ok) {
      _mode_idx = 0;
      return;
    }
    _mode_idx = std::clamp(i, 0, _num_modes - 1);
  }

  int mode_index() const { return _mode_idx; }

  /// Residual and phi variation for the active mode (call after changing mode).
  void log_current_mode_quality() const {
    if (!_spectrum_ok || _num_modes <= 0)
      return;
    bontecou::log_laplace_eigen_stats(_L_eig, _mode_idx, _evals[_mode_idx],
                                      _evecs.col(_mode_idx));
  }

  /// Per-vertex colors for the active mode: magenta-to-cyan mix with auto-scaling
  /// (min-max / z-score / 2-98% percentile). Delegates to `scalar_field_to_mesh_colors`.
  std::vector<vec4> get_mesh_colors() const {
    const int nv = __M->vert_count();
    std::vector<vec4> colors(static_cast<size_t>(nv), vec4(1.0, 0.0, 0.0, 1.0));

    if (!_spectrum_ok || _evecs.size() == 0)
      return colors;

    const int col = _mode_idx;
    const int nrows = static_cast<int>(_evecs.rows());
    const int n = std::min(nv, nrows);
    if (n <= 0)
      return colors;

    const vec4 col_a(1.0, 0.0, 1.0, 1.0);
    const vec4 col_b(0.0, 1.0, 1.0, 1.0);

    Eigen::VectorXd phi = _evecs.col(col).head(n);
    std::vector<vec4> mapped = scalar_field_to_mesh_colors(phi, col_a, col_b);
    for (int k = 0; k < static_cast<int>(mapped.size()) && k < nv; ++k)
      colors[static_cast<size_t>(k)] = mapped[static_cast<size_t>(k)];
    return colors;
  }

  shell::shell::ptr __M;

private:
  double _epsilon;
  double _mid_band_shift = -1.0;
  /// Shift–invert auto σ when `_mid_band_shift < 0`: see `create` / ctor docs.
  double _mid_slider_t = -2.0;
  double _sigma_used = 0.0;
  laplace_modes_band _band = laplace_modes_band::low_frequency;
  int _requested_modes = 32;
  int _num_modes = 0;
  int _mode_idx = 0;
  bool _spectrum_ok = false;
  /// \f$-L_{\mathrm{sym}}+\varepsilon I\f$ (PSD when \f$L\f$ is NS-definite); Spectra target.
  Eigen::SparseMatrix<double> _L_eig;
  Eigen::VectorXd _evals;
  Eigen::MatrixXd _evecs;
};

} // namespace duchamp
} // namespace gaudi

#endif
