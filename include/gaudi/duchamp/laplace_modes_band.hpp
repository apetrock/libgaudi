#ifndef GAUDI_DUCHAMP_LAPLACE_MODES_BAND_HPP
#define GAUDI_DUCHAMP_LAPLACE_MODES_BAND_HPP

namespace gaudi {
namespace duchamp {

/// Which part of the **PSD** operator \f$A=-L_{\mathrm{sym}}+\varepsilon I\f$ to request.
enum class laplace_modes_band {
  low_frequency, ///< Smallest algebraic \f$\lambda\f$ on \f$A\f$ (global harmonics).
  largest_magnitude, ///< Largest \f$|\lambda|\f$ on \f$A\f$ (high-frequency tail).
  /// \f$k\f$ eigenvalues **closest** to \f$\sigma\f$ (shift–invert). Auto \f$\sigma\f$
  /// uses `mid_slider_t` (see spectral_projection_config / spectral_modes_demo). Explicit
  /// \f$\sigma \ge 0\f$ ignores the slider.
  shift_invert_middle,
};

} // namespace duchamp
} // namespace gaudi

#endif
