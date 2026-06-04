#ifndef GAUDI_KUSAMA_DIFFUSION_HPP
#define GAUDI_KUSAMA_DIFFUSION_HPP

/// Curated diffusion timestep entry points (second-order / coupled complex).

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/complex_laplacian.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/common.h"
#include <vector>

namespace gaudi {
namespace kusama {
namespace diffusion {

/// Optional anisotropic stiffness \p C (same layout as `kusama::build_lap`).
struct anisotropic_settings {
  const kusama::laplacian::sparmat *C = nullptr;
};

/// Implicit second-order (Crank–Nicolson / midpoint) scalar diffusion on vertex range.
inline std::vector<real> second_order(asawa::shell::shell::ptr M,
                                        const std::vector<vec3> &x,
                                        std::vector<real> f, real dt,
                                        const anisotropic_settings &aniso = {}) {
  kusama::laplacian L(M, x);
  if (aniso.C && static_cast<index_t>(aniso.C->rows()) == M->vert_count() &&
      static_cast<index_t>(aniso.C->cols()) == M->vert_count()) {
    L.set_stiffness(*aniso.C);
  }
  std::vector<real> comp = asawa::shell::compress_to_vert_range<real>(*M, f);
  comp = L.diffuse2(comp, dt);
  return asawa::shell::expand_from_vert_range<real>(*M, comp);
}

inline std::vector<real> second_order(asawa::shell::shell::ptr M,
                                        const std::vector<vec3> &x,
                                        std::vector<real> f, real dt) {
  return second_order(M, x, std::move(f), dt, {});
}

/// Same as @ref second_order with explicit anisotropic stiffness pointer.
inline std::vector<real> second_order_anisotropic(
    asawa::shell::shell::ptr M, const std::vector<vec3> &x, std::vector<real> f, real dt,
    const kusama::laplacian::sparmat *C) {
  anisotropic_settings s;
  s.C = C;
  return second_order(M, x, std::move(f), dt, s);
}

/// Coupled complex linear diffusion / dispersive Laplacian step on `(u, v)`.
struct complex_settings {
  std::vector<real> alpha_verts;
  rx::cgle::linear_config linear{};
  const kusama::laplacian::sparmat *dispersive_C = nullptr;
};

inline void complex(asawa::shell::shell::ptr M, const std::vector<vec3> &x,
                    std::vector<real> &u, std::vector<real> &v, real dt,
                    const complex_settings &cfg) {
  kusama::laplacian L(M, x);
  if (cfg.dispersive_C != nullptr &&
      static_cast<index_t>(cfg.dispersive_C->rows()) == M->vert_count() &&
      static_cast<index_t>(cfg.dispersive_C->cols()) == M->vert_count()) {
    L.set_stiffness(*cfg.dispersive_C);
  }
  rx::cgle::linear_crank_nicolson(L, cfg.alpha_verts, dt, u, v, cfg.linear);
}

} // namespace diffusion
} // namespace kusama
} // namespace gaudi

#endif
