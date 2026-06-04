#ifndef __GAUDI_DUCHAMP_RX_DIFFUSE__
#define __GAUDI_DUCHAMP_RX_DIFFUSE__

#include "gaudi/asawa/shell/shell.hpp"
#include "gaudi/kusama/laplacian.hpp"
#include "gaudi/common.h"
#include <vector>

namespace gaudi {
namespace duchamp {

// Future: a single `Diffuser`/`diffuse` abstraction could unify isotropic implicit CN here
// with the anisotropic `C` path and the CGLE block solve — different matrix structures today.

/// Implicit Crank--Nicolson on vertex scalars (same as `kusama::laplacian::diffuse2`).
/// Optional anisotropic cotan stiffness (same `nv`×`nv` layout as `kusama::build_lap`).
inline void rx_diffuse_scalar_implicit_cotan(
    asawa::shell::shell::ptr M, const std::vector<vec3> &x, std::vector<real> &f, real h,
    const kusama::laplacian::sparmat *C) {
  kusama::laplacian L(M, x);
  if (C && static_cast<index_t>(C->rows()) == M->vert_count() &&
      static_cast<index_t>(C->cols()) == M->vert_count()) {
    L.set_stiffness(*C);
  }
  std::vector<real> comp = asawa::shell::compress_to_vert_range<real>(*M, f);
  comp = L.diffuse2(comp, h);
  f = asawa::shell::expand_from_vert_range<real>(*M, comp);
}

/// Crude explicit Euler graph Laplace (no mass inverse): `u += h * (C u)` in compressed form.
/// Optional; for experiments only — stiffer CFL than implicit.
inline void rx_diffuse_scalar_explicit_cotan(
    asawa::shell::shell::ptr M, const std::vector<vec3> &x, std::vector<real> &f, real h,
    const kusama::laplacian::sparmat *C) {
  kusama::laplacian L(M, x);
  if (C && static_cast<index_t>(C->rows()) == M->vert_count() &&
      static_cast<index_t>(C->cols()) == M->vert_count()) {
    L.set_stiffness(*C);
  }
  std::vector<real> comp = asawa::shell::compress_to_vert_range<real>(*M, f);
  std::vector<real> Cu = L.multC(comp);
  for (size_t i = 0; i < comp.size() && i < Cu.size(); ++i)
    comp[i] += h * Cu[i];
  f = asawa::shell::expand_from_vert_range<real>(*M, comp);
}

} // namespace duchamp
} // namespace gaudi

#endif
