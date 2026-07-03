#ifndef __GAUDI_DUCHAMP_ROD_CONSTRAINTS_DYNAMICS__
#define __GAUDI_DUCHAMP_ROD_CONSTRAINTS_DYNAMICS__

#include <vector>

#include "gaudi/asawa/rod/dynamic.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/calder/rod_integrators.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"
#include "gaudi/duchamp/utils/sdf.hpp"

namespace gaudi {
namespace duchamp {

inline std::vector<vec3> compute_boundary_gradients(const asawa::rod::rod &rod,
                                                    const sdf_base &sdf) {
  const std::vector<real> dists = sdf.distance(rod.__x);
  const std::vector<vec3> gdists = sdf.grad_distance(rod.__x);
  std::vector<vec3> f(rod.__x.size(), vec3::Zero());
  for (int i = 0; i < static_cast<int>(rod.__x.size()); ++i) {
    if (dists[i] > 0.0) {
      f[i] = -dists[i] * gdists[i];
    }
  }
  return f;
}

inline std::vector<vec3>
compute_tangent_point_gradient(asawa::rod::rod &rod,
                               const asawa::rod::dynamic &dynamic) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> l = rod.l0();
  const std::vector<vec3> T = rod.N2c();
  auto forces = calder::tangent_point_gradient(rod, x, l, T, 1.0 * eps, 6.0);
  for (auto &f : forces)
    f *= -1.0;
  return forces;
}

// Self-induced vortex field: integrates -w*kappa*pow(sin^2,q)*(dp x T) over the rod.
// q=0 disables the angular (phasor) filter; q>0 penalizes alignment with the search ray dp.
inline std::vector<vec3> compute_vortex_force(asawa::rod::rod &rod,
                                              const asawa::rod::dynamic &dynamic,
                                              real p = 4.0, real q = 1.0) {
  const real eps = dynamic._Cc;
  const std::vector<vec3> &x = rod.x();
  const std::vector<real> phi = rod.l0();
  return calder::vortex_force(rod, x, phi, eps, p, q);
}

inline sdf_base::ptr select_rod_sdf(int frame, const sdf_base::ptr &sdf0,
                                    const sdf_base::ptr &sdf1) {
return sdf0;
/*
  if ((frame / 400) % 2 == 0) {
    return sdf0;
  }
  return sdf1;
  */
}

inline void grow_rod_rest_lengths(asawa::rod::rod &rod, real factor) {
  std::vector<real> &l0 = rod.l0();
  for (auto &li : l0) {
    li *= factor;
  }
}

} // namespace duchamp
} // namespace gaudi

#endif
