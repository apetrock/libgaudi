#ifndef GAUDI_ALBERS_DARBOUX_MEDIAL_GEOMETRY_HPP
#define GAUDI_ALBERS_DARBOUX_MEDIAL_GEOMETRY_HPP

#include "gaudi/albers/darboux_cyclide.hpp"
#include "darboux_generated.hpp"
#include "medial_kernels.h"

namespace gaudi {
namespace albers {

struct darboux_geometry_bundle {
  vec14 Q;

  explicit darboux_geometry_bundle(vec14 q = vec14::Zero()) : Q(std::move(q)) {}

  vec3 G(const vec3 &x) const {
    return medial_generated::darboux_grad_generated(Q, x);
  }

  mat3 H(const vec3 &x) const {
    return medial_generated::darboux_hessian_generated(Q, x);
  }

  mat3 H_directional(const vec3 &x, const vec3 &fp) const {
    return medial_generated::darboux_H_directional_generated(Q, x, fp);
  }

  mat3 H_directional_rate(const vec3 &x, const vec3 &fp,
                          const vec3 &fpp) const {
    return medial_generated::darboux_H_directional_rate_generated(Q, x, fp,
                                                                  fpp);
  }
};

struct ray_line_bundle {
  vec3 foot = vec3::Zero();
  vec3 dir = vec3::UnitZ();

  ray_line_bundle() = default;
  ray_line_bundle(vec3 foot_in, vec3 dir_in)
      : foot(std::move(foot_in)), dir(std::move(dir_in)) {}

  vec3 f(real t) const { return foot + t * dir; }
  vec3 fp(real) const { return dir; }
  vec3 fpp(real) const { return vec3::Zero(); }
};

inline void eval_medial_energy_at_t(const darboux_geometry_bundle &geom,
                                    const ray_line_bundle &line, real t, real &E,
                                    real &E_prime, real &E_pp,
                                    real eps = 1e-12) {
  medial_generated::eval_medial_energy_at_t(geom, line, t, E, E_prime, E_pp,
                                            eps);
}

inline mat3 shape_operator_at(const vec14 &Q, const vec3 &x) {
  const vec3 g = medial_generated::darboux_grad_generated(Q, x);
  const mat3 H = medial_generated::darboux_hessian_generated(Q, x);
  mat3 W;
  medial_generated::shape_operator_from_GH(g, H, W);
  return W;
}

inline vec3 aligned_inward_ray_dir(const vec14 &Q, const vec3 &mesh_normal) {
  const vec3 g0 = medial_generated::darboux_grad_generated(Q, vec3::Zero());
  if (!g0.allFinite() || g0.norm() < 1e-12) {
    return vec3::Zero();
  }
  // March opposite the surface gradient (descent of D). With the convex kernel
  // g0 points outward (along +N), so -g0 points inward. Guarantee the returned
  // direction is inward (dir.N < 0) so the medial ridge sits at t > 0 and the
  // outward quartic artifact is excluded by the t >= 0 clamp in the search.
  vec3 dir = -g0.normalized();
  const vec3 n = mesh_normal.normalized();
  if (dir.dot(n) > 0.0) {
    dir *= -1.0;
  }
  return dir;
}

} // namespace albers
} // namespace gaudi

#endif // GAUDI_ALBERS_DARBOUX_MEDIAL_GEOMETRY_HPP
