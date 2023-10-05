
#ifndef __DUCHAMP_VORTEX_MODULE__
#define __DUCHAMP_VORTEX_MODULE__

#include "Eigen/src/Geometry/AngleAxis.h"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/calder/shell_integrators.hpp"
#include "gaudi/calder/tangent_point_integrators.hpp"

#include "gaudi/common.h"
#include "gaudi/logger.hpp"
#include "module_base_shell.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace gaudi {

// probably doesn't need to be a module, could go straight in test class...
/////////////////////////////////////////////////////////////

namespace duchamp {

class tangent_point : public module_base_shell {
public:
  DEFINE_CREATE_FUNC(tangent_point)
  tangent_point(asawa::shell::shell::ptr M) : module_base_shell(M) {
    const std::vector<vec3> &x = asawa::const_get_vec_data(*M, 0);
    _eps = asawa::shell::avg_length(*M, x);
  };

  virtual ~tangent_point(){};

  std::vector<vec3> get_frame_grad() {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    std::vector<real> w = asawa::shell::vertex_areas(M, x);
    std::vector<vec3> Gv =
        calder::tangent_point_gradient(M, x, w, Nv, 0.1 * _eps, 3.0);
    std::vector<mat3> Fv =
        calder::tangent_point_gradient_frame(M, x, w, Nv, 0.1 * _eps, 3.0);
    std::vector<vec3> Gvf(Fv.size(), vec3::Zero());
    for (int i = 0; i < Fv.size(); i++) {
      mat3 F = Fv[i];
      vec3 g = Gv[i];
      // Eigen decomp of F
      Eigen::SelfAdjointEigenSolver<mat3> es(F);
      vec3 d = es.eigenvalues();
      mat3 V = es.eigenvectors();
      // std::cout << d.transpose() << std::endl;
      vec3 v0 = va::sgn(g.dot(V.col(0))) * V.col(0);
      vec3 v1 = va::sgn(g.dot(V.col(1))) * V.col(1);
      vec3 v2 = va::sgn(g.dot(V.col(2))) * V.col(2);

      Eigen::AngleAxis<real> aa(0.1 * M_PI, v2);
      v0 = aa * v0;
      v1 = aa * v1;
      vec3 fsub = (0.6 * d[0] - 0.1 * d[2]) * v0      //
                  + (-0.75 * d[1] - 0.45 * d[2]) * v1 //
                  + d[2] * v2;
      Gvf[i] = 0.1 * g + fsub;
    }
    std::vector<vec3> Gf = asawa::shell::vert_to_face<vec3>(M, x, Gvf);
    std::vector<vec3> Gs = calder::mls_avg<vec3>(M, Gf, x, 4.0 * _eps, 2.0);

    return Gs;
  }

  std::vector<vec3> get_grad() {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    std::vector<real> w = asawa::shell::vertex_areas(M, x);
    std::vector<vec3> Gv =
        calder::tangent_point_gradient(M, x, w, Nv, 0.1 * _eps, 6.0);

    std::vector<vec3> Gf = asawa::shell::vert_to_face<vec3>(M, x, Gv);
    std::vector<vec3> Gs = calder::mls_avg<vec3>(M, Gf, x, 2.0 * _eps, 2.0);

    return Gs;
  }

  std::vector<vec3> get_smoothd_grad() {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    std::vector<vec3> Nv = asawa::shell::vertex_normals(M, x);
    std::vector<real> w = asawa::shell::vertex_areas(M, x);
    real p0 = 6.0;
    real p1 = 2.0;
    std::vector<vec3> Gv =
        calder::tangent_point_gradient(M, x, w, Nv, 0.1 * _eps, p0);

    std::vector<vec3> Gf = asawa::shell::vert_to_face<vec3>(M, x, Gv);
    std::vector<vec3> Gs = calder::smoothed_gradient(M, x, Gf, 4.0 * _eps, p1);

    std::vector<real> Kv =
        calder::tangent_point_energy(M, x, w, Nv, 0.1 * _eps, p0);
    std::vector<real> Kf = asawa::shell::vert_to_face<real>(M, x, Kv);
    std::vector<vec3> Ks = calder::gradient_scalar(M, x, Kf, 4.0 * _eps, p1);
    std::vector<vec3> G(Ks.size(), vec3::Zero());
    for (int i = 0; i < G.size(); i++) {
      // G[i] = Gs[i] - Ks[i];
      // logger::line(x[i], x[i] + 1e-8 * Gs[i], vec4(0.0, 1.0, 1.0, 1.0));
      // logger::line(x[i], x[i] - 1e-8 * Ks[i], vec4(1.0, 0.0, 1.0, 1.0));
      G[i] = Gs[i] - Ks[i];
    }
    return G;
  }

  virtual void step(real h) {
    asawa::shell::shell &M = *_M;
    std::vector<vec3> &x = asawa::get_vec_data(M, 0);
    // std::vector<vec3> G = get_grad();
    std::vector<vec3> G = get_frame_grad();
    // std::vector<vec3> G = get_smoothd_grad();

    real C = 0.0;
    for (int i = 0; i < G.size(); i++) {
      C += G[i].norm();
      C = std::max(C, G[i].norm());
    }
    // C /= real(G.size());
    //  real C = std::min()
    //  max = 16.0 * std::max(max, 50.0);
    C *= 0.005;
    for (int i = 0; i < G.size(); i++) {
      logger::line(x[i], x[i] + 0.25 * G[i] / C, vec4(1.0, 0.0, 0.0, 1.0));
      x[i] -= 1.0 * h * G[i] / C;
      // t += 0.001 * h;
    }
  }
  real t = 0.0;
  real _eps;
};

} // namespace duchamp
} // namespace gaudi

#endif