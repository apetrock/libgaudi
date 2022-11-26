#ifndef __SURF_EXPERIMENT__
#define __SURF_EXPERIMENT__

#include "GaudiGraphics/geometry_logger.h"
#include "experiment.hpp"

namespace duchamp {

////////////////////////////////////////////////////////////////////////////
// potential
////////////////////////////////////////////////////////////////////////////

template <typename SPACE>
vector<typename SPACE::vec3>
cosineGradient(asawa::surf<SPACE> *mesh,
               const std::vector<typename SPACE::vec3> &vertVals,
               const std::vector<typename SPACE::coordinate_type> &evalPoints,
               typename SPACE::real regLength = 0.5) {
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  using Avg_Integrator =
      calder::Geometry_Integrator<SPACE, vec3, triangle_type, vec3>;

  using ATree = typename Avg_Integrator::Tree;
  using ANode = typename Avg_Integrator::Node;

  std::vector<vec3> faceVals = ci::verts_to_faces<SPACE>(vertVals, mesh);

  if (evalPoints.empty())
    return std::vector<vec3>();

  auto pre = [faceVals](const vector<triangle_type> &tris, ANode &node,
                        ATree &tree, vec3 &netCharge, coordinate_type &avgPoint,
                        coordinate_type &avgNormal) -> void {
    avgPoint = coordinate_type(0, 0, 0);
    netCharge = z::zero<vec3>();

    T netWeight = 0;

    for (int i = node.begin; i < node.begin + node.size; i++) {
      int ii = tree.permutation[i];
      triangle_type tri = tris[ii];
      T w = tri.area();
      avgPoint += w * tri.center();
      avgNormal += w * tri.normal();
      netCharge += w * faceVals[ii];
      netWeight += w;
    }

    avgPoint /= netWeight;
  };

  auto computeK = [](T dist, T C) {
#if 0
      T d3 = dist * dist * dist;
      T l3 = C * C * C;
      T kappa = (1.0 - exp(-d3 / l3)) / d3;
      return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
    T dist3 = dist * dist * dist;
    T l3 = C * C * C;
    T kappa = 1.0 / (dist3 + l3);
    return kappa / pow(4.0 * M_PI, 1.5);
#elif 1
    T dist3 = C / dist / dist / dist;
    T dist6 = dist * dist * dist * dist * dist * dist;

    T l2 = C * C;
    T kappa = 1.0 / (dist3 + dist6);
    return kappa / 4.0 / M_PI;
#elif 1
    T kappa = 1.0 / (dist + C);
    return kappa / 4.0 / M_PI;
#elif 0
    T dist2 = dist * dist;
    T dt = 0.5;
    T kappa = exp(-dist2 / 4.0 / dt);
    return kappa / pow(4.0 * M_PI * dt, 1.5);
#endif
  };
  coordinate_array vertex_normals = ci::get_vertex_normals<SPACE>(mesh);

  auto compute = [&vertex_normals, &vertVals, faceVals, regLength,
                  computeK](int i_c, const vec3 &wq, const coordinate_type &pc,
                            const coordinate_type &pe, const coordinate_type &N,
                            const vector<triangle_type> &tris, ANode &node,
                            ATree &tree) -> vec3 {
    vec3 out = z::zero<vec3>();

    coordinate_type Nv = vertex_normals[i_c];
    coordinate_type qc = vertVals[i_c];
    real qcmag = qc.norm();

    auto computeKelvin = [](coordinate_type r, coordinate_type f, real eps) {
      real mu = 10.0;
      real v = 2e-1;
      real a = 1.0 / 4.0 / mu;
      real b = a / (4.0 * (1.0 - v));
      real eps2 = 1e-10 * eps * eps;
      real re = 1.0 / sqrt(r.dot(r) + eps2);
      real re3 = re * re * re;

      mat3 R = r * r.transpose();
      mat3 I;
      I.setIdentity();

      coordinate_type u =
          ((a - b) / re * I + b / re3 * R + 0.5 * a * eps2 / re3 * I) * f;

      return u;
    };

    auto compute_grad0 = [](coordinate_type X, coordinate_type N,
                            real p) -> coordinate_type {
      real x = X.norm();
      real xp = pow(x, p);
      coordinate_type Xn = X / x;
      real Nx = N.dot(Xn);
      return (N - p * Nx * Xn) / xp;
    };

    auto compute_grad = [](coordinate_type X, coordinate_type N,
                           real p) -> coordinate_type {
      real x = X.norm();
      real xp = pow(x, p);
      coordinate_type Xn = X / x;
      real Nx = N.dot(Xn);
      real Nx1 = 1.0 + Nx;
      return 0.5 * N * Nx1 / xp - p * Nx1 * Nx1 / xp / x * Xn;
    };

    if (node.isLeaf()) {
      for (int i = node.begin; i < node.begin + node.size; i++) {
        int ii = tree.permutation[i];
        vec3 qi = faceVals[ii];
        auto tri = tris[ii];
        auto w = tri.area();
        auto c = tri.center();
        auto Nf = tri.normal();
        coordinate_type dp = c - pe;

        coordinate_type dkv = compute_grad(dp, Nv, 6.0);
        coordinate_type dkf = compute_grad(-dp, Nf, 6.0);
        out += w * (dkv - dkf);
      }
    } else {
      // out += computeK(dist, regLength) * va::cross(q, dp);
      real w = N.norm();
      coordinate_type Nf = N / w;
      coordinate_type dp = pc - pe;
      coordinate_type dkv = compute_grad(dp, Nv, 6.0);
      coordinate_type dkf = compute_grad(-dp, Nf, 6.0);
      out += w * (dkv - dkf);
    }
    return out;
  };

  vector<face_ptr> &faces = mesh->get_faces();
  vector<coordinate_type> normals;
  vector<triangle_type> triangles;
  for (int i = 0; i < faces.size(); i++) {
    if (!mesh->has_face(i))
      continue;
    if (faces[i]->size() < 3)
      continue;
    std::vector<triangle_type> tris = asawa::ci::get_tris<SPACE>(faces[i]);
    triangles.insert(triangles.end(), tris.begin(), tris.end());
  }

  for (auto t : triangles) {
    normals.push_back(t.normal());
  }
  vector<vec3> u(evalPoints.size(), z::zero<vec3>());
  Avg_Integrator integrator;
  integrator.integrate(faceVals, triangles, evalPoints, u, pre, compute);
  int i = 0;

  return u;
}

template <typename SPACE> class covariant_flow {
public:
  M2_TYPEDEFS;

  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;

  typedef std::shared_ptr<covariant_flow<SPACE>> ptr;
  static ptr create(real s) {
    return std::make_shared<covariant_flow<SPACE>>(s);
  }

  covariant_flow(real scale) : _scale(scale) {}

  coordinate_array step(surf_ptr surf, const coordinate_array &p0) {

    vertex_array &verts = surf->get_vertices();
    coordinate_array normals = asawa::ci::get_vertex_normals<SPACE>(surf);

    int i = 0;
    coordinate_array p1(p0);
    coordinate_array vg(p0.size(), coordinate_type::Zero());

    i = 0;
    calder::mesh_calculator<SPACE> calc;
    std::vector<typename SPACE::mat43> cov =
        calc.template covariance(surf, p0, 4.0 * _scale);
    // std::vector<typename SPACE::mat43> cov =
    //     calc.template curvature(surf, p0, 2.0 * _scale);

    for (auto &v : verts) {
      coordinate_type N = normals[i];
      coordinate_type x = p0[i];
      typename SPACE::mat3 U = cov[i].block(0, 0, 3, 3);
      coordinate_type u0 = U.col(0).transpose();
      coordinate_type u1 = U.col(1).transpose();
      coordinate_type u2 = U.col(2).transpose();

      // coordinate_type Np = u0 * u0.transpose() * N;
      real np0 = pow(fabs(u0.dot(N)), 1.0);
      real np1 = pow(fabs(u1.dot(N)), 1.0);
      real np2 = pow(fabs(u2.dot(N)), 1.0);

      coordinate_type s = cov[i].row(3);
      auto inf1 = [](real s0, real s1, real C, real p) {
        return pow(C * (s0 / s1 - 1.0), p) + 0.001;
      };
      real s12 = inf1(s[1], s[2], 3.0, 1.0);
      real s01 = inf1(s[0], s[1], 1.0, 1.0);
      real s20 = inf1(s[2], s[0], 1.0, 1.0);
      real m0 = pow(s[0] / (s[0] + s[1] + s[2]), 1.0);
      real m1 = pow(s[1] / (s[0] + s[1] + s[2]), 1.0);
      real m2 = pow(s[2] / (s[0] + s[1] + s[2]), 1.0);
      /*
            real sx = 1.0 * np0 * (m0 / s12) - //
                      0.1 * np1 * (m1 / s20) - //
                      0.1 * np2 * (m2 / s01);
      */
      // real r0 = -0.55, r1 = 3.6, r2 = -2.0;
      real r0 = 2.8, r1 = -1.5, r2 = -0.5;
      // real r0 = -2.0, r1 = -6.0, r2 = 5.0;

      real sx = r0 * np0 * m0 + //
                r1 * np1 * m1 + //
                r2 * np2 * m2;

      /*
      real sx = 1.0 * np0 * (pow(s[0] / s[1], 4) - 4.0);
      sx = std::min(1.0, sx);
      sx = std::max(-1.0, sx);
      */
      U.block(0, 0, 3, 1) *= s[0];
      U.block(0, 1, 3, 1) *= s[1];
      U.block(0, 2, 3, 1) *= s[2];
      // gg::geometry_logger::frame(U, x, 0.01);
      // gg::geometry_logger::line(x, x + 1e-1 * sx * N, vec4(1.0, 0.0,
      // 0.5, 1.0));

      vg[i] += 1e-2 * sx * N;
      i++;
    }
    coordinate_array v1 = calder::calcAvgVerts(surf, vg, p0, 0.75 * _scale);
    for (int i = 0; i < p0.size(); i++) {
      // gg::geometry_logger::line(p0[i], p0[i] + 10.0 * v1[i],
      //                           vec4(1.0, 0.0, 0.5, 1.0));
      p1[i] += 1.0 * v1[i];
    }
    return p1;
  }

  real _scale = 0.1;
};

template <typename SPACE>
class covariant_flow_experiment : public bend_experiment<SPACE> {
public:
  M2_TYPEDEFS;

  typedef std::shared_ptr<covariant_flow_experiment<SPACE>> ptr;
  static ptr create() {
    return std::make_shared<covariant_flow_experiment<SPACE>>();
  }

  virtual void _init() {
    asawa::affine<SPACE> aff;
    this->bend_experiment<SPACE>::_init();
  };

  virtual void _step(const int &frame) {

    if (!_flower) {
      _flower = covariant_flow<SPACE>::create(this->_integrator->_max);
    }

    std::vector<typename SPACE::vec3> p0 =
        asawa::ci::get_coordinates<SPACE>(this->_surf);
    std::vector<typename SPACE::vec3> p1 = _flower->step(this->_surf, p0);

    // this->willmore(p0);
    // this->bending(p0, p1);
    this->update(p1);
  };
  typename covariant_flow<SPACE>::ptr _flower;
};

} // namespace duchamp
#endif