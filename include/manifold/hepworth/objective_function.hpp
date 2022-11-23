

#ifndef __PENKO_OBJECTIVE__
#define __PENKO_OBJECTIVE__

#include "constraints.hpp"
#include "manifold/asawa/coordinate_interface.hpp"
#include "manifold/asawa/m2.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace hepworth {
template <typename SPACE> class objective_function {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  typedef std::shared_ptr<objective_function<SPACE>> ptr;

  static ptr create() { return std::make_shared<objective_function<SPACE>>(); }
  objective_function() {}

  virtual vecX get_x() const = 0;
  virtual void set_x(const vecX &x) = 0;
  virtual void init_rest(){};
  virtual sparmat get_inv_areas() = 0;
  virtual sparmat get_areas() = 0;

  virtual sparmat get_laplacian() = 0;

  virtual void preprocess(const vecX &v) = 0;
  virtual real calc_objective() = 0;

  virtual void visualize(const sparmat &M, const vecX &x){};
  virtual void visualize(const vecX &v, const vecX &x, vec4 color, real C){};

  virtual void limit_strain(vecX &x, real h){};
  virtual vecX calc_gradient(const vecX &x) = 0;
  virtual sparmat calc_hessian(const vecX &x) = 0;
};

template <typename SPACE>
class rosenbrock_function : public objective_function<SPACE> {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  typedef std::shared_ptr<rosenbrock_function<SPACE>> ptr;

  static ptr create() { return std::make_shared<rosenbrock_function<SPACE>>(); }
  rosenbrock_function() { X = vecX(2); }

  virtual vecX get_x() const { return X; };
  virtual void set_x(const vecX &x) { X = x; };

  virtual real calc_objective() {
    real ax = a - X[0];
    real yxx = X[1] - X[0] * X[0];

    return ax * ax + b * yxx * yxx;
  };

  virtual vecX project_gradient(const vecX &x) { return x; };
  virtual vecX calc_gradient(const vecX &x) {
    double xV = X[0];
    double yV = X[1];
    double dfdx = 2 * (-a + xV + 2 * b * xV * (xV * xV - yV));
    double dfdy = 2 * b * (yV - xV * xV);

    assert(std::isfinite(dfdx));
    assert(std::isfinite(dfdy));

    // Store it
    vecX grad(2);
    grad[0] += dfdx; //+ grad[0];
    grad[1] += dfdy; //+ grad[1]; /// ? Adds ? Or replace ? Or appends ?
    return grad;
  };

  virtual sparmat calc_hessian(const vecX &x) {
    double xV = X[0];
    double yV = X[1];

    // Empty the Tripletd ?
    std::vector<triplet> hess_triplets;

    hess_triplets.push_back(
        triplet(0, 0, -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2));
    hess_triplets.push_back(triplet(0, 1, -4 * b * xV));
    hess_triplets.push_back(triplet(1, 0, -4 * b * xV));
    hess_triplets.push_back(triplet(1, 1, 2 * b));

    sparmat H(2, 2);

    H.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return H;
  };

  real a = 1.0, b = 100.0;
  vecX X;
};

template <typename SPACE>
class constraint_set : public objective_function<SPACE> {

public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 4, 1> vec4;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 3, 3> mat33;

  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matX;
  typedef Eigen::SparseMatrix<real> sparmat;

  typedef std::shared_ptr<constraint_set<SPACE>> ptr;
  static ptr create(const surf_ptr surf) {
    return std::make_shared<constraint_set<SPACE>>(surf);
  }

  constraint_set(const surf_ptr surf) : _surf(surf) {
    _reg = 2.0 * asawa::ci::geometric_mean_length<SPACE>(_surf);
  };

  template <typename MAT> mat33 get33Block(MAT M, int i, int j) {
    real c00 = M.coeff(3 * i + 0, 3 * j + 0);
    real c01 = M.coeff(3 * i + 0, 3 * j + 1);
    real c02 = M.coeff(3 * i + 0, 3 * j + 2);

    real c10 = M.coeff(3 * i + 1, 3 * j + 0);
    real c11 = M.coeff(3 * i + 1, 3 * j + 1);
    real c12 = M.coeff(3 * i + 1, 3 * j + 2);

    real c20 = M.coeff(3 * i + 2, 3 * j + 0);
    real c21 = M.coeff(3 * i + 2, 3 * j + 1);
    real c22 = M.coeff(3 * i + 2, 3 * j + 2);
    mat33 Mb;
    Mb.row(0) << c00, c01, c02;
    Mb.row(1) << c10, c11, c12;
    Mb.row(2) << c20, c21, c22;

    return Mb;
  }

  void minmax(mat33 M, real &min, real &max) {
    coordinate_type t0 = M.block(0, 0, 3, 1);
    coordinate_type t1 = M.block(0, 1, 3, 1);
    coordinate_type t2 = M.block(0, 2, 3, 1);
    real nt0 = t0.norm();
    real nt1 = t1.norm();
    real nt2 = t2.norm();
    min = std::min(min, std::min(nt0, std::min(nt1, nt2)));
    max = std::max(max, std::max(nt0, std::max(nt1, nt2)));
  }

  template <typename MAT>
  void minmax(MAT M, int i, int j, real &min, real &max) {
    mat33 Mb = get33Block(M, i, j);
    minmax(Mb, min, max);
  }

  template <typename MAT>
  void visualizeBlock(MAT M, coordinate_type c, int i, int j, real C = 1.0) {
    mat33 Mb = get33Block(M, i, j);
    gg::geometry_logger::get_instance().frame(Mb, c, C);
  }

  virtual void visualize(const sparmat &Ms, const vecX &x) {
    vertex_array verts = _surf->get_vertices();
    edge_array edges = _surf->get_edges();
    real min = std::numeric_limits<real>::max();
    real max = 0;

    for (int k = 0; k < Ms.outerSize(); ++k)
      for (typename sparmat::InnerIterator it(Ms, k); it; ++it) {
        it.value();
        it.row();   // row index
        it.col();   // col index (here it is equal to k)
        it.index(); // inner index, here it is equal to it.row()
      }

    for (auto e : edges) {
      vec3 c = vec_interface<SPACE>::center(e, x);
      int ei = e->v1()->vertex()->position_in_set();
      int ej = e->v2()->vertex()->position_in_set();
      minmax(Ms, ei, ej, min, max);
    }
    for (auto v : verts) {
      vec3 c = vec_interface<SPACE>::coord(v, x);
      int vi = v->position_in_set();
      minmax(Ms, vi, vi, min, max);
    }
    max = std::min(100.0, max);
    real scale = 1.0 / (max - min);
    std::cout << " min max scale: " << min << " " << max << " " << scale
              << std::endl;
#if 1
    for (auto e : edges) {

      vec3 c = vec_interface<SPACE>::center(e, x);
      int ei = e->v1()->vertex()->position_in_set();
      int ej = e->v2()->vertex()->position_in_set();
      visualizeBlock(Ms, c, ei, ej, 0.1 * scale);
    }
#endif
#if 1
    for (auto v : verts) {

      vec3 c = vec_interface<SPACE>::coord(v, x);
      int vi = v->position_in_set();
      visualizeBlock(Ms, c, vi, vi, 0.1 * scale);
    }
#endif
  };

  virtual void visualize(const vecX &g, const vecX &x,
                         vec4 color = vec4(0.0, 0.6, 0.7, 1.0), real C = 1.0) {
    vertex_array verts = _surf->get_vertices();
#if 1
    for (auto v : verts) {
      int vi = v->position_in_set();
      vec3 c = vec_interface<SPACE>::coord(v, x);
      vec3 gi = vec_interface<SPACE>::coord(v, g);
      gg::geometry_logger::get_instance().line(c, c + C * gi, color);
    }
#endif
  };

  virtual void update_constraints(const vecX &x) {
    for (auto c : constraints) {
      c->update(x);
    }
  }

  virtual void set_positions(const coordinate_array &positions) {
    _N = positions.size();
    vecX x = vec_interface<SPACE>::to(positions);
    set_x(x);
  }

  virtual coordinate_array get_positions() const {
    coordinate_array positions(_N, coordinate_type(0, 0, 0));
    vec_interface<SPACE>::from(positions, _x);
    return positions;
  }

  virtual vecX get_x() const {
    std::cout << " x norm: " << _x.norm() << std::endl;
    return _x;
  };

  virtual void set_x(const vecX &x) {
    update_constraints(x);
    _x = x;
  };

  virtual real calc_objective() {
    real C = 0.0;
    for (auto c : constraints) {
      C += c->evaluate_constraint();
    }
    return C;
  }

  virtual void init_rest() {
    for (auto c : constraints) {
      c->init_rest();
    }
  }

  virtual void preprocess(const vecX &v) {
    for (auto c : constraints) {
      c->preprocess(v);
    }
  }

  virtual void limit_strain(vecX &vel, real h) {
    vertex_array verts = _surf->get_vertices();
    int i = 0;
    for (auto v : verts) {
      mat3 M = mat3::Zero();
      coordinate_type p0 = ci::get_coordinate<SPACE>(v);
      for_each_vertex<SPACE>(v, [&M, p0](face_vertex_ptr fv) {
        coordinate_type p1 = ci::get_coordinate<SPACE>(fv->next()->vertex());
        coordinate_type dp = p1 - p0;
        M += dp * dp.transpose();
      });

      Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU);
      // const mat3 U = svd.matrixU();
      const vec3 S = svd.singularValues();
      coordinate_type vl = vec_interface<SPACE>::from(vel, i);
      // real Sm = sqrt(S[1]) / h;
      real Sm = sqrt(S[0]);
      real vn = vl.norm();
      if (vn > Sm)
        vl *= Sm / vn;
      vec_interface<SPACE>::to(vl, vel, i);
      i++;
    }
  }

  virtual vecX calc_gradient(const vecX &x) {
    int N = _N;
    vecX G(3 * N);
    G.setZero();
    for (auto c : constraints) {
      c->fill_gradient(G);
    }
    // visualize(G, x);
    // std::cout << " grad sum: " << G.sum() << std::endl;
    // std::cout << " grad norm: " << G.norm() << std::endl;

    return G;
  }

  sparmat get_inv_areas() {
    int N = _N;
    vertex_array vertices = _surf->get_vertices();
    std::vector<triplet> area_triplets;
    int i = 0;
    for (auto &v : vertices) {
      real a = ci::area<SPACE>(v);
      real ia = a <= 1e-8 ? 1.0 : 1.0 / a;
      area_triplets.push_back(triplet(3 * i + 0, 3 * i + 0, ia));
      area_triplets.push_back(triplet(3 * i + 1, 3 * i + 1, ia));
      area_triplets.push_back(triplet(3 * i + 2, 3 * i + 2, ia));
      i++;
    }

    sparmat A(3 * N, 3 * N);

    A.setFromTriplets(area_triplets.begin(), area_triplets.end());
    return A;
  }

  sparmat get_areas() {
    int N = _N;
    vertex_array vertices = _surf->get_vertices();
    std::vector<triplet> area_triplets;
    int i = 0;
    for (auto &v : vertices) {
      real a = ci::area<SPACE>(v);
      area_triplets.push_back(triplet(3 * i + 0, 3 * i + 0, a));
      area_triplets.push_back(triplet(3 * i + 1, 3 * i + 1, a));
      area_triplets.push_back(triplet(3 * i + 2, 3 * i + 2, a));
      i++;
    }

    sparmat A(3 * N, 3 * N);

    A.setFromTriplets(area_triplets.begin(), area_triplets.end());
    return A;
  }

  sparmat get_laplacian() {
    int N = _N;
    auto vertices = _surf->get_vertices();
    auto edges = _surf->get_edges();

    typedef Eigen::Triplet<double> triplet;
    std::vector<triplet> tripletList;
    tripletList.reserve(vertices.size() + 2 * edges.size());

    int i = 0;
    for (auto v : _surf->get_vertices()) {

      real Km = 0.0;
      asawa::for_each_vertex<SPACE>(
          v, [&Km, i, &tripletList](face_vertex_ptr fv) {
            int j = fv->next()->vertex()->position_in_set();
            face_vertex_ptr fvp = fv->vprev()->next();
            face_vertex_ptr fvn = fv->vnext()->next();
            real cotp = asawa::ci::abs_cotan<SPACE>(fvp);
            real cotn = asawa::ci::abs_cotan<SPACE>(fvn);
            assert(!isnan(cotp));
            assert(!isnan(cotn));
            real K = (cotp + cotn);
            if (K > 0) {
              tripletList.push_back(triplet(3 * i + 0, 3 * j + 0, K));
              tripletList.push_back(triplet(3 * i + 1, 3 * j + 1, K));
              tripletList.push_back(triplet(3 * i + 2, 3 * j + 2, K));
            }
          });
      tripletList.push_back(triplet(3 * i + 0, 3 * i + 0, Km));
      tripletList.push_back(triplet(3 * i + 1, 3 * i + 1, Km));
      tripletList.push_back(triplet(3 * i + 2, 3 * i + 2, Km));
      i++;
    }

    sparmat A(3 * N, 3 * N);

    A.setFromTriplets(tripletList.begin(), tripletList.end());
    return A;
  }

  virtual sparmat calc_hessian(const vecX &x) {
    int N = _N;
    std::vector<triplet> hess_triplets;
    for (auto c : constraints) {
      c->get_hess_triplets(hess_triplets);
    }
    sparmat H(3 * N, 3 * N);

    H.setFromTriplets(hess_triplets.begin(), hess_triplets.end());

    // visualize(H, x);
    //       std::cout << H.diagonal() << std::endl;
    size_t Nr = H.rows();
    size_t Nc = H.cols();

    // H.setIdentity();
    {
      real sum = 0.0;
      for (auto t : hess_triplets) {
        sum += t.value();
      }

      std::cout << " trip sum: " << sum << std::endl;
      std::cout << "    H sum: " << H.sum() << std::endl;
      std::cout << "    H norm: " << H.norm() << std::endl;
    }
    // visualize(H, x);

    //    H = I - M * H;
    //      std::cout << "    I-H sum: " << H.sum() << std::endl;
    //      std::cout << "    I-H norm: " << H.norm() << std::endl;

    return H;
  }

  void add_constraint(typename constraint<SPACE>::ptr c) {
    this->constraints.push_back(c);
  }

  size_t _N = 0;
  std::vector<typename constraint<SPACE>::ptr> constraints;
  surf_ptr _surf;
  vecX _x;
  real _reg = 0.01;
};
} // namespace hepworth

#endif