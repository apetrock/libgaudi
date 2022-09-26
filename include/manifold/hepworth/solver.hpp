#ifndef __PENKO_SOLVER__
#define __PENKO_SOLVER__

#include "objective_function.hpp"

/*archihepworth*/
namespace hepworth {

template <typename SPACE> class solver {

  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  virtual void solve(typename objective_function<SPACE>::ptr fcn) {}
  virtual void solve(typename objective_function<SPACE>::ptr fcn, const vecX &v,
                     const real &h) {}
};

template <typename SPACE> class gradient_descent_solver : public solver<SPACE> {
public:
  M2_TYPEDEFS;
  typedef Eigen::Matrix<real, 2, 1> vec2;
  typedef Eigen::Matrix<real, 6, 1> vec6;
  typedef Eigen::Matrix<real, 3, 1> vec3;
  typedef Eigen::Matrix<real, 3, 2> mat32;
  typedef Eigen::Matrix<real, 6, 6> mat66;
  typedef Eigen::Matrix<real, Eigen::Dynamic, 1> vecX;
  typedef Eigen::SparseMatrix<real> sparmat;
  gradient_descent_solver(){

  };

  real calc_lambda(const vecX &g, typename objective_function<SPACE>::ptr fcn) {
    real C = fcn->calc_objective();
    real denom = g.transpose() * g;
    return C / denom;
  }

  virtual void solve(typename objective_function<SPACE>::ptr fcn) {
    vecX xk_init = fcn->get_x();

    vecX xk0 = fcn->get_x();

    real tol = 1;
    real alpha = 1.0;
    real beta = 0.75;
    int k = 0;
    while (tol > 1e-3) {
      if (k++ > 200)
        break;
      fcn->preprocess();
      vecX g = fcn->calc_gradient(xk0);
      //      fcn->visualize(g, xk0, vec4(0.75, 0.0, 0.75, 1.0));

      std::cout << '\r' << std::flush;
      std::cout << "=== : " << k;
      std::cout << ", gnorm: " << g.norm();

      // g = fcn->project_gradient(g);
      //   std::cout << "gnorm: " << g.norm() << std::endl;

      real lambda = calc_lambda(g, fcn);
      std::cout << ", lambda: " << lambda;

      // fcn->visualize(g, xk0, vec4(0.0, 0.75, 0.75, 1.0));

      vecX xk1 = xk0 - lambda * g;
      tol = (xk1 - xk0).norm();
      // std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << " ====>>tol: " << tol;

      fcn->set_x(xk1);
      _g = g;
      std::swap(xk0, xk1);
    }
    std::cout << std::endl;
    //    vecX dx = xk0 - xk_init;

    // fcn->visualize(dx, xk0, vec4(0.75, 0.0, 0.75, 1.0));
    //  dx = fcn->project_gradient(dx);
    //  fcn->visualize(dx, xk0, vec4(0.0, 0.75, 0.75, 1.0));
    //     fcn->set_x(xk0 + dx);
  }

  virtual vecX get_gradient() { return _g; };
  vecX _g;
};

template <typename SPACE> class newton_raphson_solver : public solver<SPACE> {
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

  newton_raphson_solver(){

  };

  virtual vecX solve(sparmat &A, vecX &b) {
#if 1
    // Eigen::SimplicialLLT<sparmat> solver;
    Eigen::SimplicialLDLT<sparmat> solver;
    // Eigen::ConjugateGradient<sparmat> solver;
    solver.compute(A);
//    std::cout << " solver det: " << solver.determinant() << std::endl;
#else
    // Eigen::SparseQR<sparmat, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SparseLU<sparmat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
#endif
    if (solver.info() != Eigen::Success) {
      // decomposition failed
      std::cout << ".....decomposition error! " << std::endl;
    }
    vecX x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      // solving failed
      std::cout << ".....solve error! " << std::endl;
    }

    return x;
  }

  real calc_lambda(const vecX &g, typename objective_function<SPACE>::ptr fcn) {
    real C = fcn->calc_objective();
    real denom = g.transpose() * g;
    return C / denom;
  }

  virtual void solve(typename objective_function<SPACE>::ptr fcn) {

    vecX xk0 = fcn->get_x();

    real tol = 1;
    real alpha = 1.0;
    real beta = 1.0;
    int k = 0;
    int max = 10;
    real sig = 1e-4;
    while (tol > sig && k < max) {
      fcn->preprocess();

      std::vector<triplet> hess_triplets;
      std::vector<triplet> grad_triplets;
      vecX g = fcn->calc_gradient(xk0);
      sparmat H = fcn->calc_hessian(xk0);

#if 1
      vecX h = solve(H, g);

      /*
            std::cout << "h min/max : " << h.minCoeff() << " " << h.maxCoeff()
                      << std::endl;
            std::cout << "h norm: " << h.norm() << std::endl;
            std::cout << "h sum: " << h.sum() << std::endl;
            std::cout << "g-h: " << (g - h).norm() << std::endl;
      */
      vecX xk1 = xk0 - h;

      // fcn->visualize(H, xk0);

      // return;
      // break;
#else
      real lambda = calc_lambda(g, fcn);
      vecX xk1 = H * xk0 - g;
      xk1 = solve(H, xk1);
      _g = g;
#endif
      // h = fcn->project_gradient(h, positions);

      // real lambda = calc_lambda(h, fcn);
      alpha *= beta;
      tol = (xk1 - xk0).norm();

      if (tol < sig || k == max - 1) {
        fcn->visualize(h, xk0, vec4(0.75, 0.0, 0.75, 1.0));
        fcn->visualize(g, xk0, vec4(0.0, 0.75, 0.75, 1.0));
      }
      // if (k > 1)
      //   break;
      //  std::cout << " xki norm: " << xki.norm() << std::endl;
      std::cout << "====>>tol: " << tol << std::endl; //<< '\r'
                                                      //      if (k > 1)
                                                      //        break;
      fcn->set_x(xk1);
      std::swap(xk0, xk1);
      k++;
    }
  }
};

template <typename SPACE> class backward_euler : public solver<SPACE> {
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

  backward_euler(){

  };

  virtual vecX solve(sparmat &A, vecX &b) {
#if 1
    // Eigen::SimplicialLLT<sparmat> solver;
    Eigen::SimplicialLDLT<sparmat> solver;
    // Eigen::ConjugateGradient<sparmat> solver;
    solver.compute(A);
//    std::cout << " solver det: " << solver.determinant() << std::endl;
#else
    // Eigen::SparseQR<sparmat, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SparseLU<sparmat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
#endif
    if (solver.info() != Eigen::Success) {
      // decomposition failed
      std::cout << ".....decomposition error! " << std::endl;
    }
    vecX x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      // solving failed
      std::cout << ".....solve error! " << std::endl;
    }

    return x;
  }

  real calc_lambda(const vecX &g, typename objective_function<SPACE>::ptr fcn) {
    real C = fcn->calc_objective();
    real denom = g.transpose() * g;
    return C / denom;
  }

  virtual void solve(typename objective_function<SPACE>::ptr fcn,
                     const vecX &v0, const real &hb) {

    vecX xk0 = fcn->get_x();
    vecX xk1 = xk0;
    vecX v1 = v0;
    fcn->set_x(xk1);

    real tol = 1;
    int k = 0;
    int max = 4;
    real sig = 1e-6;
    real hs = hb;
    while (tol > sig) {
      if (k > max)
        break;
      std::vector<triplet> hess_triplets;
      std::vector<triplet> grad_triplets;
      sparmat M = fcn->get_areas();

      real ih = 1.0 / hs;
      fcn->preprocess(ih * M * v1);
      vecX g = fcn->calc_gradient(xk1);
      sparmat K = fcn->calc_hessian(xk1);

      sparmat I(K.rows(), K.cols());
      I.setIdentity();

      /*
      std::cout << "H: " << H.rows() << " " << H.cols() << std::endl;
      std::cout << "v: " << v.rows() << " " << v.cols() << std::endl;
      std::cout << "g: " << g.rows() << " " << g.cols() << std::endl;
      std::cout << "W: " << W.rows() << " " << W.cols() << std::endl;
      */
      //   xk0 = xk1;

      real gtg = g.transpose() * g;
      real gtMv = g.transpose() * M * v1;
      real hmax = gtMv / gtg;

      real h = std::min(hs, fabs(hmax));
      hs = std::min(2.0 * h, hs);
      std::cout << "h vs hmax: " << h << " " << hmax << std::endl;

      sparmat H = M - h * h * K;
      vecX r = M * v1 - h * g;
      // vecX r = K * vp + h * W * g;

      std::cout << "g: " << g.norm() << std::endl;
      std::cout << "v0: " << v0.norm() << std::endl;
      std::cout << "v1: " << v1.norm() << std::endl;

      std::cout << "r: " << r.norm() << std::endl;

      vecX dv = solve(H, r);

      xk1 += h * dv;
      // v1 = (xk1 - xk0) / h;
      v1 = v1 -= dv;
      /*
            fcn->visualize(v0, xk0, vec4(0.5, 0.5, 0.0, 1.0), 10.0);
            fcn->visualize(v1, xk0, vec4(0.75, 0.75, 0.0, 1.0), 10.0);
            fcn->visualize(dv, xk0, vec4(1.0, 1.0, 0.0, 1.0), 10.0);
            fcn->visualize(g, xk0, vec4(0.75, 0.0, 0.0, 1.0), 10.0);
            fcn->visualize(h * K * v0, xk0, vec4(0.0, 0.75, 0.0, 1.0), 10.0);
            fcn->visualize(r, xk0, vec4(0.75, 0.0, 0.75, 1.0), 10.0);
      */
      tol = r.norm();

      std::cout << "====>>tol: " << tol << std::endl; //<< '\r'

      fcn->set_x(xk1);
      k++;
    }
  }
};

} // namespace hepworth
#endif