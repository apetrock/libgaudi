/*
 *  conj_grad.hpp
 *  Manifold
 *
 *  Created by John Delaney on 5/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __ASAWA_LAPLACE_REF_MAT__
#define __ASAWA_LAPLACE_REF_MAT__

//#include "Eigen/src/SparseCore/SparseMatrix.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "manifold/asawa/datums.hpp"
#include "manifold/asawa/m2_refactor.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <math.h>

namespace bontecou_r {
class laplacian_base {
public:
  using index_t = int;
  using real = double;
  using triplet = Eigen::Triplet<real>;

  laplacian_base(asawa::manifold::ptr surf) : __M(surf) { this->init(); }

  ~laplacian_base() {}

  virtual void init() { std::cout << "in base init" << std::endl; }

  Eigen::SparseMatrix<real>
  build(std::function<real(asawa::manifold &M, index_t c, real &Km)> func_ij,
        std::function<real(asawa::manifold &M, real &Km)> func_ii) {

    std::vector<index_t> verts = __M->get_vert_range();
    std::vector<index_t> edges = __M->get_edge_vert_ids();

    std::vector<triplet> tripletList;
    tripletList.reserve(verts.size() + edges.size());

    int i = 0;
    for (auto v : __M->get_vert_range()) {

      real Km = 0.0;
      __M->for_each_vertex(
          v, [&Km, i, &tripletList, func_ij](index_t c, asawa::manifold &M) {
            index_t j = M.vert(M.next(c));
            real K = func_ij(M, c, Km);
            if (K > 0)
              tripletList.push_back(triplet(i, j, K));
          });
      Km = func_ii(*__M, Km);
      tripletList.push_back(triplet(i, i, Km));
      i++;
    }

    Eigen::SparseMatrix<real> mat(verts.size(), verts.size());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
  }

  Eigen::SparseMatrix<real>
  build(std::function<real(asawa::manifold &M, index_t c, real &Km)> func_ij) {
    return this->build(func_ij,
                       [](asawa::manifold &M, real &Km) -> real { return Km; });
  }

  asawa::manifold::ptr __M;
};
} // namespace bontecou_r

#endif
