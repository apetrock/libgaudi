#ifndef GAUDI_ENGMAN_KKT_HPP
#define GAUDI_ENGMAN_KKT_HPP

/// The Engman KKT solve (`scripts/sympy/cr_poisson/MATH.md` section 17.3):
///
///     [ H   J^T     ] [ dq ]   [ a ]
///     [ J  -I / rho ] [ dl ] = [ b ],
///
/// factored once and solved for several right-hand sides.
///
/// - `full`: one `SimplicialLDLT` of the whole quasi-definite matrix.
/// - `schur`: the variable blocks `y` paired one-to-one with 3-row
///   constraint blocks (the coupling's store columns, `J_y = -c I`) are
///   eliminated with their multipliers, per block, when `H` holds them only
///   block-diagonally (`H_yy` 3 x 3 per block, no `H_xy`) and no other row
///   reaches them. With `D = B H_yy^-1 B^T + I / rho` (3 x 3 per block) the
///   rest is
///
///       [ H_xx + J_xy^T D^-1 J_xy   J_xx^T  ] [ dx    ]
///       [ J_xx                     -I / rho ] [ dl_x  ] = ...,
///
///   positions-only and SPD when every constraint row is paired (the
///   `weingarten_gn` factor's pattern), quasi-definite with the unpaired
///   rows otherwise. `dl_y` and `dy` are recovered per block.

#include "gaudi/common.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <chrono>
#include <utility>
#include <vector>

namespace gaudi {
namespace engman {

using matS = Eigen::SparseMatrix<real>;
using trip = Eigen::Triplet<real>;
using mat3 = Eigen::Matrix<real, 3, 3>;

enum class kkt_mode { full, schur, automatic };

class kkt_system {
public:
  /// `H` (n x n, positive definite), `J` (m x n), and the pairs
  /// (constraint row of a 3-row block, first column of its 3-column `y`
  /// block). `mode` `schur` falls back to `full` when the structure does not
  /// hold; `schur()` says which was built.
  bool factor(const matS &H, const matS &J, real rho, const std::vector<std::pair<index_t, index_t>> &pairs,
              kkt_mode mode) {
    const auto t0 = std::chrono::steady_clock::now();
    _n = H.rows();
    _m = J.rows();
    _rho = rho;
    _schur = mode != kkt_mode::full && !pairs.empty() && build_schur(H, J, pairs);
    if (!_schur)
      build_full(H, J);
    _seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return ok();
  }

  bool ok() const { return _K.info() == Eigen::Success; }
  bool schur() const { return _schur; }
  /// Dimension and `L` nonzeros of the last factored matrix.
  index_t factored_size() const { return _K.rows(); }
  index_t factor_nonzeros() const { return _K.matrixL().nestedExpression().nonZeros(); }
  double seconds() const { return _seconds; }

  void solve(const vecX &a, const vecX &b, vecX &dq, vecX &dl) const {
    if (!_schur) {
      vecX rhs(_n + _m);
      rhs.head(_n) = a;
      rhs.tail(_m) = b;
      const vecX s = _K.solve(rhs);
      dq = s.head(_n);
      dl = s.tail(_m);
      return;
    }
    const index_t nx = static_cast<index_t>(_xcols.size()), mx = static_cast<index_t>(_xrows.size());
    // e = B H_yy^-1 a_y - b_y per block.
    std::vector<vec3> e(_blocks.size());
    for (size_t k = 0; k < _blocks.size(); ++k) {
      const block &bk = _blocks[k];
      e[k] = bk.B * (bk.Hinv * vec3(a.segment<3>(bk.col))) - vec3(b.segment<3>(bk.row));
    }
    vecX rhs(nx + mx);
    for (index_t i = 0; i < nx; ++i)
      rhs[i] = a[_xcols[static_cast<size_t>(i)]];
    vecX De(3 * static_cast<index_t>(_blocks.size()));
    for (size_t k = 0; k < _blocks.size(); ++k)
      De.segment<3>(3 * static_cast<index_t>(k)) = _blocks[k].Dinv * e[k];
    rhs.head(nx) -= _Jxy.transpose() * De;
    for (index_t r = 0; r < mx; ++r)
      rhs[nx + r] = b[_xrows[static_cast<size_t>(r)]];
    const vecX s = _K.solve(rhs);
    dq = vecX::Zero(_n);
    dl = vecX::Zero(_m);
    for (index_t i = 0; i < nx; ++i)
      dq[_xcols[static_cast<size_t>(i)]] = s[i];
    for (index_t r = 0; r < mx; ++r)
      dl[_xrows[static_cast<size_t>(r)]] = s[nx + r];
    const vecX Jdx = _Jxy * s.head(nx);
    for (size_t k = 0; k < _blocks.size(); ++k) {
      const block &bk = _blocks[k];
      const vec3 lk = bk.Dinv * (vec3(Jdx.segment<3>(3 * static_cast<index_t>(k))) + e[k]);
      dl.segment<3>(bk.row) = lk;
      dq.segment<3>(bk.col) = bk.Hinv * (vec3(a.segment<3>(bk.col)) - bk.B.transpose() * lk);
    }
  }

private:
  struct block {
    index_t row = 0, col = 0;
    mat3 B, Hinv, Dinv;
  };

  void build_full(const matS &H, const matS &J) {
    std::vector<trip> kt;
    kt.reserve(static_cast<size_t>(H.nonZeros() + 2 * J.nonZeros() + _m));
    for (int k = 0; k < H.outerSize(); ++k)
      for (matS::InnerIterator it(H, k); it; ++it)
        kt.emplace_back(it.row(), it.col(), it.value());
    for (int k = 0; k < J.outerSize(); ++k)
      for (matS::InnerIterator it(J, k); it; ++it) {
        kt.emplace_back(_n + it.row(), it.col(), it.value());
        kt.emplace_back(it.col(), _n + it.row(), it.value());
      }
    for (index_t r = 0; r < _m; ++r)
      kt.emplace_back(_n + r, _n + r, -real(1) / _rho);
    matS K(_n + _m, _n + _m);
    K.setFromTriplets(kt.begin(), kt.end());
    _K.compute(K);
  }

  /// False (nothing built) when the paired blocks are not eliminable.
  bool build_schur(const matS &H, const matS &J, const std::vector<std::pair<index_t, index_t>> &pairs) {
    std::vector<int> yblock(static_cast<size_t>(_n), -1), rblock(static_cast<size_t>(_m), -1);
    _blocks.assign(pairs.size(), block{});
    for (size_t k = 0; k < pairs.size(); ++k) {
      const auto [row, col] = pairs[k];
      if (row < 0 || row + 3 > _m || col < 0 || col + 3 > _n)
        return false;
      for (int c = 0; c < 3; ++c) {
        if (yblock[static_cast<size_t>(col + c)] >= 0 || rblock[static_cast<size_t>(row + c)] >= 0)
          return false;
        yblock[static_cast<size_t>(col + c)] = static_cast<int>(k);
        rblock[static_cast<size_t>(row + c)] = static_cast<int>(k);
      }
      _blocks[k].row = row;
      _blocks[k].col = col;
      _blocks[k].B.setZero();
    }
    // H holds y only in its own block.
    std::vector<mat3> Hyy(pairs.size(), mat3::Zero());
    for (int k = 0; k < H.outerSize(); ++k)
      for (matS::InnerIterator it(H, k); it; ++it) {
        const int bi = yblock[static_cast<size_t>(it.row())], bj = yblock[static_cast<size_t>(it.col())];
        if (bi < 0 && bj < 0)
          continue;
        if (bi != bj)
          return false;
        const index_t c0 = _blocks[static_cast<size_t>(bi)].col;
        Hyy[static_cast<size_t>(bi)](it.row() - c0, it.col() - c0) += it.value();
      }
    // A paired row reaches y only in its own block; an unpaired row not at all.
    _xcols.clear();
    std::vector<index_t> xidx(static_cast<size_t>(_n), -1);
    for (index_t i = 0; i < _n; ++i)
      if (yblock[static_cast<size_t>(i)] < 0) {
        xidx[static_cast<size_t>(i)] = static_cast<index_t>(_xcols.size());
        _xcols.push_back(i);
      }
    _xrows.clear();
    std::vector<index_t> ridx(static_cast<size_t>(_m), -1);
    for (index_t r = 0; r < _m; ++r)
      if (rblock[static_cast<size_t>(r)] < 0) {
        ridx[static_cast<size_t>(r)] = static_cast<index_t>(_xrows.size());
        _xrows.push_back(r);
      }
    std::vector<trip> jxy, jxx;
    for (int k = 0; k < J.outerSize(); ++k)
      for (matS::InnerIterator it(J, k); it; ++it) {
        const int rb = rblock[static_cast<size_t>(it.row())], cb = yblock[static_cast<size_t>(it.col())];
        if (cb >= 0) {
          if (rb != cb)
            return false;
          block &bk = _blocks[static_cast<size_t>(rb)];
          bk.B(it.row() - bk.row, it.col() - bk.col) += it.value();
        } else if (rb >= 0) {
          jxy.emplace_back(3 * rb + (it.row() - _blocks[static_cast<size_t>(rb)].row),
                           xidx[static_cast<size_t>(it.col())], it.value());
        } else {
          jxx.emplace_back(ridx[static_cast<size_t>(it.row())], xidx[static_cast<size_t>(it.col())], it.value());
        }
      }
    const index_t nx = static_cast<index_t>(_xcols.size()), mx = static_cast<index_t>(_xrows.size());
    std::vector<trip> dt;
    for (size_t k = 0; k < _blocks.size(); ++k) {
      block &bk = _blocks[k];
      Eigen::LDLT<mat3> hl(Hyy[k]);
      if (hl.info() != Eigen::Success || !(hl.vectorD().minCoeff() > real(0)))
        return false;
      bk.Hinv = hl.solve(mat3::Identity());
      const mat3 D = bk.B * bk.Hinv * bk.B.transpose() + mat3::Identity() / _rho;
      bk.Dinv = D.inverse();
      if (!bk.Dinv.allFinite())
        return false;
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          dt.emplace_back(3 * static_cast<index_t>(k) + r, 3 * static_cast<index_t>(k) + c, bk.Dinv(r, c));
    }
    _Jxy.resize(3 * static_cast<index_t>(_blocks.size()), nx);
    _Jxy.setFromTriplets(jxy.begin(), jxy.end());
    matS Dinv(3 * static_cast<index_t>(_blocks.size()), 3 * static_cast<index_t>(_blocks.size()));
    Dinv.setFromTriplets(dt.begin(), dt.end());
    std::vector<trip> ht;
    for (int k = 0; k < H.outerSize(); ++k)
      for (matS::InnerIterator it(H, k); it; ++it)
        if (yblock[static_cast<size_t>(it.row())] < 0)
          ht.emplace_back(xidx[static_cast<size_t>(it.row())], xidx[static_cast<size_t>(it.col())], it.value());
    matS Hxx(nx, nx);
    Hxx.setFromTriplets(ht.begin(), ht.end());
    const matS Hred = Hxx + matS(_Jxy.transpose() * Dinv * _Jxy);
    matS Jxx(mx, nx);
    Jxx.setFromTriplets(jxx.begin(), jxx.end());
    std::vector<trip> kt;
    for (int k = 0; k < Hred.outerSize(); ++k)
      for (matS::InnerIterator it(Hred, k); it; ++it)
        kt.emplace_back(it.row(), it.col(), it.value());
    for (int k = 0; k < Jxx.outerSize(); ++k)
      for (matS::InnerIterator it(Jxx, k); it; ++it) {
        kt.emplace_back(nx + it.row(), it.col(), it.value());
        kt.emplace_back(it.col(), nx + it.row(), it.value());
      }
    for (index_t r = 0; r < mx; ++r)
      kt.emplace_back(nx + r, nx + r, -real(1) / _rho);
    matS K(nx + mx, nx + mx);
    K.setFromTriplets(kt.begin(), kt.end());
    _K.compute(K);
    return true;
  }

  index_t _n = 0, _m = 0;
  real _rho = 1;
  bool _schur = false;
  double _seconds = 0;
  Eigen::SimplicialLDLT<matS> _K;
  std::vector<block> _blocks;
  std::vector<index_t> _xcols, _xrows;
  matS _Jxy;
};

} // namespace engman
} // namespace gaudi

#endif
