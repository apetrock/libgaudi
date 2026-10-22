#ifndef GAUDI_ASAWA_ROD_RESIDUAL_FOREST_HPP
#define GAUDI_ASAWA_ROD_RESIDUAL_FOREST_HPP

/// Residual forest as a rod chain:
///   rod0 → rod1+(adjacent_tuples,R) → rod2+(adjacent_tuples,R) → …
/// Restrict = half-rate vertex collapse along next/prev (not a loop flag).

#include "gaudi/arp/procrustes_pyramid.hpp"
#include "gaudi/asawa/datums.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <vector>

namespace gaudi {
namespace asawa {
namespace rod {

/// Fine corners involved in one coarse sample: prev, survivor (center), next.
using adjacent_tuple = std::array<index_t, 3>;

/// Crane J stencil on a fine corner: black {i0,-1}, white {i0,i1}.
using j_idx_pair = std::array<index_t, 2>;

/// Corner-indexed J index pairs on the fine rod datum pool.
struct j_idx_datum : public asawa::datum {
  using ptr = std::shared_ptr<j_idx_datum>;
  static ptr create(const std::vector<j_idx_pair> &data) {
    return std::make_shared<j_idx_datum>(data);
  }
  explicit j_idx_datum(const std::vector<j_idx_pair> &data)
      : asawa::datum(asawa::VERTEX), __data(data) {}

  std::vector<j_idx_pair> &data() { return __data; }
  const std::vector<j_idx_pair> &data() const { return __data; }

  void do_alloc(const size_t &sz) override {
    __data.resize(__data.size() + sz, j_idx_pair{index_t(-1), index_t(-1)});
  }
  void do_resize(const size_t &sz) override {
    __data.resize(sz, j_idx_pair{index_t(-1), index_t(-1)});
  }
  size_t size() override { return __data.size(); }
  void calc(const rod & /*R*/, const index_t & /*i*/, const index_t & /*prim*/,
            const real & /*C*/,
            const std::vector<index_t> & /*vals*/) override {}
  void do_map(const index_t /*i*/, const index_t /*it*/) override {}
  void do_permute(const std::vector<index_t> &permute) override {
    std::vector<j_idx_pair> n = __data;
    const size_t nmin = std::min(__data.size(), permute.size());
    for (size_t i = 0; i < nmin; ++i) {
      if (permute[i] >= 0 && static_cast<size_t>(permute[i]) < __data.size())
        n[i] = __data[static_cast<size_t>(permute[i])];
    }
    __data.swap(n);
  }
  void write(FILE * /*file*/) const override {}
  void read(FILE * /*file*/) override {}
  void clear() override { __data.clear(); }

  std::vector<j_idx_pair> __data;
};

inline mat3 point_window_procrustes(const std::vector<vec3> &rest,
                                    const std::vector<vec3> &live, index_t center,
                                    int k, bool loop) {
  const int n = static_cast<int>(rest.size());
  if (n < 2 || live.size() != rest.size())
    return mat3::Identity();
  const int half = std::max(1, k / 2);
  int i0 = static_cast<int>(center) - half;
  int i1 = static_cast<int>(center) + half;
  if (!loop) {
    i0 = std::max(0, i0);
    i1 = std::min(n - 1, i1);
  }
  vec3 cr = vec3::Zero();
  vec3 cl = vec3::Zero();
  int count = 0;
  auto wrap = [n, loop](int i) -> int {
    if (loop) {
      i %= n;
      if (i < 0)
        i += n;
      return i;
    }
    return i;
  };
  for (int i = i0; i <= i1; ++i) {
    const int j = wrap(i);
    if (!loop && (j < 0 || j >= n))
      continue;
    cr += rest[static_cast<size_t>(j)];
    cl += live[static_cast<size_t>(j)];
    ++count;
  }
  if (count < 2)
    return mat3::Identity();
  cr /= real(count);
  cl /= real(count);
  mat3 M = mat3::Zero();
  for (int i = i0; i <= i1; ++i) {
    const int j = wrap(i);
    if (!loop && (j < 0 || j >= n))
      continue;
    M += (live[static_cast<size_t>(j)] - cl) *
         (rest[static_cast<size_t>(j)] - cr).transpose();
  }
  Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  if (svd.singularValues()[1] < real(1e-10))
    return mat3::Identity();
  return arp::procrustes_rotation(M);
}

class residual_forest {
public:
  using ptr = std::shared_ptr<residual_forest>;

  static ptr create(rod::ptr R, int k = 5, int min_coarse = -1) {
    ptr F = std::make_shared<residual_forest>();
    F->build(R, k, min_coarse);
    return F;
  }

  bool valid() const { return _valid; }
  void invalidate() { _valid = false; }

  int n_levels() const { return static_cast<int>(_rods.size()); }
  int k() const { return _k; }
  int k_window() const { return _k; }
  int min_coarse() const { return _min_coarse; }
  bool loop() const { return _loop; }

  rod::ptr level_rod(int level) const {
    assert_valid();
    assert_level(level);
    return _rods[static_cast<size_t>(level)];
  }

  const std::vector<adjacent_tuple> &adjacent(int coarse_level) const {
    assert_valid();
    assert_level(coarse_level);
    if (coarse_level < 1)
      throw std::runtime_error("residual_forest: no adjacent on L0");
    return _adjacent[static_cast<size_t>(coarse_level)];
  }

  const std::vector<mat3> &cached_R(int coarse_level) const {
    assert_valid();
    assert_level(coarse_level);
    return _R[static_cast<size_t>(coarse_level)];
  }

  const std::vector<index_t> &parent_of(int fine_level) const {
    assert_valid();
    assert_level(fine_level);
    if (fine_level + 1 >= n_levels())
      throw std::runtime_error("residual_forest: no parent at coarsest");
    return _parent[static_cast<size_t>(fine_level)];
  }

  /// Survivor fine corner id for each coarse vert (adjacent[·][1]).
  const std::vector<index_t> &keep(int fine_level) const {
    assert_valid();
    assert_level(fine_level);
    if (fine_level + 1 >= n_levels())
      throw std::runtime_error("residual_forest: no keep at coarsest");
    return _keep[static_cast<size_t>(fine_level)];
  }

  /// Datum slot on fine rod for J index pairs; -1 if none (coarsest).
  index_t j_idx_slot(int fine_level) const {
    assert_valid();
    assert_level(fine_level);
    if (fine_level < 0 ||
        fine_level >= static_cast<int>(_j_slot.size()))
      return index_t(-1);
    return _j_slot[static_cast<size_t>(fine_level)];
  }

  const std::vector<j_idx_pair> &j_pairs(int fine_level) const {
    assert_valid();
    assert_level(fine_level);
    const index_t slot = j_idx_slot(fine_level);
    if (slot < 0)
      throw std::runtime_error("residual_forest: no J pairs at level");
    auto d = std::dynamic_pointer_cast<j_idx_datum>(
        _rods[static_cast<size_t>(fine_level)]->const_get_datum(slot));
    if (!d)
      throw std::runtime_error("residual_forest: bad J datum");
    return d->data();
  }

  /// Crane prolong: fine = J * coarse.
  std::vector<vec3> prolong_J(int fine_level,
                              const std::vector<vec3> &coarse) const {
    assert_valid();
    assert_level(fine_level);
    const auto &pairs = j_pairs(fine_level);
    const size_t nf = pairs.size();
    std::vector<vec3> fine(nf, vec3::Zero());
    for (size_t i = 0; i < nf; ++i) {
      const index_t i0 = pairs[i][0];
      const index_t i1 = pairs[i][1];
      if (i0 < 0)
        continue;
      if (i1 < 0) {
        if (static_cast<size_t>(i0) < coarse.size())
          fine[i] = coarse[static_cast<size_t>(i0)];
      } else if (static_cast<size_t>(i0) < coarse.size() &&
                 static_cast<size_t>(i1) < coarse.size()) {
        fine[i] = real(0.5) * (coarse[static_cast<size_t>(i0)] +
                               coarse[static_cast<size_t>(i1)]);
      }
    }
    return fine;
  }

  /// Crane restrict: coarse = Jᵀ * fine.
  std::vector<vec3> restrict_JT(int fine_level,
                                const std::vector<vec3> &fine) const {
    assert_valid();
    assert_level(fine_level);
    if (fine_level + 1 >= n_levels())
      return {};
    const auto &pairs = j_pairs(fine_level);
    const size_t nc = _rods[static_cast<size_t>(fine_level) + 1]->corner_count();
    std::vector<vec3> coarse(nc, vec3::Zero());
    const size_t nf = std::min(pairs.size(), fine.size());
    for (size_t i = 0; i < nf; ++i) {
      const index_t i0 = pairs[i][0];
      const index_t i1 = pairs[i][1];
      if (i0 < 0)
        continue;
      if (i1 < 0) {
        if (static_cast<size_t>(i0) < nc)
          coarse[static_cast<size_t>(i0)] += fine[i];
      } else {
        if (static_cast<size_t>(i0) < nc)
          coarse[static_cast<size_t>(i0)] += real(0.5) * fine[i];
        if (static_cast<size_t>(i1) < nc)
          coarse[static_cast<size_t>(i1)] += real(0.5) * fine[i];
      }
    }
    return coarse;
  }

  const std::vector<vec3> &live(int level) const {
    assert_valid();
    assert_level(level);
    return _rods[static_cast<size_t>(level)]->x();
  }

  const std::vector<vec3> &rest(int level) const {
    assert_valid();
    assert_level(level);
    return _rest[static_cast<size_t>(level)];
  }

  /// Rest residual vs parent at build (fine levels only). Empty at coarsest.
  const std::vector<vec3> &residual(int fine_level) const {
    assert_valid();
    assert_level(fine_level);
    return _delta[static_cast<size_t>(fine_level)];
  }

  std::vector<vec3> pickle(int level) const {
    return live(level);
  }

  rod::ptr pickle_rod(int level) const {
    assert_valid();
    assert_level(level);
    if (level == 0)
      return _rods[0];
    // Copy so callers can mutate without aliasing the forest level rod.
    return rod::create(_rods[static_cast<size_t>(level)]->x(), _loop);
  }

  /// Write absolute positions into level rod (same indexing as that rod's x).
  void write_live(int level, const std::vector<vec3> &x_abs) {
    assert_valid();
    assert_level(level);
    auto &x = _rods[static_cast<size_t>(level)]->x();
    if (x_abs.size() != x.size())
      return;
    x = x_abs;
    _rods[static_cast<size_t>(level)]->_update_frames();
  }

  /// Half-rate collapse restrict: seed rod[fine+1].x from fine survivors.
  /// Coarse positions from the prolongation normal equations, per coordinate.
  static bool pinv_restrict_positions(const std::vector<j_idx_pair> &pairs,
                                      const std::vector<vec3> &fine_x,
                                      std::vector<vec3> &coarse_x) {
    const int nf = static_cast<int>(pairs.size());
    const int nc = static_cast<int>(coarse_x.size());
    if (nf <= 0 || nc <= 0)
      return false;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(nf) * 2);
    for (int i = 0; i < nf; ++i) {
      const index_t i0 = pairs[static_cast<size_t>(i)][0];
      const index_t i1 = pairs[static_cast<size_t>(i)][1];
      if (i0 < 0 || i0 >= nc)
        continue;
      if (i1 < 0) {
        trips.emplace_back(i, static_cast<int>(i0), 1.0);
      } else if (i1 < nc) {
        trips.emplace_back(i, static_cast<int>(i0), 0.5);
        trips.emplace_back(i, static_cast<int>(i1), 0.5);
      }
    }
    if (trips.empty())
      return false;
    Eigen::SparseMatrix<double> J(nf, nc);
    J.setFromTriplets(trips.begin(), trips.end());
    const Eigen::SparseMatrix<double> JTJ = J.transpose() * J;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
    lu.compute(JTJ);
    if (lu.info() != Eigen::Success)
      return false;
    for (int c = 0; c < 3; ++c) {
      Eigen::VectorXd xf(nf);
      for (int i = 0; i < nf; ++i)
        xf(i) = (static_cast<size_t>(i) < fine_x.size())
                    ? fine_x[static_cast<size_t>(i)][c]
                    : 0.0;
      const Eigen::VectorXd xc = lu.solve(J.transpose() * xf);
      if (lu.info() != Eigen::Success)
        return false;
      for (int j = 0; j < nc; ++j)
        coarse_x[static_cast<size_t>(j)][c] = xc(j);
    }
    return true;
  }

  void restrict_up(int fine_level) {
    assert_valid();
    assert_level(fine_level);
    if (fine_level + 1 >= n_levels())
      return;
    const size_t F = static_cast<size_t>(fine_level);
    const size_t C = F + 1;
    const auto &fine_x = _rods[F]->x();
    auto &coarse_x = _rods[C]->x();
    // Reference Coarsen places the coarse curve at (JᵀJ)⁻¹ Jᵀ x.
    if (!pinv_restrict_positions(j_pairs(fine_level), fine_x, coarse_x)) {
      const auto &adj = _adjacent[C];
      for (size_t j = 0; j < adj.size(); ++j) {
        const index_t cen = adj[j][1];
        if (cen < 0 || static_cast<size_t>(cen) >= fine_x.size())
          continue;
        coarse_x[j] = fine_x[static_cast<size_t>(cen)];
      }
    }
    _rods[C]->_update_frames();
  }

  void cache_kabsch(int level) {
    assert_valid();
    assert_level(level);
    cache_R_level(level);
  }

  /// One-hop inject: fine.x[i] += R[parent]ᵀ * residual[parent].
  void inject_residual_down(int coarse_level,
                            const std::vector<vec3> &residual) {
    assert_valid();
    assert_level(coarse_level);
    if (coarse_level < 1)
      return;
    const int fine = coarse_level - 1;
    const size_t C = static_cast<size_t>(coarse_level);
    const size_t F = static_cast<size_t>(fine);
    auto &x_fine = _rods[F]->x();
    if (residual.size() != _rods[C]->x().size())
      return;
    const auto &parent = _parent[F];
    if (parent.size() != x_fine.size())
      return;
    for (size_t i = 0; i < x_fine.size(); ++i) {
      const index_t p = parent[i];
      if (p < 0)
        continue;
      const mat3 &Rm = _R[C][static_cast<size_t>(p)];
      x_fine[i] += Rm.transpose() * residual[static_cast<size_t>(p)];
    }
    _rods[F]->_update_frames();
  }

  void commit_l0() {
    assert_valid();
    if (_rods.empty() || !_rods[0])
      return;
    _rods[0]->_update_frames();
  }

  /// Absolute set: write live, cache R; re-encode δ vs parent if any.
  void set(int level, const std::vector<vec3> &x_abs) {
    write_live(level, x_abs);
    cache_R_level(level);
    if (level + 1 >= n_levels())
      return;
    const size_t L = static_cast<size_t>(level);
    const size_t C = L + 1;
    auto &delta = _delta[L];
    if (delta.size() != _rods[L]->x().size())
      delta.assign(_rods[L]->x().size(), vec3::Zero());
    const auto &parent = _parent[L];
    const auto &fx = _rods[L]->x();
    const auto &cx = _rods[C]->x();
    cache_R_level(level + 1);
    for (size_t i = 0; i < fx.size(); ++i) {
      const index_t p = parent[i];
      if (p < 0)
        continue;
      const mat3 &Rm = _R[C][static_cast<size_t>(p)];
      delta[i] = Rm.transpose() * (fx[i] - cx[static_cast<size_t>(p)]);
    }
  }

  /// Reconstruct finer from coarser: x_f = c + R δ_rest. Tests / debug only.
  void prolong() {
    assert_valid();
    for (int ell = n_levels() - 1; ell >= 1; --ell) {
      const size_t C = static_cast<size_t>(ell);
      const size_t F = C - 1;
      cache_R_level(ell);
      auto &fx = _rods[F]->x();
      const auto &cx = _rods[C]->x();
      const auto &parent = _parent[F];
      const auto &delta = _delta[F];
      for (size_t i = 0; i < fx.size(); ++i) {
        const index_t p = parent[i];
        if (p < 0)
          continue;
        const mat3 &Rm = _R[C][static_cast<size_t>(p)];
        fx[i] = cx[static_cast<size_t>(p)] + Rm * delta[i];
      }
      _rods[F]->_update_frames();
    }
    commit_l0();
  }

private:
  bool _valid = false;
  bool _loop = true;
  int _k = 5;
  int _min_coarse = 10;

  std::vector<rod::ptr> _rods;
  std::vector<std::vector<vec3>> _rest;
  std::vector<std::vector<vec3>> _delta; ///< rest residual vs parent (fine)
  std::vector<std::vector<adjacent_tuple>> _adjacent; ///< [0] unused
  std::vector<std::vector<mat3>> _R;
  std::vector<std::vector<index_t>> _parent; ///< fine→coarse parent
  std::vector<std::vector<index_t>> _keep;   ///< fine survivor corner ids
  std::vector<index_t> _j_slot;              ///< fine-level rod datum slot for J pairs

  void assert_valid() const {
    if (!_valid)
      throw std::runtime_error("residual_forest: invalid");
  }
  void assert_level(int level) const {
    if (level < 0 || level >= n_levels())
      throw std::runtime_error("residual_forest: bad level");
  }

  static std::vector<CornerId> ordered_live_corners(const rod &R) {
    return const_cast<rod &>(R).get_ordered_verts();
  }

  struct strand_walk {
    std::vector<index_t> corner;
    bool returned = false;
  };

  /// One walk per component. `returned` is next() coming back to the seed.
  std::vector<strand_walk> collect_strands(const rod &R) const {
    const int n = static_cast<int>(R.corner_count());
    std::vector<char> seen(static_cast<size_t>(std::max(n, 0)), 0);
    std::vector<strand_walk> walks;
    for (int i = 0; i < n; ++i) {
      if (seen[static_cast<size_t>(i)])
        continue;
      CornerId s = corner_id(i);
      if (R.next(s) < corner_id(0) && R.prev(s) < corner_id(0))
        continue;
      {
        CornerId p = s;
        int guard = n + 1;
        while (guard-- > 0) {
          const CornerId q = R.prev(p);
          if (q < corner_id(0) || q == s)
            break;
          p = q;
        }
        s = p;
      }
      const int si = static_cast<int>(s);
      if (si < 0 || si >= n || seen[static_cast<size_t>(si)])
        continue;

      strand_walk w;
      const CornerId seed = s;
      CornerId j = s;
      int guard = n + 1;
      while (!(j < corner_id(0)) && guard-- > 0) {
        const int ji = static_cast<int>(j);
        if (ji < 0 || ji >= n || seen[static_cast<size_t>(ji)])
          break;
        seen[static_cast<size_t>(ji)] = 1;
        w.corner.push_back(static_cast<index_t>(ji));
        const CornerId nxt = R.next(j);
        if (nxt == seed) {
          w.returned = true;
          break;
        }
        j = nxt;
      }
      if (w.corner.size() >= 2)
        walks.push_back(std::move(w));
    }
    return walks;
  }

  int k_use_at(int level, index_t center) const {
    int k_use = _k;
    if (!_loop) {
      const int n = static_cast<int>(_rest[static_cast<size_t>(level)].size());
      const int c = static_cast<int>(center);
      if (c <= 1 || c >= n - 2)
        k_use = std::min(_k, 3);
    }
    return k_use;
  }

  void cache_R_level(int level) {
    const size_t L = static_cast<size_t>(level);
    const auto &live = _rods[L]->x();
    const auto &rest = _rest[L];
    const int n = static_cast<int>(live.size());
    _R[L].assign(static_cast<size_t>(n), mat3::Identity());
    if (rest.size() != live.size())
      return;
    // Kabsch in coarse-rod index space (contiguous 0..n-1 for ℓ>0; L0 may be
    // pool-ordered — use ordered walk positions mirrored into rest/live copy).
    if (level == 0) {
      // Build dense ordered buffers for window Kabsch.
      const auto ord = ordered_live_corners(*_rods[0]);
      const int no = static_cast<int>(ord.size());
      std::vector<vec3> rest_o(static_cast<size_t>(no));
      std::vector<vec3> live_o(static_cast<size_t>(no));
      for (int i = 0; i < no; ++i) {
        const size_t ci = static_cast<size_t>(ord[static_cast<size_t>(i)]);
        rest_o[static_cast<size_t>(i)] = rest[ci];
        live_o[static_cast<size_t>(i)] = live[ci];
      }
      for (int i = 0; i < no; ++i) {
        const mat3 Rm = point_window_procrustes(
            rest_o, live_o, static_cast<index_t>(i),
            k_use_at(0, static_cast<index_t>(i)), _loop);
        _R[0][static_cast<size_t>(ord[static_cast<size_t>(i)])] = Rm;
      }
      return;
    }
    for (int i = 0; i < n; ++i) {
      _R[L][static_cast<size_t>(i)] = point_window_procrustes(
          rest, live, static_cast<index_t>(i),
          k_use_at(level, static_cast<index_t>(i)), _loop);
    }
  }

  /// Collapse-restrict: every other corner along next/prev. A white corner
  /// averages the previous and next black corners on that same walk.
  void append_collapsed_level(rod::ptr fine) {
    const std::vector<strand_walk> walks = collect_strands(*fine);
    int n_live = 0;
    for (const strand_walk &w : walks)
      n_live += static_cast<int>(w.corner.size());
    if (n_live / 2 <= _min_coarse)
      return;

    const int ncorn = static_cast<int>(fine->corner_count());
    std::vector<char> is_black(static_cast<size_t>(std::max(ncorn, 0)), 0);
    for (const strand_walk &w : walks) {
      const int m = static_cast<int>(w.corner.size());
      if (m < 2)
        continue;
      if (w.returned) {
        for (int t = 0; t < m; t += 2)
          is_black[static_cast<size_t>(w.corner[static_cast<size_t>(t)])] = 1;
      } else {
        is_black[static_cast<size_t>(w.corner.front())] = 1;
        is_black[static_cast<size_t>(w.corner.back())] = 1;
        for (int t = 2; t + 1 < m; t += 2)
          is_black[static_cast<size_t>(w.corner[static_cast<size_t>(t)])] = 1;
      }
    }

    int nc = 0;
    for (char b : is_black)
      if (b)
        ++nc;
    if (nc < 2 || nc <= _min_coarse)
      return;

    const auto &fx = fine->x();
    std::vector<index_t> keep;
    std::vector<vec3> coarse_x;
    std::vector<adjacent_tuple> adj;
    std::vector<index_t> coarse_of(static_cast<size_t>(ncorn), index_t(-1));
    keep.reserve(static_cast<size_t>(nc));
    coarse_x.reserve(static_cast<size_t>(nc));
    adj.reserve(static_cast<size_t>(nc));
    for (const strand_walk &w : walks) {
      for (index_t c : w.corner) {
        if (c < 0 || static_cast<size_t>(c) >= is_black.size() ||
            !is_black[static_cast<size_t>(c)])
          continue;
        coarse_of[static_cast<size_t>(c)] = static_cast<index_t>(keep.size());
        keep.push_back(c);
        coarse_x.push_back(fx[static_cast<size_t>(c)]);
        const CornerId cen = corner_id(static_cast<int>(c));
        const CornerId cp = fine->prev(cen);
        const CornerId cn = fine->next(cen);
        const index_t ip = cp < corner_id(0) ? index_t(-1) : static_cast<index_t>(cp);
        const index_t in = cn < corner_id(0) ? index_t(-1) : static_cast<index_t>(cn);
        adj.push_back({ip, c, in});
      }
    }
    nc = static_cast<int>(keep.size());

    auto walk_black = [&](index_t start, bool forward) -> index_t {
      if (start < 0 || start >= ncorn)
        return index_t(-1);
      if (is_black[static_cast<size_t>(start)])
        return start;
      CornerId a = forward ? fine->next(corner_id(static_cast<int>(start)))
                           : fine->prev(corner_id(static_cast<int>(start)));
      int guard = ncorn;
      while (!(a < corner_id(0)) && guard-- > 0) {
        const int ai = static_cast<int>(a);
        if (ai == static_cast<int>(start))
          return index_t(-1);
        if (ai >= 0 && ai < ncorn && is_black[static_cast<size_t>(ai)])
          return static_cast<index_t>(ai);
        a = forward ? fine->next(a) : fine->prev(a);
      }
      return index_t(-1);
    };

    std::vector<index_t> parent(fine->corner_count(), index_t(-1));
    std::vector<vec3> delta(fine->corner_count(), vec3::Zero());
    std::vector<j_idx_pair> pairs(fine->corner_count(),
                                  j_idx_pair{index_t(-1), index_t(-1)});

    for (const strand_walk &w : walks) {
      for (index_t c : w.corner) {
        if (c < 0 || c >= ncorn)
          continue;
        const size_t ci = static_cast<size_t>(c);
        const bool black = is_black[ci] != 0;
        index_t b_prev = index_t(-1);
        index_t b_next = index_t(-1);
        int hops_prev = ncorn + 1;
        int hops_next = ncorn + 1;
        if (black) {
          b_prev = c;
          hops_prev = 0;
        } else {
          CornerId a = fine->prev(corner_id(static_cast<int>(c)));
          int guard = ncorn;
          int hops = 0;
          while (!(a < corner_id(0)) && guard-- > 0) {
            const int ai = static_cast<int>(a);
            if (ai == static_cast<int>(c))
              break;
            ++hops;
            if (ai >= 0 && ai < ncorn && is_black[static_cast<size_t>(ai)]) {
              b_prev = static_cast<index_t>(ai);
              hops_prev = hops;
              break;
            }
            a = fine->prev(a);
          }
          a = fine->next(corner_id(static_cast<int>(c)));
          guard = ncorn;
          hops = 0;
          while (!(a < corner_id(0)) && guard-- > 0) {
            const int ai = static_cast<int>(a);
            if (ai == static_cast<int>(c))
              break;
            ++hops;
            if (ai >= 0 && ai < ncorn && is_black[static_cast<size_t>(ai)]) {
              b_next = static_cast<index_t>(ai);
              hops_next = hops;
              break;
            }
            a = fine->next(a);
          }
        }

        const index_t nearest =
            (hops_prev <= hops_next) ? (b_prev >= 0 ? b_prev : b_next)
                                     : (b_next >= 0 ? b_next : b_prev);
        if (nearest >= 0 &&
            static_cast<size_t>(nearest) < coarse_of.size()) {
          const index_t pc = coarse_of[static_cast<size_t>(nearest)];
          parent[ci] = pc;
          if (pc >= 0 && static_cast<size_t>(pc) < coarse_x.size())
            delta[ci] = fx[ci] - coarse_x[static_cast<size_t>(pc)];
        }

        if (black) {
          pairs[ci] = {coarse_of[ci], index_t(-1)};
        } else {
          if (b_prev < 0)
            b_prev = walk_black(c, false);
          if (b_next < 0)
            b_next = walk_black(c, true);
          const index_t cL =
              (b_prev >= 0 && static_cast<size_t>(b_prev) < coarse_of.size())
                  ? coarse_of[static_cast<size_t>(b_prev)]
                  : index_t(-1);
          const index_t cR =
              (b_next >= 0 && static_cast<size_t>(b_next) < coarse_of.size())
                  ? coarse_of[static_cast<size_t>(b_next)]
                  : index_t(-1);
          if (cL < 0)
            pairs[ci] = {cR, index_t(-1)};
          else if (cR < 0)
            pairs[ci] = {cL, index_t(-1)};
          else
            pairs[ci] = {cL, cR};
        }
      }
    }

    // Reference PolyCurveNetwork::Coarsen: X_c = (JᵀJ)⁻¹ Jᵀ X_f.
    pinv_restrict_positions(pairs, fx, coarse_x);
    for (size_t i = 0; i < parent.size(); ++i) {
      const index_t p = parent[i];
      if (p < 0 || static_cast<size_t>(p) >= coarse_x.size())
        continue;
      if (i < fx.size())
        delta[i] = fx[i] - coarse_x[static_cast<size_t>(p)];
    }

    rod::ptr coarse;
    for (const strand_walk &w : walks) {
      std::vector<vec3> xs;
      for (index_t c : w.corner) {
        if (c >= 0 && static_cast<size_t>(c) < is_black.size() &&
            is_black[static_cast<size_t>(c)]) {
          const index_t cj = coarse_of[static_cast<size_t>(c)];
          xs.push_back(cj >= 0 && static_cast<size_t>(cj) < coarse_x.size()
                           ? coarse_x[static_cast<size_t>(cj)]
                           : fx[static_cast<size_t>(c)]);
        }
      }
      if (xs.size() < 2)
        continue;
      if (!coarse)
        coarse = rod::create(xs, w.returned);
      else
        coarse->insert_strand(xs, w.returned);
    }
    if (!coarse)
      return;
    const size_t fine_level = _rods.size() - 1;
    _keep[fine_level] = std::move(keep);
    _parent[fine_level] = std::move(parent);
    _delta[fine_level] = std::move(delta);

    const index_t jslot = [&]() -> index_t {
      // Reuse existing J datum on this rod (avoids leaking slots on live L0).
      for (auto &dp : fine->get_data()) {
        if (auto j = std::dynamic_pointer_cast<j_idx_datum>(dp)) {
          j->data() = std::move(pairs);
          return j->slot_id;
        }
      }
      return fine->insert_datum(j_idx_datum::create(std::move(pairs)));
    }();
    if (fine_level >= _j_slot.size())
      _j_slot.resize(fine_level + 1, index_t(-1));
    _j_slot[fine_level] = jslot;

    _rods.push_back(coarse);
    _rest.push_back(coarse_x);
    _delta.emplace_back(); // coarse has no parent residual yet
    _adjacent.push_back(std::move(adj));
    _R.push_back(std::vector<mat3>(static_cast<size_t>(nc), mat3::Identity()));
    _parent.emplace_back();
    _keep.emplace_back();
    _j_slot.push_back(index_t(-1)); // coarsest (for now) has no child J
  }

  void build(rod::ptr R, int k, int min_coarse) {
    _k = std::max(3, k | 1);
    _min_coarse = (min_coarse < 0) ? (2 * _k) : min_coarse;
    _rods.clear();
    _rest.clear();
    _delta.clear();
    _adjacent.clear();
    _R.clear();
    _parent.clear();
    _keep.clear();
    _j_slot.clear();

    if (!R) {
      _valid = false;
      return;
    }

    const auto strands = collect_strands(*R);
    _loop = strands.size() == 1 && strands[0].returned;
    _rods.push_back(R);
    _rest.push_back(R->x());
    _delta.emplace_back();
    _adjacent.emplace_back();
    _R.push_back(std::vector<mat3>(R->corner_count(), mat3::Identity()));
    _parent.emplace_back();
    _keep.emplace_back();
    _j_slot.push_back(index_t(-1));

    while (true) {
      const int n_ord =
          static_cast<int>(ordered_live_corners(*_rods.back()).size());
      if (n_ord / 2 <= _min_coarse)
        break;
      const size_t before = _rods.size();
      append_collapsed_level(_rods.back());
      if (_rods.size() == before)
        break;
    }

    _valid = n_levels() >= 1;
  }
};

} // namespace rod
} // namespace asawa
} // namespace gaudi

#endif // GAUDI_ASAWA_ROD_RESIDUAL_FOREST_HPP
