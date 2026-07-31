#ifndef __HEP_SYMMETRIC_BEZIER_CONSTRAINT__
#define __HEP_SYMMETRIC_BEZIER_CONSTRAINT__

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gaudi/albers/symmetric_bezier.hpp"
#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/hepworth/block/block_constraint.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

// Soft absolute projector onto a screw-symmetrized cubic Bezier.
// Discrete π-screw pairing (central inversion through the station) — coil /
// flip symmetry, not planar mirror. Default stencil: ±4 → 9 points.
inline constexpr int k_symmetric_bezier_half_width = 4;

class symmetric_bezier : public block_constraint {
public:
  typedef std::shared_ptr<symmetric_bezier> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<symmetric_bezier>(ids, w, blocks);
  }

  symmetric_bezier(const std::vector<index_t> &ids, const real &w,
                   std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {}

  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    const int n = static_cast<int>(_ids.size());
    if (n < 5 || _blocks.empty())
      return;

    std::vector<vec3> pts(n);
    for (int k = 0; k < n; ++k)
      pts[k] = _blocks[0]->get_vec3(_ids[k], q);

    const int mid = n / 2;
    vec3 T = pts[std::min(mid + 1, n - 1)] - pts[std::max(mid - 1, 0)];
    if (T.norm() < 1e-12)
      T = pts[n - 1] - pts[0];
    if (T.norm() > 1e-12)
      _T_cached = T.normalized();

    const albers::symmetric_bezier_fit fit =
        albers::fit_symmetric_cubic_bezier(pts, _T_cached);
    if (!fit.ok)
      return;
    _T_cached = fit.frame.e_T;

    for (int k = 0; k < n; ++k) {
      const vec3 target = fit.eval_at_point(pts[k]);
      p.block(_id0 + 3 * k, 0, 3, 1) = _w * target;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    const int n = static_cast<int>(_ids.size());
    for (int k = 0; k < n; ++k) {
      const index_t ik = _blocks[0]->get_offset_idx(_ids[k]);
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(_id0 + 3 * k + ax, ik + ax, _w));
    }
    id0 += 3 * n;
  }

private:
  vec3 _T_cached = vec3::UnitX();
};

inline void init_symmetric_bezier(
    const asawa::rod::rod &rod,
    std::vector<projection_constraint::ptr> &constraints, const real &w,
    std::vector<sim_block::ptr> blocks,
    int half_width = k_symmetric_bezier_half_width) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;

  if (w <= 0.0 || half_width < 2)
    return;

  for (int i = 0; i < rod.corner_count(); i++) {
    CornerId ci = corner_id(i);
    std::vector<index_t> ids;
    ids.reserve(2 * half_width + 1);

    std::vector<CornerId> left;
    left.reserve(half_width);
    CornerId c = ci;
    bool ok = true;
    for (int k = 0; k < half_width; ++k) {
      c = rod.prev(c);
      if (c < 0) {
        ok = false;
        break;
      }
      left.push_back(c);
    }
    if (!ok)
      continue;

    for (int k = static_cast<int>(left.size()) - 1; k >= 0; --k)
      ids.push_back(left[k]);
    ids.push_back(ci);

    c = ci;
    for (int k = 0; k < half_width; ++k) {
      c = rod.next(c);
      if (c < 0) {
        ok = false;
        break;
      }
      ids.push_back(c);
    }
    if (!ok)
      continue;

    constraints.push_back(symmetric_bezier::create(ids, w, blocks));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
