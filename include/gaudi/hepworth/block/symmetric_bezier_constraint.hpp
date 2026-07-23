#ifndef __HEP_SYMMETRIC_BEZIER_CONSTRAINT__
#define __HEP_SYMMETRIC_BEZIER_CONSTRAINT__

#include <array>
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

// Soft 5-point absolute projector onto a data-mirrored quartic Bezier.
// ids: [i-2, i-1, i, i+1, i+2] (rod corner / vert indices in the position block).
class symmetric_bezier : public block_constraint {
public:
  typedef std::shared_ptr<symmetric_bezier> ptr;

  static ptr create(const std::vector<index_t> &ids, const real &w,
                    std::vector<sim_block::ptr> blocks) {
    return std::make_shared<symmetric_bezier>(ids, w, blocks);
  }

  symmetric_bezier(const std::vector<index_t> &ids, const real &w,
                   std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks) {
    if (_ids.size() >= 5) {
      // Seed tangent from rest ordering; updated each successful project.
      _T_cached = vec3::UnitX();
    }
  }

  virtual std::string name() { return typeid(*this).name(); }

  virtual void project(const vecX &q, vecX &p) {
    if (_ids.size() < 5 || _blocks.empty())
      return;

    std::vector<vec3> pts(5);
    for (int k = 0; k < 5; ++k)
      pts[k] = _blocks[0]->get_vec3(_ids[k], q);

    // Prefer geometric tangent from neighbors when available.
    vec3 T = pts[3] - pts[1];
    if (T.norm() < 1e-12)
      T = pts[4] - pts[0];
    if (T.norm() > 1e-12)
      _T_cached = T.normalized();

    const albers::symmetric_bezier_fit fit =
        albers::fit_symmetric_quartic_bezier(pts, _T_cached);
    if (!fit.ok) {
      // No-op: leave RHS zero contribution for this constraint.
      return;
    }
    _T_cached = fit.frame.e_T;

    for (int k = 0; k < 5; ++k) {
      const vec3 target = fit.eval_at_point(pts[k]);
      p.block(_id0 + 3 * k, 0, 3, 1) = _w * target;
    }
  }

  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    for (int k = 0; k < 5; ++k) {
      const index_t ik = _blocks[0]->get_offset_idx(_ids[k]);
      for (int ax = 0; ax < 3; ++ax)
        triplets.push_back(trip(_id0 + 3 * k + ax, ik + ax, _w));
    }
    id0 += 15;
  }

private:
  vec3 _T_cached = vec3::UnitX();
};

inline void init_symmetric_bezier(
    const asawa::rod::rod &rod,
    std::vector<projection_constraint::ptr> &constraints, const real &w,
    std::vector<sim_block::ptr> blocks) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;

  if (w <= 0.0)
    return;

  for (int i = 0; i < rod.corner_count(); i++) {
    CornerId ci = corner_id(i);
    CornerId im1 = rod.prev(ci);
    if (im1 < 0)
      continue;
    CornerId im2 = rod.prev(im1);
    if (im2 < 0)
      continue;
    CornerId ip1 = rod.next(ci);
    if (ip1 < 0)
      continue;
    CornerId ip2 = rod.next(ip1);
    if (ip2 < 0)
      continue;

    constraints.push_back(
        symmetric_bezier::create({im2, im1, ci, ip1, ip2}, w, blocks));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif
