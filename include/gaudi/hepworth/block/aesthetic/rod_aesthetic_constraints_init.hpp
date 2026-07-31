#ifndef __HEP_ROD_AESTHETIC_CONSTRAINTS_INIT__
#define __HEP_ROD_AESTHETIC_CONSTRAINTS_INIT__

#include <vector>

#include "gaudi/asawa/rod/rod.hpp"
#include "gaudi/common.h"
#include "gaudi/hepworth/block/aesthetic/rod_aesthetic_constraints.hpp"
#include "gaudi/hepworth/block/sim_block.hpp"
#include "gaudi/hepworth/projection_constraint.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

void init_helicity(const asawa::rod::rod &rod,
                   std::vector<projection_constraint::ptr> &constraints,
                   const real &w, std::vector<sim_block::ptr> blocks) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;

  for (int i = 0; i < rod.corner_count(); i++) {
    CornerId ci = corner_id(i);
    CornerId ip0 = rod.prev(ci);
    if (ip0 < 0)
      continue;
    CornerId ip1 = rod.prev(ip0);
    if (ip1 < 0)
      continue;
    CornerId ip2 = rod.prev(ip1);
    if (ip2 < 0)
      continue;
    CornerId ip3 = rod.prev(ip2);
    if (ip3 < 0)
      continue;
    CornerId in0 = rod.next(ci);
    if (in0 < 0)
      continue;
    CornerId in1 = rod.next(in0);
    if (in1 < 0)
      continue;
    CornerId in2 = rod.next(in1);
    if (in2 < 0)
      continue;
    CornerId in3 = rod.next(in2);
    if (in3 < 0)
      continue;
    constraints.push_back(helicitiy::create(
        {ci, ip3, ip2, ip1, ip0, ci, in0, in1, in2, in3}, w, blocks));
  }
}

void init_uniform_step(const asawa::rod::rod &R,
                       std::vector<projection_constraint::ptr> &constraints,
                       const real &w, std::vector<sim_block::ptr> blocks,
                       real induced_twist = 0.0) {
  using asawa::rod::corner_id;
  using asawa::rod::CornerId;
  if (w <= 0.0)
    return;
  for (int i = 0; i < R.corner_count(); ++i) {
    CornerId c = corner_id(i);
    CornerId im1 = R.prev(c);
    if (im1 == corner_id(-1))
      continue;
    CornerId im2 = R.prev(im1);
    if (im2 == corner_id(-1))
      continue;
    CornerId ip1 = R.next(c);
    if (ip1 == corner_id(-1))
      continue;
    CornerId ip2 = R.next(ip1);
    if (ip2 == corner_id(-1))
      continue;
    constraints.push_back(uniform_step::create({im2, im1, c, ip1, ip2}, w,
                                               blocks, induced_twist));
  }
}

} // namespace block
} // namespace hepworth
} // namespace gaudi

#endif // __HEP_ROD_AESTHETIC_CONSTRAINTS_INIT__
