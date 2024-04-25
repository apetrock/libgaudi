//
#ifndef __HEP_RIGID_BODY_CONSTRAINTS__
#define __HEP_RIGID_BODY_CONSTRAINTS__

#include "Eigen/src/Geometry/AngleAxis.h"
#include "block_constraint.hpp"
#include "gaudi/common.h"
#include "shell_constraints.hpp"
#include "sim_block.hpp"

namespace gaudi {
namespace hepworth {
namespace block {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class rigid_coupline_constraint : public block_constraint {
public:
  DEFINE_CREATE_FUNC(rigid_coupline_constraint)

  rigid_coupline_constraint(const std::vector<index_t> &ids, const vec3 &p, const real &w,
         std::vector<sim_block::ptr> blocks)
      : block_constraint(ids, w, blocks), _p(p) {}
  virtual std::string name() { return typeid(*this).name(); }
  virtual void project(const vecX &q, vecX &p) {
    
    index_t i = this->_ids[0];
    index_t j = this->_ids[1];
    
    vec3 q0 = _blocks[0]->get_vec3(i, q);
    p.block(_id0, 0, 3, 1) = _w * _p;
  }
  virtual void fill_A(index_t &id0, std::vector<trip> &triplets) {
    _id0 = id0;
    index_t i = _blocks[0]->get_offset_idx(this->_ids[0]);
    for (int ax = 0; ax < 3; ax++)
      triplets.push_back(trip(_id0 + ax, i + ax, _w));
    id0 += 3;
  }
  vec3 _p;
};


} // namespace block
} // namespace hepworth
} // namespace gaudi
#endif