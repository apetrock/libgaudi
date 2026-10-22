#ifndef GAUDI_HEPWORTH_CACHE_COMPUTE_HPP
#define GAUDI_HEPWORTH_CACHE_COMPUTE_HPP

#include <memory>

#include "gaudi/common.h"

namespace gaudi {
namespace hepworth {

/// Writes a cache from the current `q`. It has no columns and no constraint rows.
class cache_compute {
public:
  using ptr = std::shared_ptr<cache_compute>;
  virtual ~cache_compute() = default;
  virtual void compute(vecX &q) = 0;
};

} // namespace hepworth
} // namespace gaudi

#endif
