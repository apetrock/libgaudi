#ifndef __DUCHAMP_TEST_BASE__
#define __DUCHAMP_TEST_BASE__

#include "gaudi/define_create_func.hpp"

namespace gaudi {
namespace duchamp {

class test_base {
public:
  DEFINE_CREATE_FUNC(test_base)
  test_base(){};
  virtual ~test_base(){};
};

} // namespace duchamp
} // namespace gaudi
#endif