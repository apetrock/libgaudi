#include "gaudi/test/weingarten_tests.hpp"
#include "gaudi/test/modal_corrector_tests.hpp"
#include "gaudi/test/material_direction_tests.hpp"
#include "gaudi/test/stripes_tests.hpp"
#include "gaudi/test/test.hpp"

int main() {
  const auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
