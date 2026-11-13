#include "gaudi/test/amap_tests.hpp"
#include "gaudi/test/chiral_tests.hpp"
#include "gaudi/test/dupin_circle_tests.hpp"
#include "gaudi/test/dupin_line_tests.hpp"
#include "gaudi/test/growth_tests.hpp"
#include "gaudi/test/test.hpp"

int main() {
  const auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
