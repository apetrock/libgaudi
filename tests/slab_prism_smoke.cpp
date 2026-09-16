#include "gaudi/test/test.hpp"
#include "gaudi/test/slab_prism_tests.hpp"

int main() {
  auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
