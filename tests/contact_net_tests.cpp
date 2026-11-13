#include "gaudi/test/contact_net_tests.hpp"
#include "gaudi/test/test.hpp"

int main() {
  const auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
