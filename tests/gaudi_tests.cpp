#include "gaudi/test/test.hpp"
#include "gaudi/test/view_tests.hpp"
#include "gaudi/test/simplex_algorithm_tests.hpp"
#include "gaudi/test/morton_simplex_tests.hpp"
#include "gaudi/test/bvh_tests.hpp"
#include "gaudi/test/pyramid_datum_tests.hpp"
#include "gaudi/test/calder_winding_tests.hpp"
#include "gaudi/calder/shell_area_conservation_test.hpp"
#include "gaudi/calder/rod_length_conservation_test.hpp"
#include "gaudi/test/rod_dynamic_tests.hpp"
#include "gaudi/test/shell_dynamic_tests.hpp"

int main() {
  auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
