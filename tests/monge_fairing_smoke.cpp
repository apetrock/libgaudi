#include <Eigen/Geometry>

#include "gaudi/test/albers_fairing_tests.hpp"
#include "gaudi/test/cr_dirichlet_tests.hpp"
#include "gaudi/test/principal_curvature_tests.hpp"
#include "gaudi/test/gauss_newton_tests.hpp"
#include "gaudi/test/hopf_sweep_tests.hpp"
#include "gaudi/test/monge_fairing_tests.hpp"
#include "gaudi/test/shape_block_tests.hpp"
#include "gaudi/test/shape_match_tests.hpp"
#include "gaudi/test/test.hpp"

int main() {
  const auto result = gaudi::test::Registry::run_all();
  return result.failed == 0 ? 0 : 1;
}
