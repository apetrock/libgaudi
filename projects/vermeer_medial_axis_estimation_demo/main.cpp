#include "gaudi/duchamp/medial_axis_estimation_adapter.hpp"
#include "gaudi/vermeer/medial_duchamp_host.hpp"

int main() {
  auto demo = gaudi::duchamp::medial_axis_estimation_adapter::create();
  gaudi::vermeer::medial_duchamp_host host(std::move(demo));
  return host.run();
}
