#include "gaudi/duchamp/growth_study_adapter.hpp"
#include "gaudi/vermeer/medial_duchamp_host.hpp"

int main() {
  auto demo = gaudi::duchamp::growth_study_adapter::create();
  gaudi::vermeer::medial_duchamp_host host(std::move(demo));
  return host.run();
}
