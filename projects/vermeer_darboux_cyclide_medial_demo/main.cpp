#include "gaudi/duchamp/darboux_cyclide_medial_demo_adapter.hpp"
#include "gaudi/vermeer/medial_duchamp_host.hpp"

int main() {
  auto demo = gaudi::duchamp::darboux_cyclide_medial_demo_adapter::create();
  gaudi::vermeer::medial_duchamp_host host(std::move(demo));
  return host.run();
}
