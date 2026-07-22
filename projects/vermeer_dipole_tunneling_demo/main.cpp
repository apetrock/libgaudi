#include "gaudi/duchamp/dipole_tunneling_demo_adapter.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(gaudi::duchamp::dipole_tunneling_demo_adapter::create());
}
