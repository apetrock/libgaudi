#include "gaudi/duchamp/slab_prism_bending_adapter.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::slab_prism_bending_adapter::create());
}
