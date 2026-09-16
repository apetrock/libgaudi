#include "gaudi/duchamp/cubic_stylization_adapter.hpp"
#include "gaudi/vermeer/vermeer.hpp"

int main() {
  return gaudi::vermeer::vermeer(
      gaudi::duchamp::cubic_stylization_adapter::create());
}
