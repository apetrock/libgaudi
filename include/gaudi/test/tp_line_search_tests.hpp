#ifndef GAUDI_TP_LINE_SEARCH_TESTS_HPP
#define GAUDI_TP_LINE_SEARCH_TESTS_HPP

#include "gaudi/duchamp/modules/tp_line_search.hpp"
#include "gaudi/test/test.hpp"

namespace gaudi {
namespace test {

/// Bowl centered at 0.05, between shrink samples 0.035 and 0.07.
/// Refine should land on the bowl, not on a halving-grid point.
GAUDI_TEST(tp_line_search_picks_lowest_energy) {
  duchamp::tp_line_search_config cfg;
  cfg.shrink = 0.5;
  cfg.max_backtracks = 12;
  cfg.max_refine = 8;
  const real E0 = 10.0;
  auto eval = [](real alpha) {
    const real d = alpha - real(0.05);
    return real(10.0) + d * d - real(0.02);
  };
  const auto ls = duchamp::tp_best_energy_line_search(E0, 0.28, eval, cfg);
  GAUDI_EXPECT(std::abs(ls.alpha - real(0.05)) < real(1e-3));
  GAUDI_EXPECT(ls.E_accepted < E0);
  GAUDI_EXPECT(ls.E_accepted < eval(real(0.035)));
  GAUDI_EXPECT(!ls.used_fallback);
}

GAUDI_TEST(tp_line_search_holds_when_energy_rises) {
  duchamp::tp_line_search_config cfg;
  cfg.shrink = 0.5;
  cfg.max_backtracks = 4;
  cfg.max_refine = 0;
  auto eval = [](real alpha) { return real(10.0) + alpha; };
  const auto ls = duchamp::tp_best_energy_line_search(10.0, 1.0, eval, cfg);
  GAUDI_EXPECT(ls.alpha == real(0.0));
  GAUDI_EXPECT(ls.used_fallback);
}

} // namespace test
} // namespace gaudi

#endif
