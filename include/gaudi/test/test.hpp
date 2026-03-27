#ifndef __GAUDI_TEST_HPP__
#define __GAUDI_TEST_HPP__

#include "gaudi/console_logger.hpp"
#include <exception>
#include <string>
#include <vector>

namespace gaudi {
namespace test {

struct test_failure : public std::exception {};

struct TestContext {
  int failures = 0;

  void add_failure(const char *file, int line, const char *expr,
                   bool fatal) {
    failures++;
    console_logger::error << "TEST FAIL: " << file << ":" << line << " ("
                          << expr << ")" << std::endl;
    if (fatal) {
      throw test_failure();
    }
  }
};

inline thread_local TestContext *current_context = nullptr;

inline void set_current_context(TestContext *ctx) { current_context = ctx; }

inline void record_failure(const char *file, int line, const char *expr,
                           bool fatal) {
  if (current_context) {
    current_context->add_failure(file, line, expr, fatal);
  } else {
    console_logger::error << "TEST FAIL: " << file << ":" << line << " ("
                          << expr << ")" << std::endl;
    if (fatal) {
      throw test_failure();
    }
  }
}

struct TestCase {
  const char *name;
  void (*fn)();
};

struct Result {
  int total = 0;
  int passed = 0;
  int failed = 0;
  int assertions_failed = 0;
};

class Registry {
public:
  static void add(const char *name, void (*fn)()) {
    cases().push_back({name, fn});
  }

  static Result run_all(bool stop_on_fail = false) {
    Result result;
    for (const auto &test_case : cases()) {
      result.total++;
      TestContext ctx;
      set_current_context(&ctx);
      try {
        console_logger::info << "TEST START: " << test_case.name << std::endl;
        test_case.fn();
      } catch (const test_failure &) {
        // Failure already recorded.
      } catch (const std::exception &e) {
        ctx.add_failure(__FILE__, __LINE__, e.what(), false);
      } catch (...) {
        ctx.add_failure(__FILE__, __LINE__, "unknown exception", false);
      }

      result.assertions_failed += ctx.failures;
      if (ctx.failures == 0) {
        result.passed++;
        console_logger::info << "TEST PASS: " << test_case.name << std::endl;
      } else {
        result.failed++;
        console_logger::error << "TEST FAIL: " << test_case.name << std::endl;
        if (stop_on_fail) {
          break;
        }
      }
    }
    set_current_context(nullptr);
    return result;
  }

private:
  static std::vector<TestCase> &cases() {
    static std::vector<TestCase> test_cases;
    return test_cases;
  }
};

} // namespace test
} // namespace gaudi

#define GAUDI_TEST(name)                                                      \
  static void name();                                                         \
  namespace {                                                                 \
  struct gaudi_test_registrar_##name {                                        \
    gaudi_test_registrar_##name() { gaudi::test::Registry::add(#name, &name); } \
  };                                                                          \
  static gaudi_test_registrar_##name gaudi_test_registrar_instance_##name;    \
  }                                                                           \
  static void name()

// Variadic macros to handle commas in template arguments
// e.g., GAUDI_EXPECT(foo<A, B>::value == 2) works correctly
#define GAUDI_ASSERT(...)                                                     \
  do {                                                                        \
    if (!(__VA_ARGS__)) {                                                     \
      gaudi::test::record_failure(__FILE__, __LINE__, #__VA_ARGS__, true);    \
    }                                                                         \
  } while (0)

#define GAUDI_EXPECT(...)                                                     \
  do {                                                                        \
    if (!(__VA_ARGS__)) {                                                     \
      gaudi::test::record_failure(__FILE__, __LINE__, #__VA_ARGS__, false);   \
    }                                                                         \
  } while (0)

#endif // __GAUDI_TEST_HPP__
