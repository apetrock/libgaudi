#include <emscripten/bind.h>
#include <string>

using namespace emscripten;

int add(int a, int b) {
    return a + b;
}

std::string greet(std::string name) {
    return "Hello, " + name + "!";
}

EMSCRIPTEN_BINDINGS(bar_demo) {
    function("add", &add);
    function("greet", &greet);
}
