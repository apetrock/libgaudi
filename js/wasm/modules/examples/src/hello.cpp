#include <emscripten/bind.h>
#include <iostream>

class HelloWorld {
public:
    HelloWorld() {
        std::cout << "Hello from WebAssembly!" << std::endl;
    }
    
    std::string getMessage() const {
        return "Hello from WebAssembly C++!";
    }
    
    int add(int a, int b) const {
        return a + b;
    }
};

EMSCRIPTEN_BINDINGS(hello_world) {
    emscripten::class_<HelloWorld>("HelloWorld")
        .constructor<>()
        .function("getMessage", &HelloWorld::getMessage)
        .function("add", &HelloWorld::add);
}
