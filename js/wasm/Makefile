# Makefile for building WebAssembly "Hello World" module
# Phase 1: Basic WASM integration before Gaudi geometry_logger

# Emscripten compiler
CXX = emcc

# Source files
SOURCES = hello.cpp

# Output files
OUTPUT_JS = hello.js
OUTPUT_WASM = hello.wasm

# Compiler flags
CXXFLAGS = -std=c++17 \
		   -O2 \
		   -s WASM=1 \
		   -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
		   -s EXPORT_ES6=1 \
		   -s MODULARIZE=1 \
		   -s EXPORT_NAME="HelloModule" \
		   -s ALLOW_MEMORY_GROWTH=1 \
		   -s ENVIRONMENT='web,node' \
		   --bind

# Default target
all: $(OUTPUT_JS)

# Build the WASM module
$(OUTPUT_JS): $(SOURCES)
	$(CXX) $(SOURCES) -o $(OUTPUT_JS) $(CXXFLAGS)

# Clean build artifacts
clean:
	rm -f $(OUTPUT_JS) $(OUTPUT_WASM)

# Development build with debug info
debug: CXXFLAGS += -g -s ASSERTIONS=1 -s SAFE_HEAP=1
debug: $(OUTPUT_JS)

# Check if emscripten is available
check-emcc:
	@which emcc > /dev/null || (echo "Error: emcc not found. Please install Emscripten SDK." && exit 1)
	@echo "Emscripten found: $$(emcc --version | head -1)"

# Install instructions
install-emscripten:
	@echo "To install Emscripten SDK:"
	@echo "1. git clone https://github.com/emscripten-core/emsdk.git"
	@echo "2. cd emsdk"
	@echo "3. ./emsdk install latest"
	@echo "4. ./emsdk activate latest"
	@echo "5. source ./emsdk_env.sh"

.PHONY: all clean debug check-emcc install-emscripten
