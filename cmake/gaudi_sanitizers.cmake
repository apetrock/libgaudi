# cmake/gaudi_sanitizers.cmake
# Compiler-aware sanitizer flags for ASan, TSan, LSan, MSan, UBSan.
# Include after project() so CMAKE_CXX_COMPILER_ID is available.

# ---------------------------------------------------------------------------
# AddressSanitizer  (heap overflow, use-after-free, use-after-scope, ...)
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # MSVC cl.exe — /fsanitize=address requires VS 2019 16.9+.
    # /RTC* and /INCREMENTAL are incompatible with ASan.
    # /Z7: debug info in .obj (avoids C1041 on shared vc*.pdb; see cmake/gaudi_msvc_z7.cmake)
    set(CMAKE_C_FLAGS_ASAN   "/fsanitize=address /Z7 /Od"
        CACHE STRING "C flags for AddressSanitizer (MSVC)."   FORCE)
    set(CMAKE_CXX_FLAGS_ASAN "/fsanitize=address /Z7 /Od"
        CACHE STRING "C++ flags for AddressSanitizer (MSVC)." FORCE)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_C_FLAGS_ASAN
        "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O0"
        CACHE STRING "C flags for AddressSanitizer (Clang/GCC)."   FORCE)
    set(CMAKE_CXX_FLAGS_ASAN
        "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g -O0"
        CACHE STRING "C++ flags for AddressSanitizer (Clang/GCC)." FORCE)
endif()

if(MSVC)
    # MSVC-compatible linker (cl.exe and clang-cl both invoke link.exe)
    set(CMAKE_EXE_LINKER_FLAGS_ASAN    "/DEBUG /INCREMENTAL:NO"
        CACHE STRING "Exe linker flags for ASan (MSVC)."    FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_ASAN "/DEBUG /INCREMENTAL:NO"
        CACHE STRING "Shared linker flags for ASan (MSVC)." FORCE)
    set(CMAKE_MODULE_LINKER_FLAGS_ASAN "/DEBUG /INCREMENTAL:NO"
        CACHE STRING "Module linker flags for ASan (MSVC)." FORCE)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_EXE_LINKER_FLAGS_ASAN    "-fsanitize=address"
        CACHE STRING "Exe linker flags for ASan (Clang/GCC)."    FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_ASAN "-fsanitize=address"
        CACHE STRING "Shared linker flags for ASan (Clang/GCC)." FORCE)
    set(CMAKE_MODULE_LINKER_FLAGS_ASAN "-fsanitize=address"
        CACHE STRING "Module linker flags for ASan (Clang/GCC)." FORCE)
endif()

# ---------------------------------------------------------------------------
# ThreadSanitizer  (Clang/GCC only; no MSVC support)
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_C_FLAGS_TSAN   "-fsanitize=thread -g -O1"
        CACHE STRING "C flags for ThreadSanitizer."   FORCE)
    set(CMAKE_CXX_FLAGS_TSAN "-fsanitize=thread -g -O1"
        CACHE STRING "C++ flags for ThreadSanitizer." FORCE)
    set(CMAKE_EXE_LINKER_FLAGS_TSAN    "-fsanitize=thread"
        CACHE STRING "Exe linker flags for TSan."     FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_TSAN "-fsanitize=thread"
        CACHE STRING "Shared linker flags for TSan."  FORCE)
endif()

# ---------------------------------------------------------------------------
# LeakSanitizer  (Clang/GCC, Linux primarily)
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_C_FLAGS_LSAN   "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
        CACHE STRING "C flags for LeakSanitizer."   FORCE)
    set(CMAKE_CXX_FLAGS_LSAN "-fsanitize=leak -fno-omit-frame-pointer -g -O1"
        CACHE STRING "C++ flags for LeakSanitizer." FORCE)
    set(CMAKE_EXE_LINKER_FLAGS_LSAN    "-fsanitize=leak"
        CACHE STRING "Exe linker flags for LSan."   FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_LSAN "-fsanitize=leak"
        CACHE STRING "Shared linker flags for LSan." FORCE)
endif()

# ---------------------------------------------------------------------------
# MemorySanitizer  (Clang only, Linux only in practice)
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_C_FLAGS_MSAN
        "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
        CACHE STRING "C flags for MemorySanitizer."   FORCE)
    set(CMAKE_CXX_FLAGS_MSAN
        "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer -g -O2"
        CACHE STRING "C++ flags for MemorySanitizer." FORCE)
    set(CMAKE_EXE_LINKER_FLAGS_MSAN    "-fsanitize=memory"
        CACHE STRING "Exe linker flags for MSan."     FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_MSAN "-fsanitize=memory"
        CACHE STRING "Shared linker flags for MSan."  FORCE)
endif()

# ---------------------------------------------------------------------------
# UndefinedBehaviorSanitizer
# ---------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    set(CMAKE_C_FLAGS_UBSAN   "-fsanitize=undefined -g"
        CACHE STRING "C flags for UBSan."   FORCE)
    set(CMAKE_CXX_FLAGS_UBSAN "-fsanitize=undefined -g"
        CACHE STRING "C++ flags for UBSan." FORCE)
    set(CMAKE_EXE_LINKER_FLAGS_UBSAN    "-fsanitize=undefined"
        CACHE STRING "Exe linker flags for UBSan."    FORCE)
    set(CMAKE_SHARED_LINKER_FLAGS_UBSAN "-fsanitize=undefined"
        CACHE STRING "Shared linker flags for UBSan." FORCE)
endif()

# ---------------------------------------------------------------------------
# Multi-config generator support  (Visual Studio, Xcode, Ninja Multi-Config)
# ---------------------------------------------------------------------------
if(CMAKE_CONFIGURATION_TYPES)
    list(APPEND CMAKE_CONFIGURATION_TYPES Asan)
    list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}"
        CACHE STRING "Available build configurations" FORCE)
endif()

# ---------------------------------------------------------------------------
# Symbolizer / runtime hints  (printed at configure time)
# ---------------------------------------------------------------------------
string(TOUPPER "${CMAKE_BUILD_TYPE}" _GAUDI_BT_UPPER)

set(_GAUDI_ASAN_ACTIVE FALSE)
if(_GAUDI_BT_UPPER STREQUAL "ASAN")
    set(_GAUDI_ASAN_ACTIVE TRUE)
elseif(CMAKE_CONFIGURATION_TYPES AND "Asan" IN_LIST CMAKE_CONFIGURATION_TYPES)
    set(_GAUDI_ASAN_ACTIVE TRUE)
endif()

if(_GAUDI_ASAN_ACTIVE)
    find_program(_GAUDI_LLVM_SYMBOLIZER llvm-symbolizer
        HINTS "C:/Program Files/LLVM/bin"
              "/usr/bin" "/usr/local/bin")
    if(_GAUDI_LLVM_SYMBOLIZER)
        message(STATUS "ASan: llvm-symbolizer found at ${_GAUDI_LLVM_SYMBOLIZER}")
    else()
        message(WARNING
            "ASan: llvm-symbolizer not found — stack traces will show raw addresses.\n"
            "  Install LLVM or set ASAN_SYMBOLIZER_PATH at runtime.")
    endif()
    message(STATUS "ASan: Useful runtime environment variables:")
    message(STATUS "  ASAN_OPTIONS=symbolize=1:detect_stack_use_after_return=1:halt_on_error=0")
    message(STATUS "  ASAN_SYMBOLIZER_PATH=<path-to-llvm-symbolizer>")
endif()

unset(_GAUDI_BT_UPPER)
unset(_GAUDI_ASAN_ACTIVE)
