// Compatibility shim for legacy headers that still include <cxxabi.h>.
//
// The current codebase does not appear to use any cxxabi symbols directly in
// the native viewer targets, but MSVC does not provide this GCC/libc++ header.
// Keeping a local stub lets us preserve portability while cleaning up stale
// includes incrementally.
