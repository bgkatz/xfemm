# Building xfemm (cfemm) on Windows

Notes for compiling the `cfemm` C++ solvers on Windows with MSVC. Stock xfemm
does not build cleanly under MSVC on 64-bit Windows; the changes below are
required.

Tested with CMake 4.3 and the Visual Studio 2026 BuildTools (MSVC v145, toolset
14.51). The same recipe should work on any recent MSVC (2019/2022) — adjust the
`-G` argument to match your installed Visual Studio.

## Configure & build

From a regular PowerShell or cmd prompt (no `vcvars64.bat` needed — CMake's
Visual Studio generator handles it):

```pwsh
cd cfemm
mkdir build
cd build
cmake .. -G "Visual Studio 18 2026" -A x64 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build . --config Release
```

- `-G "Visual Studio 18 2026"` — use whatever generator matches your VS install
  (`Visual Studio 17 2022`, `Visual Studio 16 2019`, etc.).
- `-A x64` — pick an explicit architecture; the default for some generators is
  Win32, which compounds the pointer-truncation problem fixed below.
- `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` — the project's `cmake_minimum_required`
  predates CMake 4, which otherwise errors out.
- `--config Release` — essential for solver performance (`Debug` builds at
  `-Od`, which makes FEM solves 10–50× slower).

Binaries land in `cfemm/bin/Release/`.

## Verifying the build

```pwsh
cd cfemm/build
ctest -C Release
```

All 33 enabled tests should pass.

## Required source changes

### 1. GCC warning flags break MSVC
`cfemm/CMakeLists.txt` unconditionally appends `-Wall -Wextra -Wpedantic` to
`CMAKE_CXX_FLAGS`. MSVC interprets `-Wextra` as `/Wextra`, which doesn't exist,
and aborts with `error D8021: invalid numeric argument '/Wextra'`. Guard the
flags for non-MSVC compilers:

```cmake
if(MSVC)
    set(GCC_WARNING_FLAGS "")
else()
    set(GCC_WARNING_FLAGS "-Wall -Wextra -Wpedantic")
endif()
```

### 2. Triangle pointer truncation on Win64
`cfemm/fmesher/triangle/triangle.c` encodes a 2-bit orientation in the low bits
of mesh-element pointers using `unsigned long` casts, and stores intermediate
pointer arithmetic in `unsigned long` locals (e.g. `alignptr`). On Windows
(LLP64) `unsigned long` is 32-bit but pointers are 64-bit, so every round-trip
silently truncates pointers to their low 32 bits. The fmesher then dereferences
garbage and segfaults — which cascades into every downstream solver and
postprocessor test.

Two-step fix in `triangle.c`:

```c
// Add to the include block near the top:
#include <stdint.h>
```

Then replace `unsigned long` with `uintptr_t` everywhere it's used for pointer
arithmetic. Concretely, in the file as shipped:

- Every `(unsigned long)` cast → `(uintptr_t)` (52 occurrences, all in the
  decode/encode mesh-primitive macros and the debug-print helpers).
- Every `unsigned long alignptr;` local → `uintptr_t alignptr;` (7 occurrences,
  in the pool allocator's alignment math).

Leave the `unsigned long randomseed;` global alone — that one is genuinely a
counter, not a pointer.

This change is a no-op on Linux and macOS: both are LP64, so `unsigned long`
and `uintptr_t` are the same 64-bit type and the compiler emits identical code.

## Notes on the macOS-specific issues
The portability issues called out in [`BUILD_MACOS.md`](BUILD_MACOS.md) —
`sincos`, `<malloc.h>`, `std::not1` / `std::ptr_fun` — do not apply on Windows:

- MSVC ships `<malloc.h>` natively.
- The `sincos` branch in `femmcomplex.cpp` is gated behind `_GNU_SOURCE` /
  `__GLIBC__` and falls through to plain `sin`/`cos` on MSVC.
- The default MSVC C++ standard (C++14, set by the project) still has the
  deprecated `std::ptr_fun` / `std::not1`. If you bump the project to C++17+
  you'll need the lambda-based rewrite shown in the macOS guide.

## Warnings you can safely ignore
The Release build emits many MSVC warnings, all harmless on a clean checkout:

- `C4244` lossy conversions in Triangle (e.g. `uintptr_t` → `unsigned long`
  in `randomseed` updates, `double` → `int` in geometric predicates). Triangle
  intentionally truncates here.
- `C4996` `strcpy`/`fopen` deprecation notices from the CRT. Adding
  `_CRT_SECURE_NO_WARNINGS` silences them if desired.
- `C4267`/`C4305` size_t/double narrowing across the libfemm code.
