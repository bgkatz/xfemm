# Building xfemm (cfemm) on macOS

Notes for compiling the `cfemm` C++ solvers on macOS (tested on Apple Silicon,
Apple Clang). Stock xfemm does not build cleanly under modern Clang/macOS; the
changes below are required.

## Configure & build

```sh
cd cfemm
cmake . -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release
make
```

- `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` — the project's `cmake_minimum_required`
  predates CMake 4, which otherwise errors out.
- `-DCMAKE_BUILD_TYPE=Release` — **essential**. With no build type the project
  compiles at `-O0`, making the FEM solvers 10–50× slower (a single solve can
  take minutes instead of a fraction of a second).

Binaries land in `cfemm/bin/`.

> `make clean` deletes committed test-input mesh files
> (`fsolver/test/Temp.{node,ele,edge,pbc}`); restore with `git checkout --` if
> tests later fail with "Error copying file … No such file or directory".

## Required source changes

### 1. `sincos` is not available on macOS
`cfemm/libfemm/liblua/femmcomplex.cpp`, in `exp(const CComplex&)`. `sincos` is a
GNU extension; macOS provides `__sincos`. Make it portable:

```cpp
#if defined(__APPLE__)
    __sincos(x.im, &sin_x, &cos_x);
#elif defined(_GNU_SOURCE) || defined(__GLIBC__)
    sincos(x.im, &sin_x, &cos_x);
#else
    sin_x = sin(x.im);
    cos_x = cos(x.im);
#endif
```

### 2. `<malloc.h>` does not exist on macOS
macOS has no `<malloc.h>` (the allocation routines live in `<stdlib.h>`). Guard
every include. Two forms appear in the tree:

For `#include "malloc.h"` (in `libfemm/cuthill.cpp`, `libfemm/fullmatrix.cpp`):
```cpp
#ifdef __APPLE__
#include <stdlib.h>
#else
#include "malloc.h"
#endif
```

For `#include <malloc.h>` — just exclude it on Apple:
```cpp
#ifndef __APPLE__
#include <malloc.h>
#endif
```

Files needing the `<malloc.h>` guard:
`hsolver/hsolver.cpp`, `fmesher/writepoly.cpp`, `esolver/esolver.cpp`,
`fsolver/fsolver.cpp`, `fsolver/static2d.cpp`, `fsolver/staticaxi.cpp`,
`fsolver/harmonic2d.cpp`, `fsolver/harmonicaxi.cpp`.

### 3. `std::not1` / `std::ptr_fun` removed in C++17
`cfemm/libfemm/stringTools.h`, in `ltrim()` and `rtrim()`. Replace the removed
adaptors with lambdas:

```cpp
// ltrim:
s.erase(s.begin(), std::find_if(s.begin(), s.end(),
        [](int ch) { return !std::isspace(ch); }));

// rtrim:
s.erase(std::find_if(s.rbegin(), s.rend(),
        [](int ch) { return !std::isspace(ch); }).base(), s.end());
```

## Recommended (not strictly required to compile)

### Disable FP contraction for the bundled Triangle
`cfemm/fmesher/CMakeLists.txt`, on the builtin `triangle` target:

```cmake
if (NOT MSVC)
    target_compile_options(triangle PRIVATE -ffp-contract=off)
endif()
```

Triangle uses Shewchuk's adaptive-precision geometric predicates, which require
each multiply and add to be rounded separately. Clang contracts multiply-adds
into FMAs by default on arm64, which can silently break the predicates.

## Unrelated runtime fix
Compiling successfully is not sufficient to mesh air-gap-element (motor) models
on Apple Silicon — those hang due to a separate NaN bug. See
[`AIRGAP_MESH_FIX.md`](AIRGAP_MESH_FIX.md).

## Verifying the build
```sh
cd cfemm
ctest                 # full suite
ctest -E TorqueBenchmark   # skip the two heavier AGE benchmarks
```
All tests should pass (33/33 with the air-gap fix applied).
