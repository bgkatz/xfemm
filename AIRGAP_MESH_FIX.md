# Fix: xfemm hangs meshing air-gap (AGE) models on Apple Silicon

## Symptom
On arm64 macOS, `mi_createmesh()` hangs (effectively forever) for any model
using an **Air Gap Element (AGE)** — e.g. the `femmcli_TorqueBenchmark` tests
and any rotating-machine torque model. Non-AGE models mesh fine. FEMM (Windows)
meshes the same model in seconds.

## Root cause
`FemmProblem::getCircle()` computes an arc's circle center as:

```cpp
R = d / (2*sin(tta/2));            // d = chord length, tta = swept angle
c = a0 + (d/2 + I*sqrt(R*R - d*d/4)) * t;
```

For a **~180° arc**, `R == d/2`, so `R*R - d*d/4` is mathematically `0`.
Under **FMA contraction** (clang's default on arm64) the compiler fuses
`R*R - d*d/4` into `fma(R, R, -(d*d/4))`: `R*R` is computed to infinite
precision while `d*d/4` is separately rounded, so the two no longer cancel and
the radicand becomes a tiny **negative** (~`-4e-17`). `sqrt(negative)` → `NaN`.

The NaN center propagates to every discretized point on that arc (~140 `NaN`
node coordinates), which are handed to Triangle. NaN poisons every
incircle/orient predicate (all NaN comparisons are false), sending Triangle's
divide-and-conquer Delaunay into a pathological, effectively non-terminating
incircle blowup — the "hang".

**Why FEMM is unaffected:** FEMM is built with MSVC (`/fp:precise`, no FMA
contraction), so the radicand evaluates to exactly `0`. This is a
compiler/FMA difference, **not** an ARM-vs-x86 `sin()` difference (`sin(M_PI/2)`
is exactly `1.0` on both).

## The fix
Clamp the radicand at `0` before `sqrt()`. The true value there is `0`, so this
is locally exact and robust regardless of compiler, FMA, or `sin()` rounding.

### 1. `cfemm/libfemm/FemmProblem.cpp` — `getCircle()` (~line 1550)
```cpp
double radicand = R*R - d*d / 4.;
if (radicand < 0) radicand = 0;
c = a0 + (d/2. + I * sqrt(radicand)) * t;
```

### 2. `cfemm/libfemm/FemmProblem.cpp` — `GetCircle()` (~line 2334)
Same clamp (duplicate function used by the AGE path):
```cpp
double radicand = R*R-d*d/4.;
if (radicand < 0) radicand = 0;
c = a0 + (d/2. + I*sqrt(radicand))*t;
```

## Secondary fix (separate real bug, found nearby)
`cfemm/fmesher/writepoly.cpp` — AGE discretization (~line 1727).
The port of FEMM's `calloc` array to `std::vector` used `reserve()` (capacity
only; `size()` stays 0) and then wrote via `operator[]` — out-of-bounds writes.
Changed to `resize()`:
```cpp
agelst[n]->nodeNums.resize(myVector.size()+1);   // was reserve()
```
Not the cause of the hang, but undefined behaviour that should be fixed.

## Optional hardening (not required by the fix above)
`cfemm/fmesher/CMakeLists.txt` (~line 25) builds the bundled Triangle with
`-ffp-contract=off`, so Shewchuk's robust geometric predicates aren't broken by
FMA contraction. Defensive correctness only; the clamp is what resolves the AGE
hang.

## Verification
- `femmcli_TorqueBenchmark` and `femmcli_antiperiodicBC_AGE_TorqueBenchmark`:
  previously hung >11 min, now pass in <1 s with all torque values within
  tolerance.
- Full suite: `ctest` → 33/33 pass.
- Confirm input is clean: `bin/fmesher --write-poly model.fem` then
  `grep -c nan model.raw.poly` → `0`.

## Build note (macOS)
Configure with: `cmake . -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release`
(`Release` is essential — the default empty build type compiles at `-O0`).
