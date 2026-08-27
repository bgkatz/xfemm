# xfemm Performance Optimization Log

Running log of solver/toolchain speed optimizations: what's planned, what's
implemented, and measured results. Started 2026-08-27.

Machine for all measurements: Windows 11, MSVC Release build (x64).
Benchmark harness: `bench/bench.py` (untracked). Each model's mesh files are
backed up in `bench/meshes/` because fsolver deletes them after loading;
the script restores them before every run. Timing = full fsolver process
(mesh load + assembly + solve + write .ans), best of 3 unless noted.

## Benchmark models

| Model | Nodes | Elements | Notes |
|---|---|---|---|
| `temp` | 15k | — | fsolver test model, nonlinear (9-pt BH curve), 3 Newton iters |
| `tq` | 4.9k | — | femm.info TorqueBenchmark, linear, periodic BC |
| `age` | 4.9k | — | same geometry, antiperiodic air-gap element |
| `temp_dense` | 83k | — | `temp` with mesh size / 3 |
| `tq_dense` | 148k | 295k | `tq` with mesh size / 8 |
| `age_dense` | 148k | 295k | `age` with mesh size / 8 |
| `motor` | 26k | 51k | Ben's `50x14_custom.fem` (motoropt), nonlinear (13-pt + 47-pt BH), 12 Newton iters |

Correctness check for every change: relative L2 difference of the `.ans`
solution vector vs. pre-change baseline, plus identical PCG iteration counts
where the change is not supposed to alter numerics. Full `ctest -C Release`
suite must show no new failures (`femmcli_fpproc.lua` fails pre-existing on
this machine — near-zero-field H/Mu edge case in the postprocessor test,
present at baseline commit 8ed1297).

---

## 1. CSR conversion of the magnetostatic sparse solver — IMPLEMENTED (2026-08-27)

**Problem:** `CBigLinProb` (`cfemm/libfemm/spars.cpp`) stores each matrix row
as a linked list of individually heap-allocated `CEntry` nodes. The PCG hot
kernels — `MultA` (mat-vec) and `MultPC` (SSOR preconditioner, two triangular
sweeps) — traversed these pointer chains ~3x per PCG iteration. Sparse
mat-vec is memory-bound; scattered heap nodes vs. contiguous arrays is the
dominant cost.

**Change:** added `FlattenMatrix()` which copies the linked-list matrix into
flat arrays (`csrDiag[]` for the diagonal; `csrRowStart[]`/`csrCol[]`/
`csrVal[]` CSR arrays for the strictly-upper-triangle entries). Called at the
top of every `PCGSolve()` (O(nnz), negligible) because the nonlinear solvers
modify matrix values between Newton iterations. `MultA`/`MultPC` rewritten to
iterate the flat arrays. Assembly path (`Put`/`AddTo`/`SetValue`/
`Periodicity`) unchanged — still linked lists. esolver, hsolver, and fpproc's
mask builder use the same class and get the speedup for free.

Files: `cfemm/libfemm/spars.h`, `cfemm/libfemm/spars.cpp`.

**Results:**

End-to-end fsolver process time:

| Model | Baseline | CSR | Speedup |
|---|---|---|---|
| temp | 0.333 s | 0.279 s | 1.19x |
| tq / age | 0.085 s | 0.077 s | 1.1x |
| temp_dense | 3.29 s | 2.31 s | 1.43x |
| tq_dense | 5.09–5.59 s | 3.26 s | 1.56–1.7x |
| age_dense | 4.90–5.42 s | 3.26 s | 1.5–1.66x |
| motor | 3.15 s | 1.95 s | **1.61x** |

PCG-kernel-only time (temporary instrumentation, since removed):

| Model | Baseline PCG | CSR PCG | Kernel speedup |
|---|---|---|---|
| tq_dense (525 iters) | 2.99 s | 1.27 s | 2.36x |
| age_dense (525 iters) | 2.95 s | 1.26 s | 2.34x |
| temp_dense (424+386+1 iters) | 2.07 s | 1.07 s | 1.95x |
| motor (12 solves, 432→2 iters) | 2.39 s | 1.07 s | 2.22x |

Correctness: PCG iteration counts identical per Newton step on all models;
solution rel-L2 diff vs. baseline 1e-14 to 9e-12 (floating-point
reassociation only, far below solver Precision).

Status: implemented, benchmarked, **not yet committed** (working tree on
`linear-planar-age`).

---

## 2. Port CSR conversion to the harmonic solver — PLANNED

`cfemm/libfemm/cspars.cpp` (`CBigComplexLinProb`, used for AC problems) has
the same linked-list structure, and stores up to four matrices (`M`, `Mh`,
`Ma`, `Ms`) for Newton AC problems. Mechanical port of item 1. Needs an AC
benchmark model (e.g. a frequency≠0 variant of one of the existing models).

## 3. Single-walk `AddTo` — PLANNED

`spars.cpp`: `AddTo(v,p,q)` is `Put(Get(p,q)+v,p,q)` — two linked-list row
walks per assembly insertion. One combined walk (add in place, insert if
absent) halves assembly search cost. Matters more now: on the nonlinear
`motor` model, the 12x reassembly + BC application is ~0.9 s of the remaining
1.95 s runtime. Same treatment applies to `cspars.cpp`.

## 4. `Periodicity`/`AntiPeriodicity` KLUDGE full scans — PLANNED

`spars.cpp:366-474` (pre-change line numbers): `#define KLUDGE` forces
`bdw=0`, so each periodic node pair scans all n rows with two `Get()` calls
per row, re-applied every Newton iteration. Cost is NumPBCs × NumNodes row
walks per Newton step — hits sliding-band / (anti)periodic motor models
specifically. Fix: restore the banded scan (understand what the kludge
papered over first) or maintain a column-adjacency structure. `SetValue`
(Dirichlet BCs) has a smaller variant of the same pattern.

## 5. Fast mesh-file parsing — PLANNED

`fsolver.cpp LoadMesh()`: per-token `fscanf` (4 calls/node, 5+/element; MSVC
fscanf is slow and locks the stream per call). At 148k nodes / 295k elements
this plus assembly is ~2 s of the 3.26 s remaining on `tq_dense` — now the
majority of end-to-end time for linear models. Fix: read whole file, parse
with `strtol`/`strtod`. Same pattern in `.ans` writing (per-value fprintf)
and in `cuthill.cpp` (reads the .edge file with fscanf, twice).

## 6. Build system: no default CMAKE_BUILD_TYPE — PLANNED

Top-level `cfemm/CMakeLists.txt` sets no optimization flags and no default
build type: single-config generators (Linux/macOS Makefiles, Ninja, MinGW)
produce -O0 binaries on a plain `cmake .. && make` — easily 5-15x slower.
Fix: default to Release when unset on single-config generators. (MSVC
multi-config builds are unaffected.)

## 7. Parallelism — IDEAS / NOT SCOPED

- OpenMP `MultA`/`Dot` now feasible on the CSR arrays (needs full-row storage
  or per-thread accumulators to handle the symmetric scatter).
- SSOR preconditioner is inherently sequential; threaded alternative is
  (block-)Jacobi, which trades iterations for cores — benchmark, don't assume.
- Alternative worth measuring: Eigen `SimplicialLDLT` direct solve; for
  <~500k nodes a factorization is often faster than PCG outright, and Newton
  steps could reuse the symbolic factorization.

## 8. Minor items — NOTED

- Hoist `cos(t*PI/180)`/`sin(t*PI/180)` out of the per-edge magnetization
  loop (`static2d.cpp:610`, computed 3x per element per Newton iteration).
- `GetFillFactor` is O(labels × elements) (AC problems only).
- `PCGSolve` has no max-iteration cap (robustness; the complex solver has
  MAXITER).
- B-H `GetH`/`GetdHdB` do a linear scan of the curve per element per Newton
  iteration — small tables, low priority.
