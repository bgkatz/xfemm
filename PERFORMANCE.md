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
| `motor` | 26k | 51k | real outrunner motor model (`50x14_custom.fem`), nonlinear (13-pt + 47-pt BH), 12 Newton iters |

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

Status: implemented, benchmarked, committed.

---

## 2. Port CSR conversion to the harmonic solver — PLANNED

`cfemm/libfemm/cspars.cpp` (`CBigComplexLinProb`, used for AC problems) has
the same linked-list structure, and stores up to four matrices (`M`, `Mh`,
`Ma`, `Ms`) for Newton AC problems. Mechanical port of item 1. Items 3/4
apply to `cspars.cpp` too. Needs an AC benchmark model (e.g. a frequency≠0
variant of one of the existing models).

## 3. Single-walk `AddTo` — IMPLEMENTED (2026-08-27, with item 4)

**Problem:** `spars.cpp`: `AddTo(v,p,q)` was `Put(Get(p,q)+v,p,q)` — two
linked-list row walks per assembly insertion (~9 AddTo calls per element per
Newton iteration, 55 per air-gap quad element).

**Change:** single combined walk — add in place if the entry exists, insert
if absent. Bit-identical arithmetic.

## 4. Column adjacency for `SetValue`/`Periodicity`/`AntiPeriodicity` — IMPLEMENTED (2026-08-27)

**Problem:** `(Anti)Periodicity` carried a `#define KLUDGE` forcing `bdw=0`,
so each periodic node pair scanned **all n rows** with two `Get()` row walks
each, re-applied every Newton iteration. (The kludge existed because
periodicity itself creates entries outside the a-priori bandwidth, which the
original banded scan would then miss on later calls.) Cost: NumPBCs ×
NumNodes row walks per Newton step — motor: 158 pairs × 26k rows × 12 iters,
tq_dense/age_dense: 720 pairs × 148k rows. `SetValue` (Dirichlet BCs) had a
banded variant of the same scan.

**Change:** added a lazily-built column adjacency (`colRows`: for each
column, the (row, entry-pointer) pairs above the diagonal). Once built it is
maintained incrementally by the insertion paths (`Put`/`AddTo`), so it never
goes stale; entry pointers are stable because entries are never deleted.
`SetValue` now touches only the O(degree) structural entries of the node's
column/row; `(Anti)Periodicity` visits only rows holding an entry in either
paired column (collected, sorted, deduped — ~16 rows instead of n). The
KLUDGE and the banded-scan logic are gone.

Files: `cfemm/libfemm/spars.h`, `cfemm/libfemm/spars.cpp`.

**Results (items 3+4 together, end-to-end fsolver process time):**

| Model | Baseline | CSR (item 1) | + items 3/4 | vs CSR | vs baseline |
|---|---|---|---|---|---|
| temp | 0.333 s | 0.279 s | 0.211 s | 1.32x | 1.58x |
| tq | 0.085 s | 0.076 s | 0.066 s | 1.15x | 1.29x |
| age | 0.083 s | 0.077 s | 0.063 s | 1.22x | 1.32x |
| temp_dense | 3.29 s | 2.28 s | 1.73 s | 1.32x | **1.90x** |
| tq_dense | 5.09 s | 3.23 s | 2.37 s | 1.36x | **2.14x** |
| age_dense | 4.90 s | 3.25 s | 2.36 s | 1.38x | **2.07x** |
| motor | 3.15 s | 1.95 s | 1.66 s | 1.18x | **1.90x** |

Correctness: solution `.ans` files are bit-identical to the item-1 build on
all seven models (rel-L2 diffs vs. original baseline unchanged to the last
digit); full ctest suite 33/33 passed. Note: `femmcli_fpproc.lua` now passes
and failed at baseline — this test is environmental/flaky (it fails/passes
consistently within a build tree state, its bad values look like stale
data, and `makemask.cpp` doesn't use any of the changed code paths);
tracked as unrelated.

## 5. Fast mesh-file parsing — IMPLEMENTED (2026-08-27)

**Problem:** `fsolver.cpp LoadMesh()` parsed .node/.ele/.edge with per-token
`fscanf` (4 calls/node, 5/element, 4/edge; MSVC fscanf locks the stream and
re-parses the format string per call), and `cuthill.cpp` parsed the .edge
file **twice** the same way. Measured phase breakdown before the fix
(tq_dense, 148k nodes / 294k elements / 441k edges): loadmesh 0.36 s,
cuthill 0.47 s — together ~35% of the remaining runtime.

**Change:** new header-only `cfemm/libfemm/FileTokenizer.h` — reads the
whole file into memory and hands out whitespace-separated tokens via
`strtol`/`strtod`. Used for the .node/.ele/.edge loops in `LoadMesh` (the
small .pbc/AGE section keeps its line-oriented parser) and in `Cuthill()`,
which now parses .edge once into an in-memory edge list instead of two
fscanf passes. Token-stream semantics identical to the fscanf sequence it
replaces (with added truncated-file error checks). `Cuthill()` is a
FEASolver template method, so esolver/hsolver inherit the fix.

Files: `cfemm/libfemm/FileTokenizer.h` (new), `cfemm/fsolver/fsolver.cpp`,
`cfemm/libfemm/cuthill.cpp`.

**Results:** phase times on tq_dense: loadmesh 0.36 s → 0.09 s (3.8x),
cuthill 0.47 s → 0.22 s (2.1x; the remainder is the renumbering algorithm
itself, not I/O). End-to-end:

| Model | Baseline | + items 1/3/4 | + item 5 | vs baseline |
|---|---|---|---|---|
| temp | 0.333 s | 0.211 s | 0.165 s | 2.02x |
| tq | 0.085 s | 0.066 s | 0.049 s | 1.73x |
| age | 0.083 s | 0.063 s | 0.049 s | 1.69x |
| motor | 3.15 s | 1.66 s | 1.52 s | **2.07x** |
| temp_dense | 3.29 s | 1.73 s | 1.43 s | **2.30x** |
| tq_dense | 5.09 s | 2.37 s | 1.82 s | **2.80x** |
| age_dense | 4.90 s | 2.36 s | 1.83 s | **2.68x** |

Correctness: `.ans` solution diffs vs. original baseline unchanged on all
seven models (bit-identical to the item-1 build). ctest: no new failures;
`femmcli_fpproc.lua` confirmed genuinely flaky — with a single fixed binary
it passes 3/4 and fails 1/4 consecutive runs (nondeterministic test or
uninitialized read in fpproc, pre-existing; solver outputs are verified
deterministic).

Remaining end-to-end breakdown after item 5 (tq_dense): assemble+solve
~1.36 s (mostly PCG), cuthill algorithm ~0.22 s, .ans write ~0.15 s,
loadmesh ~0.09 s.

## 6. Build system: no default CMAKE_BUILD_TYPE — PLANNED

Top-level `cfemm/CMakeLists.txt` sets no optimization flags and no default
build type: single-config generators (Linux/macOS Makefiles, Ninja, MinGW)
produce -O0 binaries on a plain `cmake .. && make` — easily 5-15x slower.
Fix: default to Release when unset on single-config generators. (MSVC
multi-config builds are unaffected.)

## 7. Eigen LDLT direct solver — PROTOTYPED (2026-08-27)

**Change:** vendored Eigen 3.4.0 under `cfemm/external/eigen-3.4.0`
(auto-detected by `libfemm/CMakeLists.txt`, which then defines
`XFEMM_HAVE_EIGEN`; builds without it are unchanged). `CBigLinProb` gains
`SolveDirect()`: the symmetric upper-triangle CSR arrays from item 1 map
directly to Eigen's column-major lower triangle (zero-transform), feeding
`SimplicialLDLT`. Runtime opt-in via environment variable:
`XFEMM_DIRECT=1` = factorize+solve every call (symbolic analysis reused);
`XFEMM_DIRECT=2` = factorize once, then CG preconditioned with the stale
factorization, refactorizing on slow convergence; unset = existing PCG.
Falls back to PCG automatically if factorization fails.

**Results:** mode 1 wins everywhere — factorization is so cheap in 2D
(8–9 ms at 26k nodes, ~0.4 s at 148k) that mode 2's stale-preconditioner CG
is never worth it at these sizes.

Solve-phase comparison (all Newton steps summed):

| Model | PCG (post item 1) | LDLT mode 1 | Solve speedup |
|---|---|---|---|
| motor (12 solves) | 1.07 s | ~0.11 s | **~10x** |
| tq_dense (1 solve, 525 iters) | 1.27 s | 0.42 s | 3.0x |

End-to-end (best of 3), XFEMM_DIRECT=1:

| Model | Baseline | After items 1/3/4/5 (PCG) | LDLT | vs baseline |
|---|---|---|---|---|
| temp | 0.333 s | 0.165 s | 0.113 s | 2.9x |
| motor | 3.15 s | 1.52 s | **0.61 s** | **5.1x** |
| temp_dense | 3.29 s | 1.43 s | 1.03 s | 3.2x |
| tq_dense | 5.09 s | 1.82 s | 1.03 s | **4.9x** |
| age_dense | 4.90 s | 1.83 s | 1.04 s | 4.7x |

Correctness: linear models agree with baseline to ~1e-11 (direct solve is
exact; the diff is the baseline PCG's own convergence error). Nonlinear
models differ by 2e-7 (motor) to 4e-6 (temp) relative — expected: an exact
linear solve changes the Newton trajectory within the nonlinear stopping
tolerance (res < 100×Precision); Newton iteration counts are unchanged
(motor: 12). Default mode (env unset) is bit-identical to the item-5 build;
full ctest 33/33.

**Made the default (2026-08-27):** the direct solver is now used by default
whenever the build finds Eigen; `XFEMM_DIRECT=0` forces the legacy PCG
solver (bit-identical to the pre-Eigen build), `XFEMM_DIRECT=2` selects the
stale-factorization hybrid. Applies to every `CBigLinProb` user (fsolver,
esolver, hsolver, fpproc mask). Automatic PCG fallback if a factorization
fails. Full ctest 33/33 with direct as default.

**Licensing:** Eigen 3.4.0 is MPL2. The historically-LGPL sparse Cholesky
code (Tim Davis's LDL) was relicensed to MPL2 via an agreement with Google
(see the header of `Eigen/src/SparseCholesky/SimplicialCholesky_impl.h`),
so our entire code path is MPL2. The build defines `EIGEN_MPL2_ONLY`,
which turns any accidental inclusion of remaining LGPL-licensed files into
a compile error — a machine-checked guarantee. MPL2 is file-level
copyleft: shipping the unmodified source tree with its COPYING.* files
(kept in `cfemm/external/eigen-3.4.0/`) satisfies it, and it composes fine
with xfemm's Aladdin/FEMM-derived licensing and downstream commercial use.
The vendored tree is trimmed to headers + licenses (6.8 MB).

Remaining decision: vendored copy vs. git submodule vs. find_package with
graceful fallback (currently: vendored, auto-detected, optional). Mode-2
hybrid unmeasured on 500k+ meshes.

## 7c. LDLT silent inaccuracy on indefinite Newton Jacobians — BUG FIXED (2026-08-27)

Found via spurious Monte-Carlo results in the downstream motor-design
toolbox (jagged flux waveforms, phases not 120° apart, 8x torque ripple on
occasional designs). Root cause: `SimplicialLDLT` does no pivoting and has no
convergence criterion; on some designs the nonlinear Newton Jacobian goes
(numerically) indefinite once B-H updates kick in (cubic B-H spline
overshoot can produce locally negative differential permeability), and the
factorization then reports Success while losing ~11 digits (measured
relative residual ~3e-3 vs ~4e-14 on healthy matrices; iteration 1 with
linear init is always fine, degradation starts at iteration 2). PCG never
had this failure mode because its convergence criterion IS a residual
check.

Fix in `SolveDirect`: solve into scratch (Z) so the Newton warm start in V
survives a failed attempt; verify an explicit relative residual
||b−Ax||/||b|| against the solver Precision; run up to 5 passes of
iterative refinement (one mat-vec + one triangular solve each) for
merely-ill-conditioned cases; fall back to PCG if refinement stalls
(genuinely indefinite). Cost: one mat-vec per healthy solve (~2 ms at
148k nodes). Hostile designs run at pure-PCG speed (measured identical,
1.91 s vs 1.91 s on the repro design) — no speedup there, but correct.
Verified: the failing design's sweep now matches PCG exactly; healthy
models unchanged; ctest 33/33.

Ops note: any femmcli/fsolver binaries deployed elsewhere must be rebuilt
to pick this up, and **MSBuild's up-to-date check twice failed to relink
executables after libfemm changes** in this session — verify with
`grep -c <new string> <exe>` or delete `bin/Release/*.exe` and rebuild.
MC results produced between enabling the direct solver and this fix are
suspect for a minority of designs (garbage is obvious: huge torque
ripple / non-sinusoidal flux); re-simulate or filter those.

## 7b. Uninitialized `MuMax` in CMMaterialProp — BUG FIXED (2026-08-27)

Found while making the direct solver the default: `femmcli_fpproc.lua`
flipped from flaky to consistently failing. Root cause (pre-existing,
latent): `CMMaterialProp::MuMax` was never initialized in the default
constructor, and the copy constructor failed to copy `MuMax`, `mu_fdx`,
`mu_fdy`, and `Frequency` (all left uninitialized — the copy ctor bypasses
the default init list). `GetMu()` uses `MuMax>0` as the
incremental-permeability flag, so for linear materials (air, magnets —
`GetSlopes` never runs for them) the flag read heap garbage. When the
garbage was positive, point values of H and Mu in linear regions were
silently corrupted for DC problems (`muinc = mu_x/B` → observed Mu1 ≈
9.7e12 at B ≈ 1e-13 T). Heap-layout dependence explains the historical
flakiness; Eigen's allocations merely made the garbage deterministic.
Fixed by initializing `MuMax(0.)` and copying all four members. The fpproc
test now passes 6/6 consecutive runs.

Files: `cfemm/libfemm/CMaterialProp.cpp`.

## 8. Parallelism — IDEAS / NOT SCOPED

- OpenMP `MultA`/`Dot` now feasible on the CSR arrays (needs full-row storage
  or per-thread accumulators to handle the symmetric scatter).
- SSOR preconditioner is inherently sequential; threaded alternative is
  (block-)Jacobi, which trades iterations for cores — benchmark, don't assume.
- Eigen factorization could also be swapped for CHOLMOD/Pardiso (threaded)
  if the direct path becomes default and larger models demand it.

## 9. Minor items — NOTED

- Hoist `cos(t*PI/180)`/`sin(t*PI/180)` out of the per-edge magnetization
  loop (`static2d.cpp:610`, computed 3x per element per Newton iteration).
- `GetFillFactor` is O(labels × elements) (AC problems only).
- `PCGSolve` has no max-iteration cap (robustness; the complex solver has
  MAXITER).
- B-H `GetH`/`GetdHdB` do a linear scan of the curve per element per Newton
  iteration — small tables, low priority.
## 10. Downstream-toolbox work log

- 2026-08-27: tier-1 extraction batching implemented in the motor-design
  toolbox that drives xfemm (backends/base.py `mo_extract_batch`,
  LuaBackend single-script override, batched `extract_flux_linkage`,
  `extract_torque_and_flux` used by simulate_motor + sweep sliding-band
  path). 6 → 3 femmcli processes/design, 1.41 → 1.08 s serial,
  MC 2m47s → 1m30s.
- 2026-08-28: production validation — 200k-design overnight MC (6s/7p)
  at a sustained 540 motors/min with zero observed spurious results,
  confirming the §7c residual guard holds at scale. Baseline throughput
  before this effort was ~154 motors/min (6m30s / 1000 designs) on a
  lighter config.

## 11. Remaining time budget & next targets (as of 2026-08-27)

Context: a downstream Monte-Carlo motor-design pipeline (1000 designs)
went 6m30s → 1m30s via items 1–7 plus tier-1 extraction batching in the
toolbox (mo_extract_batch: 6 → 3 femmcli processes/design, results
digit-identical, toolbox tests 163/163). Serial per-design profile after
batching (~1.1 s + worker spawn; profiler: `femm_temp/profile_backend.py`
in the toolbox repo):

| Phase | Time | Share |
|---|---|---|
| femmcli analyze (mesh + solve) | 0.66 s | ~55% |
| Python (geometry calc + Lua emission) | 0.29 s | ~24% |
| extract_batch process | 0.10 s | ~8% |
| setup process (draw + save .fem) | 0.03 s | ~3% |
| worker python spawn (outside profile) | ~0.15–0.3 s | — |

Candidate fixes, est. serial savings per design:

- **Newton reassembly (xfemm)** — analyze's ~0.66 s is roughly 0.12 s
  meshing + ~0.35 s of 12× full matrix reassembly + ~0.11 s LDLT solves +
  I/O. The solver wipes and reassembles every element each Newton
  iteration, but only nonlinear-material element contributions change:
  cache the linear part (split b and M into fixed + nonlinear-delta),
  reassemble only nonlinear elements per iteration. Est. −0.15 to −0.25 s.
  Moderate effort, touches static2d/staticaxi assembly loops. A cheaper
  partial: hoist the per-element sin/cos and Lua-magdir work (item 9).
- **Tier-2 single femmcli process per design (toolbox)** — append
  mi_analyze + mi_loadsolution + extraction prints to the setup script:
  3 processes → 1, and the .fem/.ans re-parses disappear. Est. −0.10 to
  −0.15 s serial, likely more under parallel contention (tier 1
  over-delivered for exactly that reason).
- **Persistent workers (toolbox)** — runner spawns `python -m ..._worker`
  per design; a worker that loops over designs amortizes interpreter+numpy
  startup. Est. −0.15 to −0.25 s per design.
- **Python emission profile (toolbox)** — 0.29 s to generate a 2457-line
  Lua string looks high; cProfile it (suspects: per-call string formatting,
  geometry recalc, YAML/spec copies). Unknown until measured; maybe half.
- **fpproc .ans loader (xfemm)** — port FileTokenizer (item 5) to fpproc's
  fgets/sscanf solution parser; shrinks every extraction/femmcli
  mi_loadsolution. Mostly subsumed by tier 2, worthwhile independently for
  interactive post-processing and the losses sweep (36 steps × loads).
- **Mesh density (model-side)** — meshing is ~0.12 s and solve cost scales
  with node count; if MC fidelity allows a coarser mesh in the far field,
  everything downstream shrinks. Zero code.
- **/arch:AVX2 for MSVC Release (xfemm)** — cheap global flag; maybe
  5–15% on solve kernels. Verify no behavior change.

Realistic stack (reassembly + tier 2 + persistent workers): ~1.1 s →
~0.5–0.65 s serial per design, i.e. MC ≈ 45–60 s. Diminishing returns
beyond that without threading the solve or batching designs per process.

