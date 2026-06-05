# Scoping: a planar (linear) Air Gap Element for xfemm

## 1. Why

Linear-motor simulations in this toolbox re-mesh the entire model at every
rotor position (`use_sliding_band = not full_motor and not spec.is_linear` in
`femm/sweep.py`; the per-position `setup_motor_problem` redraws the geometry).
Each position therefore solves on a *different* unstructured triangulation, and
the Maxwell stress-tensor thrust — a local, derivative quantity — picks up the
resulting per-step variation as **spurious force ripple**.

This was characterized empirically (see the thrust-ripple convergence and
contour-averaging studies). Key findings:

- A **sinusoidal-field motor (polar-anisotropic ring)**, which physically has
  ~0 thrust ripple, still shows **28%** ripple in a 1-pole cell, converging
  toward ~2% only at a 16-pole cell. The ripple is therefore almost entirely
  **numerical**.
- Refining the **airgap element size** 4× barely moved it (26→23%): not a
  resolution problem.
- Averaging the force over several **airgap contours** roughly **halves** it,
  then **plateaus** (26→11→11% at 1/4/16 contours), and the plateau **scales
  with cell size** (≈11% at 1 pole, ≈5% at 2 poles).

Interpretation: there are two noise sources. (1) Single-contour sampling jitter
— removed by contour averaging (now the default, `n_contours=4`). (2) A
**residual floor from the per-step re-meshing perturbing the field solution
itself** — *not* removable by any post-processing, and the only known cure is to
**stop re-meshing**: share a fixed airgap mesh across rotor positions and let
the two sides slide. That is exactly what FEMM's rotary Air Gap Element (AGE)
does, and what we need a **planar/periodic analogue** of.

Rotary sector sims already enjoy this (the "draw once, modify band angle + currents
per step" sliding-band fast path). Linear has no equivalent because FEMM's AGE
is natively circular.

## 2. Background: how the rotary AGE works (FEMM 4.2 / xfemm)

The AGE replaces the meshed annular air gap between two concentric arcs with a
**closed-form coupling**. Laplace's equation in the gap (no current, no
material) has, per angular harmonic `n`, solutions `A_n(r,θ) = (a_n r^n +
b_n r^-n) e^{i n θ}`. The field on the two arc boundaries is expanded in these
harmonics; matching the analytic gap field to the FE boundary DOFs yields a
**dense coupling block** between the inner-arc and outer-arc DOFs that is added
to the global FE matrix. Relative rotation by angle `δ` enters purely as a
**per-harmonic phase shift `e^{i n δ}`** in that coupling — so the rotor angle
is changed by editing one boundary property, with **no re-meshing**.

In the xfemm/FEMM API it is a boundary type assigned to the two airgap arcs
(this toolbox already uses it: `airgap_band` BC, set via `mi_addboundprop`,
rotor angle updated with `mi_modifyboundprop('airgap_band', 10, angle)`). The
relevant solver pieces in cfemm are the air-gap-element assembly and the
harmonic machinery (e.g. `FemmProblem::getCircle`/`GetCircle`, the
`CAirGapElement` data, and the AGE contribution in the magnetics solver
assembly — `cfemm/fsolver/`). The `femmcli_*TorqueBenchmark*` tests exercise it.

## 3. The planar analogue

A linear motor's air gap is a straight rectangular strip of height `g` between
two parallel lines (stator face at `y=0`, rotor face at `y=-g`), periodic in
`x` over the simulated cell of length `L = sim_slots * slot_pitch`.

The math is the rotary AGE with the angular coordinate replaced by the periodic
linear coordinate:

| rotary AGE | planar AGE |
|---|---|
| annulus between radii `r_i`,`r_o` | strip between `y=0`, `y=-g` |
| angular harmonics `n` over `2π` | spatial harmonics `k_n = 2πn/L` |
| `r^{±n}` radial dependence | `e^{±k_n y}` transverse dependence |
| rotation `δ` → phase `e^{i n δ}` | translation `d` → phase `e^{i k_n d}` |
| full circle (periodic) | periodic **or anti-periodic** cell |

Per harmonic, `A_n(x,y) = (a_n e^{k_n y} + b_n e^{-k_n y}) e^{i k_n x}`. The
coupling between the two boundary lines is built exactly as in the rotary case,
and the rotor displacement `d` enters as the phase `e^{i k_n d}`. **No
re-meshing**; displacement is a single parameter update — the linear version of
the rotary fast path.

Anti-periodic cells (odd pole count, our `BdryFormat=5` case) use the
half-period harmonics `k_n = (2n+1)π/L` instead of `2πn/L`. This is the one
genuinely new wrinkle vs. the rotary AGE (which is always fully periodic over
`2π`); it must be threaded through the harmonic set.

## 4. Implementation plan (cfemm)

1. **Locate the rotary AGE assembly** in `cfemm/` (libfemm + fsolver): the code
   that, given the two arc boundaries and the rotor angle, builds the harmonic
   coupling block and adds it to the global matrix. Identify the harmonic count,
   how boundary DOFs are projected onto harmonics, and where rotation phase is
   applied.
2. **Add a planar AGE variant** parameterized by `(g, L, periodic|antiperiodic,
   d)`:
   - harmonic set `k_n` (periodic vs anti-periodic),
   - transverse functions `e^{±k_n y}` (replacing `r^{±n}`),
   - boundary-DOF→harmonic projection along `x` (replacing along `θ`),
   - displacement phase `e^{i k_n d}`.
   Most of the assembly structure (dense coupling block, matrix insertion)
   should be reusable; only the basis functions and the harmonic set change.
3. **Mesher (`cfemm/fmesher`)**: the two airgap boundary segments must exist as
   selectable edges with the AGE boundary assigned. As in the rotary AGE, the
   two sides are coupled by harmonic projection, so **matching node counts are
   not required** — but confirm the projection handles unequal/shifted node
   distributions on the two straight edges. The bulk meshing on either side is
   unchanged.
4. **Boundary/parameter plumbing**: a new boundary type (or extend the existing
   AGE type with a `planar` flag) carrying `g`, `L`, anti-periodic flag, and the
   displacement `d`; `d` must be updatable via `mi_modifyboundprop` so a sweep
   stays on the no-re-mesh fast path.
5. **Force/torque extraction**: the rotary AGE provides a clean gap-integral
   torque from the harmonic coefficients. Provide the planar analogue —
   **thrust directly from the AGE harmonics** (`F_x = (depth/μ0) Σ_n …`), which
   sidesteps the stress-tensor/contour machinery entirely and should be ripple-
   free by construction.

## 5. API / Lua surface

- Define once: `mi_addboundprop('airgap_band', …, <planar AGE format>, g, L,
  antiperiodic)` and assign to the two airgap segments (mirrors today's rotary
  `airgap_band`).
- Per step: `mi_modifyboundprop('airgap_band', <disp index>, d)` then
  `mi_analyze()` — no redraw, no re-mesh.
- Extraction: `mo_gapintegral('airgap_band', 0)` returns thrust (as the rotary
  AGE returns torque today; `femm/simulate.py:extract_force` already has the
  sector/`mo_gapintegral` branch to hook into).

## 6. Toolbox integration (this repo)

Once xfemm exposes it, the changes here are small and mirror the rotary path:

- `femm/problem_setup.py:_setup_linear_motor` / `_draw_linear_domain_periodic`:
  replace the sliding-rotor-domain + bridge construction with two airgap edges
  carrying the planar AGE BC; stop translating the rotor geometry per step.
- `femm/sweep.py`: drop `and not spec.is_linear` from `use_sliding_band` so
  linear gets the draw-once / modify-band fast path (currently rotary-only).
- `femm/simulate.py:extract_force`: add an AGE branch using
  `mo_gapintegral('airgap_band', …)` (the contour-averaging path stays as the
  fallback for the non-AGE / finite-length case).

## 7. Validation

- **Ring magnet → ~0 ripple at 1 pole.** The decisive test: with the AGE, a
  1-pole sinusoidal-field cell should give flat thrust (today it's 28%). 
- **Cell-size independence.** 1-pole, 2-pole, 4-pole cells should agree on mean
  *and* ripple (today ripple shrinks with cell size — the artifact).
- **Cross-check vs. the converged re-meshed result** (large cell, contour-
  averaged) for mean thrust and back-EMF.
- **Reuse the existing harness:** `examples/plot_thrust_vs_position.py` already
  sweeps cells and reports ripple; rerun it post-AGE.

## 8. Effort & risks

- **Effort:** moderate cfemm work. The hard math/assembly already exists for the
  rotary AGE; the planar version is largely a change of basis functions plus the
  anti-periodic harmonic set. Estimate the bulk of the effort in (a) replicating
  the AGE assembly with the planar basis and (b) the anti-periodic case.
- **Risks / open questions to resolve against the source:**
  - How the existing AGE projects boundary DOFs onto harmonics (does it assume a
    closed `2π` domain? the planar/periodic and anti-periodic cases need a
    different harmonic set and normalization).
  - Whether the mesher needs the two airgap edges pre-split or the projection
    tolerates arbitrary node layouts.
  - Convergence vs. number of retained harmonics (rotary AGE truncates; pick a
    default and expose it).
  - Anti-periodic sign bookkeeping at the cell boundary (the `BdryFormat=5`
    half-period harmonics).

## 9. Interim mitigation (already in place)

Until the planar AGE lands, `extract_force` averages 4 airgap contours by
default (one batched extraction), which removes the *sampling* half of the
ripple. The remaining floor is the re-mesh; for ripple-sensitive work, simulate
a larger cell (means are correct at every cell size — only ripple is
under-resolved).
