# Analytical Bipolar AMRVAC Test: Debug Process Report

## Scope

This document records the full test and debug process for the analytical bipolar AMRVAC case in:

- [analytic_bipolar_1to20_stretched](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched)

It is based on the full interaction history of this debugging session plus the local result summaries:

- [debug_summary.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly/debug_summary.txt)
- [physics_checks_summary.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly/physics_checks_summary.txt)

The goal was to construct a clean AMRVAC analytical bipolar-field test, initialize it on a spherical shell from `1` to `20 Rsun`, run lightweight local checks on a PC-safe configuration, and understand why the first evolution tests failed or produced suspicious behavior.

## Initial Goal

The original task was to construct a new AMRVAC test using:

- the instructions in [`agent.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/agent.md)
- the local AMRVAC-related code in [`amrvac_polytropic/polytropic_bipolar`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/polytropic_bipolar)
- AMRVAC version `amrvac-amrvac3.1`

The requested physics setup was:

- analytical bipolar magnetic field
- radial domain from `1` to `20 Rsun`
- Parker-like polytropic solar-wind background
- full 3D spherical geometry

The implementation was built as a new sibling test case rather than by modifying the downloaded reference case in place.

## New Case Construction

The new working directory became:

- [analytic_bipolar_1to20](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20)

and later the dedicated local stability/debug directory became:

- [analytic_bipolar_1to20_stretched](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched)

The basic structure followed the downloaded AMRVAC case template:

- `makefile`
- `amrvac.par`
- `mod_usr.t`
- `main.job`

The important simplifications relative to `polytropic_bipolar` were:

- keep only the analytical bipolar field
- remove the flux rope and electric driving logic
- preserve CT initialization, Parker wind, spherical shell geometry, and gravity

## Early Setup Corrections

Several setup issues were identified and corrected before the physics debugging started.

### 1. Domain resolution and base density

The downloaded `polytropic_bipolar` reference case used a base domain of:

- `72 x 72 x 72`

The new test kept that initial base resolution for first debugging passes.

The bottom density was reduced by a factor of 10 relative to the downloaded reference setup:

- from `5.0d9/100.0d0/unit_numberdensity`
- to `5.0d9/1000.0d0/unit_numberdensity`

### 2. Correct AMRVAC tree

The case was explicitly switched to build against:

- `/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1`

### 3. CPU-safe local testing

The user requested PC-safe runs, so all local evolution tests were kept to:

- one process
- `OMP_NUM_THREADS=1`

and only short smoke tests or initialization checks were run.

### 4. Full-sphere theta range issue

An early mistake was identified in the theta extent:

- the case originally used pole boundaries while theta did not properly span the full spherical interval

This was corrected so that the full-pole spherical test used:

- `xprobmin2=0.0d0`
- `xprobmax2=0.5d0`

which corresponds to `theta in [0, pi]` in this AMRVAC normalized setup.

## Radial Stretching Clarification

There was a brief misunderstanding about whether the radial mesh was already stretched.

The final clarification was:

- `stretch_dim(1)='uni'` in AMRVAC does **not** mean uniform spacing
- it means a stretched unidirectional mesh
- the run output confirmed this with messages such as `Stretched dimension 1`

So the radial mesh was already stretched and did not require further correction for that purpose.

## Bipole Burial Depth Discussion

After inspecting ParaView output, the user asked whether the bipole was unrealistically hanging above the surface.

Inspection of [`mod_usr.t`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20/mod_usr.t) showed:

- the magnetic charges were actually buried well beneath `r=1`
- the visual impression came from the geometry of the buried-charge field, not from sources being placed above the surface

The source depth was then made even deeper by changing the runtime controls:

- `f_d`
- `f_L`

This was done through the parameter files, not by changing the analytical field formula itself.

## Local Stability Ladder

To test whether the case could evolve on a local machine, a PC-safe staged run structure was introduced in:

- [analytic_bipolar_1to20_stretched](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched)

The staged parameter files included:

- [pc_init.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_init.par)
- [pc_smoke.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke.par)
- [pc_medium.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_medium.par)
- [pc_long.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_long.par)

The user explicitly requested that small-value fixing should **not** be enabled.

This led to an important clarification.

### `fix_small_values` clarification

The user correctly pointed out that the `.par` namelist logic had initially been misunderstood.

The corrected understanding was:

- each `&list ... /` block is self-contained
- parameters written after `/` are not part of that namelist

It was then verified from the AMRVAC source that the default behavior is already:

- `fix_small_values = .false.`
- `small_values_method = 'error'`

So the smoke tests were intentionally exposing true failures rather than masking them.

## First Physics Failure: Full-Pole Run

The first real evolution test was the full-pole smoke run:

- [pc_smoke.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke.par)

This test failed immediately.

### Observed behavior

- initialization completed
- the code started time stepping
- at `Iteration 0`, before completing the first full time step, it encountered negative gas pressure and aborted

The diagnostic files were:

- [pc_smoke.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_smoke.log)
- [pc_smoke.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_smoke.stdout)

### Why internal-energy evolution did not prevent this

It was confirmed from the AMRVAC source and docs that:

- `mhd_internal_e=.true.` avoids pressure recovery through total-energy subtraction
- but it does **not** guarantee positivity of pressure
- the evolved internal energy can still be driven negative by the numerical update and source terms

So the negative pressure in this case was not a total-energy reconstruction problem. It was a true internal-energy failure in the first step.

## Hypotheses for the First-Step Crash

Inspection of [`mod_usr.t`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t) suggested three immediate suspects.

### 1. Inner radial boundary mismatch

The initial state used Parker outflow, but the inner boundary had been setting radial momentum to zero.

This was a strong inconsistency:

- interior: outward Parker wind
- ghost cells: static inner boundary

### 2. Magnetic boundary inconsistency under CT

The initial magnetic field was analytical, but the inner radial boundary CT treatment was extrapolating face fields numerically and then reconstructing the centered field, rather than directly enforcing the same analytical state.

### 3. Polar-region sensitivity

Thetheta `case(3)` and `case(4)` were largely empty, so the run relied on AMRVAC’s generic full-pole treatment.

Given the analytical field geometry, the south-pole ring was the first place where any inconsistency was likely to blow up.

## First Boundary Fixes

Two local fixes were made to test whether the immediate crash was primarily due to the inner radial boundary.

### Fix A. Inner radial momentum made Parker-consistent

In the inner boundary branch of [`mod_usr.t`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t):

- `mom(1)` was changed from zero
- to `vout * rho`

### Fix B. Inner ghost-cell centered magnetic field reset to the analytical bipole

After the CT face update and `mhd_face_to_center`, the ghost-cell centered magnetic field was overwritten with the same analytical bipolar field used in the interior initialization.

### Result of these two fixes

The full-pole smoke test still failed immediately.

This was the key clue that the dominant instability source was not just the inner radial boundary.

## Off-Pole Diagnostic Run

To separate pole handling from the rest of the model, a trimmed theta test was created:

- [pc_smoke_offpole.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke_offpole.par)

with:

- `xprobmin2=0.005d0`
- `xprobmax2=0.495d0`

This intentionally removed the exact north and south poles from the domain.

### Important boundary interpretation

The user asked whether off-pole runs could still meaningfully use `pole` boundaries.

The final conclusion from AMRVAC source inspection was:

- no, true pole treatment is only activated if theta exactly begins at `0` and/or ends at `pi`
- the message `There is no northpole!` / `There is no southpole!` comes from geometry recognition, not from the boundary-type switch itself
- for off-pole runs, theta should be handled explicitly through `specialbound_usr`, not conceptually treated as a real pole boundary

## Off-Pole Theta Boundary: First Working Version

The first off-pole working version used explicit `special` theta boundaries and implemented `case(3)` and `case(4)` with a pragmatic mixed strategy:

- hydro quantities: zero-gradient in theta
- CT magnetic faces: extrapolated tangential components and a divergence-based reconstruction for the normal component
- centered magnetic field: overwritten with the analytical bipole

This was numerically stable.

### Numerical result

The run reached `it=50` without crashing.

### Physics interpretation

This was **not** a clean solar-wind MHD relaxation.

The main quantitative signs were:

- density changed only mildly
- magnetic-field changes were small
- strong transverse velocities appeared near the off-pole theta cuts
- some cells developed negative radial velocity near the inner boundary

This suggested:

- the off-pole boundary was numerically usable
- but it was still likely injecting non-physical transverse disturbances

## Off-Pole Theta Boundary: Analytical Continuation Attempt

To make the off-pole theta boundary more self-consistent, a second version was introduced while preserving the old version as comments.

The new logic added a helper routine:

- [`set_analytic_parker_bipole`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t#L338)

and changed theta `case(3)` / `case(4)` to:

- fill ghost cells with Parker wind quantities directly
- fill centered magnetic field with the analytical bipole directly
- map this state into the staggered magnetic arrays for CT and then recenter

The user specifically requested that the older working implementation should not be deleted, so the old `case(3)` and `case(4)` logic was commented out and preserved above the new implementation.

## Result of the Analytical Continuation Attempt

This second off-pole version also ran to `it=50`.

So in terms of pure stability, it was still acceptable.

However, the physics quality became worse.

### Main new symptom

Pressure changed dramatically:

- the old off-pole working version had `p_max ~ 8.98e-3`
- the analytical-continuation version had `p_max ~ 4.46e-1`

This was a very large pressure enhancement compared to the previous variant.

### Localization of the pressure anomaly

Detailed cell-location checks showed that the pressure enhancement was:

- not domain-wide
- not centered on the interior bipolar structure
- concentrated near the truncated theta boundaries
- located on the first radial shell above the inner boundary

The strongest cells were found at approximately:

- `r ~ 1.109688`
- `theta ~ 0.050639`
- `theta ~ 3.086728`

which are exactly the north and south theta cut locations.

The pressure maximum was therefore interpreted as a **boundary-generated pressure ring** rather than a genuine global MHD relaxation feature.

## Meaning of `debug_summary.txt`

The file:

- [debug_summary.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly/debug_summary.txt)

contains a minimal quick check:

- number of cells
- final maximum pressure

It is not a full narrative summary, but it is a useful marker that the final VTU could be read and that the pressure peak was:

- `pmax=0.4456518292427063`

which matches the later detailed diagnosis.

## Brief Comparison of the Main Test Results

### Test A. Full-pole smoke test

- configuration: full spherical theta with poles included
- result: failed immediately at `Iteration 0`
- main issue: negative gas pressure near the south-pole ring
- conclusion: unusable in current form

### Test B. Off-pole test with first working theta boundary

- configuration: theta trimmed away from poles; mixed extrapolation-based special theta boundary
- result: numerically stable through 50 steps
- main issue: strong `Vt`/`Vp` near theta cuts and some negative `Vr` near the inner boundary
- conclusion: good diagnostic configuration, but not yet physically clean

### Test C. Off-pole test with Parker+bipole analytical continuation theta boundary

- configuration: theta trimmed away from poles; ghost cells filled from a direct Parker+bipole continuation
- result: numerically stable through 50 steps
- main issue: much stronger pressure enhancement localized at the theta cuts and first radial shell
- conclusion: physically worse than Test B despite comparable stability

## Current Best Interpretation

The evidence collected so far suggests:

- the original immediate failure is strongly tied to full-pole handling for this configuration
- trimming away the poles stabilizes the time stepping
- but the off-pole theta boundary remains the dominant source of non-physical behavior
- simply making the ghost-cell centered state more analytical is not sufficient

The most likely unresolved issue is:

- the treatment of the **staggered CT magnetic faces** at the theta-cut boundaries

In other words, the main remaining problem is probably not the cell-centered hydro state alone, but the geometric consistency of the face-centered magnetic boundary state under CT.

## Practical Outcome

At this stage, the case can be used as:

- a local diagnostic off-pole shell test

but not yet as:

- a clean, physically self-consistent solar-wind MHD relaxation benchmark

The most defensible current working baseline is the earlier off-pole special-theta version, because:

- it is numerically stable
- its artifacts are milder
- it does not create the large pressure ring seen in the later analytical-continuation attempt

## Suggested Next Step

The next debugging step should focus on:

- reworking the theta-cut **CT face magnetic boundary**

rather than on:

- more long integrations of the current analytical-continuation version
- or further interpretation of the current pressure ring as physical relaxation

Until that is done, the abnormal pressure ring should be treated as a boundary artifact.
