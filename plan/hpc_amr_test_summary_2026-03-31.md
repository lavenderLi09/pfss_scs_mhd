# HPC AMR Test Summary

Date: 2026-03-31

## Scope

This note summarizes the AMR-related HPC tests for
`/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test`.

There are two layers:

- Historical baseline tests: 3 groups
- Current large-mesh debug campaign: 6 groups

Total counted here: 9 groups

## Historical Baselines

### 1. `d2_50`

- Purpose: earliest failing AMR baseline
- Main setup:
  - active regrid
  - no later `jb_eta`-based AMR strengthening
- Outcome:
  - NaN started around `it ~ 16`
- Meaning:
  - established the original failure timescale

### 2. `rgfreeze`

- Purpose: isolate whether regrid itself is the trigger
- Main setup:
  - effectively froze AMR/regrid
- Outcome:
  - run passed the early danger region
- Meaning:
  - strongly suggested the initial condition and basic time stepping were not the primary problem
  - regrid/AMR became the main suspect

### 3. `jb_50`

- Purpose: test the `jb_eta` AMR route
- Main setup:
  - introduced `jb_eta_max + jb_eta_mean` driven refinement logic
- Outcome:
  - earlier smaller-grid case could pass 50 steps
- Meaning:
  - `jb_eta` was a useful direction
  - but this did not prove stability for the newer large mesh

## Current Large-Mesh Campaign

Common large-mesh parameters used across most tests:

- `domain_nx = 400, 96, 128`
- `block_nx = 20, 12, 8`
- `refine_max_level = 2`
- `refine_threshold = 0.18d0, 0.24d0, 18*0.18d0`
- `derefine_ratio = 20*0.08d0`
- `xprobmin2 = 0.005d0`
- `xprobmax2 = 0.495d0`
- `fix_small_values = .false.`

### 1. `data_large/regrid4`

- Main changes:
  - active AMR
  - `ditregrid = 4`
  - `mod_usr` used the early "globalized `jb_eta` + theta-adjacent no-coarsen" version
- Key files:
  - `/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test/data_large/regrid4/amrvac.par`
  - `/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test/data_large/regrid4/mod_usr.t`
- Outcome:
  - could run substantially longer than the earliest failures
  - later entered NaN
- Meaning:
  - improvement was real
  - but AMR-triggered bad states still existed

### 2. `data_large/regrid8`

- Main changes:
  - same large mesh
  - real run used `ditregrid = 8`
  - actual run file:
    - `/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test/data_large/regrid8/amrvac_amr.par`
- Important note:
  - the directory also contains an `amrvac.par` with `ditregrid = 4`, but that was not the real submitted parameter file
- Outcome:
  - NaN still appeared
- Meaning:
  - the problem was not just "regrid too frequent"
  - some regridded mesh/state pattern itself remained unsafe

### 3. `data_large/amr_theta`

- Main changes:
  - switched to `inteldebug_nantrap`
  - tested stronger theta-cut AMR protection
- Outcome:
  - did not reach normal evolution
  - Intel runtime first reported `error (63)` in output formatting
  - later also exposed `error (75)` floating-point failure
- Meaning:
  - this run was useful mainly as a trap-based diagnostic stage
  - it helped expose an initialization-time AMR bug

### 4. `data_large/amr_theta1`

- Main changes:
  - fixed the `special_refine_grid` bug where `get_current(...)` results were used over `ixI^S`
  - restricted derived AMR diagnostics to `ixO^S`:
    - `bmag`
    - `current_mag`
    - `delta_min`
    - `jb_eta`
  - kept stronger theta protection
- Outcome:
  - initialization-phase floating-point failure disappeared
  - run progressed into normal evolution
  - later still developed NaN
- Meaning:
  - confirmed the initialization AMR bug was real and fixed
  - remaining problem moved to evolution/reconstruction stage

### 5. `data_large/amr_debug`

- Main changes:
  - switched to `debug` arch
  - real submitted parameter file was the root:
    - `/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test/amrvac_amr.par`
  - `ditregrid = 8`
  - used the stronger theta-freeze logic plus the `ixO^S` AMR fix
- Outcome:
  - no need to wait for full-field NaN
  - `debug` caught the earlier failure directly in
    `mhd_get_cbounds()` when evaluating `sqrt(wLp(rho_))`
  - `bipo_test.log` stayed finite through the last written step
- Meaning:
  - the first visible failure in this configuration is a bad reconstructed interface density
  - this is earlier and more useful than the later silent NaN pattern

### 6. `data_large/amr_debug1`

- Main changes:
  - still uses `debug` arch
  - server root `mod_usr.t` was patched into a user-side probe version only
  - no AMRVAC core source was modified
- Added user-side probes:
  - after `specialbound_usr`, check ghost-strip state for bad `rho / pth / mom / mag`
  - at the beginning of `special_refine_grid`, check incoming block state
  - print `AMRDBG` information for:
    - `jb_eta_max`
    - `jb_eta_mean`
    - theta protection distance
    - `force_inner`
    - `force_sheet`
    - `refine`
    - `coarsen`
- Backup made before patching:
  - `/share/home/guoyang/liyihua/outer_corona_relaxation/bipolar_hpc_test/mod_usr.t.bak_20260331_usrdbg`
- Current status:
  - in progress / probe stage
  - early output exists, but this test is being used mainly to catch the first upstream corruption layer

## What Actually Changed Across the Campaign

The important modifications across these tests fall into three categories.

### A. Regrid cadence

- `ditregrid = 4`
- `ditregrid = 8`
- effectively frozen regrid

### B. Theta-cut protection strength

- theta-adjacent blocks only prevented from coarsening
- later stronger local freeze near theta-cut:
  - `refine = -1`
  - `coarsen = -1`

### C. `special_refine_grid` AMR diagnostic implementation

- older bug:
  - `get_current(...)` was followed by derived calculations over `ixI^S`
- fixed version:
  - restricted those derived calculations to `ixO^S`

## Current Best Interpretation

At this stage the evidence supports the following picture:

- the earliest initialization-time AMR bug in `special_refine_grid` has been fixed
- stronger theta-cut protection helps, but does not fully eliminate the failure
- under `debug`, the earliest currently observed evolution-time failure is not "field already NaN everywhere"
- instead, a reconstructed interface state reaches invalid density first, and this appears in `wLp(rho_)` before the usual later NaN cascade

## Most Useful Current Probe

The most useful current test is `data_large/amr_debug1`, because it combines:

- `debug` arch
- the fixed `ixO^S` AMR diagnostic logic
- stronger theta protection
- user-side probes in `mod_usr.t` only

This is the cleanest setup for deciding whether the next bad state first appears:

- in theta-cut ghost filling
- in AMR-input block state before regridding
- or only later inside finite-volume reconstruction
