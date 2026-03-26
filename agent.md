# PFSS + SCS to AMRVAC: Project Context

## Purpose

This project aims to build a practical pipeline from a synoptic magnetogram to an AMRVAC-ready initial magnetic field for MHD relaxation in the low corona and inner heliosphere.

The intended end-to-end workflow is:

`magnetogram -> PFSS coefficients -> SCS coefficients -> field evaluation on the AMRVAC target mesh -> AMRVAC MHD relaxation`

This document serves two roles:

- a compact project brief for ongoing work
- a context file for future coding agents

It records the project direction and the debugging conclusions that have already been established. Open design choices are listed explicitly rather than hidden.

## Confirmed Technical Direction

### 1. The target is AMRVAC initialization, not a standalone extrapolation product

The PFSS and SCS solutions are not being developed as independent archive products. Their main purpose here is to prepare an initial magnetic field for AMRVAC.

That means implementation choices should be judged mainly by:

- consistency with the AMRVAC mesh
- numerical cleanliness of the initial condition
- usefulness for subsequent MHD relaxation

### 2. Direct evaluation on the target AMRVAC mesh is preferred

The current preferred workflow is to evaluate PFSS and SCS directly on the final AMRVAC mesh whenever possible.

This is preferred over:

`magnetogram -> PFSS/SCS on a uniform storage grid -> interpolation to AMRVAC`

for three reasons:

- it avoids one extra interpolation layer
- it preserves more of the benefit of the analytical or harmonic representation
- it reduces the chance of introducing extra divergence-related noise into the MHD initial condition

### 3. PFSS and SCS should share one target mesh

The intended split is:

- below the PFSS-to-SCS transition radius: evaluate the PFSS solution
- above the transition radius: evaluate the SCS solution

But both should be evaluated on the same final target mesh, especially the same stretched radial mesh used by AMRVAC.

The important idea is:

- one target mesh
- two radial regimes
- no separate archived grid products as the main path

### 4. SCS should support harmonic-style direct evaluation

The SCS implementation should not be limited to interpolating from a precomputed grid.

The desired capability is conceptually similar to harmonic PFSS evaluation:

- given SCS coefficients such as `glm/hlm`
- evaluate the magnetic field directly at arbitrary `(r, theta, phi)` locations

This file does not lock down a final API signature, but it does record the behavioral requirement:

- SCS direct pointwise evaluation is a project goal

### 5. Use the project-local modified PFSS/SCS code, not an older installed copy

If multiple PFSS-related packages exist on the system, project work must use the local project version rather than an outdated conda-installed copy.

This is important for:

- notebooks
- scripts
- future packaging
- debugging imports

## Current Repository Structure

The current working repository is:

- [pfss_scs_mhd](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd)

The most relevant top-level items are:

- [agent.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/agent.md)
  - this context file
- [knowledge_base.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/knowledge_base.md)
  - additional project notes
- [cursor_pfss_and_scs_model_grid_discussi.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/cursor_pfss_and_scs_model_grid_discussi.md)
  - related design discussion notes
- [references](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/references)
  - papers and background references
- [hmi_synoptic_maps](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/hmi_synoptic_maps)
  - synoptic map data and plotting helpers
- [2026_scs_relaxation](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation)
  - PFSS/SCS experiments, notebooks, stretched-grid products, and local PFSS code
- [amrvac_polytropic](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic)
  - AMRVAC-side test problems and debug cases

### `2026_scs_relaxation`

This directory contains the current PFSS/SCS experimentation layer.

Important contents include:

- [initial_magnetic.ipynb](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation/initial_magnetic.ipynb)
  - notebook for building and inspecting the initial magnetic configuration
- [codes/pfss](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation/codes/pfss)
  - project-local PFSS-related source tree
- [Brtp_data](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation/Brtp_data)
  - saved PFSS/SCS field samples
- [rtp_stretched.npz](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation/rtp_stretched.npz)
  - stretched spherical mesh data
- [vts](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/2026_scs_relaxation/vts)
  - VTK outputs for field inspection

This directory is the main PFSS/SCS construction side of the project.

### `amrvac_polytropic`

This directory contains AMRVAC test cases derived from or compared against a downloaded reference setup.

Important subdirectories are:

- [polytropic_bipolar](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/polytropic_bipolar)
  - downloaded reference AMRVAC case used to inspect coding patterns
- [analytic_bipolar_1to20](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20)
  - first clean analytical bipolar test case from `1` to `20 Rsun`
- [analytic_bipolar_1to20_stretched](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched)
  - main active AMRVAC debug case with stretched radial grid and staged local tests

### `analytic_bipolar_1to20_stretched`

This is the current AMRVAC debugging focus.

Important files:

- [mod_usr.t](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t)
  - user physics, initialization, boundaries, and auxiliary output
- [amrvac.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/amrvac.par)
  - main default run file
- [pc_init.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_init.par)
  - initialization-only check
- [pc_smoke.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke.par)
  - short full-pole smoke test
- [pc_smoke_offpole.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke_offpole.par)
  - short off-pole smoke test
- [run_pc_stability.sh](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/run_pc_stability.sh)
  - single-core local runner
- [analysis/pressure_anomaly](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly)
  - summaries and debug records for the pressure anomaly investigation

## Current AMRVAC Debugging Conclusions

### 1. The radial grid is already stretched

The AMRVAC setup already uses a stretched radial grid.

In this project:

- `stretch_dim(1)='uni'` means stretched radial spacing, not uniform radial spacing

This has already been verified from AMRVAC behavior and run output.

### 2. Full-pole runs are currently unstable

The full-pole smoke test:

- [pc_smoke.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke.par)

fails immediately at the first evolution step.

Observed behavior:

- initialization succeeds
- at `Iteration 0`, before completing the first full step, negative gas pressure appears
- the failure is concentrated near the south-pole ring

So the current full-pole configuration is not yet usable for physical runs.

### 3. Off-pole runs are numerically stable but still boundary-sensitive

The off-pole test:

- [pc_smoke_offpole.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/pc_smoke_offpole.par)

trims theta away from the exact poles and uses `special` theta boundaries.

This version can run stably for at least 50 steps on one core.

However, stability does not yet imply physical cleanliness.

### 4. Off-pole theta boundaries are the current main difficulty

Two off-pole theta boundary variants have been tested:

- older working variant
  - zero-gradient hydro
  - extrapolated CT magnetic field
  - centered analytical magnetic overwrite
- newer analytical-continuation variant
  - Parker-wind continuation in ghost cells
  - analytical bipolar continuation in ghost cells

The newer variant did **not** solve the physical-quality problem. Instead, it produced a stronger pressure anomaly.

### 5. Current best interpretation of the pressure anomaly

The pressure anomaly investigation is summarized in:

- [physics_checks_summary.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly/physics_checks_summary.txt)
- [debug_process_report.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly/debug_process_report.md)

The main conclusion is:

- the large pressure increase is localized near the truncated theta boundaries and the first radial shell above the inner boundary
- it should currently be treated as a boundary artifact
- the most likely unresolved issue is the CT staggered magnetic-face treatment at the theta-cut boundaries

In short:

- full-pole case: crashes
- off-pole case: runs
- off-pole physics: still not clean enough to call it a trustworthy solar-wind relaxation baseline

## Practical Guidance For Future Agents

If future work continues from this repository, the default interpretation should be:

- this is not just a PFSS/SCS extrapolation project
- it is a PFSS/SCS-to-AMRVAC initialization project
- direct evaluation on the final mesh remains the preferred direction
- the current active AMRVAC debugging bottleneck is theta-boundary consistency in the off-pole bipolar test

Future work should be careful not to assume:

- that the current off-pole AMRVAC case is physically validated
- that full-pole stability has been solved
- that the analytical continuation attempt fixed the theta boundary problem

## Open Technical Questions

The following are still open and should not be treated as finalized:

- the final stretched radial mesh design for production PFSS/SCS-to-AMRVAC runs
- the final stretching law and resolution
- the production angular resolution
- the final PFSS/SCS `lmax`
- the detailed PFSS/SCS transition treatment
- whether the production AMRVAC initialization should use direct `B` or be rebuilt through a vector potential `A`
- what additional divergence-control strategy should be used for the final production initial condition
- how to construct a geometrically consistent CT face magnetic boundary for off-pole AMRVAC tests

## One-Sentence Summary

The confirmed project direction is to construct PFSS + SCS magnetic fields and evaluate them directly on the final stretched AMRVAC mesh for MHD relaxation, while the current active debugging focus is the boundary consistency of the analytical bipolar AMRVAC test, especially the off-pole theta treatment under CT.
