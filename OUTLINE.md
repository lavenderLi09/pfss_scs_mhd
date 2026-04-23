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

This directory is the AMRVAC-side workspace. It now contains several generations of cases:

- inherited or downloaded reference cases
- reduced analytical bipolar debug cases
- isolated HPC-only variants
- PFSS-driven production-like cases
- comparison / archive / planning material that records why later cases were created

The current working map is:

- [polytropic_bipolar](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/polytropic_bipolar)
  - original downloaded/reference polytropic bipolar AMRVAC case
  - mainly used to inspect baseline coding patterns and retain a known runnable reference
  - recorded test status:
    - [data/out.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/polytropic_bipolar/data/out.txt) shows a completed 24-core HPC relaxation job on 2025-02-27

- [analytic_bipolar_1to20](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20)
  - first stripped-down analytical bipolar precursor on `1..20 Rsun`
  - used to remove flux-rope and data-driven complications and test a Parker-plus-bipole setup
  - useful files:
    - [init_check.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20/init_check.par)
    - [data/init_check_deeper.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20/data/init_check_deeper.log)
  - recorded test status:
    - `init_check_deeper0000` was produced as an initialization snapshot
    - older logs such as [data/init_check.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20/data/init_check.log) and [data/init_check_theta.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20/data/init_check_theta.log) contain `NaN` output at `t=0`
    - this directory should therefore be treated as an early precursor, not the current validated local baseline

- [analytic_bipolar_1to20_stretched](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched)
  - current local analytical-bipolar debugging focus
  - keeps the reduced physics of `analytic_bipolar_1to20`, but adds the stretched radial mesh and a staged PC stability ladder
  - important files:
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
      - single-core local runner for staged tests
    - [analysis/pressure_anomaly](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/analysis/pressure_anomaly)
      - summaries and debug records for the off-pole pressure anomaly
  - recorded test status:
    - [data/pc_init.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_init.stdout) finishes cleanly and writes `pc_init0000`
    - [data/pc_smoke.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_smoke.stdout) fails at iteration `0` with negative gas pressure and `MPI_ABORT`
    - [data/pc_smoke_offpole.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_smoke_offpole.stdout) reaches `it=50` and finishes
    - [data/pc_smoke_offpole_oldbc.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/data/pc_smoke_offpole_oldbc.stdout) also finishes, giving a comparison point for boundary-condition variants
    - the pressure-anomaly reports conclude that the remaining off-pole issue is likely a theta-boundary / CT face-treatment artifact rather than the same catastrophic full-pole failure

- [bipolar_hpc_test](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test)
  - isolated HPC-only sibling of `analytic_bipolar_1to20_stretched`
  - used when mesh size, AMR rules, and `mod_usr.t` need to diverge from the PC ladder without polluting the parent local case
  - key subdirectories:
    - [data](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data)
      - local preflight outputs
    - [data_jb_50](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_jb_50)
      - dynamic-AMR test using `|J|*Delta/|B|` sheet forcing
    - [data_jb_50_clean](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_jb_50_clean)
      - clean-build rerun of the same `jb_50` variant
    - [data_rgfreeze](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_rgfreeze)
      - regridding-frozen isolation test
  - recorded test status:
    - [data/hpc_pilot_init_offpole.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data/hpc_pilot_init_offpole.stdout) finishes cleanly
    - [data/hpc_pilot_offpole.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data/hpc_pilot_offpole.stdout) finishes and produces snapshots through `00010`
    - [data_rgfreeze/hpc_pilot_offpole_rgfreeze.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_rgfreeze/hpc_pilot_offpole_rgfreeze.stdout) finishes, isolating the no-regrid variant
    - [data_jb_50/hpc_pilot_offpole_jb_50.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_jb_50/hpc_pilot_offpole_jb_50.stdout) finishes the dynamic `J/B` AMR test
    - [data_jb_50_clean/hpc_pilot_offpole_jb_50_clean.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/bipolar_hpc_test/data_jb_50_clean/hpc_pilot_offpole_jb_50_clean.stdout) reaches `it=50` after a strict clean rebuild and is the clearest archived rerun of that branch

- [off_2270](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270)
  - the main PFSS-driven production-like AMRVAC case in this repository
  - loads external initial magnetic-field data from [initial/OFF_combined_lmax10_q1.008_nr600.bin](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/initial/OFF_combined_lmax10_q1.008_nr600.bin)
  - [analysis](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis) stores post-run diagnostics for timestep bottlenecks, AMR coverage, and physical evolution
  - recorded test / analysis status:
    - the run log [data/off.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/data/off.log) shows a long ongoing/stable baseline with `dt ~ 1.08e-06` at late sampled times
    - [analysis/keyframe_dt_stability_report_2026-04-06.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/keyframe_dt_stability_report_2026-04-06.md) analyzes `off0000`, `off0016`, `off0032`, `off0064`, `off0096`, and `off0112`
    - [analysis/off_2270_summary_2026-04-07.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/off_2270_summary_2026-04-07.md) concludes that the dominant timestep limit is the small-theta inner-boundary `phi` geometry, not a global physical blow-up
    - the scripts [analysis/amr_level_coverage_stats.py](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/amr_level_coverage_stats.py), [analysis/analyze_snapshot_dt_bottleneck.py](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/analyze_snapshot_dt_bottleneck.py), and [analysis/test_snapshot_dt_outputs.py](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/test_snapshot_dt_outputs.py) exist to reproduce or extend those diagnostics

- [off_2270_initwind](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind)
  - local init-only reproduction / sanity-check case derived from `off_2270`
  - used to verify whether the newly rebuilt initial atmosphere is numerically reasonable before a full production relaxation
  - important files:
    - [init_check.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/init_check.par)
    - [run_init_check.sh](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/run_init_check.sh)
    - [analysis/check_init_atmosphere.py](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/check_init_atmosphere.py)
  - recorded test status:
    - [data/init_check.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/data/init_check.stdout) finishes cleanly
    - [analysis/init_atmosphere_summary.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/init_atmosphere_summary.md) reports PASS for conservative checks, primitive checks, and radial-profile checks
    - [data/amr_probe.stdout](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/data/amr_probe.stdout) finishes and provides an AMR-probe comparison output
    - [analysis/hpc_off0000_vs_initcheck_report_2026-04-08.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/hpc_off0000_vs_initcheck_report_2026-04-08.md) compares the local init-only snapshot against a downloaded HPC reference and concludes that the local init is reasonable as an initialization even though it is not identical to the restart-context HPC `off0000`

- [off_2270_lts_32_linde](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_lts_32_linde)
  - separate AMRVAC 3.2 experiment for `local_timestep=.true.` with `typedivbfix='linde'`
  - exists to test LTS feasibility without modifying the production `off_2270` directory
  - recorded test status:
    - [smoke.out](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_lts_32_linde/smoke.out) shows a short smoke run that finishes 5 steps from the `off0000` restart
    - [data/lts_smoke.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_lts_32_linde/data/lts_smoke.log) confirms the short run remains finite
    - [data/lts_off.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_lts_32_linde/data/lts_off.log) shows a longer run collapsing around iterations `30-31`, with `dt` dropping to near zero and state variables blowing up
    - current interpretation: the LTS experiment is not yet a production-safe acceleration path

- [ring_average_design](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design)
  - dedicated development area for the polar ring-average idea based on Zhang et al. (JCP 2019)
  - contains paper notes, implementation plans, a small documentation site, and one active AMRVAC case cloned from the `off_2270` lineage
  - important top-level files:
    - [2026-04-08-amrvac32-ring-average-plan.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/2026-04-08-amrvac32-ring-average-plan.md)
    - [2026-04-08-ring-average-paper-notes.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/2026-04-08-ring-average-paper-notes.md)
    - [git_workflow_ring_average.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/git_workflow_ring_average.md)
    - [paper.pdf](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/paper.pdf)

- [ring_average_design/off_2270_ringavg_32](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32)
  - active ring-average test case, derived from the `off_2270` family but used specifically for AMRVAC 3.2 ring-average implementation and verification
  - this directory keeps both case inputs and a local consolidated context file:
    - [agent.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/agent.md)
    - [RING_AVERAGE_METHOD_SUMMARY.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_METHOD_SUMMARY.md)
    - [RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md)
    - [RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md)
    - [ringavg_full_on.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_full_on.par)
    - [ringavg_full_off.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_full_off.par)
    - [ringavg_nodat_1step.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_nodat_1step.par)
    - [ringavg_nodat_4steps.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_nodat_4steps.par)
    - [ringavg_poleamr2_test.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_poleamr2_test.par)
  - [analysis](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis) reuses and extends the `off_2270`-style diagnostics, while `out_*` files store ring-average on/off and smoke-test logs
  - recorded test / implementation status:
    - the ring-average code path is reported as active in both the CFL path and the reconstruction/evolution path, as summarized in [RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md)
    - on/off smoke tests in [out_ringavg_smoke_off](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_ringavg_smoke_off) and [out_ringavg_smoke_on](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_ringavg_smoke_on) show no immediate global `dt` increase
    - 4-step full-implementation comparison in [out_fullimpl_off_4steps](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_off_4steps) and [out_fullimpl_on_4steps](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_on_4steps) shows the phi-direction CFL term is strongly reduced when ring-average is on, but the global limiter remains radial-dominated so `dt` stays nearly unchanged
    - the pole-AMR2 comparison documented in [RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md) reaches the same practical conclusion: ring-average clearly activates and suppresses the phi term, but the snapshot `dt` is still controlled elsewhere
    - current implementation level is best described as `cell-centered averaging-reconstruction + CFL effective scale + optional post-ring clean`; it is not yet the full CT/staggered magnetic-field version of the paper

- [real_bipo_test](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test)
  - larger and more production-like bipolar/HPC experiments kept mainly as artifacts and alternative parameter variants
  - [data_large/regrid4](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test/data_large/regrid4) and [data_large/regrid8](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test/data_large/regrid8) retain large HPC runs with different regridding settings
  - [data_large/local](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test/data_large/local) is a local debugging copy with `debug_run.err/out` and `lldb_run.log`
  - recorded test status:
    - `regrid4` and `regrid8` both contain HPC job records and generated `.dat` outputs
    - [data_large/regrid4/out.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test/data_large/regrid4/out.txt) and [data_large/regrid8/out.txt](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/real_bipo_test/data_large/regrid8/out.txt) show the runs were terminated by owner after producing early snapshots, so these are exploratory artifacts rather than a settled baseline

- [archive](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/archive)
  - stores archived test inputs for reproducibility
  - each subdirectory keeps the exact `mod_usr.t`, `.par`, and parameter-change note used for one test branch
  - current archived branches include:
    - `2026-03-27_bipolar_hpc_test_jb_50`
    - `2026-03-27_bipolar_hpc_test_rgfreeze`
    - `2026-03-30_bipolar_hpc_test_jb_50_clean`

- [plan](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/plan)
  - planning notes and helper checks for future AMRVAC work
  - currently focused on `off_2270` acceleration options and the separate LTS experiment path

- [context_localtime.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/context_localtime.md)
  - source-inspection note on AMRVAC `local_timestep` support in versions `3.2` and `3.3`
  - this is not a case directory, but it is an important context file for interpreting why `off_2270_lts_32_linde` was created

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

### 6. Ring-average development is active, but it is not yet a proven timestep-speedup for `off_2270`

The ring-average development branch is tracked in:

- [off_2270_ringavg_32/agent.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/agent.md)

What is currently established:

- the new ring-average logic does activate in the intended code paths
- it can strongly reduce the polar `phi` CFL contribution in on/off comparisons
- but several preserved tests still show the same global `dt`, because the actual limiter in those snapshots remains dominated by another contribution, usually the radial term

So the correct current interpretation is:

- ring-average is a live implementation branch with real activation evidence
- it is **not** yet evidence that the production-like `off_2270` case will immediately run faster

### 7. The current ring-average implementation is not the full paper-grade CT version

The current AMRVAC 3.2 ring-average work appears to have reached:

- cell-centered averaging-reconstruction
- effective CFL rescaling
- optional post-ring cleaning hook

But it has **not** yet reached the full staggered / CT magnetic update path described in the paper.

So future work should not assume:

- that the present ring-average branch already preserves CT structure in the same way as the full published method
- that ring-average is production-ready for the main `off_2270` workflow

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

Future work must also preserve generated run artifacts unless the user explicitly approves otherwise:

- do **not** delete, discard, or clean `data/` directories, logs, `.dat`, `.vtu`, `.vts`, or other simulation outputs without explicit user permission
- do **not** overwrite existing simulation outputs (`.dat`, `.vtu`, `.vts`, logs) by reusing the same `base_filename` or convert target names
- for every new run or convert task, use a new output prefix (recommended: append a timestamp, e.g. `*_YYYYMMDD_HHMMSS`) so prior artifacts remain intact
- this no-overwrite rule applies to every thread and every AMRVAC case in this project unless the user explicitly asks to overwrite
- do **not** assume worktree cleanup or folder reorganization implies permission to remove generated results
- when reorganizing cases, move or copy requested outputs first and confirm before removing the original location

Future work must also archive test inputs for every AMRVAC run across this whole project:

- before each test run, create or use an `archive/` directory at the same directory level as the current AMRVAC case directory
- copy the exact post-edit `mod_usr.t` used for that test into the archive
- copy the exact `.par` file used for that test into the archive
- save a short accompanying record of the modified parameters for that test, so the archived inputs make clear what changed relative to the previous or baseline setup
- treat this archiving step as required for every test, not as an optional cleanup task after the run
- if `mod_usr.t` was modified before a test, then the build step is mandatory and must be run in this order: `make clean` -> `setup.pl -d=3 -arch=<arch>` -> `make`
- this rule applies to every thread and every AMRVAC case in this project; do not skip `setup.pl` even when source files already exist
- do not assume pre-existing `makefile`/`amrvac` binaries in a case directory; the case build structure is expected to be materialized by `setup.pl` followed by `make`

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
