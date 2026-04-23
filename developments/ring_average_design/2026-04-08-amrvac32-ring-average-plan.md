# Ring-Average Sandbox Plan on AMRVAC 3.2

## Summary
- Work in an isolated AMRVAC source copy based on the current local `amrvac3.2` branch.
- Implement only **V1 CFL-only ring-average** in that copy.
- Create one new test case folder based on the current `off_2270` setup, because it is the cleanest true-pole `phi`-geometry bottleneck case in this project.
- Download the original Ring Average paper PDF into that new test folder if the publisher PDF is accessible; otherwise stop at the DOI landing page and report the access block.

## Workspace Layout
- Source-tree copy:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`
  - Create it as a full filesystem copy of `/Users/zhaoyan/Documents/codes/amrvac3.2` at the current checked-out `amrvac3.2` branch state.
  - Keep the copied repo on a new local working branch named `ring-average-v1`, created from that copied checkout's current `amrvac3.2` HEAD.
- New working test folder:
  - `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_ringavg_32`
  - Create it by copying the current production-like case `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270`.
- Paper location inside the new test folder:
  - `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_ringavg_32/references/zhang_2019_ring_average_jcp.pdf`

## Implementation Changes
- Apply the code change only in `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`.
- Add four new runtime parameters in the copied AMRVAC source:
  - `polar_ring_average = .false.`
  - `polar_ring_theta_cap_north = 0.d0`
  - `polar_ring_theta_cap_south = 0.d0`
  - `polar_ring_max_chunk = 8`
- Implement the helper as a CFL-only utility used from [`mod_dt.t`](/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_dt.t):
  - active only for `coordinate==spherical`, `ndim==3`, and `idim==phi_`
  - north-cap mask: `theta <= polar_ring_theta_cap_north`
  - south-cap mask: `theta >= pi - polar_ring_theta_cap_south`
  - chunk rule: `nchunk = min(polar_ring_max_chunk, largest power of two <= max(1, ds_theta/ds_phi))`
  - effective width: `ds_phi_eff = nchunk * ds_phi`
- Use `ds_phi_eff` only in the CFL denominator inside `getdt_courant`.
- Do not modify:
  - stored `block%ds`, `block%surfaceC`, `block%dvolume`
  - `mod_finite_volume`
  - `mod_finite_difference`
  - `mhd/mod_mhd_phys`
  - `mod_fix_conserve`
  - `mod_ghostcells_update`

## Test Case Setup
- Base the new case on `off_2270`, not on the off-pole bipolar test.
- Keep the copied case numerically as close as possible to current tests:
  - retain the same physics, grid, AMR behavior, and divB method as the current `off_2270`
  - use the same `amrvac.par` as baseline input
- Add one new short-run input file in the new folder:
  - `amrvac_ringavg_smoke.par`
- In that smoke file:
  - enable `polar_ring_average = .true.`
  - set an initial north-cap-only test, because the current `off_2270` limiter is on the small-`theta` side
  - start with `polar_ring_theta_cap_north = 0.10`
  - set `polar_ring_theta_cap_south = 0.d0`
  - set `polar_ring_max_chunk = 8`
  - keep all other settings matched to the copied baseline
- Add one short runner script in the new folder, mirroring the style of existing local test scripts:
  - build against `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`
  - run the smoke file
  - capture log output for timestep comparison

## Validation
- Primary acceptance check in `off_2270_ringavg_32`:
  - verify the worst-cell limiter is materially relaxed relative to the copied baseline
  - verify no new inner-boundary or polar artifact appears in the short run
- Secondary guardrail:
  - reuse the existing full-pole bipolar case as a separate regression only after the `off_2270` smoke works
  - do not make a second new folder for that in this first pass
- Non-goal:
  - do not expect this to fix off-pole theta-cut artifacts

## Paper Acquisition
- Download target: the original publisher paper for Zhang et al. 2019, JCP 376, 276-294.
- Canonical source:
  - [ScienceDirect article / DOI landing page](https://www.sciencedirect.com/science/article/abs/pii/S0021999118305436)
  - [DOI metadata summary at NCAR/UCAR](https://impacts.ucar.edu/en/publications/conservative-averaging-reconstruction-techniques-ring-average-for/)
- Download rule:
  - first attempt the publisher PDF from the ScienceDirect/DOI route and save it as `zhang_2019_ring_average_jcp.pdf`
  - if access is blocked, do not substitute a non-authoritative mirror; instead leave the folder with no PDF and report that the original PDF requires access
- This keeps the folder aligned with your request for the original paper rather than a secondary copy.

## Assumptions
- "same branch as amrvac3.2" is interpreted as "same current code state as the local `amrvac3.2` branch," but edits happen on a new working branch inside the copied repo so the original source tree stays untouched.
- One new test folder is sufficient for the first pass, and it should be the `off_2270`-based true-pole case.
- The first implementation target remains CFL-only ring-average, not full conservative ring-chunk flux evolution.
