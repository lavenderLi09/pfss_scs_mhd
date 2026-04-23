# off_2270 Acceleration Plan, Keeping `threestep`

## Summary
- Keep `time_stepper='threestep'` unchanged in [amrvac.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/amrvac.par#L46).
- Treat `courantpar` increase as the first solver-parameter speed test. Start from the current `0.5`, then test `0.55`, `0.60`, and only then `0.65` if diagnostics stay clean.
- Do **not** enable `local_timestep` directly in the current `off_2270` setup. Even if `linde` were substituted for `glm`, the case still has gravity, and AMRVAC’s MHD source code explicitly says local timestep support is not added for gravity-like nonideal/source terms:
  - `3.1`: `/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1/src/mhd/mod_mhd_phys.t:4323`
  - `3.2`: `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mhd/mod_mhd_phys.t:4435`
  - `3.3`: [amrvac3.3 `src/mhd/mod_mhd_phys.t`](https://github.com/amrvac/amrvac/blob/amrvac3.3/src/mhd/mod_mhd_phys.t) around lines `4973-4975`

## Why `linde` Helps LTS More Than `glm`
- Yes, `linde` is generally the more diffusive divB option than `glm`; AMRVAC documents it as a diffusive/parabolic fix in [methods.md](/Users/zhaoyan/Documents/codes/amrvac3.2/doc/methods.md#L185).
- The reason LTS is said to work with `linde` is **not** that `linde` is universally better, but that its divB source update does not explicitly use the global `qdt` in the way `glm` does.
- In `3.2`, `glm` uses `qdt` directly in the psi damping and momentum/energy updates in `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mhd/mod_mhd_phys.t:5301`.
- In `3.2`, `linde` updates `B` through a diffusive `grad(divB)` term without an explicit `qdt` factor in `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mhd/mod_mhd_phys.t:5428`.
- So switching `glm -> linde` only removes the **divB-cleaning** incompatibility for LTS. It does **not** remove the gravity incompatibility for `off_2270`.

## Key Changes
- **Production-safe path for `off_2270`**
  - Fix the AMR-window bug first in [mod_usr.t](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/mod_usr.t#L621): `inwin(ixO^S)` is used but never assigned because the window logic is commented out.
  - Keep `glm`, gravity, and `threestep` unchanged.
  - Run a controlled `courantpar` ladder: `0.50 -> 0.55 -> 0.60 -> 0.65`.
  - If stable, test less frequent regridding by increasing `ditregrid` from `8` to `16`, then `32`.
  - Reduce monitoring overhead in test copies only: larger `ditsave_log`, and disable `autoconvert` during long relaxations unless outputs are needed immediately.

- **Separate LTS experiment**
  - Put this in a separate case cloned from `off_2270`.
  - Change only:
    - `local_timestep=.true.`
    - keep `typecourant='maxsum'`
    - switch `typedivbfix='glm'` to `typedivbfix='linde'`
  - Use this experiment only to test whether stock LTS can run numerically with the current geometry and boundaries.
  - Do **not** treat success there as validation for production use while gravity remains active.

- **If LTS must be made truly usable for off_2270**
  - The needed next step is AMRVAC source work so gravity-related source updates use cell-local `block%dt` consistently.
  - That is a code-modification project, not a parameter-only change.

## Test Plan
- Baseline from the same restart file and same runtime window:
  - measure wall time per `1000` iterations
  - record average `dt`, min `dt`, active block counts, and limiting block indices from [off.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/data/off.log)
- Safe acceleration sequence:
  1. AMR-window fix only
  2. `courantpar=0.55`
  3. `courantpar=0.60`
  4. `ditregrid=16`
  5. reduced log/convert overhead
- LTS experiment sequence:
  1. `glm -> linde` without LTS
  2. `linde + local_timestep`
- Compare at matched simulation time:
  - `rho`, thermal pressure, `Vr`
  - `Br/Bt/Bp` near the inner boundary
  - `divB`
  - new boundary artifacts or stronger smoothing
- Reject a variant if it gains speed but noticeably worsens field structure, inner-boundary behavior, or `divB`.

## Assumptions
- Your priority is faster relaxation with the existing physics, not changing the production numerics unless clearly justified.
- `linde` is acceptable only in a separate test because it may be more diffusive than `glm`.
- Recommendation: first pursue `courantpar` and AMR/regridding cleanup; keep LTS as a separate diagnostic experiment unless AMRVAC source code is extended for gravity-aware local timesteps.
