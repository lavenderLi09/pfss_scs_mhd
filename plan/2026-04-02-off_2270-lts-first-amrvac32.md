# off_2270 LTS-First Plan on AMRVAC 3.2

This case does not use a theta-phi AMR window; refinement is not window-restricted.

## Summary
- Clarify the plan to state explicitly: there is **no theta-phi refinement window** for `off_2270`; the commented window code in [mod_usr.t](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/mod_usr.t#L621) is dead code and should not be restored.
- The first implementation task is a **separate LTS experiment** on **AMRVAC 3.2**, not a direct change to the production `off_2270` case.
- Keep the production case on `glm`, gravity, and `threestep` for now. Do not change `time_stepper` in either the production case or the first experiment.

## Implementation Changes
- **Create a separate experiment case**
  - Clone `off_2270` to:
    - `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_lts_32_linde`
- **In the experiment case**
  - keep `time_stepper='threestep'`
  - keep `typecourant='maxsum'`
  - set `local_timestep=.true.` in `&methodlist`
  - change `typedivbfix='glm'` to `typedivbfix='linde'`
  - leave gravity enabled so the test reflects the real case, but treat it as a numerical experiment only, not production validation
- **Fix the AMR logic in the experiment copy only**
  - remove the unused `inwin` filter path from `special_refine_grid`
  - make the custom `J/B` refinement criterion global
  - do not introduce any theta-phi window
- Do not modify the original `off_2270` case during this first experiment.

## Test Plan
- **Baseline**
  - use the current `off_2270` run and [off.log](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/data/off.log) as the baseline
  - record wall time per `1000` iterations, average `dt`, active block counts, and limiting block indices
- **Experiment sequence**
  1. clone case and make only the AMR global-criterion fix
  2. switch `glm -> linde`
  3. enable `local_timestep=.true.`
  4. build against AMRVAC 3.2
  5. run a short smoke test from the same restart
  6. if stable, run a longer matched-time test
- **Compare against baseline at the same simulation time**
  - `rho`
  - thermal pressure
  - `Vr`
  - `Br/Bt/Bp` near the inner boundary
  - `divB`
  - any stronger smoothing or boundary artifacts
- **Acceptance**
  - the experiment is acceptable only if it remains numerically clean and gives a meaningful wall-time improvement
  - even if it runs faster, it is **not** promoted to production use unless gravity-coupled LTS behavior also looks physically acceptable

## Follow-Up After LTS Experiment
- If the LTS experiment fails or is too diffusive, revert to the production-safe path:
  - keep `glm`
  - keep `threestep`
  - test `courantpar = 0.55 -> 0.60 -> 0.65`
  - test larger `ditregrid`
- If the LTS experiment is promising but still questionable, the next project is AMRVAC source work so gravity-related source terms use cell-local `block%dt` consistently.

## Assumptions
- The first implementation target is AMRVAC 3.2 because it is available locally and matches the current case workflow.
- `off_2270` intentionally does not use a theta-phi refinement window.
- The first experiment is allowed to use `linde`, but only in the separate test case.
