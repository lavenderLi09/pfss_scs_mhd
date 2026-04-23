# Ring-Average Full-Path Smoke Test (2026-04-10)

## Scope
- Target code: `/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average` (branch `ra/test-github-push`)
- Case: `off_2270_ringavg_32`
- Objective: verify that ring-average is active in both:
  - CFL effective `ds_phi` path (`mod_dt.t`)
  - reconstruction/evolution path (`mod_finite_volume.t`, `polar_ring_evolve`)

## Run Setup
- 2 MPI ranks
- Base input: `amrvac.par`
- No-data 4-step override: `ringavg_nodat_4steps.par`

### Baseline (full off)
```bash
mpirun -np 2 ./amrvac -i amrvac.par ringavg_nodat_4steps.par ringavg_full_off.par > out_fullimpl_off_4steps 2>&1
```

### Ring full on
```bash
mpirun -np 2 ./amrvac -i amrvac.par ringavg_nodat_4steps.par ringavg_full_on.par > out_fullimpl_on_4steps 2>&1
```

## Activation Evidence
- CFL ring active:
  - `RING_DEBUG active=T max_nchunk=32 ...` in `out_fullimpl_on_4steps`
- Evolution-path ring active:
  - `RING_EVOLVE_DEBUG ...` in `out_fullimpl_on_4steps`

So the new reconstruction-path switch (`polar_ring_evolve`) is being executed.

## Follow-up Fix
- In reconstruction-path chunk sizing, use `abs(ds_theta)` and `abs(ds_phi)` to avoid sign-driven chunk inflation.
- Quick rerun (`out_fullimpl_on_1step_absfix`) confirms:
  - `RING_EVOLVE_DEBUG nchunk=16 ds_phi=+...`
  - activation still works and run remains stable.

## dt / Limiter Outcome
- `it=0..4` timestep is effectively unchanged between `full_off` and `full_on`.
- Limiter remains the same block:
  - `DT_CFL_LIMITER igrid=209`
- Dominant directional contribution in limiter print remains radial (`dirmax(1)` largest).

Interpretation:
- Ring-average logic is active.
- For this specific case state, global minimum dt is still controlled by non-phi contribution (radial bottleneck), so enabling ring-average does not increase dt here.
