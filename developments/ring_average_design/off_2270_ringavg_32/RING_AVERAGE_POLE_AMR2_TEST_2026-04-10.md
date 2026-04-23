# Ring-Average Test Under Pole AMR Level-2 Condition (2026-04-10)

## Goal
Check whether ring-average is still active/effective when polar AMR is allowed to reach level 2.

## Test Setup
- Case: `off_2270_ringavg_32`
- MPI: 2 ranks
- Base input: `amrvac.par`
- Overrides:
  - `ringavg_poleamr2_test.par`
  - `ringavg_nodat_1step.par`
  - `ringavg_full_off.par` or `ringavg_full_on.par`

## Important Note
For this test only, `mod_usr.t:special_refine_grid` was temporarily modified during the run to:
- disable pole guard exclusion
- force refinement in near-pole cells up to `refine_max_level`
- print `POLE_AMR_DEBUG ...` once as evidence

After collecting logs, `mod_usr.t` was restored.

## Commands Used
```bash
mpirun -np 2 ./amrvac -i amrvac.par ringavg_poleamr2_test.par ringavg_nodat_1step.par ringavg_full_off.par > out_poleamr2_fulloff_1step 2>&1
mpirun -np 2 ./amrvac -i amrvac.par ringavg_poleamr2_test.par ringavg_nodat_1step.par ringavg_full_on.par  > out_poleamr2_fullon_1step 2>&1
```

## Evidence That Pole AMR Became Finer
- Both runs printed:
  - `POLE_AMR_DEBUG forced refine level=1 ...`
- Ring-on run printed:
  - `RING_DEBUG active=T max_nchunk=32 min_ds_phi=1.3763E-04 ...`
- Previous non-forced baseline had:
  - `min_ds_phi=5.3578E-04`

So polar `ds_phi` became about `3.89x` smaller, consistent with stronger/finer polar AMR.

## Ring On/Off Comparison At it=0
- Off (`out_poleamr2_fulloff_1step`):
  - `CFL_DEBUG ring=F ... dirmax_phi=4.0855E+04 ...`
  - `dt = 6.0612E-06`
- On (`out_poleamr2_fullon_1step`):
  - `RING_DEBUG active=T ...`
  - `CFL_DEBUG ring=T ... dirmax_phi=2.6419E+03 ...`
  - `dt = 6.0612E-06`

Interpretation:
- Ring-average clearly **activated** and strongly reduced the phi-direction CFL term.
- But global minimum dt remained unchanged because limiter block/term stayed the same (`DT_CFL_LIMITER igrid=2543`, radial-dominated contribution), so timestep did not increase in this snapshot.

