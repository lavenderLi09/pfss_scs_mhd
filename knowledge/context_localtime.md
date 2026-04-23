# AMRVAC 3.2 / 3.3 local time stepping context

## What we checked

- Official repo: `amrvac/amrvac`
- Relevant branches on GitHub:
  - `amrvac3.2`
  - `amrvac3.3`
- Local clone used for source inspection:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2`

## Version conclusion

- `local time stepping` was added before `3.2` and `3.3`.
- The GitHub release notes for `v3.1` explicitly mention:
  - `Addition of local time stepping option`
- Therefore:
  - `3.2` supports it
  - `3.3` supports it
  - earliest confirmed official version is `3.1`

## Whether 3.2 has a test/example using local time stepping

- I searched the `amrvac3.2` branch for:
  - `local_timestep`
  - `local_timestep = .true.`
  - related `.par` files under `tests/`
- Result:
  - `local_timestep` appears in source code
  - I did **not** find a shipped test/example `.par` file that actually enables it
- So the code supports it, but the repository does not appear to include a ready-made test case for it in `3.2`.

## Where local time stepping is defined

- Global switch definition:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_global_parameters.t`
  - `logical :: local_timestep = .false.`
- Input namelist location:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/io/mod_input_output.t`
  - `local_timestep` is part of `&methodlist`
- Per-cell storage:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_physicaldata.t`
  - `block%dt` stores cell-local timesteps
- Allocation/deallocation:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/amr/mod_amr_solution_node.t`

## How local time stepping is computed in 3.2

- Main computation happens in:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_dt.t`
- The local timestep is computed from a CFL-like expression using `typecourant='maxsum'`.
- For each cell, the code builds:
  - `cmaxtot = sum_over_dimensions(cmax / ds)`
- Then sets:
  - `block%dt = courantpar / cmaxtot`
- In code this is the key line:
  - `block%dt(hxO^S) = courantpar/cmaxtot(hxO^S)`

## Important compatibility restriction

- In `3.2`, `local_timestep` is only implemented with:
  - `typecourant='maxsum'`
- The code explicitly aborts for:
  - `summax`
  - `minimum`

## Where local dt is actually used

- Flux updates use `block%dt` instead of a single global `qdt` when `local_timestep` is enabled.
- Main files:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_finite_volume.t`
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mod_finite_difference.t`
- Some source terms also branch on `local_timestep`, for example:
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/physics/mod_rotating_frame.t`
  - `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mhd/mod_mhd_phys.t` in `add_pe0_divv()`

## Important limitation in MHD

The key comment found in:

- `/Users/zhaoyan/Documents/codes/amrvac3.2/src/mhd/mod_mhd_phys.t`

states:

- `local_timestep support is only added for splitting`
- `but not for other nonideal terms such gravity, RC, viscosity,..`
- `it will also only work for divbfix 'linde'`

## What that comment means

This is the practical interpretation:

- `local_timestep` is **not** fully wired through all MHD source-term paths.
- It is supported in the main flux/CFL update path and in some split source/update paths.
- But many additional MHD physics modules still use a single global `qdt` rather than per-cell `block%dt`.
- The code specifically warns about nonideal/source terms such as:
  - gravity
  - radiative cooling (`RC`)
  - viscosity
- For divB cleaning, the comment says it is only safe with:
  - `divbfix='linde'`
  because that update path does not require the same dt-dependent modification.

## Practical takeaway

- In `AMRVAC 3.2`, `local_timestep` looks like a partial/limited feature rather than something universally safe for arbitrary MHD setups.
- It may be reasonable for simpler/cleaner MHD cases.
- It is riskier when combining MHD with extra physics or source terms, especially:
  - gravity
  - radiative cooling
  - viscosity
  - non-`linde` divB fixes
- If trying it in `3.2`, safest baseline assumptions are:
  - use `typecourant='maxsum'`
  - avoid extra nonideal/source terms unless verified carefully
  - treat it as a feature that may need case-by-case validation

## Suggested minimal parameter sketch

```fortran
&methodlist
  local_timestep = .true.
  typecourant    = 'maxsum'
/
```

## Good next step in a new chat

If continuing elsewhere, the next useful task is:

- inspect a specific `.par` file / physics setup
- decide whether that exact MHD case is suitable for `local_timestep` in `3.2`
- if needed, build a minimal safe test case before using it in production
