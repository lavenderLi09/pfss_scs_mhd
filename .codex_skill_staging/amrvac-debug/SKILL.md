---

## name: amrvac-debug
description: Use when building, running, or debugging spherical MPI-AMRVAC cases locally, especially for initialization failures, negative pressure, timestep collapse, par-file confusion, or boundary-versus-physics diagnosis in heliosphere and coronal setups.

# AMRVAC Debug

## Overview

Use this skill for safe local debugging of **spherical** AMRVAC cases, especially solar-wind and coronal-field problems.

Focus on evidence, not guesses: confirm the geometry, boundary logic, namelist placement, runtime log, and output fields before claiming a fix.

This skill is intentionally **not** a general Cartesian AMRVAC skill. Its heuristics assume spherical coordinates, radial shells, theta handling, and heliosphere-style interpretation.

## When To Use

- AMRVAC crashes at `Iteration 0` or the first step
- a case initializes but does not evolve physically
- `fix_small_values`, `small_values_method`, or other `.par` settings may not actually be taking effect
- pole or off-pole spherical domains behave differently than expected
- you need a PC-safe debug ladder instead of a large MPI run
- VTU or DAT outputs must be checked against logs and initialization assumptions

Do not use this skill for PFSS/SCS coefficient generation itself; use it once the work reaches an AMRVAC case or output.

## Scope Classification

### Primary scope: spherical AMRVAC

Use this skill directly when the case involves one or more of these:

- spherical coordinates
- inner and outer radial boundaries
- theta pole logic or trimmed-theta domains
- radial mesh stretching
- Parker-wind or heliosphere-like interpretation
- spherical components such as `Vr/Vt/Vp` or `Br/Bt/Bp`

### Limited overlap: Cartesian AMRVAC

Some generic checks still carry over:

- verify `.par` namelist placement
- use staged init and smoke tests
- compare log evidence with actual outputs
- distinguish numerical stability from physical plausibility

But do **not** apply this skill's spherical heuristics to Cartesian cases. In Cartesian setups, do not reason in terms of:

- poles or off-pole domains
- theta cuts
- radial shells
- Parker background consistency
- spherical velocity or magnetic components

For Cartesian AMRVAC debugging, use this skill only for the generic run-discipline pieces above, and rely on general debugging skills for the physics interpretation.

## Workflow

1. Read the active `amrvac.par`, the user module, and any case-specific runner files before changing anything.
2. Confirm which domain is actually being run:
  - radial range
  - theta range
  - phi range
  - base resolution
  - radial stretching settings
3. Check `.par` namelist placement carefully. In AMRVAC, parameters only belong to the active namelist block before its closing `/`.
4. Start with the lightest safe run:
  - one core by default
  - initialization-only if structure is still uncertain
  - then a short smoke test
5. Compare runtime evidence:
  - stdout/log messages
  - whether the code reached `Iteration 1`
  - whether VTU/DAT files were written
6. Separate numerical stability from physical quality:
  - “did not crash” is not the same as “physically self-consistent”
7. After each change, record what was changed, what was rerun, and what improved or regressed.

## Safe Local Run Pattern

- Prefer one core unless the user explicitly asks for more.
- Use staged cases or staged parameter files:
  - init
  - smoke
  - medium
  - long
- Stop the ladder when the first hard failure appears.
- Keep output prefixes distinct so failed and successful runs are easy to compare.

## AMRVAC Checks

### Par-file checks

- Verify every edited parameter is inside the correct namelist.
- If a result contradicts expectations, suspect namelist placement before suspecting the solver.

### Geometry checks

- Confirm whether theta includes the true poles.
- `stretch_dim(1)='uni'` in AMRVAC means stretched unidirectional mesh, not uniform radial cells.
- If AMRVAC prints `There is no northpole!` or `There is no southpole!`, that is a geometry statement about the theta extent.

### Runtime checks

- Check whether the failure happens:
  - during initialization
  - at `Iteration 0`
  - after several steps
- Distinguish among:
  - negative pressure or density
  - NaN or Inf
  - `dt=0` or timestep collapse
  - boundary-driven oscillation

### Physics checks

Inspect at least:

- `rho`
- `p`
- `Vr`, `Vt`, `Vp`
- `Br`, `Bt`, `Bp`

Useful interpretations:

- magnetic field almost unchanged with large `Vt/Vp` near theta cuts usually means boundary artifacts, not meaningful relaxation
- local negative `Vr` near the inner boundary often points to boundary inconsistency or impulsive adjustment
- pressure rings near the inner shell and theta cuts usually implicate the theta boundary plus CT face-field handling

## Spherical Heliosphere Heuristics

- Full-pole and off-pole cases are different numerical problems; do not treat them as interchangeable.
- If a full-pole case fails immediately but a trimmed-theta case runs, suspect polar regularity or boundary consistency before changing core physics.
- For Parker-plus-magnetic-field tests, keep inner radial boundary states consistent with initialization whenever possible.
- Do not claim “solar-wind MHD relaxation” from a short run just because the code stayed stable; check whether the flow changes are physically located and plausible.

## Reporting

Summaries should always include:

- test case name
- geometry and resolution
- whether radial stretching is on
- whether the run is full-pole or off-pole
- whether `fix_small_values` is on or off
- where the first anomaly occurs
- whether the result is numerically stable, physically plausible, both, or neither

## Common Mistakes

- editing `.par` values outside their namelist block
- assuming a stable run is already a good physical solution
- treating pole warnings as proof that `pole` boundaries are still active
- changing too many things between smoke tests
- running heavy MPI tests before the initialization and first-step behavior are understood

