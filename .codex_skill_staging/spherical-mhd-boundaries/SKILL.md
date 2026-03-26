---
name: spherical-mhd-boundaries
description: Use when implementing or debugging spherical MHD boundaries, especially full-pole versus off-pole theta domains, analytical continuation on trimmed theta cuts, CT versus non-CT magnetic handling, or boundary-driven artifacts near theta cuts and inner radial shells.
---

# Spherical MHD Boundaries

## Overview

Use this skill when boundary logic, not core physics, is likely driving the result.

The main goal is to keep geometry, ghost cells, and magnetic representation self-consistent, especially in spherical coordinates with constrained transport.

This skill is for spherical and polar-like geometries. Do not reuse it as a Cartesian boundary guide.

## When To Use

- choosing between `pole` and `special` theta boundaries
- debugging `specialbound_usr`
- comparing full-pole and off-pole spherical domains
- deciding whether a trimmed-theta domain should behave as a cut surface or a true pole
- implementing analytic continuation across spherical cut surfaces
- handling staggered face fields together with cell-centered magnetic fields
- locating pressure rings, strong `Vt/Vp`, or local negative `Vr` near theta cuts or the inner boundary

## First Decision: Full Pole Or Off Pole

### Full-pole domain

Use true pole handling only when the theta range includes the actual poles:

- north pole: `theta = 0`
- south pole: `theta = pi`

If the domain includes both, `pole` boundary logic is appropriate.

AMRVAC test anchor:

- [`tests/hd/blast_wave_spherical_stretched_3D/amrvac.par`](/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1/tests/hd/blast_wave_spherical_stretched_3D/amrvac.par) uses
  - `xprobmin2=0`
  - `xprobmax2=0.5`
  - `typeboundary_min2='pole'`
  - `typeboundary_max2='pole'`

Use this branch when you want a genuine spherical shell with active pole logic.

### Off-pole trimmed domain

If theta does not include the true poles, treat the theta ends as ordinary cut surfaces, not poles.

That means:

- do not reason about them as pi-periodic poles
- prefer explicit `special` theta boundaries if the physical continuation matters
- interpret AMRVAC messages like `There is no northpole!` as geometry diagnostics, not as proof of a boundary bug

Use this branch when theta is intentionally trimmed to avoid polar regularity issues or to isolate a local spherical patch.

## Boundary Families

### Full-pole spherical

Typical use:

- full spherical shell
- theta includes the true poles
- `pole` boundaries are part of the intended geometry

Checks:

- confirm the theta interval really includes the poles
- confirm the `.par` file still uses `pole`, not an inherited `special`
- if the run fails near a polar ring at the first step, suspect polar regularity or pole-plus-inner-boundary inconsistency before changing the bulk physics

### Off-pole spherical

Typical use:

- trimmed-theta heliosphere test
- local spherical patch
- diagnostic runs that intentionally remove the pole singularity

Checks:

- treat the theta ends as cut surfaces with their own physical meaning
- decide whether the cut is:
  - a passive outflow-like cut
  - a reflective wall
  - an analytic continuation of a larger spherical field
- do not describe the result as a pole-boundary success or failure

## Analytical Continuation Boundaries

Use analytical continuation when the computational shell is only a subdomain of a larger known analytic or semi-analytic state.

Best fit cases:

- Parker wind backgrounds
- analytic dipoles or bipolar fields
- diagnostic off-pole shells cut out of a larger spherical solution

Preferred rule:

- build ghost cells from the same underlying state used for initialization

That usually means reconstructing, in the ghost region:

- density
- pressure
- velocity
- magnetic field

Do not mix:

- zero-gradient fluid ghosts
- extrapolated magnetic ghosts
- analytic overwrite of only part of the field representation

If the boundary is supposed to represent continuation of the initialized state, all coupled quantities should come from that same continuation model.

## Boundary Design Rules

1. Decide what the boundary physically represents:
   - analytic continuation
   - zero-gradient outflow
   - reflective wall
   - fixed inflow or anchoring
2. Apply the same physical idea to all coupled quantities:
   - density
   - pressure
   - velocity
   - magnetic field
3. In CT runs, make cell-centered and face-centered magnetic fields consistent with the same boundary model.
4. If a boundary is only a diagnostic cut through a larger analytic field, prefer continuation of the analytic state over ad hoc extrapolation.

## Ghost-Cell Consistency

For spherical MHD tests, check three layers of consistency:

- fluid state in ghost cells
- cell-centered magnetic field in ghost cells
- face-centered magnetic field used by CT

If these represent different physical boundaries, the solver often responds with localized artifacts rather than clean global evolution.

## Magnetic Handling: CT Versus Non-CT

### Non-CT or cell-centered magnetic evolution

Use this branch when the run does not rely on staggered face-centered magnetic fields.

Typical guidance:

- build consistent ghost cells for the fluid state
- fill cell-centered magnetic components according to the intended boundary model
- if using analytic continuation, evaluate the analytic field directly in ghost cells rather than mixing extrapolation and overwrite

AMRVAC test anchor:

- [`tests/mhd/CAKwind_Magnetosphere_spherical_2D/mod_usr.t`](/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1/tests/mhd/CAKwind_Magnetosphere_spherical_2D/mod_usr.t) explicitly stops when `typedivbfix == 'ct'`, so its `special_bound` is a useful example of a spherical special boundary in a non-CT branch.

### CT or staggered magnetic evolution

Use this branch when magnetic fields live both in cell centers and on faces.

Key requirement:

- the face-centered magnetic boundary must represent the same physical boundary as the cell-centered magnetic field and the fluid ghosts

AMRVAC test anchor:

- [`tests/mhd/solar_atmosphere_3D/mod_usr.t`](/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1/tests/mhd/solar_atmosphere_3D/mod_usr.t) shows a staggered-grid `specialbound_usr` pattern where tangential face fields are extrapolated and the normal face field is reconstructed through `divB` control before `mhd_face_to_center`.

#### CT face-field consistency steps

1. Decide which magnetic quantity is normal to the boundary face.
2. Fill tangential face components in a way consistent with the intended physical boundary.
3. Reconstruct or constrain the normal face component so the discrete magnetic field remains compatible with the solver's `divB` handling.
4. Update the cell-centered magnetic field from the face representation.
5. Check whether the resulting ghost magnetic field still matches the intended boundary physics.

CT cautions:

- Avoid mixing:
  - fluid zero-gradient ghosts
  - magnetic extrapolation ghosts
  - analytic overwrite of only cell-centered `B`
- If possible, evaluate the analytic magnetic field at the geometrically correct face locations for staggered storage.
- Filling face fields from cell-centered values is only an approximation and can create localized magnetic-pressure errors at the boundary.
- A boundary that is numerically stable may still be physically inconsistent.

#### Non-CT consistency steps

1. Build fluid ghosts from the intended boundary model.
2. Fill cell-centered magnetic components from the same boundary model.
3. Check whether pressure, velocity, and magnetic geometry all imply the same physical boundary.
4. Compare initialization and ghost states directly before advancing the run.

## Diagnostic Patterns

### Pattern: crash at the first step near a polar ring

Likely causes:

- true-pole regularity issue
- inner radial boundary and pole treatment inconsistent with initialization
- magnetic field and fluid states not compatible near the pole

### Pattern: off-pole case runs, but `Vt` and `Vp` grow mainly at theta cuts

Likely cause:

- theta boundary is injecting nonphysical tangential motion

### Pattern: pressure ring at the first radial shell and both theta cuts

Likely cause:

- boundary inconsistency involving theta cuts plus magnetic-boundary treatment

Interpretation split:

- in CT runs, first suspect face-field inconsistency
- in non-CT runs, first suspect ghost-cell model mismatch or mixed analytic-versus-extrapolated `B`

### Pattern: magnetic field barely changes but flow develops strong local artifacts

Likely cause:

- the boundary is driving a local transient rather than a meaningful MHD relaxation

## Recommended Debug Order

1. Confirm the theta extent and whether the poles are actually included.
2. Confirm the boundary types in the `.par` file.
3. Classify the run:
   - full-pole or off-pole
   - CT or non-CT
   - analytic continuation or non-analytic boundary
4. Read the exact `specialbound_usr` cases used by the active boundaries.
5. Compare the boundary ghost state with the initialization state.
6. Check whether anomalies cluster:
   - near theta cuts
   - near the inner radial shell
   - near the pole ring
7. Only after that change the boundary implementation.

## What Good Looks Like

A better spherical boundary treatment should reduce:

- concentrated `Vt/Vp` at theta cuts
- negative `Vr` patches adjacent to the inner boundary
- narrow pressure rings at cut boundaries

without needing small-value repair to hide the problem.

## Common Mistakes

- using `pole` reasoning on a trimmed-theta domain
- treating AMRVAC pole-geometry messages as proof that the wrong boundary type is active
- fixing only cell-centered `B` while leaving CT face fields inconsistent
- judging success only by “the run finished”
- changing the theta boundary and inner radial boundary simultaneously, making the source of improvement unclear
