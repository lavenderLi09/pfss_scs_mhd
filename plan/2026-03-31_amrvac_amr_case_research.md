# AMRVAC AMR Case Research for the Bipolar Spherical CT Run

## Scope

This note compares three official AMRVAC test cases against the current
`amrvac_polytropic/real_bipo_test` spherical bipolar run, with the specific goal
of understanding why the present `jb_eta + theta-cut` AMR logic still develops
NaNs when AMR is active.

The three reference cases are:

1. `tests/mhd/icarus`
2. `tests/mhd/potential_field_source_surface_3D`
3. `tests/hd/blast_wave_spherical_stretched_3D`

The most relevant questions are:

1. What AMR strategy does each case use?
2. How much of the AMR decision is left to default Lohner refinement, and how
   much is fully user-controlled?
3. Do the closest official spherical large-domain cases use `typedivbfix='ct'`?
4. What does this imply for the current bipolar spherical CT test?

## AMRVAC AMR Semantics Relevant Here

From the local AMRVAC 3.1 documentation:

- `refine_criterion=3` means standard Lohner-type error estimation based on a
  weighted subset of conserved or auxiliary variables, followed by an optional
  user hook `usr_refine_grid` for extra refinement or derefinement.
- `ditregrid` reconstructs the AMR grid only once every `ditregrid` iterations.
- `usr_refine_grid` is called after the default estimator, so user logic can
  augment or suppress the default refine/coarsen decision.

Relevant local references:

- `doc/par.md`
- `doc/amrstructure.md`
- `src/amr/mod_errest.t`

Important implementation detail from `src/amr/mod_errest.t`:

- `my_coarsen = -1` suppresses coarsening.
- `my_refine = -1` suppresses refinement.
- `my_refine = 1` forces refinement if the current block is below
  `refine_max_level`.

This matters directly for the current theta-cut protection logic, because
`refine=-1, coarsen=-1` really does freeze AMR level changes in the protected
zone.

### A Note on `refine_criterion=0` vs `1`

The local 3.1 tree contains an inconsistency:

- `src/io/mod_input_output.t` prints `refine_criterion=0` as `user defined`.
- `src/amr/mod_errest.t` comments that `case (1)` is the all-user-refinement
  path.
- `doc/par.md` documents `refine_criterion=1` as the user-only case.

So for any production change that intends to switch fully away from default
Lohner AMR, the actual executable should be checked at runtime via the startup
message `Refine estimation: ...`, rather than relying only on the docs.

## Case 1: Icarus

### Geometry and Numerics

`icarus` is the closest official case to the present run in terms of geometry
and scale:

- `spherical_3D`
- large radial range
- radial stretching
- MHD
- active AMR

Key references:

- `tests/mhd/icarus/mod_usr.t`
- `tests/mhd/icarus/amrvac.par`

Important settings:

- `typedivbfix = 'linde'`
- `stretch_dim(1)='uni'`
- `refine_criterion=0`
- `refine_max_level=1`
- `xprobmin1=21.5`, `xprobmax1=432.5`
- `xprobmin2=8.146...e-2`, `xprobmax2=4.185...e-1`

So this is a large, trimmed-theta, stretched spherical shell MHD run, but it
does **not** use CT.

### AMR Logic

`icarus` does not rely on default variable-weighted Lohner refinement. Instead,
it wires in a custom `usr_refine_grid => specialrefine_grid` and then chooses
the refinement behavior through `amr_criterion` in `&icarus_list`.

The default example uses:

- `amr_criterion = 'tracing'`

The custom refinement logic then does:

- refine when the tracer exceeds a threshold in the block
- forbid coarsening in that case
- otherwise explicitly disable refine and request coarsen

In the `shock` mode it uses `div(V)` instead, again through a fully custom
block-level condition.

### Interpretation

The key design choice is not the tracer itself. The key design choice is this:

`icarus` uses a large stretched spherical shell with custom boundaries, but it
does **not** mix that setup with generic Lohner refinement. It uses
problem-aware AMR logic.

That is the strongest official precedent relevant to the current bipolar run.

## Case 2: PFSS (`potential_field_source_surface_3D`)

### Geometry and Numerics

This is the closest official case to a global PFSS-like spherical field setup.

Key references:

- `tests/mhd/potential_field_source_surface_3D/mod_usr.t`
- `tests/mhd/potential_field_source_surface_3D/amrvac.par`

Important settings:

- `set_coordinate_system('spherical_3D')`
- `typedivbfix='linde'`
- `refine_criterion=3`
- `refine_max_level=1`
- `ditregrid=2`
- `xprobmin1=1.0`, `xprobmax1=2.5`
- full pole handling in theta: `typeboundary_min2='pole'`, `typeboundary_max2='pole'`

Notable differences from the current bipolar run:

- no CT
- no radial stretching
- much smaller radial span
- only one AMR level
- pole geometry, not off-pole theta cuts

### AMR Logic

PFSS uses standard Lohner AMR, not a custom geometric/physics hook:

- `refine_criterion=3`
- `refine_threshold=20*0.2`
- `derefine_ratio=20*0.1`
- `w_refine_weight(5)=0.4`
- `w_refine_weight(6)=0.3`
- `w_refine_weight(7)=0.3`

Inference:

Because this is MHD with `mhd_energy=.false.`, the weighted components are very
likely the magnetic field variables. So PFSS is using the default estimator to
track magnetic structure, not a hand-written AMR map.

Also important:

- no `usr_refine_grid`
- no `specialthreshold`
- no explicit no-coarsen or freeze region

### Interpretation

PFSS shows that the official codebase is comfortable using default Lohner AMR
for spherical magnetic problems when the setup is relatively controlled:

- no CT
- no stretched radial shell
- no trimmed-theta cut surfaces
- shallow AMR

This is not the same risk profile as the current bipolar off-pole CT run.

## Case 3: `blast_wave_spherical_stretched_3D`

### Geometry and Numerics

This case gives an official example of spherical coordinates plus radial
stretching plus AMR, but it is HD, not MHD.

Key references:

- `tests/hd/blast_wave_spherical_stretched_3D/amrvac.par`
- `tests/hd/blast_wave_spherical_stretched_3D/mod_usr.t`

Important settings:

- spherical geometry
- `stretch_dim(1)='uni'`
- `refine_criterion=3`
- `refine_max_level=3`
- no MHD
- therefore no `typedivbfix='ct'`

### AMR Logic

This is a plain standard-Lohner case:

- no `usr_refine_grid`
- no `specialthreshold`
- no custom no-coarsen region
- no boundary-aware AMR handling

### Interpretation

This case proves that spherical stretching and default AMR can coexist
numerically in a relatively clean HD setup. It does **not** provide evidence
that the same approach should remain robust once CT magnetic fields and
sensitive theta boundaries are introduced.

## Cross-Case Comparison

### What official cases actually do

1. Large spherical MHD with stretching and special boundaries:
   `icarus`
   - uses `linde`
   - uses user-controlled AMR
   - does not mix the hardest geometry with generic Lohner AMR

2. Spherical PFSS-like magnetic configuration:
   `potential_field_source_surface_3D`
   - uses `linde`
   - uses shallow default AMR
   - no CT, no radial stretching, no off-pole theta cut

3. Spherical stretched AMR example:
   `blast_wave_spherical_stretched_3D`
   - uses default Lohner AMR
   - but only for HD

### What official cases do not show

There is no close official template in this local tree for:

- `spherical_3D`
- radial stretching
- off-pole theta cuts
- `typedivbfix='ct'`
- active AMR

all at once.

That absence matters. It means the current bipolar run is already in a numerically
harder combination than the official examples directly cover.

## Comparison to the Current Bipolar Run

Current AMR structure in the bipolar test:

- base AMR path is still `refine_criterion=3`, i.e. default Lohner AMR
- user hook adds `force_inner` and `force_sheet`
- `force_sheet` is based on blockwise `jb_eta_max` and `jb_eta_mean`
- theta-cut region is now partially or fully protected by:
  - `coarsen=-1` in the first guarded version
  - `refine=-1, coarsen=-1` in the stronger freeze version
- `specialthreshold` still relaxes the outer radial region

This means the current design is still a hybrid:

default Lohner AMR + custom sheet forcing + boundary freeze patches

That is different from both official spherical precedents:

- unlike `PFSS`, the geometry and divergence treatment are harder
- unlike `icarus`, AMR is not fully problem-defined

## What the Existing Failure Pattern Suggests

From the recent runs:

- `rgfreeze` survives
- `regrid4` crashes at about `it=77`
- the `regrid8` run crashes earlier, around `it=67`
- NaNs appear in all core variables almost simultaneously
- in `regrid4`, NaN appears immediately after a block-state switch
- in `regrid8`, NaN appears a few steps after a switch, not exactly at the same
  iteration

Interpretation:

1. The root problem is still AMR/regrid related, not basic time stepping.
2. The failure is not only "regrid too often".
3. A dangerous AMR topology can be created by a regrid event, then amplified by
   subsequent steps even if the next step does not regrid again.
4. The hardest zone is still the theta-cut neighborhood combined with CT and
   spherical geometry.

## Recommended Improvements for the Current Run

### Priority 1: Move away from hybrid AMR toward fully user-controlled AMR

This is the most important lesson from `icarus`.

For the present bipolar case, the official precedent argues for reducing the
role of generic Lohner AMR, not increasing it. A cleaner design would be:

1. Use only user-controlled AMR for this case
2. Let `jb_eta` decide refinement globally
3. Let geometry rules decide where coarsening or any level change is forbidden

In other words:

- sheet/current structure should be handled by `jb_eta`
- boundary-hazard regions should be handled by geometry
- default Lohner refinement should not still be competing in the background

Because of the local 3.1 inconsistency around `refine_criterion=0/1`, this
should be enabled only after checking what the executable prints at startup.

### Priority 2: Keep the theta-cut zone topology fixed, not just hard to coarsen

The current stronger test version with:

- `refine=-1`
- `coarsen=-1`

inside a theta protection band is the correct next experiment.

Reason:

the failure pattern suggests that the dangerous quantity is not only local
resolution, but local AMR topology changes near a CT-sensitive cut surface.

For the current case, "local rgfreeze" near theta cuts is more physically
justified than global `rgfreeze`, because it isolates the numerically fragile
boundary while preserving AMR elsewhere.

### Priority 3: Add hysteresis to the `jb_eta` logic

The present `force_sheet` uses:

- one threshold for `jb_eta_max`
- one threshold for `jb_eta_mean`

This can still create block-map oscillation. A more stable strategy would use
different refine and derefine conditions, for example:

- refine when `jb_eta_max >= A` and `jb_eta_mean >= B`
- allow derefine only when `jb_eta_max < A_low` and `jb_eta_mean < B_low`

with `A_low < A` and `B_low < B`

This follows the same logic as AMRVAC's own `refine_threshold` /
`derefine_ratio` hysteresis, but applied directly to the user physics
criterion.

### Priority 4: Protect the inner-shell plus theta-cut corner more explicitly

The most fragile region is likely not the whole theta-cut strip equally, but
the combined region:

- near the inner radial boundary
- near the theta cuts
- where CT face fields and special spherical boundary logic are both active

So the protection mask should likely be strongest there, for example:

- larger freeze zone for `r < r_corner`
- smaller freeze zone farther out

This would preserve more freedom for AMR in the outer shell while still
protecting the geometrically hardest corner.

### Priority 5: Use standard Lohner AMR only as a secondary diagnostic, not as the main driver

If `refine_criterion=3` is kept at all, its role should be minimal:

- low weights
- no theta-cut relaxation
- no expectation that it will choose stable geometry on its own

The official cases do not provide evidence that generic Lohner AMR is robust
enough for a stretched off-pole spherical CT problem with custom boundaries.

### Priority 6: Consider non-CT divergence control only as a diagnostic branch

The closest official spherical large-domain MHD examples (`icarus`, `PFSS`)
use `linde`, not `ct`.

This does not mean the present production route should abandon CT immediately.
But it does mean:

- if CT keeps failing after AMR topology is stabilized,
- a comparison run with `linde` or `glm`

would be a valuable diagnostic to isolate whether the dominant fragility is:

1. AMR topology alone
2. CT plus AMR topology

This should be treated as a diagnostic branch, not the first-line fix.

## Recommended Debug Ladder

1. Run the stronger local-theta-freeze version:
   - theta protection band uses `refine=-1, coarsen=-1`
   - keep everything else unchanged

2. If that stabilizes the run:
   - conclude that boundary-adjacent topology change is the main trigger
   - then migrate toward fully user-controlled AMR

3. If that still fails:
   - switch the case to user-only AMR
   - keep only `jb_eta` plus geometry freeze logic

4. If that still fails:
   - widen the freeze region near the inner-shell/theta-cut corner

5. If that still fails:
   - run a diagnostic non-CT comparison branch

## Bottom Line

The official cases suggest a clear lesson:

- When AMRVAC runs a large stretched spherical MHD shell with custom geometry,
  it tends to use problem-aware AMR and non-CT divergence control.
- When AMRVAC uses default Lohner AMR, the geometry is simpler, or the physics
  is less fragile than the current bipolar spherical CT setup.

So the current bipolar AMR problem is unlikely to be solved by only tuning
`ditregrid` or slightly shifting thresholds. The stronger likely direction is:

1. freeze theta-cut-adjacent AMR topology
2. reduce or remove dependence on generic Lohner AMR
3. let `jb_eta` become the real global physics criterion
4. treat CT plus off-pole theta cuts as a special hazard zone in the AMR design

