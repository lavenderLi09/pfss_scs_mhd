# Conversation & Implementation Summary — `I1_ppbound`

This document summarizes the technical discussion and the code changes made for the **positivity-preserving (PP) inner boundary** workflow, from notebook analysis through a dedicated AMRVAC case directory **`off_2270_highstream/I1_ppbound`**.

---

## 1. Original scientific goal (notebook)

### 1.1 Inner boundary diagnostics

Using:

- Magnetic field from `OFF_combined_lmax10_q1.008_nr600.bin` (Fortran unformatted stream with `int32` shape header + `float64` data, Fortran column-major order, with theta reversal consistent with `mod_usr.t`).
- Polytropic Parker solar wind density/pressure from `mod_usr.t` normalization and solver.

Compute at the inner boundary (approximately `r = 1 R☉` in normalized coordinates):

- Alfvén speed (CGS-consistent): \(v_A = B/\sqrt{4\pi\rho}\)
- Plasma beta: \(\beta = 8\pi p / B^2\)

Implementation target notebook:

- `2026_scs_relaxation/initial_density_pressure.ipynb`

### 1.2 Angular coordinate convention (critical)

The user clarified that in `amrvac.par`, the angular coordinates are scaled such that **`1` in normalized coordinates corresponds to `2π` radians physically**.

Therefore for a global simulation:

- \(\theta\) domain maps to physical \([0,\pi]\) when `xprobmax2 = 0.5` (because \(0.5 \times 2\pi = \pi\)).
- \(\phi\) domain maps to physical \([0,2\pi]\) when `xprobmax3 = 1.0`.

This affects how `theta`/`phi` arrays are built for plotting and indexing consistency with Fortran.

### 1.3 PP strategy (Wang et al. 2025, adapted)

The user wanted a PP-style adjustment at the inner boundary to:

- Increase density where local Alfvén speed is too high (density PP).
- Increase pressure where plasma beta is too low (pressure PP).

**Important design choice:** thresholds should be derived from **this simulation’s own distributions** (percentiles of `vA` and `beta`), not fixed literature values.

Density PP (CGS-consistent cap form used in notebook):

\[
\rho_\text{new} = \Upsilon_\rho \frac{B^2}{4\pi V_{A,\max}^2} + (1-\Upsilon_\rho)\rho_o
\]

Pressure PP:

\[
p_\text{new} = \Upsilon_p \frac{\beta_{\min} B^2}{8\pi} + (1-\Upsilon_p) p_o
\]

with hyperbolic tangent blending weights \(\Upsilon_\rho,\Upsilon_p\).

### 1.4 Matplotlib `LogNorm` error (user hit)

`ValueError: Invalid vmin or vmax` occurred because `LogNorm` requires finite `vmin>0` and `vmax>vmin`. Using `vmax=array.max()` when the masked array is all `NaN` (or `vmax <= vmin`) triggers invalid norms.

Fix direction: compute `vmax` only from finite positive values, or fall back to linear plotting if no positive values exist.

### 1.5 Contour “graininess” (visualization explanation)

Single-level contours on discrete grids can appear “speckled” when:

- The threshold cuts a high-gradient field into many small islands.
- Percentile thresholds isolate tail regions into scattered cells.
- \(\beta \propto 1/B^2\) amplifies small-scale magnetic structure.

This is usually not a numerical bug; it is a discrete contouring artifact.

---

## 2. Bridge goal: export PP-adjusted state for Fortran / AMRVAC

The user wanted `rho_new` and `p_new` from the notebook to be usable in Fortran similarly to scalar `rhob` in `mod_usr.t`.

Key normalization mapping (consistent with `mod_usr.t`):

- `rhob` in Fortran is **normalized number density**: \(n / \texttt{unit\_numberdensity}\).
- Notebook mass density `rho` in g/cm³ maps to the same dimensionless `rhob` via:

\[
\texttt{rhob\_map} = \frac{\rho_\text{new}}{\texttt{unit\_density}}
\]

because `unit_density = 1.4 m_H * unit_numberdensity`.

For pressure / temperature mapping at `r=1` in the existing polytropic Parker closure used in `mod_usr.t`:

At `r=1`, the original closure yields \(p_\text{th} = T_\text{iso}\cdot \rho_b\) in normalized units (since \(\rho/\rho_b = 1\)).

Therefore a local normalized “MK temperature factor” map can be defined as:

\[
\texttt{Tiso\_map} = \frac{p_\text{new}/\texttt{unit\_pressure}}{\texttt{rhob\_map}}
\]

This reproduces the modified pressure at `r=1` while keeping the same algebraic relationship as the scalar `Tiso` case.

### 2.1 Binary file format (notebook → Fortran)

Planned/produced file:

- `./initial/pp_inner_bc_rho_Tiso.bin`

Layout:

1. `int32(2)` header: `[nth, nph]`
2. `float64(nth,nph)` `rhob_map` in Fortran column-major order
3. `float64(nth,nph)` `Tiso_map` in Fortran column-major order

---

## 3. Fortran case: `off_2270_highstream/I1_ppbound`

### 3.1 What was requested

Copy:

- `amrvac_polytropic/off_2270_initwind/amrvac.par`
- `amrvac_polytropic/off_2270_initwind/mod_usr.t`

into:

- `amrvac_polytropic/off_2270_highstream/I1_ppbound/`

Then modify **only** the copied `mod_usr.t` to:

1. Read `rhob_map` and `Tiso_map`.
2. Use them in **initialization** and **inner boundary special BC (`iB=1`)**.
3. For each \((\theta,\phi)\), solve Parker wind using **local** `Tiso_map` (user selected “map Tiso”, not scalar-only).
4. Keep Parker density radial scaling consistent with mass conservation:

\[
\rho(r,\theta,\phi)=\frac{\rho_b(\theta,\phi)\,V_\text{surface}(\theta,\phi)}{v_r(r,\theta,\phi)\,r^2}
\]

User explicitly chose the “A” option earlier in discussion: **apply 2D maps at the base and still use Parker outward extrapolation in radius**.

### 3.2 Important coding constraint from the user

The codebase historically avoided certain primitive/conservative conversion calls because they could error.

Therefore the plan explicitly required:

- **Manual conservative writes** (`rho_`, `mom`, `p_` as internal energy when `mhd_internal_e=.true.`) consistent with existing style.

### 3.3 Debug printing requirement

Add **lightweight** diagnostics:

- Only `mype==0` prints for heavy diagnostics.
- On successful read: min/max of maps; NaN checks; positivity checks.
- A few sample points for consistency at `r=1`.
- Boundary (`iB=1`) prints once (first call) to avoid log spam.

---

## 4. Implemented changes (this repository)

### 4.1 Files created/updated

Created by copying:

- `amrvac_polytropic/off_2270_highstream/I1_ppbound/amrvac.par`
- `amrvac_polytropic/off_2270_highstream/I1_ppbound/mod_usr.t`

Then edited:

- `amrvac_polytropic/off_2270_highstream/I1_ppbound/mod_usr.t`

**Not edited:** the original `off_2270_initwind` sources (only the new copy).

### 4.2 New module-level state in `I1_ppbound/mod_usr.t`

Added:

- `real(8), allocatable :: rhob_map(:,:), Tiso_map(:,:)`
- `logical, save :: pp_bc_ready = .false.`

### 4.3 New reader

Added subroutine:

- `read_pp_inner_bc(filename, nth, nph, rhob_map, Tiso_map, ok)`

Called from `initglobaldata_usr` after `nth`/`nph` are known:

- Default path: `./initial/pp_inner_bc_rho_Tiso.bin`

If shape mismatch or read fails, `pp_bc_ready` remains false and the code falls back to the original scalar `rhob`/`Tiso` behavior.

### 4.4 New local Parker background builder

Added:

- `set_background_polytropic_state_pp(r_in, rhob_loc, Tiso_loc, rho_out, vr_out, pth_out, eint_out)`

Mechanics:

- Compute local surface speed:

  - `call parker_solar_wind_polytropic(1.0d0, mhd_gamma, Tiso_loc*unit_temperature, V_surface_loc)`

- Compute local radial speed at `r_in`:

  - `call parker_solar_wind_polytropic(r_in, mhd_gamma, Tiso_loc*unit_temperature, vr_out)`

- Density:

  - `rho_out = (rhob_loc*V_surface_loc)/(max(vr_out,1d-14)*max(r_in,1d-12)^2)`

- Thermal pressure + internal energy:

  - `pth_out = Tiso_loc*rhob_loc*(rho_out/rhob_loc)**mhd_gamma`
  - `eint_out = pth_out/(mhd_gamma-1)`

This mirrors the scalar `set_background_polytropic_state` structure but parameterized by local maps.

### 4.5 Initialization: `initonegrid_usr`

For each cell:

- Compute `(irr, ith, iph)` indices (same style as existing magnetic field indexing).
- If `pp_bc_ready`:

  - `call set_background_polytropic_state_pp(x(ix^D,1), rhob_map(ith,iph), Tiso_map(ith,iph), ...)`

Else:

  - fall back to `set_background_polytropic_state`

Then manually set:

- `w(rho_) = rho_bg`
- `w(mom(1)) = rho_bg*vr_bg`
- `w(mom(2:3)) = 0`
- `w(p_) = eint_bg`

Magnetic field assignment from `B_init` remains unchanged.

### 4.6 Boundary condition: `specialbound_usr`, `case(1)` (`iB=1`)

Previously, `case(1)` set `mom` to zero after filling density/pressure-like quantities, which prevents imposing a non-zero Parker radial momentum.

Updated `case(1)` loop to:

- compute clamped `(ith1,iph1)` indices per boundary point
- call `set_background_polytropic_state_pp` when `pp_bc_ready` else scalar routine
- set:

  - `w(rho_)`
  - `w(mom(1)) = rho*vr`
  - `w(mom(2:3)) = 0`
  - `w(p_) = eint`

Added a **one-time** debug print on first `iB=1` application (min/max of `rho`, `mom(1)`, `p_`).

**Explicitly unchanged:** `case(2)` outer boundary extrapolation logic (per plan).

---

## 5. Operational notes / pitfalls

### 5.1 The `./initial/` directory must exist at runtime

`mod_usr.t` reads:

- `./initial/pp_inner_bc_rho_Tiso.bin`

At the time this summary was written, `I1_ppbound/initial/` may not exist yet in the repository tree. You must either:

- create `amrvac_polytropic/off_2270_highstream/I1_ppbound/initial/` and place the bin there **before running** from that working directory, or
- adjust the filename path in `initglobaldata_usr` to match your run directory layout.

### 5.2 Indexing / theta reversal consistency

Magnetic field `B_init` is read and then reversed along the theta index in `read_initial_magnetic_field`.

Whether `rhob_map`/`Tiso_map` need an additional theta flip depends on how the notebook exported the maps relative to Fortran indexing.

**Rule of thumb:** if maps appear rotated/flipped relative to `B` in AMRVAC, apply the same theta reversal to the maps after read (or export them already aligned).

### 5.3 Performance / MPI considerations

`set_background_polytropic_state_pp` calls `parker_solar_wind_polytropic` twice per evaluation (surface + local radius). With per-cell usage in initialization and boundary loops, this can be expensive.

If this becomes a bottleneck, future optimization could precompute `V_surface_loc` on the `(nth,nph)` grid once at startup and reuse.

---

## 6. Related files elsewhere in the project (context)

Notebook work (not necessarily committed in the same change set as Fortran):

- `2026_scs_relaxation/initial_density_pressure.ipynb`

Original reference case:

- `amrvac_polytropic/off_2270_initwind/mod_usr.t`
- `amrvac_polytropic/off_2270_initwind/amrvac.par`

---

## 7. Suggested verification checklist

1. Confirm `pp_inner_bc_rho_Tiso.bin` exists at `./initial/pp_inner_bc_rho_Tiso.bin` **from the run cwd**.
2. Run with `mype==0` logging enabled; confirm:
   - `PP inner BC loaded ...`
   - sample lines show `rho(r=1) ≈ rhob_loc` (within floating noise)
3. Visual sanity: compare a theta–phi slice of `rho` at inner radius against notebook `rho_new` pattern (up to alignment flips).

---

## 8. What was *not* done in this pass (possible follow-ups)

- Automatic creation of `initial/` in the repo (runtime packaging choice).
- Automated compile/run regression test in CI (not requested).
- Optional periodic boundary debug prints keyed to `it` (currently only first-call boundary stats).

---

## 9. Meta

This `.md` is a narrative consolidation of the conversation goals, constraints, physics interpretation, debugging notes, and the concrete Fortran changes under:

- `amrvac_polytropic/off_2270_highstream/I1_ppbound/`

If you want this document to also embed exact code excerpts from `mod_usr.t`, say so and we can add a “Appendix: code excerpts” section with line-anchored citations.
