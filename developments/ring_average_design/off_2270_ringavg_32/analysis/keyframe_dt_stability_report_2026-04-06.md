# `off_2270` Key-Frame DT / Stability Report

## Scope

Key frames analyzed on HPC and copied back locally:

- `off0000`
- `off0016`
- `off0032`
- `off0064`
- `off0096`
- `off0112`

The corresponding CSV summary is:

- [`sample_evolution_summary.csv`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/sample_evolution_summary.csv)
- [`amr_level_overall.csv`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/amr_level_overall.csv)
- [`amr_theta_cap_summary.csv`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/amr_theta_cap_summary.csv)

## Geometry / Normalization Correction

For this case, the angular coordinates in `amrvac.par` are normalized by `2\pi`.

Therefore:

- `x2 = 0.0` corresponds to physical `theta = 0`
- `x2 = 0.5` corresponds to physical `theta = \pi`

So this run covers the full polar angle range `0 \le theta \le \pi`, not just a north-polar cap.

Also, this dataset only contains AMR levels `L1` and `L2`. In other words:

- `L1` is the base level
- `L2` is the refined level
- there is no separate `L0` in these outputs

## Core Findings

### 1. DT bottleneck stays in the same class of region

Across all six key frames, the dt bottleneck remains:

- very close to the inner boundary: `r_min ≈ 1.001`
- at very small theta: `theta_min ≈ 0.0013–0.0026`
- dominated by `phi`
- classified as `geometry`

Representative values:

- `off0000`: `dt_pred_min = 1.046e-06`, `r_min = 1.00064`, `theta_min = 0.00260`, `level = 1`
- `off0016`: `dt_pred_min = 2.536e-07`, `r_min = 1.00161`, `theta_min = 0.00130`, `level = 2`
- `off0112`: `dt_pred_min = 2.375e-07`, `r_min = 1.00161`, `theta_min = 0.00130`, `level = 2`

Interpretation:

- the global timestep bottleneck does **not** migrate to the high-speed region
- it stays tied to the small-`theta` inner-shell geometry, where `r sin(theta) dphi` is extremely small

### 2. The overall physical evolution looks stable

No evidence of runaway instability appears in the sampled frames.

Selected trends:

- `dt_log` stays in a narrow range after startup:
  - `off0016`: `1.089e-06`
  - `off0032`: `1.125e-06`
  - `off0064`: `1.107e-06`
  - `off0096`: `1.108e-06`
  - `off0112`: `1.162e-06`
- `c_s` remains very stable:
  - `c_s p50`: `1.453 -> 1.432`
  - `c_s p99`: `1.453 -> 1.488`
- `v_A` and `c_f,max` develop stronger high tails, but smoothly:
  - `v_A p99`: `14.30 -> 26.76`
  - `c_f,max p99`: `14.37 -> 26.80`

Interpretation:

- the solution is evolving physically, not blowing up numerically
- magnetic/fast-wave high tails strengthen over time, but in a controlled way

### 3. The near-surface high-speed patch is real and stable in angle

For `r <= 2`, the strongest `|v|` patch after startup stays concentrated near:

- `theta ≈ 0.2005 rad`
- `phi ≈ 0.0859–0.0911 rad`

Representative values:

- `off0016`: `|v| = 15.31`, `r = 1.377`, `rho = 7.49e-04`, `|B| = 1.35`
- `off0032`: `|v| = 12.88`, `r = 1.770`, `rho = 3.34e-04`, `|B| = 0.429`
- `off0112`: `|v| = 13.61`, `r = 1.756`, `rho = 3.11e-04`, `|B| = 0.409`

Interpretation:

- this is not a random full-domain artifact
- the strongest near-surface flow enhancement remains tied to one angular sector
- it looks more like a persistent local acceleration channel / flux-tube funnel

## Physical Interpretation

### Why the timestep bottleneck remains `phi + geometry`

The bottleneck cell sits near very small theta, so the azimuthal physical cell size

- `r sin(theta) dphi`

becomes tiny.

Even when magnetic and fast-mode speeds become large elsewhere, this small azimuthal metric factor keeps:

- `C_ph = (|v_phi| + c_f(phi)) / (r sin(theta) dphi)`

dominant.

So the timestep is still primarily controlled by spherical geometry near the inner boundary, not by the globally fastest flow.

### Why the near-surface fast patch is likely physical

The near-surface strongest-`|v|` region:

- is localized in one angular sector
- persists across multiple late frames
- is accompanied by low density and moderate/strong magnetic field

This is consistent with a local open-field / flux-tube acceleration channel rather than random noise.

## Numerical Interpretation

### AMR coverage over the full domain keeps increasing

The whole-domain AMR statistics show that `L2` expands steadily with time:

- `off0000`: `L2 = 46.41%` of all cells
- `off0016`: `L2 = 62.54%`
- `off0032`: `L2 = 65.41%`
- `off0064`: `L2 = 69.05%`
- `off0096`: `L2 = 70.78%`
- `off0112`: `L2 = 70.94%`

At the inner boundary surface, the refined fraction becomes even more dominant:

- `off0000`: `L2 inner-area fraction = 0.00%`
- `off0016`: `86.23%`
- `off0032`: `87.35%`
- `off0064`: `87.35%`
- `off0096`: `89.43%`
- `off0112`: `89.43%`

So this is not a case where only a tiny corner is refined. By late times, most of the domain is already at `L2`, and the inner boundary surface is overwhelmingly `L2`.

### Polar small-theta caps also gain stronger `L2` coverage

After correcting the angle normalization, the cap statistics should be read as physical `theta < 0.5^\circ`, `1^\circ`, `2^\circ`, and `5^\circ`.

The strongest result is at the very smallest cap:

- for `theta < 0.5^\circ`, once cells appear there, they are entirely `L2`
  - `off0016`: `100% L2`
  - `off0112`: `100% L2`

Broader caps also show systematic `L2` growth:

- `theta < 1^\circ`
  - `off0016`: `L2 = 10.26%` of cap cells
  - `off0112`: `L2 = 25.23%`
- `theta < 2^\circ`
  - `off0016`: `L2 = 18.60%`
  - `off0112`: `L2 = 40.29%`
- `theta < 5^\circ`
  - `off0000`: `L2 = 0%`
  - `off0016`: `L2 = 16.00%`
  - `off0032`: `L2 = 23.67%`
  - `off0064`: `L2 = 31.27%`
  - `off0096`: `L2 = 35.99%`
  - `off0112`: `L2 = 35.99%`

On the inner boundary area itself, the same trend appears:

- for `theta < 5^\circ`, `L2` already covers `58.15%` of cap area by `off0032`
- and stays at `58.15%` through `off0112`
- for `theta < 2^\circ`, `L2` covers `66.67%` of cap area from `off0032` onward

This means the timestep bottleneck is not just geometrically close to the pole; it is also increasingly kept on refined cells there.

### What looks numerical

- In `off0000`, angular `|v|` patchiness on shell maps should not be over-interpreted.
- Since `off0000` is a Parker-wind-like initial condition, absolute `|v|` should be nearly angle-uniform at fixed radius.
- Small patchy structure there is more likely due to shell selection / AMR binning / map projection than true physics.

### What does **not** currently look numerically dangerous

- `dt_log` does not collapse
- the bottleneck class does not hop around
- no sampled frame shows abrupt blow-up in `c_s`, `v_A`, or `c_f,max`

## Bottom Line

The six key frames support three main conclusions:

1. The timestep bottleneck is structurally stable and remains an inner-boundary, small-theta, `phi`-geometry limitation.
2. The global MHD solution appears numerically stable over the sampled interval.
3. The localized near-surface high-speed region is likely a real, persistent acceleration structure, not the same feature that limits the timestep.

## Practical AMR Recommendation

If the goal is to loosen the timestep, the most effective first lever is **not** to reduce resolution everywhere. It is to **limit `L2` coverage near the smallest-theta polar caps first**.

Why this is the best first target:

- the dt bottleneck is always found at very small `theta`
- it is always classified as `phi + geometry`
- the smallest polar caps become partially or fully `L2`
- so the numerically most dangerous region is also being kept on the finest angular spacing

In practical terms, the most targeted intervention would be:

1. first restrict or freeze `L2` inside the smallest polar caps
2. only consider broader/global derefinement if that is still insufficient

This should be more effective than a domain-wide resolution reduction, because the current bottleneck is highly localized in geometry rather than spread across the full volume.
