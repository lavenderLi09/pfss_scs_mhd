## dt Physics Diagnostics

This document is the running diagnostic note for per-snapshot timestep-control analysis in
`[off_2270](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270)`.

It is meant to be updated frame by frame as new snapshot analyses are run.

Snapshot interpretation used in this analysis:

- The restart `.dat` payload is treated as **conservative**, not primitive.
- For this case, `[amrvac.par](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/amrvac.par)` sets `mhd_internal_e=.true.`, so:
  - `m1,m2,m3` are momenta and are converted by `v_i = m_i / rho`
  - `e` is internal energy density and is converted by `p = (\gamma-1)e`
- This is the correct `off_2270` interpretation for all results below.

Primary analysis artifacts:

- `[dt_bottleneck_summary.csv](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/dt_bottleneck_summary.csv)`
- `[dt_bottleneck_topcells.csv](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/dt_bottleneck_topcells.csv)`
- `[physical_speed_summary.csv](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/physical_speed_summary.csv)`
- `[analyze_snapshot_dt_bottleneck.py](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/hao_code/analyze_snapshot_dt_bottleneck.py)`

### Parameter Glossary

`|v|`

- Definition: `sqrt(v_r^2 + v_theta^2 + v_phi^2)`.
- Meaning: bulk flow speed of the plasma.
- CFL role: contributes to the advective part of the timestep limit. Large `|v|` means material motion alone can force smaller `dt`.

`c_s`

- Definition: `sqrt(gamma * p / rho)`.
- Meaning: sound speed, the propagation speed of thermal pressure disturbances.
- CFL role: if `c_s` is large, pressure waves travel quickly and can dominate the timestep limit even when bulk flow is slow.

`v_A`

- Definition: `|B| / sqrt(rho)`.
- Meaning: Alfvén speed, the characteristic propagation speed of magnetic-tension disturbances.
- CFL role: large `v_A` usually points to strong-field and/or low-density regions, which are common magnetic bottleneck candidates in MHD.

`c_f,max`

- Definition: the maximum of the directional fast magnetosonic speeds.
- Meaning: the strongest local MHD wave speed after combining thermal and magnetic effects.
- CFL role: often closer than either `|v|` or `c_s` alone to what the Riemann solver actually uses as a characteristic speed bound.

### 2026-04-03: `off0000`

Input:

- `[off0000.dat](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/data/off0000.dat)`

Relevant outputs:

- `[off0000_cfl_breakdown.png](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/off0000_cfl_breakdown.png)`
- `[off0000_speed_distribution.png](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/off0000_speed_distribution.png)`

Bottleneck location:

- `r ≈ 1.00064`
- `theta ≈ 0.00260`
- `phi ≈ 0.32031`
- `level = 1`
- `block_id = 37`

Predicted timestep bottleneck:

- `dt_pred ≈ 1.046e-06`
- Local `dt_log` is unavailable for this frame because the local `off.log` does not include the `it=0` entry.

Directional CFL split at the worst cell:

- `C_r / C_sum ≈ 1.028%`
- `C_th / C_sum ≈ 0.257%`
- `C_ph / C_sum ≈ 98.715%`

Interpretation:

- The bottleneck is overwhelmingly `phi`-direction controlled.
- The current classifier marks it as `geometry`, because near very small `theta`, the metric length `r sin(theta) dphi` becomes extremely small.
- In other words, this frame is not currently limited by bulk flow speed, and not primarily by a pressure-wave spike either; the angular cell geometry near the lower-theta boundary appears to be the first-order limiter.

Physical speed distribution summary for this frame:

- `|v|`: `p50 ≈ 1.124`, `p99 ≈ 2.561`
- `c_s`: `p50 = p99 ≈ 1.453`
- `v_A`: `p50 ≈ 1.80`, `p99 ≈ 14.37`, `max ≈ 49.03`
- `c_f,max`: `p50 ≈ 2.305`, `p99 ≈ 14.44`, `max ≈ 49.05`

Interpretation of these distributions:

- `|v|` is no longer negligible once the conservative momenta are converted correctly; the flow is dynamically relevant, even though it still does not beat the angular geometry term at the worst cell.
- `c_s` is almost spatially constant, which is consistent with the internal-energy formulation and the intended nearly isothermal background.
- `v_A` has a broad high-value tail, so magnetic structure and/or low density create physically fast regions that could matter in later snapshots.
- `c_f,max` still tracks the magnetic high tail, but after the correction its median is much lower than in the earlier, incorrect primitive interpretation.

Current takeaway:

- For `off0000`, geometry still dominates before the strongest magnetic wave-speed regions do.
- This does not rule out Alfvén/fast-speed control in later frames; it only says the initial snapshot is still controlled mainly by the near-boundary angular metric.

Next intended additions:

- `off0016`
- `off0015/off0016/off0017` continuity check
- comparison between bottleneck-cell physical speeds and whole-frame distributions

### 2026-04-03: `off0000` top-`|v|` location check

Relevant outputs:

- [`top_vabs_cells.csv`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/top_vabs_cells.csv)
- [`top_vabs_inner_cells.csv`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/top_vabs_inner_cells.csv)
- [`off0000_top_vabs_br_overlay.png`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/off0000_top_vabs_br_overlay.png)
- [`off0000_top_vabs_inner_br_overlay.png`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270/analysis/off0000_top_vabs_inner_br_overlay.png)

Global top-`|v|` result:

- If the whole domain is searched without any radial filter, the strongest `|v|` cells all come from the outermost shell:
  - `r ≈ 19.955`
  - `shell_index = 599`
  - `level = 2`
  - `block_id = 2106`
- This is expected for a solar-wind solution because the radial flow keeps accelerating outward.
- Therefore, the global top-`|v|` view is useful for the full-domain flow maximum, but it is **not** the right diagnostic for the small, near-surface hotspots seen in the rendered volume view.

Near-surface top-`|v|` result (`r <= 2`):

- The filtered near-surface top-`|v|` cells all collapse to:
  - `r ≈ 1.994`
  - `shell_index = 247`
  - `level = 1`
  - `block_id = 456`
- Their angular spread is compact:
  - `theta ≈ 0.0234 - 0.0807 rad`
  - `phi ≈ 0.0026 - 0.0807 rad`
- Typical values in this cluster are:
  - `|v| ≈ 0.4709`
  - `rho ≈ 2.90e-03`
  - `|B| ≈ 0.15 - 0.20`

Interpretation:

- The near-surface high-speed cells are not scattered randomly over the shell.
- They form a compact patch on a single inner shell, which is consistent with a localized magnetic-footpoint/funnel interpretation.
- In other words, there is already evidence that the small near-surface speed enhancements are tied to a specific region on the `Br(theta,phi)` map, rather than being domain-wide numerical speckle.

Current takeaway:

- Use the **global** top-`|v|` output to understand where the wind reaches its absolute maximum.
- Use the **near-surface** top-`|v|` output to diagnose the compact speed patches close to the inner boundary.
- For this case, the second diagnostic is the one that matches the visual features near the Sun.
