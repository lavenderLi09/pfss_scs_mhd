# Rho Bottleneck Diagnosis Outputs

`run_rho_bottleneck_diagnosis.py` writes three files into the chosen output directory and
reserves one optional plot path:

## 1. `rho_bottleneck_diagnostic_per_frame.csv`

One row per analyzed frame. The current minimal pipeline includes:

- `frame`, `it`, `time`, `dt_pred_min`
- `rho`, `p`, `v_abs`, `b_abs`, `v_A`
- `rho_ratio_to_initial`
- `shell_percentile`
- `flux_jump`
- `persistence_frames`
- `score`, `S1`, `S2`, `S3`, `S4`
- `gate_failed`, `final_label`

When `--frames` is provided, only the emitted rows are filtered. The normalization terms
(`rho_ratio_to_initial`, `shell_percentile`, `persistence_frames`, and the frame-to-frame
`flux_jump`) are still computed from the full merged timeline before filtering.

## 2. `rho_bottleneck_gate_report.csv`

One row per analyzed frame with the numerical gate fields returned by
`evaluate_numerical_gate`:

- `neg_rho_or_p`
- `nan_inf`
- `divb_flag`
- `flux_residual_flag`
- `gate_failed`
- `divb_median`
- `divb_p999`

## 3. `rho_bottleneck_diagnosis_summary.md`

A short markdown summary with:

- the number of analyzed frames
- the final-label distribution
- the latest frame id, iteration, time, and final label
- the deterministic assumptions used by the minimal pipeline

## 4. `rho_bottleneck_trend.png` (optional)

The current Task 5 pipeline does not generate plots, but `output_paths()` exposes this path
for callers that want a stable location for a future trend figure.

## Current deterministic assumptions

- Rank selection uses only `rank == 1` rows from `dt_cfl_topcells.csv`.
- `flux_jump` is the frame-to-frame relative change in rank-1 `b_abs` over the full merged
  timeline.
- `shell_percentile` is `100 * shell_index / max(shell_index)` over the full merged timeline.
- If merged rows do not contain explicit `divB` or flux-residual diagnostics, the pipeline
  passes non-finite values into `evaluate_numerical_gate` as a conservative safe default.
- That safe default intentionally forces a numerical gate failure, so missing diagnostics
  cannot yield a physical or mixed final label.
