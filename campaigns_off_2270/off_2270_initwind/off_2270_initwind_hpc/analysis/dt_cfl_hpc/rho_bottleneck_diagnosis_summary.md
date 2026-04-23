# Rho Bottleneck Diagnosis Summary

- Frames analyzed: 28
- Assumptions:
  - flux_jump is the frame-to-frame relative change in rank-1 b_abs over the full merged timeline.
  - shell_percentile is 100 * shell_index / max(shell_index) over the full merged timeline.
  - divB and flux-residual are computed from DAT when available.
  - if DAT diagnostics are unavailable for a frame, non-finite values are passed to the numerical gate as a conservative safe default.

## Final Label Distribution
- MIXED_UNCERTAIN: 3
- NUMERICAL_DOMINANT: 25

## Latest Frame
- frame=amr_probe0027 it=86372 time=0.27
- Latest frame label: MIXED_UNCERTAIN
