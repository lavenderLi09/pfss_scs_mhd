> Copied to project knowledge from case: `off_2270_initwind`
> Original file: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_initwind/analysis/hpc_off0000_vs_initcheck_report_2026-04-08.md`
> Copied on: 2026-04-20# HPC `off0000.dat` vs Local `init_check0000.dat` Detailed Comparison (2026-04-08)

## 1) Data sources
- HPC reference (downloaded read-only):
  - `off0000.dat`: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/hpc_ref/off0000.dat`
  - `amrvac.par`: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/hpc_ref/amrvac.par`
- Local new initialization:
  - `init_check0000.dat`: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/data/init_check0000.dat`
  - `init_check.par`: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/init_check.par`

## 2) AMR and file-size differences
- File size:
  - HPC: `1,340,934,000` bytes (`1.249 GiB`)
  - Local: `796,375,072` bytes (`0.742 GiB`)
  - Difference: `544,558,928` bytes (`519.33 MiB`)
- AMR header:
  - HPC: `levmax=2`, `nleafs=3637`
  - Local: `levmax=1`, `nleafs=2160`
- Level distribution:
  - HPC: `L1=1949`, `L2=1688`
  - Local: `L1=2160`, `L2=0`

Why the size differs:
- Each extra leaf block (payload) is roughly:
  - `block_nx1*block_nx2*block_nx3*nw*8 + 2*ndim*4 = 368,664` bytes
- `nleaf` difference is `1477`, which explains:
  - `1477 * 368,664 = 544,516,728` bytes
- This matches almost all observed size gap; residual `42,200` bytes is metadata/tree overhead.

## 3) Why HPC `off0000` already has AMR at `it=0`
From downloaded HPC `amrvac.par`, `&filelist` contains:
- `restart_from_file='./data/off0115.dat'`
- `snapshotnext = 116`

So that `off0000` is not a pure “fresh t=0 from analytic init only” snapshot; it is consistent with a restart/reset workflow where evolved AMR structure can already exist at output index `0000`.

## 4) Initial rho / Vr / p comparison
All values below are code units (`Vr` also converted to km/s with `unit_velocity=116.4488468 km/s`).

### Global stats
- `rho`
  - HPC min/max/mean: `5.262946e-06 / 9.955690e-02 / 8.422156e-03`
  - Local min/max/mean: `1.049778e-05 / 9.964800e-02 / 1.605877e-02`
  - Mean ratio (local/HPC): `1.9067`
- `Vr`
  - HPC min/max/mean: `5.451279e-02 / 2.592926e+00 / 1.181939e+00`
  - Local min/max/mean: `1.392941e-01 / 3.337379e+00 / 1.319209e+00`
  - Mean ratio (local/HPC): `1.1161`
  - In km/s:
    - HPC mean: `137.64 km/s`
    - Local mean: `153.62 km/s`
- `p` (from conservative `e` with `p=(gamma-1)e`, `gamma=1.05`)
  - HPC min/max/mean: `1.057852e-05 / 2.001094e-01 / 1.692853e-02`
  - Local min/max/mean: `1.334592e-05 / 2.002572e-01 / 3.105753e-02`
  - Mean ratio (local/HPC): `1.8346`

### Radial sample (mean per r-bin)
- Near `r~1.05 Rs`:
  - HPC: `rho=5.8269e-02`, `Vr=9.99 km/s`, `p=1.1712e-01`
  - Local: `rho=6.3944e-02`, `Vr=22.40 km/s`, `p=1.2599e-01`
- Near `r~5 Rs`:
  - HPC: `rho=1.6637e-04`, `Vr=157.81 km/s`, `p=3.3440e-04`
  - Local: `rho=3.1209e-04`, `Vr=215.63 km/s`, `p=4.7011e-04`
- Near `r~19.5 Rs`:
  - HPC: `rho=5.6206e-06`, `Vr=299.39 km/s`, `p=1.1297e-05`
  - Local: `rho=1.1154e-05`, `Vr=385.84 km/s`, `p=1.4223e-05`

## 5) Is the current local initialization reasonable?
Numerically and physically as an initialization, yes:
- Conservative checks (`init_check0000.dat`): no nonfinite, no non-positive `rho/e`.
- Primitive checks (`init_check0000.vtu`): no nonfinite, no non-positive `rho/p`, outward `Vr` everywhere.
- Radial trends: `rho` monotonically decreases and `Vr` outward fraction is `1.0` (see `init_atmosphere_summary.md`).

Interpretation caveat:
- HPC `off0000` is a restart-context snapshot with strong AMR refinement (`L2=1688`), not the same kind of “fresh init-only 0000” as local `init_check0000`.
- Therefore, absolute differences (especially size and AMR pattern) are expected and do not imply the local init is wrong by itself.

## 6) Reproducibility artifacts
- Summary CSV: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/hpc_vs_initcheck_summary.csv`
- Local init sanity report: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind/analysis/init_atmosphere_summary.md`
