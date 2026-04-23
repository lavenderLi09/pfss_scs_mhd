# Copied to project analysis task toolbox
# Source case path: off_2270_initwind/analysis/check_init_atmosphere.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_initwind/analysis/check_init_atmosphere.py

#!/usr/bin/env python3
"""Check init_check atmosphere consistency from conservative DAT + primitive VTU."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

OUT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".mplconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = OUT_DIR.parent
HAO_CASE_DIR = CASE_DIR.parent / "off_2270"
if str(HAO_CASE_DIR) not in sys.path:
    sys.path.insert(0, str(HAO_CASE_DIR))

from hao_code.datfile_io import SIZE_INT, get_header, get_tree_info  # type: ignore  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate init_check atmosphere from DAT/VTU.")
    parser.add_argument("--case-dir", type=Path, default=CASE_DIR)
    parser.add_argument("--stem", default="init_check")
    parser.add_argument("--output-dir", default="analysis")
    return parser.parse_args()


def latest_match(case_dir: Path, pattern: str) -> Path:
    matches = sorted(case_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matches {pattern} in {case_dir}")
    return matches[-1]


def read_block_field(
    handle,
    offset: int,
    block_shape: np.ndarray,
    ndim: int,
    field_idx: int,
) -> np.ndarray:
    count = int(np.prod(block_shape))
    byte_size_field = count * 8
    handle.seek(offset + 2 * ndim * SIZE_INT + field_idx * byte_size_field)
    arr = np.fromfile(handle, dtype="=f8", count=count)
    if arr.size != count:
        raise IOError(f"Failed to read field {field_idx} at offset {offset}")
    arr = arr.reshape(tuple(block_shape[::-1]), order="C").T
    while arr.ndim < 3:
        arr = arr[..., np.newaxis]
    return arr


def scan_dat(dat_path: Path) -> dict[str, float]:
    with dat_path.open("rb") as handle:
        header = get_header(handle)
        block_lvls, _block_ixs, block_offsets = get_tree_info(handle)

        names = [name.strip() for name in header["w_names"]]
        try:
            rho_idx = names.index("rho")
            eint_idx = names.index("e")
        except ValueError as exc:
            raise RuntimeError(f"Expected conservative fields rho/e in DAT w_names={names}") from exc

        block_shape = header["block_nx"].astype(int)
        ndim = int(header["ndim"])

        rho_min = np.inf
        rho_max = -np.inf
        eint_min = np.inf
        eint_max = -np.inf
        rho_nonfinite = 0
        eint_nonfinite = 0
        rho_nonpos = 0
        eint_nonpos = 0
        n_cells = 0

        for _lvl, offset in zip(block_lvls.astype(int), block_offsets.astype(int)):
            rho = read_block_field(handle, int(offset), block_shape, ndim, rho_idx)
            eint = read_block_field(handle, int(offset), block_shape, ndim, eint_idx)

            n_cells += rho.size
            rho_finite = np.isfinite(rho)
            eint_finite = np.isfinite(eint)
            rho_nonfinite += int(rho.size - np.count_nonzero(rho_finite))
            eint_nonfinite += int(eint.size - np.count_nonzero(eint_finite))

            if np.any(rho_finite):
                rho_min = min(rho_min, float(np.min(rho[rho_finite])))
                rho_max = max(rho_max, float(np.max(rho[rho_finite])))
            if np.any(eint_finite):
                eint_min = min(eint_min, float(np.min(eint[eint_finite])))
                eint_max = max(eint_max, float(np.max(eint[eint_finite])))

            rho_nonpos += int(np.count_nonzero(rho_finite & (rho <= 0.0)))
            eint_nonpos += int(np.count_nonzero(eint_finite & (eint <= 0.0)))

    return {
        "dat_cells": float(n_cells),
        "dat_rho_min": float(rho_min),
        "dat_rho_max": float(rho_max),
        "dat_e_min": float(eint_min),
        "dat_e_max": float(eint_max),
        "dat_rho_nonfinite": float(rho_nonfinite),
        "dat_e_nonfinite": float(eint_nonfinite),
        "dat_rho_nonpos": float(rho_nonpos),
        "dat_e_nonpos": float(eint_nonpos),
    }


def read_vtu(vtu_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    grid = reader.GetOutput()

    centers = vtk.vtkCellCenters()
    centers.SetInputData(grid)
    centers.Update()
    pts = vtk_to_numpy(centers.GetOutput().GetPoints().GetData()).astype(float)
    radius = np.sqrt((pts[:, 0] ** 2) + (pts[:, 1] ** 2) + (pts[:, 2] ** 2))

    cell_data = grid.GetCellData()
    required = ["rho", "p", "Vr"]
    arrays: dict[str, np.ndarray] = {}
    for name in required:
        arr = cell_data.GetArray(name)
        if arr is None:
            raise RuntimeError(f"Missing VTU cell array '{name}' in {vtu_path}")
        arrays[name] = vtk_to_numpy(arr).astype(float)

    return radius, arrays["rho"], arrays["p"], arrays["Vr"]


def radial_median_profile(
    radius: np.ndarray,
    values: np.ndarray,
    nbins: int = 80,
) -> tuple[np.ndarray, np.ndarray]:
    rmin = float(np.nanmin(radius))
    rmax = float(np.nanmax(radius))
    edges = np.linspace(rmin, rmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(nbins, np.nan)

    ibin = np.digitize(radius, edges) - 1
    for i in range(nbins):
        mask = ibin == i
        if np.any(mask):
            sample = values[mask]
            sample = sample[np.isfinite(sample)]
            if sample.size:
                med[i] = float(np.median(sample))
    return centers, med


def profile_checks(rho_med: np.ndarray, vr_med: np.ndarray) -> tuple[float, float]:
    rho_valid = np.isfinite(rho_med)
    rho_inc_frac = np.nan
    if np.count_nonzero(rho_valid) > 2:
        drho = np.diff(rho_med[rho_valid])
        rho_inc_frac = float(np.count_nonzero(drho > 0.0) / max(drho.size, 1))

    vr_valid = np.isfinite(vr_med)
    vr_out_frac = np.nan
    if np.count_nonzero(vr_valid) > 0:
        vr_out_frac = float(np.count_nonzero(vr_med[vr_valid] > 0.0) / np.count_nonzero(vr_valid))

    return rho_inc_frac, vr_out_frac


def write_profile_csv(path: Path, r: np.ndarray, rho: np.ndarray, pth: np.ndarray, vr: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["r", "rho_median", "p_median", "Vr_median"])
        writer.writeheader()
        for i in range(r.size):
            writer.writerow(
                {
                    "r": f"{r[i]:.10e}",
                    "rho_median": f"{rho[i]:.10e}" if np.isfinite(rho[i]) else "nan",
                    "p_median": f"{pth[i]:.10e}" if np.isfinite(pth[i]) else "nan",
                    "Vr_median": f"{vr[i]:.10e}" if np.isfinite(vr[i]) else "nan",
                }
            )


def plot_profiles(path: Path, r: np.ndarray, rho: np.ndarray, pth: np.ndarray, vr: np.ndarray) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].set_title("init_check radial median profiles")
    axes[0].plot(r, rho, lw=1.4)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("rho")
    axes[0].grid(alpha=0.3)

    axes[1].plot(r, pth, lw=1.4)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("p")
    axes[1].grid(alpha=0.3)

    axes[2].plot(r, vr, lw=1.4)
    axes[2].axhline(0.0, color="k", ls="--", lw=0.8)
    axes[2].set_xlabel("r")
    axes[2].set_ylabel("Vr")
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    out_dir = (case_dir / args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    dat_path = latest_match(case_dir, f"data/{args.stem}*.dat")
    vtu_path = latest_match(case_dir, f"data/{args.stem}*.vtu")

    dat_stats = scan_dat(dat_path)

    radius, rho, pth, vr = read_vtu(vtu_path)
    vtu_stats = {
        "vtu_cells": float(rho.size),
        "vtu_rho_nonfinite": float(np.count_nonzero(~np.isfinite(rho))),
        "vtu_p_nonfinite": float(np.count_nonzero(~np.isfinite(pth))),
        "vtu_vr_nonfinite": float(np.count_nonzero(~np.isfinite(vr))),
        "vtu_rho_nonpos": float(np.count_nonzero(np.isfinite(rho) & (rho <= 0.0))),
        "vtu_p_nonpos": float(np.count_nonzero(np.isfinite(pth) & (pth <= 0.0))),
        "vtu_vr_min": float(np.nanmin(vr)),
        "vtu_vr_max": float(np.nanmax(vr)),
    }

    r_prof, rho_prof = radial_median_profile(radius, rho)
    _, p_prof = radial_median_profile(radius, pth)
    _, vr_prof = radial_median_profile(radius, vr)

    rho_inc_frac, vr_out_frac = profile_checks(rho_prof, vr_prof)

    csv_path = out_dir / "init_checks.csv"
    profile_csv_path = out_dir / "init_radial_profiles.csv"
    png_path = out_dir / "init_radial_profiles.png"
    md_path = out_dir / "init_atmosphere_summary.md"

    all_stats = {**dat_stats, **vtu_stats}
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_stats.keys()))
        writer.writeheader()
        writer.writerow(all_stats)

    write_profile_csv(profile_csv_path, r_prof, rho_prof, p_prof, vr_prof)
    plot_profiles(png_path, r_prof, rho_prof, p_prof, vr_prof)

    dat_ok = (
        dat_stats["dat_rho_nonfinite"] == 0.0
        and dat_stats["dat_e_nonfinite"] == 0.0
        and dat_stats["dat_rho_nonpos"] == 0.0
        and dat_stats["dat_e_nonpos"] == 0.0
    )
    vtu_ok = (
        vtu_stats["vtu_rho_nonfinite"] == 0.0
        and vtu_stats["vtu_p_nonfinite"] == 0.0
        and vtu_stats["vtu_vr_nonfinite"] == 0.0
        and vtu_stats["vtu_rho_nonpos"] == 0.0
        and vtu_stats["vtu_p_nonpos"] == 0.0
    )
    profile_ok = (
        np.isfinite(rho_inc_frac)
        and np.isfinite(vr_out_frac)
        and rho_inc_frac <= 0.25
        and vr_out_frac >= 0.80
    )

    summary_lines = [
        "# init_check Atmosphere Summary",
        "",
        f"- DAT file: `{dat_path}`",
        f"- VTU file: `{vtu_path}`",
        "",
        "## Conservative DAT checks",
        f"- rho min/max: {dat_stats['dat_rho_min']:.6e} / {dat_stats['dat_rho_max']:.6e}",
        f"- e min/max: {dat_stats['dat_e_min']:.6e} / {dat_stats['dat_e_max']:.6e}",
        f"- rho nonfinite/nonpos: {int(dat_stats['dat_rho_nonfinite'])} / {int(dat_stats['dat_rho_nonpos'])}",
        f"- e nonfinite/nonpos: {int(dat_stats['dat_e_nonfinite'])} / {int(dat_stats['dat_e_nonpos'])}",
        f"- DAT verdict: {'PASS' if dat_ok else 'FAIL'}",
        "",
        "## Primitive VTU checks",
        f"- Vr min/max: {vtu_stats['vtu_vr_min']:.6e} / {vtu_stats['vtu_vr_max']:.6e}",
        f"- rho nonfinite/nonpos: {int(vtu_stats['vtu_rho_nonfinite'])} / {int(vtu_stats['vtu_rho_nonpos'])}",
        f"- p nonfinite/nonpos: {int(vtu_stats['vtu_p_nonfinite'])} / {int(vtu_stats['vtu_p_nonpos'])}",
        f"- VTU verdict: {'PASS' if vtu_ok else 'FAIL'}",
        "",
        "## Radial profile checks",
        f"- rho median increasing-step fraction: {rho_inc_frac:.3f}",
        f"- Vr median outward-bin fraction: {vr_out_frac:.3f}",
        f"- Profile verdict: {'PASS' if profile_ok else 'WARN'}",
        "",
        "## Outputs",
        f"- [`init_checks.csv`]({csv_path})",
        f"- [`init_radial_profiles.csv`]({profile_csv_path})",
        f"- [`init_radial_profiles.png`]({png_path})",
    ]

    md_path.write_text("\n".join(summary_lines) + "\n")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
