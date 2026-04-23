# Copied to project analysis task toolbox
# Source case path: off_2270/analysis/amr_level_coverage_stats.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270/analysis/amr_level_coverage_stats.py

#!/usr/bin/env python3
"""Summarize AMR level coverage for key snapshots in off_2270."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np

CASE_DIR = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = Path(__file__).resolve().parent
ANGLE_SCALE = 2.0 * math.pi

if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from analyze_snapshot_dt_bottleneck import (  # noqa: E402
    block_uniform_edges,
    get_header,
    get_tree_info,
    parse_log,
    read_par_config,
)


DEFAULT_FRAMES = "off0000,off0016,off0032,off0064,off0096,off0112"
DEFAULT_THRESHOLDS_DEG = "0.1,0.5,1,2,5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Report AMR level coverage over the whole domain and polar theta caps."
    )
    parser.add_argument("--case-dir", type=Path, default=CASE_DIR)
    parser.add_argument("--par", default="amrvac.par")
    parser.add_argument("--log", default="data/off.log")
    parser.add_argument("--frames", default=DEFAULT_FRAMES)
    parser.add_argument("--theta-thresholds-deg", default=DEFAULT_THRESHOLDS_DEG)
    parser.add_argument("--output-dir", default="analysis")
    return parser.parse_args()


def parse_float_list(raw: str) -> List[float]:
    return [float(token.strip()) for token in raw.split(",") if token.strip()]


def parse_frame_list(raw: str) -> List[str]:
    return [token.strip() for token in raw.split(",") if token.strip()]


def inner_surface_cell_areas(th_edges: np.ndarray, ph_edges: np.ndarray) -> np.ndarray:
    theta_phys = ANGLE_SCALE * th_edges
    phi_phys = ANGLE_SCALE * ph_edges
    dphi = np.diff(phi_phys)[None, :]
    theta_lo = theta_phys[:-1][:, None]
    theta_hi = theta_phys[1:][:, None]
    return (np.cos(theta_lo) - np.cos(theta_hi)) * dphi


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    output_dir = (case_dir / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    par_cfg = read_par_config(case_dir / args.par)
    log_by_it = parse_log(case_dir / args.log)
    thresholds_deg = parse_float_list(args.theta_thresholds_deg)
    thresholds_rad = [math.radians(val) for val in thresholds_deg]
    frames = parse_frame_list(args.frames)

    overall_rows: List[Dict[str, object]] = []
    theta_rows: List[Dict[str, object]] = []
    md_lines: List[str] = []
    md_lines.append("# AMR Level Coverage Summary")
    md_lines.append("")
    md_lines.append(
        "Domain theta range (physical): "
        f"{ANGLE_SCALE * par_cfg['xprobmin2']:.6f} to {ANGLE_SCALE * par_cfg['xprobmax2']:.6f} rad"
    )
    md_lines.append("")

    for frame in frames:
        dat_path = case_dir / "data" / f"{frame}.dat"
        with dat_path.open("rb") as handle:
            header = get_header(handle)
            block_lvls, block_ixs, _block_offsets = get_tree_info(handle)

        block_shape = header["block_nx"].astype(int)
        domain_nx = header["domain_nx"].astype(int)
        ndim = int(header["ndim"])
        if ndim != 3:
            raise ValueError(f"{frame}: expected 3D data, got ndim={ndim}")

        xmins = [par_cfg[f"xprobmin{i}"] for i in range(1, 4)]
        xmaxs = [par_cfg[f"xprobmax{i}"] for i in range(1, 4)]

        total_cells = 0
        level_cells: Counter[int] = Counter()
        inner_area_total = 0.0
        inner_area_by_level: defaultdict[int, float] = defaultdict(float)

        cap_total_cells: Dict[float, int] = {deg: 0 for deg in thresholds_deg}
        cap_level_cells: Dict[float, Counter[int]] = {deg: Counter() for deg in thresholds_deg}
        cap_inner_area_total: Dict[float, float] = {deg: 0.0 for deg in thresholds_deg}
        cap_inner_area_by_level: Dict[float, defaultdict[int, float]] = {
            deg: defaultdict(float) for deg in thresholds_deg
        }

        for level, bix in zip(block_lvls.astype(int), block_ixs.astype(int)):
            block_ncells = int(np.prod(block_shape))
            total_cells += block_ncells
            level_cells[level] += block_ncells

            th_edges = block_uniform_edges(
                xmins[1], xmaxs[1], int(domain_nx[1]), int(level), int(bix[1]), int(block_shape[1])
            )
            ph_edges = block_uniform_edges(
                xmins[2], xmaxs[2], int(domain_nx[2]), int(level), int(bix[2]), int(block_shape[2])
            )
            th_centers = 0.5 * (th_edges[:-1] + th_edges[1:])
            th_centers_phys = ANGLE_SCALE * th_centers

            # volume-cell statistics in theta caps across all radii/phi
            for deg, thr in zip(thresholds_deg, thresholds_rad):
                n_theta = int(np.count_nonzero(th_centers_phys < thr))
                if n_theta <= 0:
                    continue
                count_here = n_theta * int(block_shape[0]) * int(block_shape[2])
                cap_total_cells[deg] += count_here
                cap_level_cells[deg][level] += count_here

            # inner-surface area statistics only for blocks touching the inner boundary
            if int(bix[0]) == 1:
                area2d = inner_surface_cell_areas(th_edges, ph_edges)
                area_sum = float(np.sum(area2d))
                inner_area_total += area_sum
                inner_area_by_level[level] += area_sum

                for deg, thr in zip(thresholds_deg, thresholds_rad):
                    theta_rows_mask = th_centers_phys < thr
                    area_cap = float(np.sum(area2d[theta_rows_mask, :]))
                    if area_cap <= 0.0:
                        continue
                    cap_inner_area_total[deg] += area_cap
                    cap_inner_area_by_level[deg][level] += area_cap

        log_info = log_by_it.get(int(header["it"]), {})
        base_level = int(np.min(block_lvls))
        max_level = int(np.max(block_lvls))
        frame_overall = {
            "frame": frame,
            "it": int(header["it"]),
            "time": float(header["time"]),
            "dt_log": float(log_info.get("dt", np.nan)),
            "nleafs": int(header["nleafs"]),
            "base_level": base_level,
            "max_level": max_level,
        }
        for level in sorted(level_cells):
            frame_overall[f"cells_L{level}"] = int(level_cells[level])
            frame_overall[f"cell_frac_L{level}"] = float(level_cells[level] / max(total_cells, 1))
            frame_overall[f"inner_area_L{level}"] = float(inner_area_by_level[level])
            frame_overall[f"inner_area_frac_L{level}"] = float(
                inner_area_by_level[level] / max(inner_area_total, 1.0e-30)
            )
        overall_rows.append(frame_overall)

        md_lines.append(f"## {frame}")
        md_lines.append("")
        md_lines.append(
            f"- base/max level: `L{base_level}` / `L{max_level}`; "
            f"`nleafs={int(header['nleafs'])}`; "
            f"`dt_log={float(log_info.get('dt', np.nan)):.3e}`"
        )
        level_bits = []
        for level in sorted(level_cells):
            level_bits.append(
                f"`L{level}` whole-domain cell frac = {level_cells[level] / max(total_cells, 1):.3%}, "
                f"inner-surface area frac = {inner_area_by_level[level] / max(inner_area_total, 1.0e-30):.3%}"
            )
        md_lines.append("- " + "; ".join(level_bits))
        md_lines.append("")

        for deg in thresholds_deg:
            row = {
                "frame": frame,
                "it": int(header["it"]),
                "time": float(header["time"]),
                "theta_deg": deg,
                "cap_total_cells": int(cap_total_cells[deg]),
                "cap_domain_cell_frac": float(cap_total_cells[deg] / max(total_cells, 1)),
                "cap_inner_area_total": float(cap_inner_area_total[deg]),
                "cap_inner_area_domain_frac": float(cap_inner_area_total[deg] / max(inner_area_total, 1.0e-30)),
            }
            for level in sorted(level_cells):
                row[f"cap_cells_L{level}"] = int(cap_level_cells[deg][level])
                row[f"cap_cell_frac_within_L{level}"] = float(
                    cap_level_cells[deg][level] / max(cap_total_cells[deg], 1)
                )
                row[f"cap_inner_area_L{level}"] = float(cap_inner_area_by_level[deg][level])
                row[f"cap_inner_area_frac_within_L{level}"] = float(
                    cap_inner_area_by_level[deg][level] / max(cap_inner_area_total[deg], 1.0e-30)
                )
            theta_rows.append(row)

            if deg in (0.1, 0.5, 1.0, 2.0, 5.0):
                md_lines.append(
                    f"- `theta < {deg:.1f} deg`: "
                    f"cell frac of domain = {cap_total_cells[deg] / max(total_cells, 1):.3%}, "
                    f"inner-area frac = {cap_inner_area_total[deg] / max(inner_area_total, 1.0e-30):.3%}"
                )
        md_lines.append("")

    write_csv(output_dir / "amr_level_overall.csv", overall_rows)
    write_csv(output_dir / "amr_theta_cap_summary.csv", theta_rows)
    (output_dir / "amr_level_coverage_summary.md").write_text("\n".join(md_lines) + "\n")

    print(f"Wrote {output_dir / 'amr_level_overall.csv'}")
    print(f"Wrote {output_dir / 'amr_theta_cap_summary.csv'}")
    print(f"Wrote {output_dir / 'amr_level_coverage_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
