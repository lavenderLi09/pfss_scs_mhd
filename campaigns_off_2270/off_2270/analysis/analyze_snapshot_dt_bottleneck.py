#!/usr/bin/env python3
"""Analyze per-snapshot CFL bottlenecks for stretched spherical AMRVAC data."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

CASE_DIR = Path(__file__).resolve().parent.parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))

from hao_code.datfile_io import SIZE_INT, get_header, get_tree_info  # noqa: E402


EPS = 1.0e-30
TOP_VABS_N = 30
TOP_VABS_INNER_N = 30
TOP_VABS_INNER_RMAX = 2.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each AMRVAC snapshot, locate the cell that would control the "
            "minimum dt and report whether geometry, flow, sound speed, or "
            "magnetosonic speed dominates the CFL constraint."
        )
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        default=CASE_DIR,
        help="Case directory containing amrvac.par and data/",
    )
    parser.add_argument(
        "--pattern",
        default="data/off*.dat",
        help="Glob pattern, relative to --case-dir, for snapshots to analyze.",
    )
    parser.add_argument(
        "--log",
        default="data/off.log",
        help="Log path relative to --case-dir.",
    )
    parser.add_argument(
        "--topn",
        type=int,
        default=10,
        help="Number of top bottleneck cells to save per frame.",
    )
    parser.add_argument(
        "--frames",
        default="",
        help="Optional comma-separated frame stems such as off0016,off0017.",
    )
    parser.add_argument(
        "--par",
        default="amrvac.par",
        help="Parameter file relative to --case-dir.",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis",
        help="Output directory relative to --case-dir.",
    )
    return parser.parse_args()


def read_par_config(par_path: Path) -> Dict[str, float]:
    text = par_path.read_text()
    keys = [
        "domain_nx1",
        "domain_nx2",
        "domain_nx3",
        "block_nx1",
        "block_nx2",
        "block_nx3",
        "xprobmin1",
        "xprobmin2",
        "xprobmin3",
        "xprobmax1",
        "xprobmax2",
        "xprobmax3",
        "qstretch_baselevel",
        "courantpar",
        "mhd_gamma",
    ]
    config: Dict[str, float] = {}
    for key in keys:
        match = re.search(
            rf"{re.escape(key)}\s*=\s*([-+0-9.dDeE]+)",
            text,
            flags=re.IGNORECASE,
        )
        if not match:
            raise ValueError(f"Could not find {key} in {par_path}")
        config[key] = float(match.group(1).replace("d", "e").replace("D", "e"))
    return config


def parse_log(log_path: Path) -> Dict[int, Dict[str, float]]:
    by_it: Dict[int, Dict[str, float]] = {}
    with log_path.open() as handle:
        for raw in handle:
            if "|" in raw:
                raw = raw.split("|", 1)[0]
            parts = raw.split()
            if len(parts) < 3:
                continue
            try:
                it = int(parts[0])
                time = float(parts[1])
                dt = float(parts[2])
            except ValueError:
                continue
            by_it[it] = {"time": time, "dt": dt}
    return by_it


def select_frames(paths: Sequence[Path], frames_arg: str) -> List[Path]:
    if not frames_arg.strip():
        return list(paths)
    wanted = {token.strip() for token in frames_arg.split(",") if token.strip()}
    return [path for path in paths if path.stem in wanted]


def build_base_radial_edges(xmin: float, xmax: float, ncells: int, q: float) -> np.ndarray:
    if abs(q - 1.0) < 1.0e-14:
        return np.linspace(xmin, xmax, ncells + 1, dtype=np.float64)
    dr0 = (xmax - xmin) * (1.0 - q) / (1.0 - q**ncells)
    widths = dr0 * q ** np.arange(ncells, dtype=np.float64)
    edges = np.empty(ncells + 1, dtype=np.float64)
    edges[0] = xmin
    edges[1:] = xmin + np.cumsum(widths)
    edges[-1] = xmax
    return edges


def block_uniform_edges(
    xmin: float,
    xmax: float,
    domain_nx: int,
    level: int,
    block_ix: int,
    block_nx: int,
) -> np.ndarray:
    refine_ratio = 2 ** (level - 1)
    dx = (xmax - xmin) / (domain_nx * refine_ratio)
    start = (block_ix - 1) * block_nx
    return xmin + np.arange(start, start + block_nx + 1, dtype=np.float64) * dx


def block_radial_edges(
    base_edges: np.ndarray,
    level: int,
    block_ix: int,
    block_nx: int,
) -> np.ndarray:
    refine_ratio = 2 ** (level - 1)
    start = (block_ix - 1) * block_nx
    fine_edges = np.empty(block_nx + 1, dtype=np.float64)
    nbase = len(base_edges) - 1
    max_fine = nbase * refine_ratio
    for local_edge, fine_index in enumerate(range(start, start + block_nx + 1)):
        if fine_index >= max_fine:
            fine_edges[local_edge] = base_edges[-1]
            continue
        parent = fine_index // refine_ratio
        sub = fine_index % refine_ratio
        width = base_edges[parent + 1] - base_edges[parent]
        fine_edges[local_edge] = base_edges[parent] + width * (sub / refine_ratio)
    return fine_edges


def read_block_fields(handle, offset: int, block_shape: np.ndarray, ndim: int, nw: int) -> np.ndarray:
    count = int(np.prod(block_shape))
    handle.seek(offset + 2 * ndim * SIZE_INT)
    raw = np.fromfile(handle, dtype="=f8", count=nw * count)
    if raw.size != nw * count:
        raise IOError(f"Failed to read block fields at offset {offset}")
    fields = np.empty((nw, *block_shape), dtype=np.float64)
    for iw in range(nw):
        field = raw[iw * count : (iw + 1) * count]
        field = field.reshape(tuple(block_shape[::-1]), order="C").T
        while field.ndim < 3:
            field = field[..., np.newaxis]
        fields[iw] = field
    return fields


def conservative_internal_to_primitive(
    rho: np.ndarray,
    m1: np.ndarray,
    m2: np.ndarray,
    m3: np.ndarray,
    eint: np.ndarray,
    gamma: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rho_safe = np.maximum(rho, EPS)
    v_r = m1 / rho_safe
    v_th = m2 / rho_safe
    v_ph = m3 / rho_safe
    pth = np.maximum((gamma - 1.0) * eint, 0.0)
    return v_r, v_th, v_ph, pth


def classify_bottleneck(
    c_r: float,
    c_th: float,
    c_ph: float,
    l_r: float,
    l_th: float,
    l_ph: float,
    v_n: float,
    c_fast_n: float,
    c_s: float,
    v_a: float,
) -> Tuple[str, str]:
    components = np.array([c_r, c_th, c_ph], dtype=np.float64)
    directions = ["radial", "theta", "phi"]
    dominant_idx = int(np.nanargmax(components))
    dominant_direction = directions[dominant_idx]

    lengths = np.array([l_r, l_th, l_ph], dtype=np.float64)
    inv_lengths = 1.0 / np.maximum(lengths, EPS)
    geom_ratio = inv_lengths[dominant_idx] / np.maximum(np.median(inv_lengths), EPS)
    speed_ratio = (v_n + c_fast_n) / np.maximum(np.median([v_n + c_fast_n, c_s + v_a, EPS]), EPS)

    if geom_ratio > max(1.5, 1.2 * speed_ratio):
        return dominant_direction, "geometry"
    if v_n >= c_fast_n:
        return dominant_direction, "flow"
    if c_s >= 0.9 * v_a:
        return dominant_direction, "sound"
    return dominant_direction, "alfven/fast"


def summarize_cell(
    frame: str,
    it: int,
    time: float,
    dt_log: float,
    level: int,
    block_id: int,
    i_cell: int,
    j_cell: int,
    k_cell: int,
    r: float,
    theta: float,
    phi: float,
    rho: float,
    pth: float,
    v_r: float,
    v_th: float,
    v_ph: float,
    b_r: float,
    b_th: float,
    b_ph: float,
    c_s: float,
    v_a: float,
    cfr: float,
    cfth: float,
    cfph: float,
    dr: float,
    ds_th: float,
    ds_ph: float,
    shell_index: int,
    theta_bin: int,
    phi_bin: int,
    c_r: float,
    c_th: float,
    c_ph: float,
    courantpar: float,
) -> Dict[str, float | str | int]:
    c_sum = c_r + c_th + c_ph
    dt_pred = courantpar / c_sum
    vabs = math.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
    babs = math.sqrt(b_r * b_r + b_th * b_th + b_ph * b_ph)
    dom_idx = int(np.argmax([c_r, c_th, c_ph]))
    dominant_direction, dominant_physics = classify_bottleneck(
        c_r,
        c_th,
        c_ph,
        dr,
        ds_th,
        ds_ph,
        abs([v_r, v_th, v_ph][dom_idx]),
        [cfr, cfth, cfph][dom_idx],
        c_s,
        v_a,
    )

    return {
        "frame": frame,
        "it": it,
        "time": time,
        "dt_log": dt_log,
        "dt_pred": dt_pred,
        "level": level,
        "block_id": block_id,
        "i_cell": i_cell,
        "j_cell": j_cell,
        "k_cell": k_cell,
        "r": r,
        "theta": theta,
        "phi": phi,
        "rho": rho,
        "p": pth,
        "v_abs": vabs,
        "b_abs": babs,
        "v_r": v_r,
        "v_theta": v_th,
        "v_phi": v_ph,
        "b_r": b_r,
        "b_theta": b_th,
        "b_phi": b_ph,
        "c_s": c_s,
        "v_A": v_a,
        "cf_r": cfr,
        "cf_theta": cfth,
        "cf_phi": cfph,
        "dr": dr,
        "r_dtheta": ds_th,
        "r_sintheta_dphi": ds_ph,
        "shell_index": shell_index,
        "theta_bin": theta_bin,
        "phi_bin": phi_bin,
        "C_r": c_r,
        "C_th": c_th,
        "C_ph": c_ph,
        "C_sum": c_sum,
        "dominant_direction": dominant_direction,
        "dominant_physics": dominant_physics,
    }


def summarize_vabs_cell(
    frame: str,
    it: int,
    time: float,
    level: int,
    block_id: int,
    r: float,
    theta: float,
    phi: float,
    shell_index: int,
    theta_bin: int,
    phi_bin: int,
    rho: float,
    pth: float,
    v_r: float,
    v_th: float,
    v_ph: float,
    b_r: float,
    b_th: float,
    b_ph: float,
) -> Dict[str, float | str | int]:
    vabs = math.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
    babs = math.sqrt(b_r * b_r + b_th * b_th + b_ph * b_ph)
    return {
        "frame": frame,
        "it": it,
        "time": time,
        "level": level,
        "block_id": block_id,
        "r": r,
        "theta": theta,
        "phi": phi,
        "shell_index": shell_index,
        "theta_bin": theta_bin,
        "phi_bin": phi_bin,
        "rho": rho,
        "p": pth,
        "v_abs": vabs,
        "b_abs": babs,
        "v_r": v_r,
        "v_theta": v_th,
        "v_phi": v_ph,
        "b_r": b_r,
        "b_theta": b_th,
        "b_phi": b_ph,
    }


def cyclic_bin_distance(a: int, b: int, period: int) -> int:
    direct = abs(a - b)
    return min(direct, period - direct)


def rows_are_neighbors(a: Dict[str, object], b: Dict[str, object], phi_bins: int) -> bool:
    return (
        abs(int(a["shell_index"]) - int(b["shell_index"])) <= 1
        and abs(int(a["theta_bin"]) - int(b["theta_bin"])) <= 1
        and cyclic_bin_distance(int(a["phi_bin"]), int(b["phi_bin"]), phi_bins) <= 1
    )


def cluster_frame_rows(frame_rows: Sequence[Dict[str, object]], phi_bins: int) -> List[Dict[str, object]]:
    if not frame_rows:
        return []
    n = len(frame_rows)
    seen = [False] * n
    clusters: List[List[int]] = []
    for i in range(n):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        members: List[int] = []
        while stack:
            cur = stack.pop()
            members.append(cur)
            for j in range(n):
                if seen[j]:
                    continue
                if rows_are_neighbors(frame_rows[cur], frame_rows[j], phi_bins):
                    seen[j] = True
                    stack.append(j)
        clusters.append(sorted(members))

    hotspot_rows: List[Dict[str, object]] = []
    for cluster_id, members in enumerate(clusters, start=1):
        rows = [frame_rows[idx] for idx in members]
        rep = min(rows, key=lambda row: float(row["dt_pred"]))
        block_counter = Counter(int(row["block_id"]) for row in rows)
        level_counter = Counter(int(row["level"]) for row in rows)
        dir_counter = Counter(str(row["dominant_direction"]) for row in rows)
        phys_counter = Counter(str(row["dominant_physics"]) for row in rows)
        hotspot_rows.append(
            {
                "frame": rep["frame"],
                "cluster_id": cluster_id,
                "n_cells": len(rows),
                "rank_min": min(int(row["rank"]) for row in rows),
                "dt_pred_min": float(rep["dt_pred"]),
                "r_rep": float(rep["r"]),
                "theta_rep": float(rep["theta"]),
                "phi_rep": float(rep["phi"]),
                "shell_index_min": min(int(row["shell_index"]) for row in rows),
                "shell_index_max": max(int(row["shell_index"]) for row in rows),
                "theta_min": min(float(row["theta"]) for row in rows),
                "theta_max": max(float(row["theta"]) for row in rows),
                "phi_min": min(float(row["phi"]) for row in rows),
                "phi_max": max(float(row["phi"]) for row in rows),
                "block_id_dom": block_counter.most_common(1)[0][0],
                "level_dom": level_counter.most_common(1)[0][0],
                "dominant_direction": dir_counter.most_common(1)[0][0],
                "dominant_physics": phys_counter.most_common(1)[0][0],
            }
        )
    hotspot_rows.sort(key=lambda row: (int(row["rank_min"]), float(row["dt_pred_min"])))
    return hotspot_rows


def build_hotspot_rows(top_rows: Sequence[Dict[str, object]], phi_bins: int) -> List[Dict[str, object]]:
    by_frame: Dict[str, List[Dict[str, object]]] = {}
    for idx, row in enumerate(sorted(top_rows, key=lambda r: (str(r["frame"]), float(r["dt_pred"]))), start=1):
        copied = dict(row)
        copied["rank"] = idx
        by_frame.setdefault(str(copied["frame"]), []).append(copied)

    hotspot_rows: List[Dict[str, object]] = []
    for frame, rows in by_frame.items():
        rows.sort(key=lambda row: float(row["dt_pred"]))
        for rank, row in enumerate(rows, start=1):
            row["rank"] = rank
        hotspot_rows.extend(cluster_frame_rows(rows, phi_bins))
    return hotspot_rows


def write_hotspot_markdown(path: Path, hotspot_rows: Sequence[Dict[str, object]], topn_by_frame: Dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_frame: Dict[str, List[Dict[str, object]]] = {}
    for row in hotspot_rows:
        by_frame.setdefault(str(row["frame"]), []).append(row)

    lines: List[str] = ["## Bottleneck Hotspots", ""]
    for frame in sorted(by_frame):
        rows = sorted(by_frame[frame], key=lambda row: int(row["rank_min"]))
        n_cells = sum(int(row["n_cells"]) for row in rows)
        topn = topn_by_frame.get(frame, n_cells)
        if len(rows) == 1:
            summary = f"All top-{topn} cells collapse into one hotspot near shell {rows[0]['shell_index_min']}."
        else:
            summary = f"Top-{topn} cells split into {len(rows)} distinct hotspots."
        lines.extend([f"### Top bottleneck hotspots for `{frame}`", "", summary, ""])
        for row in rows:
            lines.append(
                f"- `H{int(row['cluster_id'])}`: `n={int(row['n_cells'])}`, "
                f"`rank_min={int(row['rank_min'])}`, `dt_pred_min={float(row['dt_pred_min']):.3e}`, "
                f"`shell={int(row['shell_index_min'])}-{int(row['shell_index_max'])}`, "
                f"`r≈{float(row['r_rep']):.4f}`, `theta≈{float(row['theta_rep']):.4e}`, `phi≈{float(row['phi_rep']):.4f}`, "
                f"`level={int(row['level_dom'])}`, `block={int(row['block_id_dom'])}`, "
                f"`dir={row['dominant_direction']}`, `physics={row['dominant_physics']}`"
            )
        lines.append("")
    path.write_text("\n".join(lines))


def summarize_distribution(name: str, values: np.ndarray) -> Dict[str, float]:
    finite = values[np.isfinite(values) & (values > 0.0)]
    if finite.size == 0:
        return {
            f"{name}_count": 0.0,
            f"{name}_min": math.nan,
            f"{name}_p50": math.nan,
            f"{name}_p90": math.nan,
            f"{name}_p99": math.nan,
            f"{name}_max": math.nan,
            f"{name}_mean": math.nan,
        }
    return {
        f"{name}_count": float(finite.size),
        f"{name}_min": float(np.min(finite)),
        f"{name}_p50": float(np.percentile(finite, 50.0)),
        f"{name}_p90": float(np.percentile(finite, 90.0)),
        f"{name}_p99": float(np.percentile(finite, 99.0)),
        f"{name}_max": float(np.max(finite)),
        f"{name}_mean": float(np.mean(finite)),
    }


def init_shell_stats(nshell: int) -> Dict[str, np.ndarray]:
    return {
        "sum": np.zeros(nshell, dtype=np.float64),
        "count": np.zeros(nshell, dtype=np.float64),
        "max": np.full(nshell, -np.inf, dtype=np.float64),
    }


def update_shell_stats(stats: Dict[str, np.ndarray], shell_idx: np.ndarray, values: np.ndarray) -> None:
    idx = shell_idx.ravel().astype(np.int64)
    vals = values.ravel()
    finite = np.isfinite(vals) & (vals > 0.0)
    if not np.any(finite):
        return
    idx = idx[finite]
    vals = vals[finite]
    stats["sum"] += np.bincount(idx, weights=vals, minlength=stats["sum"].size)
    stats["count"] += np.bincount(idx, minlength=stats["count"].size)
    np.maximum.at(stats["max"], idx, vals)


def init_surface_stats(shape: Tuple[int, int]) -> Dict[str, np.ndarray]:
    return {
        "sum": np.zeros(shape, dtype=np.float64),
        "count": np.zeros(shape, dtype=np.float64),
        "max": np.full(shape, -np.inf, dtype=np.float64),
    }


def update_surface_stats(
    stats: Dict[str, np.ndarray], row_idx: np.ndarray, col_idx: np.ndarray, values: np.ndarray
) -> None:
    rows = row_idx.ravel().astype(np.int64)
    cols = col_idx.ravel().astype(np.int64)
    vals = values.ravel()
    finite = np.isfinite(vals) & (vals > 0.0)
    if not np.any(finite):
        return
    rows = rows[finite]
    cols = cols[finite]
    vals = vals[finite]
    ncols = stats["sum"].shape[1]
    flat = rows * ncols + cols
    stats["sum"].ravel()[:] += np.bincount(flat, weights=vals, minlength=stats["sum"].size)
    stats["count"].ravel()[:] += np.bincount(flat, minlength=stats["count"].size)
    np.maximum.at(stats["max"].ravel(), flat, vals)


def shell_mean(stats: Dict[str, np.ndarray]) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return stats["sum"] / np.maximum(stats["count"], 1.0)


def plot_surface_maps(
    output_dir: Path,
    frame: str,
    selector: str,
    maps: Dict[str, np.ndarray],
    shell_rows: Sequence[Dict[str, object]],
    theta_edges_deg: np.ndarray,
    phi_edges_deg: np.ndarray,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path = output_dir / f"{frame}_surface_maps_{selector}.png"
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    quantities = [
        ("vabs", "|v|"),
        ("cs", "c_s"),
        ("va", "v_A"),
        ("cfmax", "c_f,max"),
    ]
    shell_lookup = {str(row["quantity"]): row for row in shell_rows}
    X, Y = np.meshgrid(phi_edges_deg, theta_edges_deg)
    for ax, (key, label) in zip(axes.ravel(), quantities):
        grid = maps[key]
        finite = grid[np.isfinite(grid) & (grid > 0.0)]
        if finite.size == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
            ax.set_title(label)
            continue
        dynamic = float(np.max(finite) / max(np.min(finite), EPS))
        if dynamic > 100.0:
            plot_values = np.log10(np.where(grid > 0.0, grid, np.nan))
            cbar_label = f"log10({label})"
        else:
            plot_values = grid
            cbar_label = label
        mesh = ax.pcolormesh(X, Y, plot_values, shading="auto")
        shell_row = shell_lookup[key]
        ax.set_title(
            f"{label}: r={float(shell_row['shell_radius']):.3f}\n"
            f"{selector}={float(shell_row['selector_value']):.3e}"
        )
        ax.set_xlabel("phi [deg]")
        ax.set_ylabel("theta [deg]")
        fig.colorbar(mesh, ax=ax, pad=0.01, label=cbar_label)
    fig.suptitle(f"{frame} theta-phi surface maps ({selector}-selected shells)", fontsize=14)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def analyze_snapshot(
    dat_path: Path,
    par_cfg: Dict[str, float],
    log_by_it: Dict[int, Dict[str, float]],
    topn: int,
    base_edges: np.ndarray,
) -> Tuple[
    Dict[str, float | str | int],
    List[Dict[str, float | str | int]],
    List[Dict[str, float | str | int]],
    List[Dict[str, float | str | int]],
    Dict[str, float | str | int],
    Dict[str, np.ndarray],
    List[Dict[str, object]],
    Dict[str, Dict[str, np.ndarray]],
    Dict[str, int],
]:
    with dat_path.open("rb") as handle:
        header = get_header(handle)
        block_lvls, block_ixs, block_offsets = get_tree_info(handle)

        name_to_idx = {name: idx for idx, name in enumerate(header["w_names"])}
        required = ["rho", "m1", "m2", "m3", "e", "b1", "b2", "b3"]
        missing = [name for name in required if name not in name_to_idx]
        if missing:
            raise KeyError(f"{dat_path.name} missing required fields: {missing}")

        gamma = float(
            header["params"][header["param_names"].index("gamma")]
            if "gamma" in header["param_names"]
            else par_cfg["mhd_gamma"]
        )
        dt_log = float(log_by_it.get(int(header["it"]), {}).get("dt", math.nan))
        time_log = float(log_by_it.get(int(header["it"]), {}).get("time", header["time"]))

        top_rows: List[Dict[str, float | str | int]] = []
        top_vabs_rows: List[Dict[str, float | str | int]] = []
        top_vabs_inner_rows: List[Dict[str, float | str | int]] = []
        cs_values: List[np.ndarray] = []
        va_values: List[np.ndarray] = []
        cfmax_values: List[np.ndarray] = []
        vabs_values: List[np.ndarray] = []
        shell_stats = {
            "vabs": init_shell_stats(int(par_cfg["domain_nx1"])),
            "cs": init_shell_stats(int(par_cfg["domain_nx1"])),
            "va": init_shell_stats(int(par_cfg["domain_nx1"])),
            "cfmax": init_shell_stats(int(par_cfg["domain_nx1"])),
        }
        block_shape = header["block_nx"].astype(int)
        domain_nx = header["domain_nx"].astype(int)
        xmins = [par_cfg[f"xprobmin{i}"] for i in range(1, 4)]
        xmaxs = [par_cfg[f"xprobmax{i}"] for i in range(1, 4)]
        courantpar = par_cfg["courantpar"]

        for block_id, (level, bix, offset) in enumerate(
            zip(block_lvls, block_ixs, block_offsets), start=1
        ):
            fields = read_block_fields(
                handle,
                int(offset),
                block_shape,
                int(header["ndim"]),
                int(header["nw"]),
            )

            r_edges = block_radial_edges(base_edges, int(level), int(bix[0]), int(block_shape[0]))
            th_edges = block_uniform_edges(
                xmins[1], xmaxs[1], int(domain_nx[1]), int(level), int(bix[1]), int(block_shape[1])
            )
            ph_edges = block_uniform_edges(
                xmins[2], xmaxs[2], int(domain_nx[2]), int(level), int(bix[2]), int(block_shape[2])
            )

            r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
            th_centers = 0.5 * (th_edges[:-1] + th_edges[1:])
            ph_centers = 0.5 * (ph_edges[:-1] + ph_edges[1:])
            dr = np.diff(r_edges)
            dth = np.diff(th_edges)
            dph = np.diff(ph_edges)
            refine_ratio = 2 ** (int(level) - 1)
            start_r = (int(bix[0]) - 1) * int(block_shape[0])
            start_t = (int(bix[1]) - 1) * int(block_shape[1])
            start_p = (int(bix[2]) - 1) * int(block_shape[2])
            base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)
            base_t_idx = (np.arange(start_t, start_t + int(block_shape[1])) // refine_ratio).astype(np.int64)
            base_p_idx = (np.arange(start_p, start_p + int(block_shape[2])) // refine_ratio).astype(np.int64)

            rr, tt, pp = np.meshgrid(r_centers, th_centers, ph_centers, indexing="ij")
            dr3 = dr[:, None, None]
            ds_th = rr * dth[None, :, None]
            ds_ph = rr * np.maximum(np.sin(tt), EPS) * dph[None, None, :]

            rho = fields[name_to_idx["rho"]]
            m1 = fields[name_to_idx["m1"]]
            m2 = fields[name_to_idx["m2"]]
            m3 = fields[name_to_idx["m3"]]
            eint = fields[name_to_idx["e"]]
            b_r = fields[name_to_idx["b1"]]
            b_th = fields[name_to_idx["b2"]]
            b_ph = fields[name_to_idx["b3"]]
            v_r, v_th, v_ph, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)

            valid = (rho > 0.0) & (pth > 0.0)
            if not np.any(valid):
                continue

            b_sq = b_r * b_r + b_th * b_th + b_ph * b_ph
            c_s_sq = np.where(valid, gamma * pth / np.maximum(rho, EPS), np.nan)
            v_a_sq = np.where(valid, b_sq / np.maximum(rho, EPS), np.nan)
            c_s = np.sqrt(np.maximum(c_s_sq, 0.0))
            v_a = np.sqrt(np.maximum(v_a_sq, 0.0))

            def fast_speed(normal_b: np.ndarray) -> np.ndarray:
                b_n_sq = normal_b * normal_b
                disc = (c_s_sq + v_a_sq) ** 2 - 4.0 * c_s_sq * b_n_sq / np.maximum(rho, EPS)
                disc = np.maximum(disc, 0.0)
                cf_sq = 0.5 * ((c_s_sq + v_a_sq) + np.sqrt(disc))
                return np.sqrt(np.maximum(cf_sq, 0.0))

            cf_r = fast_speed(b_r)
            cf_th = fast_speed(b_th)
            cf_ph = fast_speed(b_ph)
            cf_max = np.maximum(np.maximum(cf_r, cf_th), cf_ph)
            vabs = np.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)

            c_r = np.where(valid, (np.abs(v_r) + cf_r) / np.maximum(dr3, EPS), np.nan)
            c_th = np.where(valid, (np.abs(v_th) + cf_th) / np.maximum(ds_th, EPS), np.nan)
            c_ph = np.where(valid, (np.abs(v_ph) + cf_ph) / np.maximum(ds_ph, EPS), np.nan)
            c_sum = c_r + c_th + c_ph

            cs_values.append(c_s[valid].ravel())
            va_values.append(v_a[valid].ravel())
            cfmax_values.append(cf_max[valid].ravel())
            vabs_values.append(vabs[valid].ravel())
            shell_idx_3d = np.broadcast_to(base_r_idx[:, None, None], rho.shape)
            update_shell_stats(shell_stats["vabs"], shell_idx_3d, vabs)
            update_shell_stats(shell_stats["cs"], shell_idx_3d, c_s)
            update_shell_stats(shell_stats["va"], shell_idx_3d, v_a)
            update_shell_stats(shell_stats["cfmax"], shell_idx_3d, cf_max)

            finite_mask = np.isfinite(c_sum) & (c_sum > 0.0)
            if not np.any(finite_mask):
                continue

            flat_indices = np.flatnonzero(finite_mask)
            values = c_sum.ravel()[flat_indices]
            n_keep = min(topn, values.size)
            chosen = flat_indices[np.argpartition(values, -n_keep)[-n_keep:]]
            chosen = chosen[np.argsort(c_sum.ravel()[chosen])[::-1]]

            for flat_idx in chosen:
                i_cell, j_cell, k_cell = np.unravel_index(int(flat_idx), c_sum.shape)
                row = summarize_cell(
                    frame=dat_path.stem,
                    it=int(header["it"]),
                    time=float(header["time"]),
                    dt_log=dt_log,
                    level=int(level),
                    block_id=block_id,
                    i_cell=int(i_cell),
                    j_cell=int(j_cell),
                    k_cell=int(k_cell),
                    r=float(rr[i_cell, j_cell, k_cell]),
                    theta=float(tt[i_cell, j_cell, k_cell]),
                    phi=float(pp[i_cell, j_cell, k_cell]),
                    rho=float(rho[i_cell, j_cell, k_cell]),
                    pth=float(pth[i_cell, j_cell, k_cell]),
                    v_r=float(v_r[i_cell, j_cell, k_cell]),
                    v_th=float(v_th[i_cell, j_cell, k_cell]),
                    v_ph=float(v_ph[i_cell, j_cell, k_cell]),
                    b_r=float(b_r[i_cell, j_cell, k_cell]),
                    b_th=float(b_th[i_cell, j_cell, k_cell]),
                    b_ph=float(b_ph[i_cell, j_cell, k_cell]),
                    c_s=float(c_s[i_cell, j_cell, k_cell]),
                    v_a=float(v_a[i_cell, j_cell, k_cell]),
                    cfr=float(cf_r[i_cell, j_cell, k_cell]),
                    cfth=float(cf_th[i_cell, j_cell, k_cell]),
                    cfph=float(cf_ph[i_cell, j_cell, k_cell]),
                    dr=float(dr3[i_cell, 0, 0]),
                    ds_th=float(ds_th[i_cell, j_cell, k_cell]),
                    ds_ph=float(ds_ph[i_cell, j_cell, k_cell]),
                    shell_index=int(base_r_idx[i_cell]),
                    theta_bin=int(base_t_idx[j_cell]),
                    phi_bin=int(base_p_idx[k_cell]),
                    c_r=float(c_r[i_cell, j_cell, k_cell]),
                    c_th=float(c_th[i_cell, j_cell, k_cell]),
                    c_ph=float(c_ph[i_cell, j_cell, k_cell]),
                    courantpar=courantpar,
                )
                top_rows.append(row)

            vabs_mask = np.isfinite(vabs) & (vabs > 0.0)
            if np.any(vabs_mask):
                flat_indices_v = np.flatnonzero(vabs_mask)
                values_v = vabs.ravel()[flat_indices_v]
                n_keep_v = min(TOP_VABS_N, values_v.size)
                chosen_v = flat_indices_v[np.argpartition(values_v, -n_keep_v)[-n_keep_v:]]
                chosen_v = chosen_v[np.argsort(vabs.ravel()[chosen_v])[::-1]]
                for flat_idx in chosen_v:
                    i_cell, j_cell, k_cell = np.unravel_index(int(flat_idx), vabs.shape)
                    top_vabs_rows.append(
                        summarize_vabs_cell(
                            frame=dat_path.stem,
                            it=int(header["it"]),
                            time=float(header["time"]),
                            level=int(level),
                            block_id=block_id,
                            r=float(rr[i_cell, j_cell, k_cell]),
                            theta=float(tt[i_cell, j_cell, k_cell]),
                            phi=float(pp[i_cell, j_cell, k_cell]),
                            shell_index=int(base_r_idx[i_cell]),
                            theta_bin=int(base_t_idx[j_cell]),
                            phi_bin=int(base_p_idx[k_cell]),
                            rho=float(rho[i_cell, j_cell, k_cell]),
                            pth=float(pth[i_cell, j_cell, k_cell]),
                            v_r=float(v_r[i_cell, j_cell, k_cell]),
                            v_th=float(v_th[i_cell, j_cell, k_cell]),
                            v_ph=float(v_ph[i_cell, j_cell, k_cell]),
                            b_r=float(b_r[i_cell, j_cell, k_cell]),
                            b_th=float(b_th[i_cell, j_cell, k_cell]),
                            b_ph=float(b_ph[i_cell, j_cell, k_cell]),
                        )
                    )
            vabs_inner_mask = vabs_mask & (rr <= TOP_VABS_INNER_RMAX)
            if np.any(vabs_inner_mask):
                flat_indices_v_inner = np.flatnonzero(vabs_inner_mask)
                values_v_inner = vabs.ravel()[flat_indices_v_inner]
                n_keep_v_inner = min(TOP_VABS_INNER_N, values_v_inner.size)
                chosen_v_inner = flat_indices_v_inner[
                    np.argpartition(values_v_inner, -n_keep_v_inner)[-n_keep_v_inner:]
                ]
                chosen_v_inner = chosen_v_inner[np.argsort(vabs.ravel()[chosen_v_inner])[::-1]]
                for flat_idx in chosen_v_inner:
                    i_cell, j_cell, k_cell = np.unravel_index(int(flat_idx), vabs.shape)
                    top_vabs_inner_rows.append(
                        summarize_vabs_cell(
                            frame=dat_path.stem,
                            it=int(header["it"]),
                            time=float(header["time"]),
                            level=int(level),
                            block_id=block_id,
                            r=float(rr[i_cell, j_cell, k_cell]),
                            theta=float(tt[i_cell, j_cell, k_cell]),
                            phi=float(pp[i_cell, j_cell, k_cell]),
                            shell_index=int(base_r_idx[i_cell]),
                            theta_bin=int(base_t_idx[j_cell]),
                            phi_bin=int(base_p_idx[k_cell]),
                            rho=float(rho[i_cell, j_cell, k_cell]),
                            pth=float(pth[i_cell, j_cell, k_cell]),
                            v_r=float(v_r[i_cell, j_cell, k_cell]),
                            v_th=float(v_th[i_cell, j_cell, k_cell]),
                            v_ph=float(v_ph[i_cell, j_cell, k_cell]),
                            b_r=float(b_r[i_cell, j_cell, k_cell]),
                            b_th=float(b_th[i_cell, j_cell, k_cell]),
                            b_ph=float(b_ph[i_cell, j_cell, k_cell]),
                        )
                    )

    top_rows.sort(key=lambda row: row["dt_pred"])
    top_rows = top_rows[:topn]
    if not top_rows:
        raise RuntimeError(f"No valid CFL cells found in {dat_path}")
    top_vabs_rows.sort(key=lambda row: float(row["v_abs"]), reverse=True)
    top_vabs_rows = top_vabs_rows[:TOP_VABS_N]
    top_vabs_inner_rows.sort(key=lambda row: float(row["v_abs"]), reverse=True)
    top_vabs_inner_rows = top_vabs_inner_rows[:TOP_VABS_INNER_N]

    cs_all = np.concatenate(cs_values) if cs_values else np.array([], dtype=np.float64)
    va_all = np.concatenate(va_values) if va_values else np.array([], dtype=np.float64)
    cfmax_all = np.concatenate(cfmax_values) if cfmax_values else np.array([], dtype=np.float64)
    vabs_all = np.concatenate(vabs_values) if vabs_values else np.array([], dtype=np.float64)

    summary = {
        "frame": dat_path.stem,
        "it": top_rows[0]["it"],
        "time": top_rows[0]["time"],
        "dt_log": top_rows[0]["dt_log"],
        "dt_pred_min": top_rows[0]["dt_pred"],
        "dominant_direction": top_rows[0]["dominant_direction"],
        "dominant_physics": top_rows[0]["dominant_physics"],
        "r_min": top_rows[0]["r"],
        "theta_min": top_rows[0]["theta"],
        "phi_min": top_rows[0]["phi"],
        "level": top_rows[0]["level"],
        "block_id": top_rows[0]["block_id"],
        "nleafs": header["nleafs"],
        "time_log": time_log,
        "gamma": gamma,
    }
    physical_summary: Dict[str, float | str | int] = {
        "frame": dat_path.stem,
        "it": int(header["it"]),
        "time": float(header["time"]),
    }
    for key, values in (
        ("vabs", vabs_all),
        ("cs", cs_all),
        ("va", va_all),
        ("cfmax", cfmax_all),
    ):
        physical_summary.update(summarize_distribution(key, values))
    physical_samples = {
        "vabs": vabs_all,
        "cs": cs_all,
        "va": va_all,
        "cfmax": cfmax_all,
    }
    shell_rows: List[Dict[str, object]] = []
    shell_selectors: Dict[str, Dict[str, int]] = {}
    for quantity, stats in shell_stats.items():
        mean_arr = shell_mean(stats)
        max_arr = stats["max"]
        mean_idx = int(np.nanargmax(np.where(stats["count"] > 0.0, mean_arr, np.nan)))
        max_idx = int(np.nanargmax(np.where(np.isfinite(max_arr), max_arr, np.nan)))
        shell_selectors[quantity] = {"mean": mean_idx, "max": max_idx}
        shell_rows.append(
            {
                "frame": dat_path.stem,
                "quantity": quantity,
                "selector": "mean",
                "shell_index": mean_idx,
                "shell_radius": float(0.5 * (base_edges[mean_idx] + base_edges[mean_idx + 1])),
                "selector_value": float(mean_arr[mean_idx]),
                "shell_max": float(max_arr[mean_idx]),
                "shell_count": float(stats["count"][mean_idx]),
            }
        )
        shell_rows.append(
            {
                "frame": dat_path.stem,
                "quantity": quantity,
                "selector": "max",
                "shell_index": max_idx,
                "shell_radius": float(0.5 * (base_edges[max_idx] + base_edges[max_idx + 1])),
                "selector_value": float(max_arr[max_idx]),
                "shell_mean": float(mean_arr[max_idx]),
                "shell_count": float(stats["count"][max_idx]),
            }
        )

    surface_maps: Dict[str, Dict[str, np.ndarray]] = {
        "mean": {
            "vabs": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "cs": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "va": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "cfmax": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
        },
        "max": {
            "vabs": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "cs": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "va": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
            "cfmax": np.full((int(domain_nx[1]), int(domain_nx[2])), np.nan, dtype=np.float64),
        },
    }
    surface_sums = {
        selector: {key: np.zeros_like(grid) for key, grid in maps.items()}
        for selector, maps in surface_maps.items()
    }
    surface_counts = {
        selector: {key: np.zeros((int(domain_nx[1]), int(domain_nx[2])), dtype=np.float64) for key in maps}
        for selector, maps in surface_maps.items()
    }

    theta_bins = int(domain_nx[1])
    phi_bins = int(domain_nx[2])
    with dat_path.open("rb") as handle:
        header = get_header(handle)
        block_lvls, block_ixs, block_offsets = get_tree_info(handle)
        name_to_idx = {name: idx for idx, name in enumerate(header["w_names"])}
        for level, bix, offset in zip(block_lvls, block_ixs, block_offsets):
            fields = read_block_fields(
                handle,
                int(offset),
                block_shape,
                int(header["ndim"]),
                int(header["nw"]),
            )
            refine_ratio = 2 ** (int(level) - 1)
            start_r = (int(bix[0]) - 1) * int(block_shape[0])
            start_t = (int(bix[1]) - 1) * int(block_shape[1])
            start_p = (int(bix[2]) - 1) * int(block_shape[2])
            base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)
            base_t_idx = (np.arange(start_t, start_t + int(block_shape[1])) // refine_ratio).astype(np.int64)
            base_p_idx = (np.arange(start_p, start_p + int(block_shape[2])) // refine_ratio).astype(np.int64)

            rho = fields[name_to_idx["rho"]]
            valid = rho > 0.0
            if not np.any(valid):
                continue
            m1 = fields[name_to_idx["m1"]]
            m2 = fields[name_to_idx["m2"]]
            m3 = fields[name_to_idx["m3"]]
            eint = fields[name_to_idx["e"]]
            b_r = fields[name_to_idx["b1"]]
            b_th = fields[name_to_idx["b2"]]
            b_ph = fields[name_to_idx["b3"]]
            v_r, v_th, v_ph, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)
            valid = valid & (pth > 0.0)
            if not np.any(valid):
                continue
            vabs = np.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
            b_sq = b_r * b_r + b_th * b_th + b_ph * b_ph
            c_s_sq = np.where(valid, gamma * pth / np.maximum(rho, EPS), np.nan)
            v_a_sq = np.where(valid, b_sq / np.maximum(rho, EPS), np.nan)
            c_s = np.sqrt(np.maximum(c_s_sq, 0.0))
            v_a = np.sqrt(np.maximum(v_a_sq, 0.0))

            def fast_speed(normal_b: np.ndarray) -> np.ndarray:
                b_n_sq = normal_b * normal_b
                disc = (c_s_sq + v_a_sq) ** 2 - 4.0 * c_s_sq * b_n_sq / np.maximum(rho, EPS)
                disc = np.maximum(disc, 0.0)
                cf_sq = 0.5 * ((c_s_sq + v_a_sq) + np.sqrt(disc))
                return np.sqrt(np.maximum(cf_sq, 0.0))

            cf_max = np.maximum(np.maximum(fast_speed(b_r), fast_speed(b_th)), fast_speed(b_ph))
            quantities = {"vabs": vabs, "cs": c_s, "va": v_a, "cfmax": cf_max}
            for selector in ("mean", "max"):
                for quantity, array in quantities.items():
                    target_shell = shell_selectors[quantity][selector]
                    for i_local, shell_idx in enumerate(base_r_idx):
                        if shell_idx != target_shell:
                            continue
                        for j_local, t_idx in enumerate(base_t_idx):
                            if t_idx < 0 or t_idx >= theta_bins:
                                continue
                            for k_local, p_idx in enumerate(base_p_idx):
                                if p_idx < 0 or p_idx >= phi_bins:
                                    continue
                                value = float(array[i_local, j_local, k_local])
                                if not (math.isfinite(value) and value > 0.0):
                                    continue
                                surface_sums[selector][quantity][t_idx, p_idx] += value
                                surface_counts[selector][quantity][t_idx, p_idx] += 1.0

    for selector in ("mean", "max"):
        for quantity in ("vabs", "cs", "va", "cfmax"):
            count = surface_counts[selector][quantity]
            with np.errstate(divide="ignore", invalid="ignore"):
                surface_maps[selector][quantity] = np.where(
                    count > 0.0,
                    surface_sums[selector][quantity] / count,
                    np.nan,
                )

    shell_lookup = {
        f"{row['quantity']}::{row['selector']}": int(row["shell_index"])
        for row in shell_rows
    }
    return (
        summary,
        top_rows,
        top_vabs_rows,
        top_vabs_inner_rows,
        physical_summary,
        physical_samples,
        shell_rows,
        surface_maps,
        shell_lookup,
    )


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        seen = set()
        for row in rows[1:]:
            for key in row.keys():
                if key not in seen and key not in fieldnames:
                    fieldnames.append(key)
                    seen.add(key)
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_frame_plot(output_dir: Path, frame_row: Dict[str, object]) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path = output_dir / f"{frame_row['frame']}_cfl_breakdown.png"
    vals = [float(frame_row["C_r"]), float(frame_row["C_th"]), float(frame_row["C_ph"])]
    labels = ["C_r", "C_th", "C_ph"]
    colors = ["#1f77b4", "#ff7f0e", "#d62728"]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, vals, color=colors)
    ax.set_ylabel("CFL contribution")
    ax.set_title(f"{frame_row['frame']} bottleneck CFL contributions")
    text = (
        f"dt_pred={float(frame_row['dt_pred']):.3e}\n"
        f"dt_log={float(frame_row['dt_log']):.3e}\n"
        f"r={float(frame_row['r']):.4f}, theta={float(frame_row['theta']):.4e}\n"
        f"phi={float(frame_row['phi']):.4f}, level={int(frame_row['level'])}, block={int(frame_row['block_id'])}\n"
        f"physics={frame_row['dominant_physics']}"
    )
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_speed_distribution_plot(output_dir: Path, frame: str, samples: Dict[str, np.ndarray]) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path = output_dir / f"{frame}_speed_distribution.png"
    series = [
        ("|v|", samples["vabs"], "#1f77b4"),
        ("c_s", samples["cs"], "#ff7f0e"),
        ("v_A", samples["va"], "#2ca02c"),
        ("c_f,max", samples["cfmax"], "#d62728"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes = axes.ravel()
    for ax, (label, values, color) in zip(axes, series):
        finite = values[np.isfinite(values) & (values > 0.0)]
        if finite.size == 0:
            ax.text(0.5, 0.5, "no finite data", ha="center", va="center")
            ax.set_title(label)
            continue
        plot_vals = np.log10(finite)
        vmin = float(np.min(plot_vals))
        vmax = float(np.max(plot_vals))
        if not math.isfinite(vmin) or not math.isfinite(vmax):
            ax.text(0.5, 0.5, "non-finite log data", ha="center", va="center")
            ax.set_title(label)
            continue
        if abs(vmax - vmin) < 1.0e-12:
            ax.axvline(vmin, color=color, linewidth=2.0)
            ax.set_xlim(vmin - 0.5, vmax + 0.5)
        else:
            ax.hist(plot_vals, bins=min(80, max(10, int(np.sqrt(plot_vals.size)))), color=color, alpha=0.85)
        ax.set_title(label)
        ax.set_xlabel("log10(value)")
        ax.set_ylabel("count")
    fig.suptitle(f"{frame} physical speed distributions", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_rtheta_maps(
    output_dir: Path,
    frame: str,
    rtheta_maps: Dict[str, Dict[str, np.ndarray]],
    r_edges: np.ndarray,
    theta_edges_deg: np.ndarray,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None

    out_path = output_dir / f"{frame}_rtheta_maps.png"
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    quantities = [
        ("vabs", "|v|"),
        ("cs", "c_s"),
        ("va", "v_A"),
        ("cfmax", "c_f,max"),
    ]
    R, T = np.meshgrid(r_edges, theta_edges_deg)
    for ax, (key, label) in zip(axes.ravel(), quantities):
        grid = rtheta_maps["max"][key]
        finite = grid[np.isfinite(grid) & (grid > 0.0)]
        if finite.size == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
            ax.set_title(label)
            continue
        dynamic = float(np.max(finite) / max(np.min(finite), EPS))
        if dynamic > 100.0:
            plot_values = np.log10(np.where(grid > 0.0, grid, np.nan))
            cbar_label = f"log10(phi-max {label})"
        else:
            plot_values = grid
            cbar_label = f"phi-max {label}"
        mesh = ax.pcolormesh(R, T, plot_values.T, shading="auto")
        ax.set_xlabel("r")
        ax.set_ylabel("theta [deg]")
        ax.set_title(f"{label} (phi-max over each r-theta bin)")
        fig.colorbar(mesh, ax=ax, pad=0.01, label=cbar_label)
    fig.suptitle(f"{frame} r-theta physical speed maps", fontsize=14)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_bottleneck_scatter(
    output_dir: Path,
    frame: str,
    top_rows: Sequence[Dict[str, object]],
    hotspot_rows: Sequence[Dict[str, object]],
    r_edges: np.ndarray,
    theta_edges_deg: np.ndarray,
    br_map: np.ndarray | None,
    phi_deg: float | None,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    if not top_rows:
        return None

    out_path = output_dir / f"{frame}_bottleneck_scatter.png"
    theta_deg = np.array([float(row["theta"]) * 180.0 / math.pi for row in top_rows])
    radius = np.array([float(row["r"]) for row in top_rows])
    color = np.array([float(row["dt_pred"]) for row in top_rows])
    levels = np.array([int(row["level"]) for row in top_rows])

    fig, ax = plt.subplots(figsize=(7, 5))
    if br_map is not None:
        R, T = np.meshgrid(r_edges, theta_edges_deg)
        finite_br = br_map[np.isfinite(br_map)]
        if finite_br.size > 0:
            vmax = float(np.nanpercentile(np.abs(finite_br), 99.0))
            vmax = max(vmax, EPS)
            mesh = ax.pcolormesh(
                R,
                T,
                br_map.T,
                shading="auto",
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
                alpha=0.82,
            )
            cbar = fig.colorbar(mesh, ax=ax, pad=0.01, label="Br background")
            cbar.ax.tick_params(labelsize=8)
    scatter = ax.scatter(
        radius,
        theta_deg,
        c=np.log10(color),
        cmap="viridis",
        s=50 + 25 * levels,
        alpha=0.95,
        edgecolors="black",
        linewidths=0.5,
    )
    ax.set_xlabel("r")
    ax.set_ylabel("theta [deg]")
    title = f"{frame} top bottleneck cells"
    if phi_deg is not None:
        title += f"\nBr background at phi={phi_deg:.2f} deg"
    ax.set_title(title)
    fig.colorbar(scatter, ax=ax, pad=0.01, label="log10(dt_pred)")
    hotspot_lookup = {
        int(row["cluster_id"]): row
        for row in hotspot_rows
    }
    for cluster_id, row in hotspot_lookup.items():
        ax.annotate(
            f"H{cluster_id}",
            (float(row["r_rep"]), float(row["theta_rep"]) * 180.0 / math.pi),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=9,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
        )
    for row, r, th in zip(top_rows[:10], radius[:10], theta_deg[:10]):
        ax.annotate(
            f"r{int(row['rank'])}",
            (r, th),
            textcoords="offset points",
            xytext=(4, -10),
            fontsize=7,
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_top_vabs_overlay(
    output_dir: Path,
    frame: str,
    top_vabs_rows: Sequence[Dict[str, object]],
    br_surface: np.ndarray | None,
    theta_edges_deg: np.ndarray,
    phi_edges_deg: np.ndarray,
    shell_index_bg: int | None,
    shell_radius_bg: float | None,
    suffix: str = "top_vabs",
    title_prefix: str = "top-|v| cells",
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return None
    if not top_vabs_rows or br_surface is None:
        return None

    out_path = output_dir / f"{frame}_{suffix}_br_overlay.png"
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    X, Y = np.meshgrid(phi_edges_deg, theta_edges_deg)
    finite_br = br_surface[np.isfinite(br_surface)]
    if finite_br.size > 0:
        vmax = float(np.nanpercentile(np.abs(finite_br), 99.0))
        vmax = max(vmax, EPS)
        mesh = ax.pcolormesh(
            X,
            Y,
            br_surface,
            shading="auto",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            alpha=0.86,
        )
        fig.colorbar(mesh, ax=ax, pad=0.01, label="Br")

    phi_deg = np.array([float(row["phi"]) * 180.0 / math.pi for row in top_vabs_rows])
    theta_deg = np.array([float(row["theta"]) * 180.0 / math.pi for row in top_vabs_rows])
    vabs = np.array([float(row["v_abs"]) for row in top_vabs_rows])
    babs = np.array([float(row["b_abs"]) for row in top_vabs_rows])
    shell_indices = np.array([int(row["shell_index"]) for row in top_vabs_rows])
    sizes = 40.0 + 180.0 * babs / max(float(np.nanmax(babs)), EPS)
    sc = ax.scatter(
        phi_deg,
        theta_deg,
        c=np.log10(np.maximum(vabs, EPS)),
        s=sizes,
        cmap="inferno",
        edgecolors="black",
        linewidths=0.5,
        alpha=0.95,
    )
    fig.colorbar(sc, ax=ax, pad=0.01, label="log10(|v|)")
    for rank, row in enumerate(top_vabs_rows[:8], start=1):
        ax.annotate(
            f"v{rank}\n$\\rho$={float(row['rho']):.2e}\n|B|={float(row['b_abs']):.2f}",
            (float(row["phi"]) * 180.0 / math.pi, float(row["theta"]) * 180.0 / math.pi),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=7,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.75),
        )

    shell_min = int(np.min(shell_indices))
    shell_max = int(np.max(shell_indices))
    title = f"{frame} {title_prefix} on Br(theta,phi)"
    if shell_index_bg is not None and shell_radius_bg is not None:
        title += f"\nBr background: shell {shell_index_bg}, r~{shell_radius_bg:.3f}; point shells {shell_min}-{shell_max}"
    ax.set_title(title)
    ax.set_xlabel("phi [deg]")
    ax.set_ylabel("theta [deg]")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def print_summary(summary_rows: Sequence[Dict[str, object]]) -> None:
    print("frame       it      time        dt_log       dt_pred      dir     physics       r        theta    level")
    for row in summary_rows:
        print(
            f"{row['frame']:<10} {int(row['it']):>7d} "
            f"{float(row['time']):>10.6f} {float(row['dt_log']):>11.3e} "
            f"{float(row['dt_pred_min']):>11.3e} {str(row['dominant_direction']):>7} "
            f"{str(row['dominant_physics']):>11} {float(row['r_min']):>8.3f} "
            f"{float(row['theta_min']):>8.4f} {int(row['level']):>5d}"
        )
        dt_log = float(row["dt_log"])
        dt_pred = float(row["dt_pred_min"])
        if math.isfinite(dt_log) and dt_log > 0.0:
            ratio = dt_pred / dt_log
            if ratio > 5.0 or ratio < 0.2:
                print(
                    f"  warning: {row['frame']} dt_pred/dt_log={ratio:.3e}, "
                    "check variable semantics or geometry reconstruction"
                )


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    par_cfg = read_par_config(case_dir / args.par)
    log_by_it = parse_log(case_dir / args.log)

    dat_paths = sorted(case_dir.glob(args.pattern))
    dat_paths = select_frames(dat_paths, args.frames)
    if not dat_paths:
        raise FileNotFoundError(f"No frames matched {args.pattern!r} in {case_dir}")

    base_edges = build_base_radial_edges(
        par_cfg["xprobmin1"],
        par_cfg["xprobmax1"],
        int(par_cfg["domain_nx1"]),
        par_cfg["qstretch_baselevel"],
    )

    summary_rows: List[Dict[str, object]] = []
    top_rows: List[Dict[str, object]] = []
    top_vabs_rows_all: List[Dict[str, object]] = []
    top_vabs_inner_rows_all: List[Dict[str, object]] = []
    physical_rows: List[Dict[str, object]] = []
    shell_rows_all: List[Dict[str, object]] = []
    physical_samples_by_frame: Dict[str, Dict[str, np.ndarray]] = {}
    surface_maps_by_frame: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    rtheta_maps_by_frame: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    br_background_by_frame: Dict[str, np.ndarray] = {}
    phi_background_deg_by_frame: Dict[str, float] = {}
    shell_lookup_by_frame: Dict[str, Dict[str, int]] = {}
    for dat_path in dat_paths:
        summary, top, top_vabs_rows, top_vabs_inner_rows, physical_summary, physical_samples, shell_rows, surface_maps, shell_lookup = analyze_snapshot(
            dat_path, par_cfg, log_by_it, args.topn, base_edges
        )
        summary_rows.append(summary)
        top_rows.extend(top)
        top_vabs_rows_all.extend(top_vabs_rows)
        top_vabs_inner_rows_all.extend(top_vabs_inner_rows)
        physical_rows.append(physical_summary)
        shell_rows_all.extend(shell_rows)
        physical_samples_by_frame[str(summary["frame"])] = physical_samples
        surface_maps_by_frame[str(summary["frame"])] = surface_maps
        shell_lookup_by_frame[str(summary["frame"])] = shell_lookup

    output_dir = case_dir / args.output_dir
    write_csv(output_dir / "dt_bottleneck_summary.csv", summary_rows)
    write_csv(output_dir / "dt_bottleneck_topcells.csv", top_rows)
    write_csv(output_dir / "top_vabs_cells.csv", top_vabs_rows_all)
    write_csv(output_dir / "top_vabs_inner_cells.csv", top_vabs_inner_rows_all)
    hotspot_rows = build_hotspot_rows(top_rows, int(par_cfg["domain_nx3"]))
    write_csv(output_dir / "bottleneck_hotspots.csv", hotspot_rows)
    topn_by_frame = {}
    for row in top_rows:
        topn_by_frame[str(row["frame"])] = topn_by_frame.get(str(row["frame"]), 0) + 1
    write_hotspot_markdown(output_dir / "bottleneck_hotspots.md", hotspot_rows, topn_by_frame)
    write_csv(output_dir / "physical_speed_summary.csv", physical_rows)
    write_csv(output_dir / "surface_shell_summary.csv", shell_rows_all)
    plot_paths: List[Path] = []
    top_by_frame = {}
    for row in top_rows:
        frame = str(row["frame"])
        if frame not in top_by_frame:
            top_by_frame[frame] = row
    top_vabs_by_frame: Dict[str, List[Dict[str, object]]] = {}
    for row in top_vabs_rows_all:
        top_vabs_by_frame.setdefault(str(row["frame"]), []).append(dict(row))
    top_vabs_inner_by_frame: Dict[str, List[Dict[str, object]]] = {}
    for row in top_vabs_inner_rows_all:
        top_vabs_inner_by_frame.setdefault(str(row["frame"]), []).append(dict(row))
    theta_edges_deg = np.linspace(
        par_cfg["xprobmin2"] * 180.0 / math.pi,
        par_cfg["xprobmax2"] * 180.0 / math.pi,
        int(par_cfg["domain_nx2"]) + 1,
    )
    phi_edges_deg = np.linspace(
        par_cfg["xprobmin3"] * 180.0 / math.pi,
        par_cfg["xprobmax3"] * 180.0 / math.pi,
        int(par_cfg["domain_nx3"]) + 1,
    )
    br_surface_vabs_by_frame: Dict[str, np.ndarray] = {}
    br_surface_vabs_inner_by_frame: Dict[str, np.ndarray] = {}
    shell_bg_radius_by_frame: Dict[str, float] = {}
    shell_bg_index_by_frame: Dict[str, int] = {}
    shell_bg_inner_radius_by_frame: Dict[str, float] = {}
    shell_bg_inner_index_by_frame: Dict[str, int] = {}
    # Build r-theta maps, Br(phi-slice), and Br(theta-phi) backgrounds from a lightweight second pass through snapshots.
    for dat_path in dat_paths:
        frame = dat_path.stem
        theta_bins = int(par_cfg["domain_nx2"])
        phi_bins = int(par_cfg["domain_nx3"])
        shape = (int(par_cfg["domain_nx1"]), theta_bins)
        rtheta_stats = {
            "max": {key: init_surface_stats(shape) for key in ("vabs", "cs", "va", "cfmax")},
        }
        br_stats = init_surface_stats(shape)
        br_surface_stats = init_surface_stats((theta_bins, phi_bins))
        br_surface_inner_stats = init_surface_stats((theta_bins, phi_bins))
        top_row = top_by_frame.get(frame)
        phi_target_deg = None
        phi_target_idx = None
        vabs_shell_index = shell_lookup_by_frame.get(frame, {}).get("vabs::max")
        if vabs_shell_index is not None:
            shell_bg_index_by_frame[frame] = int(vabs_shell_index)
            shell_bg_radius_by_frame[frame] = float(0.5 * (base_edges[vabs_shell_index] + base_edges[vabs_shell_index + 1]))
        frame_inner_rows = top_vabs_inner_by_frame.get(frame, [])
        if frame_inner_rows:
            shell_counts: Dict[int, int] = {}
            for row in frame_inner_rows:
                idx = int(row["shell_index"])
                shell_counts[idx] = shell_counts.get(idx, 0) + 1
            inner_shell_index = sorted(shell_counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
            shell_bg_inner_index_by_frame[frame] = int(inner_shell_index)
            shell_bg_inner_radius_by_frame[frame] = float(
                0.5 * (base_edges[inner_shell_index] + base_edges[inner_shell_index + 1])
            )
        if top_row is not None:
            phi_target_deg = float(top_row["phi"]) * 180.0 / math.pi
            phi_target_idx = int(np.digitize([phi_target_deg], phi_edges_deg)[0] - 1)
            phi_target_idx = max(0, min(phi_target_idx, len(phi_edges_deg) - 2))
        with dat_path.open("rb") as handle:
            header = get_header(handle)
            block_lvls, block_ixs, block_offsets = get_tree_info(handle)
            name_to_idx = {name: idx for idx, name in enumerate(header["w_names"])}
            gamma = float(
                header["params"][header["param_names"].index("gamma")]
                if "gamma" in header["param_names"]
                else par_cfg["mhd_gamma"]
            )
            block_shape = header["block_nx"].astype(int)
            for level, bix, offset in zip(block_lvls, block_ixs, block_offsets):
                fields = read_block_fields(handle, int(offset), block_shape, int(header["ndim"]), int(header["nw"]))
                refine_ratio = 2 ** (int(level) - 1)
                start_r = (int(bix[0]) - 1) * int(block_shape[0])
                start_t = (int(bix[1]) - 1) * int(block_shape[1])
                start_p = (int(bix[2]) - 1) * int(block_shape[2])
                base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)
                base_t_idx = (np.arange(start_t, start_t + int(block_shape[1])) // refine_ratio).astype(np.int64)
                base_p_idx = (np.arange(start_p, start_p + int(block_shape[2])) // refine_ratio).astype(np.int64)
                rho = fields[name_to_idx["rho"]]
                m1 = fields[name_to_idx["m1"]]
                m2 = fields[name_to_idx["m2"]]
                m3 = fields[name_to_idx["m3"]]
                eint = fields[name_to_idx["e"]]
                v_r, v_th, v_ph, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)
                valid = (rho > 0.0) & (pth > 0.0)
                if not np.any(valid):
                    continue
                b_r = fields[name_to_idx["b1"]]
                b_th = fields[name_to_idx["b2"]]
                b_ph = fields[name_to_idx["b3"]]
                vabs = np.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
                b_sq = b_r * b_r + b_th * b_th + b_ph * b_ph
                c_s_sq = np.where(valid, gamma * pth / np.maximum(rho, EPS), np.nan)
                v_a_sq = np.where(valid, b_sq / np.maximum(rho, EPS), np.nan)
                c_s = np.sqrt(np.maximum(c_s_sq, 0.0))
                v_a = np.sqrt(np.maximum(v_a_sq, 0.0))
                def fast_speed(normal_b: np.ndarray) -> np.ndarray:
                    b_n_sq = normal_b * normal_b
                    disc = (c_s_sq + v_a_sq) ** 2 - 4.0 * c_s_sq * b_n_sq / np.maximum(rho, EPS)
                    disc = np.maximum(disc, 0.0)
                    cf_sq = 0.5 * ((c_s_sq + v_a_sq) + np.sqrt(disc))
                    return np.sqrt(np.maximum(cf_sq, 0.0))
                cfmax = np.maximum(np.maximum(fast_speed(b_r), fast_speed(b_th)), fast_speed(b_ph))
                row_idx = np.broadcast_to(base_r_idx[:, None, None], rho.shape)
                col_idx = np.broadcast_to(base_t_idx[None, :, None], rho.shape)
                update_surface_stats(rtheta_stats["max"]["vabs"], row_idx, col_idx, vabs)
                update_surface_stats(rtheta_stats["max"]["cs"], row_idx, col_idx, c_s)
                update_surface_stats(rtheta_stats["max"]["va"], row_idx, col_idx, v_a)
                update_surface_stats(rtheta_stats["max"]["cfmax"], row_idx, col_idx, cfmax)
                if phi_target_idx is not None:
                    phi_mask = np.broadcast_to((base_p_idx == phi_target_idx)[None, None, :], rho.shape)
                    br_values = np.where(phi_mask, b_r, np.nan)
                    update_surface_stats(br_stats, row_idx, col_idx, br_values)
                if vabs_shell_index is not None:
                    shell_mask = np.broadcast_to((base_r_idx == int(vabs_shell_index))[:, None, None], rho.shape)
                    theta_idx = np.broadcast_to(base_t_idx[None, :, None], rho.shape)
                    phi_idx = np.broadcast_to(base_p_idx[None, None, :], rho.shape)
                    br_values_surface = np.where(shell_mask, b_r, np.nan)
                    update_surface_stats(br_surface_stats, theta_idx, phi_idx, br_values_surface)
                inner_shell_index = shell_bg_inner_index_by_frame.get(frame)
                if inner_shell_index is not None:
                    shell_mask_inner = np.broadcast_to((base_r_idx == int(inner_shell_index))[:, None, None], rho.shape)
                    theta_idx = np.broadcast_to(base_t_idx[None, :, None], rho.shape)
                    phi_idx = np.broadcast_to(base_p_idx[None, None, :], rho.shape)
                    br_values_surface_inner = np.where(shell_mask_inner, b_r, np.nan)
                    update_surface_stats(br_surface_inner_stats, theta_idx, phi_idx, br_values_surface_inner)
        rtheta_maps_by_frame[frame] = {
            "max": {
                key: np.where(stats["count"] > 0.0, stats["max"], np.nan)
                for key, stats in rtheta_stats["max"].items()
            }
        }
        br_background_by_frame[frame] = np.where(br_stats["count"] > 0.0, br_stats["sum"] / np.maximum(br_stats["count"], 1.0), np.nan)
        br_surface_vabs_by_frame[frame] = np.where(
            br_surface_stats["count"] > 0.0,
            br_surface_stats["sum"] / np.maximum(br_surface_stats["count"], 1.0),
            np.nan,
        )
        br_surface_vabs_inner_by_frame[frame] = np.where(
            br_surface_inner_stats["count"] > 0.0,
            br_surface_inner_stats["sum"] / np.maximum(br_surface_inner_stats["count"], 1.0),
            np.nan,
        )
        if phi_target_deg is not None:
            phi_background_deg_by_frame[frame] = phi_target_deg
    for frame in sorted(top_by_frame):
        plot_path = write_frame_plot(output_dir, top_by_frame[frame])
        if plot_path is not None:
            plot_paths.append(plot_path)
        speed_plot = write_speed_distribution_plot(output_dir, frame, physical_samples_by_frame[frame])
        if speed_plot is not None:
            plot_paths.append(speed_plot)
        rtheta_plot = write_rtheta_maps(output_dir, frame, rtheta_maps_by_frame[frame], base_edges, theta_edges_deg)
        if rtheta_plot is not None:
            plot_paths.append(rtheta_plot)
        frame_top_rows = [dict(row) for row in top_rows if str(row["frame"]) == frame]
        frame_top_rows.sort(key=lambda row: float(row["dt_pred"]))
        for rank, row in enumerate(frame_top_rows, start=1):
            row["rank"] = rank
        frame_hotspot_rows = [row for row in hotspot_rows if str(row["frame"]) == frame]
        scatter_plot = write_bottleneck_scatter(
            output_dir,
            frame,
            frame_top_rows,
            frame_hotspot_rows,
            base_edges,
            theta_edges_deg,
            br_background_by_frame.get(frame),
            phi_background_deg_by_frame.get(frame),
        )
        if scatter_plot is not None:
            plot_paths.append(scatter_plot)
        frame_top_vabs_rows = top_vabs_by_frame.get(frame, [])
        frame_top_vabs_rows.sort(key=lambda row: float(row["v_abs"]), reverse=True)
        vabs_overlay = write_top_vabs_overlay(
            output_dir,
            frame,
            frame_top_vabs_rows,
            br_surface_vabs_by_frame.get(frame),
            theta_edges_deg,
            phi_edges_deg,
            shell_bg_index_by_frame.get(frame),
            shell_bg_radius_by_frame.get(frame),
        )
        if vabs_overlay is not None:
            plot_paths.append(vabs_overlay)
        frame_top_vabs_inner_rows = top_vabs_inner_by_frame.get(frame, [])
        frame_top_vabs_inner_rows.sort(key=lambda row: float(row["v_abs"]), reverse=True)
        vabs_inner_overlay = write_top_vabs_overlay(
            output_dir,
            frame,
            frame_top_vabs_inner_rows,
            br_surface_vabs_inner_by_frame.get(frame),
            theta_edges_deg,
            phi_edges_deg,
            shell_bg_inner_index_by_frame.get(frame),
            shell_bg_inner_radius_by_frame.get(frame),
            suffix="top_vabs_inner",
            title_prefix=f"top-|v| cells (r<={TOP_VABS_INNER_RMAX:g})",
        )
        if vabs_inner_overlay is not None:
            plot_paths.append(vabs_inner_overlay)
        frame_shell_rows = [row for row in shell_rows_all if str(row["frame"]) == frame]
        for selector in ("mean", "max"):
            surface_plot = plot_surface_maps(
                output_dir,
                frame,
                selector,
                surface_maps_by_frame[frame][selector],
                [row for row in frame_shell_rows if str(row["selector"]) == selector],
                theta_edges_deg,
                phi_edges_deg,
            )
            if surface_plot is not None:
                plot_paths.append(surface_plot)
    print_summary(summary_rows)
    print(f"\nWrote {output_dir / 'dt_bottleneck_summary.csv'}")
    print(f"Wrote {output_dir / 'dt_bottleneck_topcells.csv'}")
    print(f"Wrote {output_dir / 'top_vabs_cells.csv'}")
    print(f"Wrote {output_dir / 'top_vabs_inner_cells.csv'}")
    print(f"Wrote {output_dir / 'bottleneck_hotspots.csv'}")
    print(f"Wrote {output_dir / 'bottleneck_hotspots.md'}")
    print(f"Wrote {output_dir / 'physical_speed_summary.csv'}")
    print(f"Wrote {output_dir / 'surface_shell_summary.csv'}")
    for plot_path in plot_paths:
        print(f"Wrote {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
