# Copied to project analysis task toolbox
# Source case path: off_2270_initwind/off_2270_initwind_hpc/analysis/dt_cfl_common.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_initwind/off_2270_initwind_hpc/analysis/dt_cfl_common.py

#!/usr/bin/env python3
"""Common helpers for DT/CFL diagnostics on MPI-AMRVAC DAT snapshots."""

from __future__ import annotations

import csv
import math
import re
import struct
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

SIZE_INT = 4
NAME_LEN = 16
ALIGN = "="
EPS = 1.0e-30


def parse_fortran_float(raw: str) -> float:
    return float(raw.replace("d", "e").replace("D", "e"))


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
    cfg: Dict[str, float] = {}
    for key in keys:
        m = re.search(rf"{re.escape(key)}\s*=\s*([-+0-9.dDeE]+)", text, flags=re.IGNORECASE)
        if not m:
            raise ValueError(f"Could not find '{key}' in {par_path}")
        cfg[key] = parse_fortran_float(m.group(1))
    return cfg


def parse_log(log_path: Path) -> Dict[int, Dict[str, float]]:
    by_it: Dict[int, Dict[str, float]] = {}
    if not log_path.exists():
        return by_it
    with log_path.open() as f:
        for raw in f:
            left = raw.split("|", 1)[0].strip()
            if not left or left.startswith("it global_time"):
                continue
            parts = left.split()
            if len(parts) < 3:
                continue
            try:
                it = int(parts[0])
                time = float(parts[1])
                dt = float(parts[2])
            except ValueError:
                continue

            row: Dict[str, float] = {
                "time": time,
                "dt": dt,
            }
            # AMRVAC logs usually place n1/n2 in columns 15/16 (0-based: 14/15),
            # but keep this optional for portability across log styles.
            if len(parts) >= 16:
                try:
                    row["n1"] = float(int(float(parts[14])))
                    row["n2"] = float(int(float(parts[15])))
                except ValueError:
                    pass
            by_it[it] = row
    return by_it


def infer_angle_mode(xprobmax2: float, angle_mode: str) -> str:
    mode = angle_mode.lower().strip()
    if mode == "auto":
        return "normalized_2pi" if xprobmax2 <= 0.5 + 1.0e-12 else "radian"
    if mode not in {"normalized_2pi", "radian"}:
        raise ValueError("angle_mode must be one of: auto, normalized_2pi, radian")
    return mode


def get_angle_scales(par_cfg: Dict[str, float], angle_mode: str) -> Tuple[str, float, float]:
    mode_used = infer_angle_mode(par_cfg["xprobmax2"], angle_mode)
    if mode_used == "normalized_2pi":
        theta_scale = 2.0 * math.pi
        phi_scale = 2.0 * math.pi if par_cfg["xprobmax3"] <= 1.0 + 1.0e-12 else 1.0
    else:
        theta_scale = 1.0
        phi_scale = 1.0
    return mode_used, theta_scale, phi_scale


def get_header(handle) -> Dict[str, object]:
    handle.seek(0)
    h: Dict[str, object] = {}
    (datver,) = struct.unpack(ALIGN + "i", handle.read(4))
    if datver < 3:
        raise OSError(f"Unsupported datfile version: {datver}")
    h["datfile_version"] = datver

    hdr = struct.unpack(ALIGN + 9 * "i" + "d", handle.read(struct.calcsize(ALIGN + 9 * "i" + "d")))
    (
        h["offset_tree"],
        h["offset_blocks"],
        h["nw"],
        h["ndir"],
        h["ndim"],
        h["levmax"],
        h["nleafs"],
        h["nparents"],
        h["it"],
        h["time"],
    ) = hdr

    ndim = int(h["ndim"])
    fmt_d = ALIGN + ndim * "d"
    h["xmin"] = np.array(struct.unpack(fmt_d, handle.read(struct.calcsize(fmt_d))))
    h["xmax"] = np.array(struct.unpack(fmt_d, handle.read(struct.calcsize(fmt_d))))

    fmt_i = ALIGN + ndim * "i"
    h["domain_nx"] = np.array(struct.unpack(fmt_i, handle.read(struct.calcsize(fmt_i))))
    h["block_nx"] = np.array(struct.unpack(fmt_i, handle.read(struct.calcsize(fmt_i))))

    if datver >= 4:
        handle.read(struct.calcsize(fmt_i))  # periodic
        handle.read(NAME_LEN)  # geometry
        handle.read(4)  # staggered

    w_names: List[str] = []
    for _ in range(int(h["nw"])):
        token = handle.read(NAME_LEN).decode(errors="ignore").replace("\x00", "")
        w_names.append(token.strip())
    h["w_names"] = w_names

    h["physics_type"] = handle.read(NAME_LEN).decode(errors="ignore").strip()
    (n_par,) = struct.unpack(ALIGN + "i", handle.read(4))
    h["n_par"] = n_par
    if n_par > 0:
        params = struct.unpack(ALIGN + n_par * "d", handle.read(8 * n_par))
        names_raw = handle.read(NAME_LEN * n_par)
        names: List[str] = []
        for i in range(n_par):
            chunk = names_raw[i * NAME_LEN : (i + 1) * NAME_LEN]
            token = chunk.decode(errors="ignore").replace("\x00", "")
            names.append(token.strip())
        h["params"] = np.array(params)
        h["param_names"] = names
    else:
        h["params"] = np.array([], dtype=np.float64)
        h["param_names"] = []

    # snapshotnext/slicenext/collapsenext
    handle.read(12)
    return h


def get_tree_info(handle) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    h = get_header(handle)
    nleafs = int(h["nleafs"])
    nparents = int(h["nparents"])
    ndim = int(h["ndim"])

    handle.seek(int(h["offset_tree"]))
    handle.read((nleafs + nparents) * SIZE_INT)  # leaf/logical array

    lv = np.array(struct.unpack(ALIGN + nleafs * "i", handle.read(nleafs * SIZE_INT)), dtype=np.int64)
    bix = np.array(
        struct.unpack(ALIGN + (nleafs * ndim) * "i", handle.read(nleafs * ndim * SIZE_INT)),
        dtype=np.int64,
    ).reshape((nleafs, ndim))
    off = np.array(struct.unpack(ALIGN + nleafs * "q", handle.read(nleafs * 8)), dtype=np.int64)
    return lv, bix, off


def read_block_fields_selected(
    handle,
    offset: int,
    block_shape: np.ndarray,
    ndim: int,
    field_map: Dict[str, int],
) -> Dict[str, np.ndarray]:
    count = int(np.prod(block_shape))
    byte_size = count * 8
    out: Dict[str, np.ndarray] = {}
    for name, idx in field_map.items():
        handle.seek(offset + 2 * ndim * SIZE_INT + idx * byte_size)
        data = np.fromfile(handle, dtype="=f8", count=count)
        if data.size != count:
            raise IOError(f"Failed to read field '{name}' at offset {offset}")
        arr = data.reshape(tuple(block_shape[::-1]), order="C").T
        while arr.ndim < 3:
            arr = arr[..., np.newaxis]
        out[name] = arr
    return out


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


def block_radial_edges(base_edges: np.ndarray, level: int, block_ix: int, block_nx: int) -> np.ndarray:
    refine_ratio = 2 ** (level - 1)
    start = (block_ix - 1) * block_nx
    out = np.empty(block_nx + 1, dtype=np.float64)
    nbase = len(base_edges) - 1
    max_fine = nbase * refine_ratio
    for local_edge, fine_idx in enumerate(range(start, start + block_nx + 1)):
        if fine_idx >= max_fine:
            out[local_edge] = base_edges[-1]
            continue
        parent = fine_idx // refine_ratio
        sub = fine_idx % refine_ratio
        width = base_edges[parent + 1] - base_edges[parent]
        out[local_edge] = base_edges[parent] + width * (sub / refine_ratio)
    return out


def conservative_to_primitive(
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
    dom_idx = int(np.nanargmax(components))
    dom_dir = directions[dom_idx]

    lengths = np.array([l_r, l_th, l_ph], dtype=np.float64)
    inv_lengths = 1.0 / np.maximum(lengths, EPS)
    geom_ratio = inv_lengths[dom_idx] / np.maximum(np.median(inv_lengths), EPS)
    speed_ratio = (v_n + c_fast_n) / np.maximum(np.median([v_n + c_fast_n, c_s + v_a, EPS]), EPS)

    if geom_ratio > max(1.5, 1.2 * speed_ratio):
        return dom_dir, "geometry"
    if v_n >= c_fast_n:
        return dom_dir, "flow"
    if c_s >= 0.9 * v_a:
        return dom_dir, "sound"
    return dom_dir, "alfven/fast"


def safe_float(x: object) -> float:
    if x is None:
        return math.nan
    try:
        return float(x)
    except Exception:
        return math.nan


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: List[str] = list(rows[0].keys())
    seen = set(fieldnames)
    for row in rows[1:]:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
