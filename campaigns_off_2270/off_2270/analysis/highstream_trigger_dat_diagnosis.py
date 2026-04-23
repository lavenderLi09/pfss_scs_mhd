#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

CASE_DIR = Path(__file__).resolve().parent.parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))

from hao_code.datfile_io import get_header, get_tree_info  # noqa: E402

EPS = 1.0e-30


@dataclass
class FrameMetrics:
    frame: str
    it: int
    time: float
    high_any: int
    high_frac_global: float
    comp_size: int
    comp_rmin: float
    comp_r95: float
    anchored_inner115: int


class PrioritySampler:
    """Keep an unbiased approximate sample using top-k random priorities."""

    def __init__(self, max_size: int, rng: np.random.Generator):
        self.max_size = max(1, int(max_size))
        self.rng = rng
        self._vals: List[np.ndarray] = []
        self._pri: List[np.ndarray] = []
        self._n_buffer = 0

    def add(self, values: np.ndarray) -> None:
        if values.size == 0:
            return
        vals = values[np.isfinite(values)]
        if vals.size == 0:
            return
        pri = self.rng.random(vals.size)
        if vals.size > self.max_size:
            keep = np.argpartition(pri, -self.max_size)[-self.max_size :]
            vals = vals[keep]
            pri = pri[keep]
        self._vals.append(vals)
        self._pri.append(pri)
        self._n_buffer += vals.size
        if self._n_buffer > 4 * self.max_size:
            self._compress()

    def _compress(self) -> None:
        if not self._vals:
            return
        vals = np.concatenate(self._vals)
        pri = np.concatenate(self._pri)
        if vals.size > self.max_size:
            keep = np.argpartition(pri, -self.max_size)[-self.max_size :]
            vals = vals[keep]
            pri = pri[keep]
        self._vals = [vals]
        self._pri = [pri]
        self._n_buffer = vals.size

    def values(self) -> np.ndarray:
        self._compress()
        if not self._vals:
            return np.empty(0, dtype=np.float64)
        return self._vals[0]


def _to_float(token: str) -> float:
    return float(token.replace("d", "e").replace("D", "e"))


def read_par_config(par_path: Path) -> Dict[str, float]:
    text = par_path.read_text()

    def req(key: str) -> float:
        m = re.search(rf"{re.escape(key)}\s*=\s*([-+0-9.dDeE]+)", text, flags=re.IGNORECASE)
        if not m:
            raise ValueError(f"Could not find {key} in {par_path}")
        return _to_float(m.group(1))

    def opt(key: str, default: float) -> float:
        m = re.search(rf"{re.escape(key)}\s*=\s*([-+0-9.dDeE]+)", text, flags=re.IGNORECASE)
        return _to_float(m.group(1)) if m else default

    cfg: Dict[str, float] = {
        "domain_nx1": req("domain_nx1"),
        "domain_nx2": req("domain_nx2"),
        "domain_nx3": req("domain_nx3"),
        "block_nx1": req("block_nx1"),
        "block_nx2": req("block_nx2"),
        "block_nx3": req("block_nx3"),
        "xprobmin1": req("xprobmin1"),
        "xprobmin2": req("xprobmin2"),
        "xprobmin3": req("xprobmin3"),
        "xprobmax1": req("xprobmax1"),
        "xprobmax2": req("xprobmax2"),
        "xprobmax3": req("xprobmax3"),
        "qstretch_baselevel": req("qstretch_baselevel"),
        "courantpar": req("courantpar"),
        "mhd_gamma": req("mhd_gamma"),
        "small_density": opt("small_density", math.nan),
        "small_pressure": opt("small_pressure", math.nan),
    }
    return cfg


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
    from hao_code.datfile_io import SIZE_INT

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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Diagnose high-speed stream generation/accumulation from dat only.")
    p.add_argument("--case-dir", type=Path, default=CASE_DIR)
    p.add_argument("--par", default="amrvac.par")
    p.add_argument("--mod-usr", default="mod_usr.t")
    p.add_argument("--pattern", default="data/off*.dat")
    p.add_argument("--start", type=int, default=None)
    p.add_argument("--end", type=int, default=None)
    p.add_argument("--stride", type=int, default=1)

    p.add_argument("--vfix", type=float, default=6.0, help="Fixed speed threshold in code unit (~700 km/s if unit~116.5km/s).")
    p.add_argument("--frac-thresh", type=float, default=0.005, help="Per-shell onset threshold of high-speed fraction.")
    p.add_argument("--inner-anchor-r", type=float, default=1.15)
    p.add_argument("--near-inner-r", type=float, default=1.25)
    p.add_argument("--outer-check-r", type=float, default=1.60)
    p.add_argument("--growth-frames", type=int, default=20)
    p.add_argument("--min-comp-cells", type=int, default=30)

    p.add_argument("--root-rmin", type=float, default=1.03)
    p.add_argument("--root-rmax", type=float, default=1.20)
    p.add_argument("--gravity-g0", type=float, default=-14.072, help="g(r)=g0*(SRadius^2/r^2), in code unit.")
    p.add_argument("--sradius", type=float, default=1.0)

    p.add_argument("--sample-max", type=int, default=50000, help="Max sample size per metric/group/frame for quantile estimation.")
    p.add_argument("--floor-factor", type=float, default=10.0, help="Near-floor threshold factor: floor*factor.")
    p.add_argument("--persist-frames", type=int, default=3, help="Persistence required for first-anomaly ranking.")
    p.add_argument("--make-plots", action="store_true", default=True)
    p.add_argument("--no-make-plots", dest="make_plots", action="store_false")

    p.add_argument("--output-prefix", default="analysis/highstream_trigger_diagnosis")
    p.add_argument("--seed", type=int, default=20260413)
    return p.parse_args()


def frame_id(stem: str) -> int | None:
    m = re.search(r"(\d+)$", stem)
    return int(m.group(1)) if m else None


def select_paths(paths: Sequence[Path], start: int | None, end: int | None, stride: int) -> List[Path]:
    out: List[Path] = []
    base = start if start is not None else 0
    for p in sorted(paths):
        n = frame_id(p.stem)
        if n is None:
            continue
        if start is not None and n < start:
            continue
        if end is not None and n > end:
            continue
        if stride > 1 and ((n - base) % stride != 0):
            continue
        out.append(p)
    return out


def _linear_idx(r: int, t: int, p: int, nt: int, np_: int) -> int:
    return (r * nt + t) * np_ + p


def _unpack_idx(idx: int, nt: int, np_: int) -> Tuple[int, int, int]:
    r = idx // (nt * np_)
    x = idx % (nt * np_)
    t = x // np_
    p = x % np_
    return r, t, p


def largest_component(mask: np.ndarray) -> np.ndarray:
    coords = np.argwhere(mask)
    if coords.size == 0:
        return np.zeros((0, 3), dtype=np.int32)

    nr, nt, np_ = mask.shape
    active = {_linear_idx(int(r), int(t), int(p), nt, np_) for r, t, p in coords}
    visited = set()
    best: List[int] = []

    while active:
        seed = active.pop()
        q = deque([seed])
        comp = [seed]
        visited.add(seed)

        while q:
            cur = q.popleft()
            r, t, p = _unpack_idx(cur, nt, np_)
            neigh = []
            if r > 0:
                neigh.append((r - 1, t, p))
            if r < nr - 1:
                neigh.append((r + 1, t, p))
            if t > 0:
                neigh.append((r, t - 1, p))
            if t < nt - 1:
                neigh.append((r, t + 1, p))
            neigh.append((r, t, (p - 1) % np_))
            neigh.append((r, t, (p + 1) % np_))

            for rr, tt, pp in neigh:
                nid = _linear_idx(rr, tt, pp, nt, np_)
                if nid in visited:
                    continue
                if nid in active:
                    active.remove(nid)
                    visited.add(nid)
                    q.append(nid)
                    comp.append(nid)

        if len(comp) > len(best):
            best = comp

    out = np.empty((len(best), 3), dtype=np.int32)
    for i, idx in enumerate(best):
        out[i] = _unpack_idx(idx, nt, np_)
    return out


def rankdata(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=np.float64)
    i = 0
    n = len(x)
    while i < n:
        j = i
        while j + 1 < n and x[order[j + 1]] == x[order[i]]:
            j += 1
        rank = 0.5 * (i + j) + 1.0
        ranks[order[i : j + 1]] = rank
        i = j + 1
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 3 or y.size < 3:
        return math.nan
    rx = rankdata(x)
    ry = rankdata(y)
    xc = rx - rx.mean()
    yc = ry - ry.mean()
    den = math.sqrt(float((xc * xc).sum()) * float((yc * yc).sum()))
    if den <= 0:
        return math.nan
    return float((xc * yc).sum() / den)


def _ensure_sampler(d: Dict[str, PrioritySampler], key: str, max_size: int, rng: np.random.Generator) -> PrioritySampler:
    if key not in d:
        d[key] = PrioritySampler(max_size=max_size, rng=rng)
    return d[key]


def _qstats(vals: np.ndarray) -> Dict[str, float]:
    if vals.size == 0:
        return {
            "count": 0.0,
            "p05": math.nan,
            "p50": math.nan,
            "p95": math.nan,
            "mean": math.nan,
            "min": math.nan,
            "max": math.nan,
        }
    q05, q50, q95 = np.quantile(vals, [0.05, 0.5, 0.95])
    return {
        "count": float(vals.size),
        "p05": float(q05),
        "p50": float(q50),
        "p95": float(q95),
        "mean": float(np.mean(vals)),
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
    }


def _add_group_stats(row: Dict[str, float | str | int], prefix: str, stats: Dict[str, float]) -> None:
    row[f"{prefix}_count"] = stats["count"]
    row[f"{prefix}_p05"] = stats["p05"]
    row[f"{prefix}_p50"] = stats["p50"]
    row[f"{prefix}_p95"] = stats["p95"]
    row[f"{prefix}_mean"] = stats["mean"]
    row[f"{prefix}_min"] = stats["min"]
    row[f"{prefix}_max"] = stats["max"]


def _gradient_axis(a: np.ndarray, x: np.ndarray, axis: int) -> np.ndarray:
    if a.shape[axis] < 2:
        return np.zeros_like(a)
    return np.gradient(a, x, axis=axis, edge_order=1)


def _safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den):
        return math.nan
    return num / den if abs(den) > EPS else math.nan


def _first_persistent(cond: np.ndarray, persist: int) -> int:
    n = cond.size
    if n == 0:
        return -1
    p = max(1, int(persist))
    for i in range(0, n - p + 1):
        if np.all(cond[i : i + p]):
            return i
    return -1


def analyze_frame(
    dat_path: Path,
    par_cfg: Dict[str, float],
    base_edges: np.ndarray,
    vfix: float,
    inner_anchor_r: float,
    min_comp_cells: int,
    root_rmin: float,
    root_rmax: float,
    gravity_g0: float,
    sradius: float,
    sample_max: int,
    floor_factor: float,
    rng: np.random.Generator,
) -> Tuple[FrameMetrics, np.ndarray, np.ndarray, Dict[str, float | str | int], Dict[str, float | str | int], Dict[str, float | str | int]]:
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

        domain_nx = header["domain_nx"].astype(int)
        block_shape = header["block_nx"].astype(int)
        nr, nt, np_ = int(domain_nx[0]), int(domain_nx[1]), int(domain_nx[2])

        shell_count = np.zeros(nr, dtype=np.float64)
        shell_sum = np.zeros(nr, dtype=np.float64)
        shell_sumsq = np.zeros(nr, dtype=np.float64)

        xmins = [par_cfg[f"xprobmin{i}"] for i in range(1, 4)]
        xmaxs = [par_cfg[f"xprobmax{i}"] for i in range(1, 4)]

        # pass 1: shell speed stats
        for level, bix, offset in zip(block_lvls, block_ixs, block_offsets):
            fields = read_block_fields(handle, int(offset), block_shape, int(header["ndim"]), int(header["nw"]))
            refine_ratio = 2 ** (int(level) - 1)
            start_r = (int(bix[0]) - 1) * int(block_shape[0])
            base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)

            rho = fields[name_to_idx["rho"]]
            m1 = fields[name_to_idx["m1"]]
            m2 = fields[name_to_idx["m2"]]
            m3 = fields[name_to_idx["m3"]]
            eint = fields[name_to_idx["e"]]
            v_r, v_th, v_ph, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)
            valid = (rho > 0.0) & (pth > 0.0)
            if not np.any(valid):
                continue

            vabs = np.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
            for i_local, ridx in enumerate(base_r_idx):
                vals = vabs[i_local][valid[i_local]]
                if vals.size == 0:
                    continue
                shell_count[ridx] += float(vals.size)
                shell_sum[ridx] += float(vals.sum())
                shell_sumsq[ridx] += float((vals * vals).sum())

        shell_mean = np.divide(shell_sum, np.maximum(shell_count, 1.0))
        shell_var = np.divide(shell_sumsq, np.maximum(shell_count, 1.0)) - shell_mean * shell_mean
        shell_std = np.sqrt(np.maximum(shell_var, 0.0))
        shell_thr = np.maximum(vfix, shell_mean + 3.0 * shell_std)

        # pass 2: high-speed mask on base grid
        shell_high = np.zeros(nr, dtype=np.float64)
        base_mask = np.zeros((nr, nt, np_), dtype=np.bool_)

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
            vabs = np.sqrt(v_r * v_r + v_th * v_th + v_ph * v_ph)
            valid = (rho > 0.0) & (pth > 0.0)

            thr_3d = shell_thr[base_r_idx][:, None, None]
            high = valid & np.isfinite(vabs) & (vabs > thr_3d) & (v_r > 0.0)

            for i_local, ridx in enumerate(base_r_idx):
                shell_high[ridx] += float(np.count_nonzero(high[i_local]))

            if np.any(high):
                ridx_3d = np.broadcast_to(base_r_idx[:, None, None], high.shape)
                tidx_3d = np.broadcast_to(base_t_idx[None, :, None], high.shape)
                pidx_3d = np.broadcast_to(base_p_idx[None, None, :], high.shape)
                rr = ridx_3d[high]
                tt = tidx_3d[high]
                pp = pidx_3d[high]
                base_mask[rr, tt, pp] = True

        shell_frac = np.divide(shell_high, np.maximum(shell_count, 1.0))

        comp = largest_component(base_mask)
        r_centers = 0.5 * (base_edges[:-1] + base_edges[1:])
        comp_base_mask = np.zeros((nr, nt, np_), dtype=np.bool_)
        if comp.size > 0:
            comp_base_mask[comp[:, 0], comp[:, 1], comp[:, 2]] = True

        if comp.shape[0] >= min_comp_cells:
            rvals = r_centers[comp[:, 0]]
            comp_rmin = float(np.min(rvals))
            comp_r95 = float(np.quantile(rvals, 0.95))
            comp_size = int(comp.shape[0])
        else:
            comp_rmin = math.nan
            comp_r95 = math.nan
            comp_size = int(comp.shape[0])

        metrics = FrameMetrics(
            frame=dat_path.stem,
            it=int(header["it"]),
            time=float(header["time"]),
            high_any=int(np.any(base_mask)),
            high_frac_global=float(np.mean(base_mask)),
            comp_size=comp_size,
            comp_rmin=comp_rmin,
            comp_r95=comp_r95,
            anchored_inner115=int(comp_size >= min_comp_cells and comp_rmin <= inner_anchor_r),
        )

        # pass 3: thermo chain, force budget, boundary health
        thermo_q = ["e", "pth", "rho", "beta", "cs", "va", "babs", "vr", "vabs"]
        force_q = ["f_pres_r", "jxb_r", "rho_g_r", "f_res_r", "divv"]
        samplers: Dict[str, PrioritySampler] = {}

        rho_floor = par_cfg.get("small_density", math.nan)
        p_floor = par_cfg.get("small_pressure", math.nan)
        rho_floor_thr = rho_floor * floor_factor if math.isfinite(rho_floor) else math.nan
        p_floor_thr = p_floor * floor_factor if math.isfinite(p_floor) else math.nan

        seg_bounds = {
            "inner": (0.0, 1.15),
            "mid": (1.15, 1.8),
            "outer": (1.8, 1.0e9),
        }
        seg_counts = {k: 0 for k in seg_bounds}
        seg_floor_rho = {k: 0 for k in seg_bounds}
        seg_floor_p = {k: 0 for k in seg_bounds}

        n_valid_global = 0
        n_floor_rho_global = 0
        n_floor_p_global = 0
        rho_min_global = math.inf
        p_min_global = math.inf

        for level, bix, offset in zip(block_lvls, block_ixs, block_offsets):
            fields = read_block_fields(handle, int(offset), block_shape, int(header["ndim"]), int(header["nw"]))
            r_edges = block_radial_edges(base_edges, int(level), int(bix[0]), int(block_shape[0]))
            th_edges = block_uniform_edges(
                xmins[1], xmaxs[1], int(domain_nx[1]), int(level), int(bix[1]), int(block_shape[1])
            )
            ph_edges = block_uniform_edges(
                xmins[2], xmaxs[2], int(domain_nx[2]), int(level), int(bix[2]), int(block_shape[2])
            )
            r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
            th_cent = 0.5 * (th_edges[:-1] + th_edges[1:])
            ph_cent = 0.5 * (ph_edges[:-1] + ph_edges[1:])
            rr, tt, _pp = np.meshgrid(r_cent, th_cent, ph_cent, indexing="ij")

            refine_ratio = 2 ** (int(level) - 1)
            start_r = (int(bix[0]) - 1) * int(block_shape[0])
            start_t = (int(bix[1]) - 1) * int(block_shape[1])
            start_p = (int(bix[2]) - 1) * int(block_shape[2])
            base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)
            base_t_idx = (np.arange(start_t, start_t + int(block_shape[1])) // refine_ratio).astype(np.int64)
            base_p_idx = (np.arange(start_p, start_p + int(block_shape[2])) // refine_ratio).astype(np.int64)

            in_comp = comp_base_mask[np.ix_(base_r_idx, base_t_idx, base_p_idx)]

            rho = fields[name_to_idx["rho"]]
            m1 = fields[name_to_idx["m1"]]
            m2 = fields[name_to_idx["m2"]]
            m3 = fields[name_to_idx["m3"]]
            eint = fields[name_to_idx["e"]]
            b_r = fields[name_to_idx["b1"]]
            b_t = fields[name_to_idx["b2"]]
            b_p = fields[name_to_idx["b3"]]

            v_r, v_t, v_p, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)
            valid = (rho > 0.0) & (pth > 0.0) & np.isfinite(v_r)

            b2 = b_r * b_r + b_t * b_t + b_p * b_p
            babs = np.sqrt(np.maximum(b2, 0.0))
            cs = np.sqrt(np.maximum(gamma * pth / np.maximum(rho, EPS), 0.0))
            va = np.sqrt(np.maximum(b2 / np.maximum(rho, EPS), 0.0))
            beta = 2.0 * pth / np.maximum(b2, EPS)
            vabs = np.sqrt(v_r * v_r + v_t * v_t + v_p * v_p)

            # anomaly component vs same-r envelope background
            if math.isfinite(comp_rmin) and math.isfinite(comp_r95):
                env = (rr >= comp_rmin) & (rr <= comp_r95)
            else:
                env = (rr >= root_rmin) & (rr <= root_rmax)

            mask_comp = valid & in_comp
            mask_bg = valid & (~in_comp) & env

            thermo_arrays = {
                "e": eint,
                "pth": pth,
                "rho": rho,
                "beta": beta,
                "cs": cs,
                "va": va,
                "babs": babs,
                "vr": v_r,
                "vabs": vabs,
            }
            for qn, qa in thermo_arrays.items():
                _ensure_sampler(samplers, f"comp.{qn}", sample_max, rng).add(qa[mask_comp])
                _ensure_sampler(samplers, f"bg.{qn}", sample_max, rng).add(qa[mask_bg])

            # force budget in root region
            sin_t = np.maximum(np.sin(tt), 1.0e-8)
            dp_dr = _gradient_axis(pth, r_cent, axis=0)

            d_r2vr_dr = _gradient_axis((rr * rr) * v_r, r_cent, axis=0)
            d_sinvth_dth = _gradient_axis(sin_t * v_t, th_cent, axis=1)
            d_vp_dph = _gradient_axis(v_p, ph_cent, axis=2)
            divv = d_r2vr_dr / np.maximum(rr * rr, EPS) + d_sinvth_dth / np.maximum(rr * sin_t, EPS) + d_vp_dph / np.maximum(rr * sin_t, EPS)

            d_sinbp_dth = _gradient_axis(sin_t * b_p, th_cent, axis=1)
            d_bt_dph = _gradient_axis(b_t, ph_cent, axis=2)
            curl_r = (d_sinbp_dth - d_bt_dph) / np.maximum(rr * sin_t, EPS)

            d_br_dph = _gradient_axis(b_r, ph_cent, axis=2)
            d_rbp_dr = _gradient_axis(rr * b_p, r_cent, axis=0)
            curl_t = (d_br_dph / np.maximum(sin_t, EPS) - d_rbp_dr) / np.maximum(rr, EPS)

            d_rbt_dr = _gradient_axis(rr * b_t, r_cent, axis=0)
            d_br_dth = _gradient_axis(b_r, th_cent, axis=1)
            curl_p = (d_rbt_dr - d_br_dth) / np.maximum(rr, EPS)

            _ = curl_r  # keep for consistency
            jxb_r = curl_t * b_p - curl_p * b_t
            g_r = gravity_g0 * (sradius**2) / np.maximum(rr * rr, EPS)
            rho_g_r = rho * g_r
            f_pres_r = -dp_dr
            f_res_r = f_pres_r + jxb_r + rho_g_r

            mask_root = valid & (rr >= root_rmin) & (rr <= root_rmax)
            mask_root_comp = mask_root & in_comp
            mask_root_bg = mask_root & (~in_comp)

            force_arrays = {
                "f_pres_r": f_pres_r,
                "jxb_r": jxb_r,
                "rho_g_r": rho_g_r,
                "f_res_r": f_res_r,
                "divv": divv,
            }
            for qn, qa in force_arrays.items():
                _ensure_sampler(samplers, f"root_comp.{qn}", sample_max, rng).add(qa[mask_root_comp])
                _ensure_sampler(samplers, f"root_bg.{qn}", sample_max, rng).add(qa[mask_root_bg])

            # boundary health stats
            if np.any(valid):
                rho_min_global = min(rho_min_global, float(np.min(rho[valid])))
                p_min_global = min(p_min_global, float(np.min(pth[valid])))
                n_valid_global += int(np.count_nonzero(valid))
                if math.isfinite(rho_floor_thr):
                    n_floor_rho_global += int(np.count_nonzero(valid & (rho <= rho_floor_thr)))
                if math.isfinite(p_floor_thr):
                    n_floor_p_global += int(np.count_nonzero(valid & (pth <= p_floor_thr)))

            for seg_name, (r0, r1) in seg_bounds.items():
                seg_mask = valid & (rr >= r0) & (rr < r1)
                seg_counts[seg_name] += int(np.count_nonzero(seg_mask))
                if math.isfinite(rho_floor_thr):
                    seg_floor_rho[seg_name] += int(np.count_nonzero(seg_mask & (rho <= rho_floor_thr)))
                if math.isfinite(p_floor_thr):
                    seg_floor_p[seg_name] += int(np.count_nonzero(seg_mask & (pth <= p_floor_thr)))
                _ensure_sampler(samplers, f"{seg_name}.rho", sample_max, rng).add(rho[seg_mask])
                _ensure_sampler(samplers, f"{seg_name}.pth", sample_max, rng).add(pth[seg_mask])
                _ensure_sampler(samplers, f"{seg_name}.vr", sample_max, rng).add(v_r[seg_mask])
                _ensure_sampler(samplers, f"{seg_name}.b1", sample_max, rng).add(b_r[seg_mask])

            sh0 = valid & (np.broadcast_to((base_r_idx == 0)[:, None, None], valid.shape))
            sh1 = valid & (np.broadcast_to((base_r_idx == 1)[:, None, None], valid.shape))
            _ensure_sampler(samplers, "shell0.rho", sample_max, rng).add(rho[sh0])
            _ensure_sampler(samplers, "shell1.rho", sample_max, rng).add(rho[sh1])
            _ensure_sampler(samplers, "shell0.pth", sample_max, rng).add(pth[sh0])
            _ensure_sampler(samplers, "shell1.pth", sample_max, rng).add(pth[sh1])
            _ensure_sampler(samplers, "shell0.vr", sample_max, rng).add(v_r[sh0])
            _ensure_sampler(samplers, "shell1.vr", sample_max, rng).add(v_r[sh1])
            _ensure_sampler(samplers, "shell0.b1", sample_max, rng).add(b_r[sh0])
            _ensure_sampler(samplers, "shell1.b1", sample_max, rng).add(b_r[sh1])

        # collect thermo row
        thermo_row: Dict[str, float | str | int] = {
            "frame": dat_path.stem,
            "it": int(header["it"]),
            "time": float(header["time"]),
            "comp_size": comp_size,
            "comp_rmin": comp_rmin,
            "comp_r95": comp_r95,
        }
        for grp in ["comp", "bg"]:
            for qn in thermo_q:
                st = _qstats(samplers.get(f"{grp}.{qn}", PrioritySampler(1, rng)).values() if f"{grp}.{qn}" in samplers else np.empty(0))
                _add_group_stats(thermo_row, f"{grp}_{qn}", st)

        for qn in ["e", "pth", "rho", "beta", "cs", "va", "babs", "vr", "vabs"]:
            c = float(thermo_row.get(f"comp_{qn}_p50", math.nan))
            b = float(thermo_row.get(f"bg_{qn}_p50", math.nan))
            thermo_row[f"ratio_comp_bg_{qn}_p50"] = _safe_ratio(c, b)
        c95 = float(thermo_row.get("comp_vr_p95", math.nan))
        b95 = float(thermo_row.get("bg_vr_p95", math.nan))
        thermo_row["ratio_comp_bg_vr_p95"] = _safe_ratio(c95, b95)

        # collect force row
        force_row: Dict[str, float | str | int] = {
            "frame": dat_path.stem,
            "it": int(header["it"]),
            "time": float(header["time"]),
            "root_rmin": root_rmin,
            "root_rmax": root_rmax,
        }
        for grp in ["root_comp", "root_bg"]:
            for qn in force_q:
                st = _qstats(samplers.get(f"{grp}.{qn}", PrioritySampler(1, rng)).values() if f"{grp}.{qn}" in samplers else np.empty(0))
                _add_group_stats(force_row, f"{grp}_{qn}", st)

        # residual imbalance indicator in root component
        fres = abs(float(force_row.get("root_comp_f_res_r_p50", math.nan)))
        fb = (
            abs(float(force_row.get("root_comp_f_pres_r_p50", math.nan)))
            + abs(float(force_row.get("root_comp_jxb_r_p50", math.nan)))
            + abs(float(force_row.get("root_comp_rho_g_r_p50", math.nan)))
        )
        force_row["root_comp_residual_imbalance"] = _safe_ratio(fres, fb)

        # collect boundary row
        bound_row: Dict[str, float | str | int] = {
            "frame": dat_path.stem,
            "it": int(header["it"]),
            "time": float(header["time"]),
            "rho_min_global": rho_min_global if math.isfinite(rho_min_global) else math.nan,
            "pth_min_global": p_min_global if math.isfinite(p_min_global) else math.nan,
            "frac_floor_rho_global": _safe_ratio(float(n_floor_rho_global), float(n_valid_global)),
            "frac_floor_p_global": _safe_ratio(float(n_floor_p_global), float(n_valid_global)),
        }
        for seg_name in seg_bounds:
            bound_row[f"{seg_name}_n_valid"] = float(seg_counts[seg_name])
            bound_row[f"{seg_name}_frac_floor_rho"] = _safe_ratio(float(seg_floor_rho[seg_name]), float(seg_counts[seg_name]))
            bound_row[f"{seg_name}_frac_floor_p"] = _safe_ratio(float(seg_floor_p[seg_name]), float(seg_counts[seg_name]))
            for qn in ["rho", "pth", "vr", "b1"]:
                st = _qstats(samplers.get(f"{seg_name}.{qn}", PrioritySampler(1, rng)).values() if f"{seg_name}.{qn}" in samplers else np.empty(0))
                _add_group_stats(bound_row, f"{seg_name}_{qn}", st)

        for qn in ["rho", "pth", "vr", "b1"]:
            s0 = _qstats(samplers.get(f"shell0.{qn}", PrioritySampler(1, rng)).values() if f"shell0.{qn}" in samplers else np.empty(0))
            s1 = _qstats(samplers.get(f"shell1.{qn}", PrioritySampler(1, rng)).values() if f"shell1.{qn}" in samplers else np.empty(0))
            bound_row[f"shell0_{qn}_p50"] = s0["p50"]
            bound_row[f"shell1_{qn}_p50"] = s1["p50"]
            bound_row[f"jump_{qn}_01"] = _safe_ratio(s1["p50"] - s0["p50"], abs(s1["p50"]) + EPS)

        return metrics, shell_frac, r_centers, thermo_row, force_row, bound_row


def write_csv(path: Path, rows: List[Dict[str, float | str | int]]) -> None:
    if not rows:
        return
    keys: List[str] = []
    for row in rows:
        for k in row.keys():
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def build_anomaly_ranking(
    per_frame: List[FrameMetrics],
    thermo_rows: List[Dict[str, float | str | int]],
    force_rows: List[Dict[str, float | str | int]],
    bound_rows: List[Dict[str, float | str | int]],
    persist_frames: int,
) -> List[Dict[str, float | str | int]]:
    n = len(per_frame)
    if n == 0:
        return []

    def arr(rows: List[Dict[str, float | str | int]], key: str) -> np.ndarray:
        out = np.full(n, np.nan, dtype=np.float64)
        for i, r in enumerate(rows):
            v = r.get(key, math.nan)
            try:
                out[i] = float(v)
            except (TypeError, ValueError):
                out[i] = math.nan
        return out

    ratio_pth = arr(thermo_rows, "ratio_comp_bg_pth_p50")
    ratio_beta = arr(thermo_rows, "ratio_comp_bg_beta_p50")
    comp_beta = arr(thermo_rows, "comp_beta_p50")
    ratio_va = arr(thermo_rows, "ratio_comp_bg_va_p50")
    ratio_cs = arr(thermo_rows, "ratio_comp_bg_cs_p50")
    ratio_vr95 = arr(thermo_rows, "ratio_comp_bg_vr_p95")

    divv_comp = arr(force_rows, "root_comp_divv_p50")
    divv_bg = arr(force_rows, "root_bg_divv_p50")
    res_imb = arr(force_rows, "root_comp_residual_imbalance")

    frac_floor_p_inner = arr(bound_rows, "inner_frac_floor_p")
    frac_floor_rho_inner = arr(bound_rows, "inner_frac_floor_rho")
    jump_p = np.abs(arr(bound_rows, "jump_pth_01"))
    jump_rho = np.abs(arr(bound_rows, "jump_rho_01"))

    indicators: List[Tuple[str, np.ndarray, str]] = [
        ("pth_drop_vs_bg", np.isfinite(ratio_pth) & (ratio_pth < 0.70), "ratio_comp_bg_pth_p50 < 0.70"),
        (
            "beta_collapse",
            np.isfinite(comp_beta) & np.isfinite(ratio_beta) & (comp_beta < 1.0e-3) & (ratio_beta < 0.70),
            "comp_beta_p50 < 1e-3 and ratio_comp_bg_beta_p50 < 0.70",
        ),
        ("va_rise_vs_bg", np.isfinite(ratio_va) & (ratio_va > 1.30), "ratio_comp_bg_va_p50 > 1.30"),
        ("cs_drop_vs_bg", np.isfinite(ratio_cs) & (ratio_cs < 0.80), "ratio_comp_bg_cs_p50 < 0.80"),
        (
            "divv_positive_root",
            np.isfinite(divv_comp) & (divv_comp > 0.0) & np.isfinite(divv_bg) & (divv_comp > 1.5 * divv_bg),
            "root_comp_divv_p50 > 0 and > 1.5*root_bg_divv_p50",
        ),
        (
            "force_residual_imbalance",
            np.isfinite(res_imb) & (res_imb > 0.15),
            "root_comp_residual_imbalance > 0.15",
        ),
        (
            "inner_floor_p_rise",
            np.isfinite(frac_floor_p_inner) & (frac_floor_p_inner > 1.0e-4),
            "inner_frac_floor_p > 1e-4",
        ),
        (
            "inner_floor_rho_rise",
            np.isfinite(frac_floor_rho_inner) & (frac_floor_rho_inner > 1.0e-4),
            "inner_frac_floor_rho > 1e-4",
        ),
        ("jump_pth_01_large", np.isfinite(jump_p) & (jump_p > 0.20), "abs(jump_pth_01) > 0.20"),
        ("jump_rho_01_large", np.isfinite(jump_rho) & (jump_rho > 0.20), "abs(jump_rho_01) > 0.20"),
        ("vr95_boost_vs_bg", np.isfinite(ratio_vr95) & (ratio_vr95 > 1.50), "ratio_comp_bg_vr_p95 > 1.50"),
    ]

    rows: List[Dict[str, float | str | int]] = []
    for name, cond, expr in indicators:
        idx = _first_persistent(cond, persist_frames)
        if idx >= 0:
            rows.append(
                {
                    "indicator": name,
                    "first_idx": idx,
                    "frame": per_frame[idx].frame,
                    "time": per_frame[idx].time,
                    "persist_frames": persist_frames,
                    "condition": expr,
                }
            )
        else:
            rows.append(
                {
                    "indicator": name,
                    "first_idx": -1,
                    "frame": "NA",
                    "time": math.nan,
                    "persist_frames": persist_frames,
                    "condition": expr,
                }
            )

    rows.sort(key=lambda r: (int(r["first_idx"]) if int(r["first_idx"]) >= 0 else 10**9, str(r["indicator"])))
    for rank, r in enumerate(rows, start=1):
        r["rank"] = rank
    return rows


def maybe_make_plots(
    out_prefix: Path,
    per_frame: List[FrameMetrics],
    thermo_rows: List[Dict[str, float | str | int]],
    force_rows: List[Dict[str, float | str | int]],
    bound_rows: List[Dict[str, float | str | int]],
) -> List[Path]:
    paths: List[Path] = []
    if not per_frame:
        return paths

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return paths

    t = np.array([m.time for m in per_frame], dtype=np.float64)

    def s(rows: List[Dict[str, float | str | int]], key: str) -> np.ndarray:
        out = np.full(len(rows), np.nan, dtype=np.float64)
        for i, r in enumerate(rows):
            try:
                out[i] = float(r.get(key, math.nan))
            except (TypeError, ValueError):
                out[i] = math.nan
        return out

    # thermo chain
    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(t, s(thermo_rows, "comp_pth_p50"), label="comp pth P50")
    axes[0].plot(t, s(thermo_rows, "bg_pth_p50"), label="bg pth P50")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("pth")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].plot(t, s(thermo_rows, "comp_beta_p50"), label="comp beta P50")
    axes[1].plot(t, s(thermo_rows, "bg_beta_p50"), label="bg beta P50")
    axes[1].axhline(1.0e-3, color="k", linestyle="--", linewidth=1, label="beta=1e-3")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("beta")
    axes[1].legend(loc="best", fontsize=8)

    axes[2].plot(t, s(thermo_rows, "comp_va_p50"), label="comp vA P50")
    axes[2].plot(t, s(thermo_rows, "bg_va_p50"), label="bg vA P50")
    axes[2].plot(t, s(thermo_rows, "comp_cs_p50"), label="comp cs P50")
    axes[2].plot(t, s(thermo_rows, "bg_cs_p50"), label="bg cs P50")
    axes[2].set_ylabel("speed")
    axes[2].legend(loc="best", fontsize=8)

    axes[3].plot(t, s(thermo_rows, "comp_vr_p95"), label="comp vr P95")
    axes[3].plot(t, s(thermo_rows, "bg_vr_p95"), label="bg vr P95")
    axes[3].set_ylabel("vr")
    axes[3].set_xlabel("time")
    axes[3].legend(loc="best", fontsize=8)

    fig.suptitle("Thermodynamics and Low-beta Chain")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    p1 = Path(str(out_prefix) + "_thermo_chain.png")
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    paths.append(p1)

    # force budget
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(t, s(force_rows, "root_comp_f_pres_r_p50"), label="-(dp/dr) P50")
    axes[0].plot(t, s(force_rows, "root_comp_jxb_r_p50"), label="(JxB)_r P50")
    axes[0].plot(t, s(force_rows, "root_comp_rho_g_r_p50"), label="rho g_r P50")
    axes[0].plot(t, s(force_rows, "root_comp_f_res_r_p50"), label="F_res,r P50", linewidth=2)
    axes[0].set_ylabel("force")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].plot(t, s(force_rows, "root_comp_divv_p50"), label="comp div(v) P50")
    axes[1].plot(t, s(force_rows, "root_bg_divv_p50"), label="bg div(v) P50")
    axes[1].axhline(0.0, color="k", linestyle="--", linewidth=1)
    axes[1].set_ylabel("div(v)")
    axes[1].set_xlabel("time")
    axes[1].legend(loc="best", fontsize=8)

    fig.suptitle("Root-Region Force Budget")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    p2 = Path(str(out_prefix) + "_force_budget.png")
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    paths.append(p2)

    # boundary health
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    axes[0].plot(t, s(bound_rows, "inner_frac_floor_rho"), label="inner frac floor rho")
    axes[0].plot(t, s(bound_rows, "inner_frac_floor_p"), label="inner frac floor p")
    axes[0].set_ylabel("floor frac")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].plot(t, np.abs(s(bound_rows, "jump_rho_01")), label="|jump rho 0->1|")
    axes[1].plot(t, np.abs(s(bound_rows, "jump_pth_01")), label="|jump pth 0->1|")
    axes[1].plot(t, np.abs(s(bound_rows, "jump_vr_01")), label="|jump vr 0->1|")
    axes[1].set_ylabel("jump")
    axes[1].legend(loc="best", fontsize=8)

    axes[2].plot(t, s(bound_rows, "inner_rho_p50"), label="inner rho P50")
    axes[2].plot(t, s(bound_rows, "mid_rho_p50"), label="mid rho P50")
    axes[2].plot(t, s(bound_rows, "outer_rho_p50"), label="outer rho P50")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("rho")
    axes[2].set_xlabel("time")
    axes[2].legend(loc="best", fontsize=8)

    fig.suptitle("Boundary Coupling and Numerical Health")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    p3 = Path(str(out_prefix) + "_boundary_health.png")
    fig.savefig(p3, dpi=150)
    plt.close(fig)
    paths.append(p3)

    return paths


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    par_cfg = read_par_config(case_dir / args.par)
    base_edges = build_base_radial_edges(
        par_cfg["xprobmin1"],
        par_cfg["xprobmax1"],
        int(par_cfg["domain_nx1"]),
        par_cfg["qstretch_baselevel"],
    )

    paths = select_paths(case_dir.glob(args.pattern), args.start, args.end, args.stride)
    if not paths:
        raise FileNotFoundError("No dat snapshots selected.")

    per_frame: List[FrameMetrics] = []
    shell_frac_by_frame: Dict[str, np.ndarray] = {}
    thermo_rows: List[Dict[str, float | str | int]] = []
    force_rows: List[Dict[str, float | str | int]] = []
    bound_rows: List[Dict[str, float | str | int]] = []
    r_centers_ref: np.ndarray | None = None

    rng = np.random.default_rng(args.seed)

    for i, dat_path in enumerate(paths, start=1):
        m, shell_frac, r_centers, thermo_row, force_row, bound_row = analyze_frame(
            dat_path=dat_path,
            par_cfg=par_cfg,
            base_edges=base_edges,
            vfix=args.vfix,
            inner_anchor_r=args.inner_anchor_r,
            min_comp_cells=args.min_comp_cells,
            root_rmin=args.root_rmin,
            root_rmax=args.root_rmax,
            gravity_g0=args.gravity_g0,
            sradius=args.sradius,
            sample_max=args.sample_max,
            floor_factor=args.floor_factor,
            rng=rng,
        )
        per_frame.append(m)
        shell_frac_by_frame[m.frame] = shell_frac
        thermo_rows.append(thermo_row)
        force_rows.append(force_row)
        bound_rows.append(bound_row)
        r_centers_ref = r_centers
        print(
            f"[{i}/{len(paths)}] {m.frame} high={m.high_any} comp={m.comp_size} rmin={m.comp_rmin:.3f} r95={m.comp_r95:.3f}",
            flush=True,
        )

    if r_centers_ref is None:
        raise RuntimeError("No r-centers computed.")

    # Original 3 hard criteria
    n_shell = r_centers_ref.size
    t_on = np.full(n_shell, np.nan, dtype=np.float64)
    f_on = np.full(n_shell, -1, dtype=np.int64)
    for fi, m in enumerate(per_frame):
        frac = shell_frac_by_frame[m.frame]
        hit = (frac > args.frac_thresh) & np.isnan(t_on)
        t_on[hit] = m.time
        f_on[hit] = fi

    valid_on = np.isfinite(t_on)
    near_inner = valid_on & (r_centers_ref < args.near_inner_r)
    outer_zone = valid_on & (r_centers_ref >= args.outer_check_r)

    earliest_inner = int(np.nanmin(f_on[near_inner])) if np.any(near_inner) else -1
    earliest_outer = int(np.nanmin(f_on[outer_zone])) if np.any(outer_zone) else -1
    outer_leads = int(earliest_outer >= 0 and earliest_inner >= 0 and earliest_outer + 2 <= earliest_inner)

    r_on = r_centers_ref[valid_on]
    t_on_valid = t_on[valid_on]
    c1_spearman = spearman(r_on, t_on_valid)
    c1_near_inner_first = int(np.any(near_inner) and earliest_inner == int(np.nanmin(f_on[valid_on])))
    c1_pass = int((not outer_leads) and c1_near_inner_first == 1 and (math.isfinite(c1_spearman) and c1_spearman >= 0.8))

    onset_frame_idx = next((i for i, m in enumerate(per_frame) if m.high_any == 1), None)
    if onset_frame_idx is None:
        onset_frame_idx = 0
    end_idx = min(len(per_frame), onset_frame_idx + max(1, args.growth_frames))
    win = per_frame[onset_frame_idx:end_idx]
    anchored_ratio = float(np.mean([m.anchored_inner115 for m in win])) if win else math.nan
    c2_pass = int(math.isfinite(anchored_ratio) and anchored_ratio >= 0.8)

    fit_t: List[float] = []
    fit_r95: List[float] = []
    t0 = per_frame[onset_frame_idx].time if per_frame else 0.0
    for m in win:
        if m.comp_size >= args.min_comp_cells and math.isfinite(m.comp_r95):
            fit_t.append(m.time - t0)
            fit_r95.append(m.comp_r95)

    if len(fit_t) >= 4 and (max(fit_t) - min(fit_t)) > 0:
        coef = np.polyfit(np.array(fit_t), np.array(fit_r95), 1)
        uf = float(coef[0])
        r0_fit = float(coef[1])
        c3_pass = int((uf > 0.0) and (1.0 <= r0_fit <= args.near_inner_r))
    else:
        uf = math.nan
        r0_fit = math.nan
        c3_pass = 0

    pass_count = c1_pass + c2_pass + c3_pass
    if pass_count >= 2:
        final = "boundary_trigger_dominant"
    elif pass_count <= 1 and ((1 - c1_pass) + (1 - c2_pass) + (1 - c3_pass) >= 2):
        final = "volume_trigger_possible_or_mixed"
    else:
        final = "inconclusive"

    out_prefix = case_dir / args.output_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # original outputs
    frame_path = Path(str(out_prefix) + "_framewise.dat")
    with frame_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "frame",
            "it",
            "time",
            "high_any",
            "high_frac_global",
            "comp_size",
            "comp_rmin",
            "comp_r95",
            "anchored_inner115",
        ])
        for m in per_frame:
            writer.writerow(
                [
                    m.frame,
                    m.it,
                    f"{m.time:.10e}",
                    m.high_any,
                    f"{m.high_frac_global:.6e}",
                    m.comp_size,
                    f"{m.comp_rmin:.6f}",
                    f"{m.comp_r95:.6f}",
                    m.anchored_inner115,
                ]
            )

    summary_path = Path(str(out_prefix) + "_summary.dat")
    with summary_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["n_frames", len(per_frame)])
        writer.writerow(["onset_frame_index", onset_frame_idx])
        writer.writerow(["criterion1_pass", c1_pass])
        writer.writerow(["criterion1_spearman_r_vs_ton", f"{c1_spearman:.6f}"])
        writer.writerow(["criterion1_near_inner_first", c1_near_inner_first])
        writer.writerow(["criterion1_outer_leads_by_2frames", outer_leads])
        writer.writerow(["criterion2_pass", c2_pass])
        writer.writerow(["criterion2_anchor_ratio_growth_window", f"{anchored_ratio:.6f}"])
        writer.writerow(["criterion3_pass", c3_pass])
        writer.writerow(["criterion3_front_speed", f"{uf:.6e}"])
        writer.writerow(["criterion3_r0_fit", f"{r0_fit:.6f}"])
        writer.writerow(["pass_count", pass_count])
        writer.writerow(["final_judgement", final])

    # new outputs
    thermo_path = Path(str(out_prefix) + "_thermo_timeseries.dat")
    force_path = Path(str(out_prefix) + "_force_root_timeseries.dat")
    bound_path = Path(str(out_prefix) + "_boundary_health_timeseries.dat")
    ranking_path = Path(str(out_prefix) + "_first_anomaly_ranking.dat")

    write_csv(thermo_path, thermo_rows)
    write_csv(force_path, force_rows)
    write_csv(bound_path, bound_rows)

    ranking_rows = build_anomaly_ranking(
        per_frame=per_frame,
        thermo_rows=thermo_rows,
        force_rows=force_rows,
        bound_rows=bound_rows,
        persist_frames=args.persist_frames,
    )
    write_csv(ranking_path, ranking_rows)

    plot_paths: List[Path] = []
    if args.make_plots:
        plot_paths = maybe_make_plots(out_prefix, per_frame, thermo_rows, force_rows, bound_rows)

    print(f"Wrote {frame_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {thermo_path}")
    print(f"Wrote {force_path}")
    print(f"Wrote {bound_path}")
    print(f"Wrote {ranking_path}")
    for pp in plot_paths:
        print(f"Wrote {pp}")

    first_hit = next((r for r in ranking_rows if int(r.get("first_idx", -1)) >= 0), None)
    if first_hit is not None:
        print(
            f"FIRST_ANOMALY: indicator={first_hit['indicator']} frame={first_hit['frame']} "
            f"time={first_hit['time']:.6e} idx={first_hit['first_idx']}"
        )

    print(f"FINAL: {final} (pass={pass_count}/3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
