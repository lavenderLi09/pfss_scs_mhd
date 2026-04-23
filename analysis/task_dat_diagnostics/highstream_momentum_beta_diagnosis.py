# Copied to project analysis task toolbox
# Source case path: off_2270/analysis/highstream_momentum_beta_diagnosis.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270/analysis/highstream_momentum_beta_diagnosis.py

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

CASE_DIR = Path(__file__).resolve().parent.parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))

from hao_code.datfile_io import get_header, get_tree_info  # noqa: E402
from highstream_trigger_dat_diagnosis import (  # noqa: E402
    EPS,
    PrioritySampler,
    _add_group_stats,
    _ensure_sampler,
    _gradient_axis,
    _qstats,
    block_radial_edges,
    block_uniform_edges,
    build_base_radial_edges,
    conservative_internal_to_primitive,
    frame_id,
    largest_component,
    read_block_fields,
    read_par_config,
    select_paths,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Diagnose anomalous high-speed component evolution from dat only, "
            "with emphasis on plasma beta, Alfven speed, momentum force budget, "
            "residual and acceleration direction."
        )
    )
    p.add_argument("--case-dir", type=Path, default=CASE_DIR)
    p.add_argument("--par", default="amrvac.par")
    p.add_argument("--pattern", default="data/off*.dat")
    p.add_argument("--start", type=int, default=None)
    p.add_argument("--end", type=int, default=None)
    p.add_argument("--stride", type=int, default=1)

    p.add_argument(
        "--vfix",
        type=float,
        default=6.0,
        help="Fixed speed threshold in code unit (~700 km/s if unit~116.5km/s).",
    )
    p.add_argument(
        "--adaptive-threshold",
        action="store_true",
        default=False,
        help="If enabled, use max(vfix, shell_mean+3sigma); otherwise use fixed vfix only (faster).",
    )
    p.add_argument("--min-comp-cells", type=int, default=30)
    p.add_argument("--root-rmin", type=float, default=1.03)
    p.add_argument("--root-rmax", type=float, default=1.20)
    p.add_argument("--gravity-g0", type=float, default=-14.072)
    p.add_argument("--sradius", type=float, default=1.0)
    p.add_argument("--sample-max", type=int, default=50000)
    p.add_argument("--seed", type=int, default=20260417)

    p.add_argument("--make-plots", action="store_true", default=True)
    p.add_argument("--no-make-plots", dest="make_plots", action="store_false")
    p.add_argument("--output-prefix", default="analysis/highstream_momentum_beta")
    return p.parse_args()


def _safe_ratio(num: float, den: float) -> float:
    if not math.isfinite(num) or not math.isfinite(den):
        return math.nan
    if abs(den) <= EPS:
        return math.nan
    return num / den


def _init_sign_counter() -> Dict[str, float]:
    return {"pos": 0.0, "neg": 0.0, "nz": 0.0}


def _update_sign_counter(counter: Dict[str, float], arr: np.ndarray) -> None:
    finite = np.isfinite(arr)
    if not np.any(finite):
        return
    x = arr[finite]
    nz = np.abs(x) > 1.0e-20
    if not np.any(nz):
        return
    xx = x[nz]
    counter["nz"] += float(xx.size)
    counter["pos"] += float(np.count_nonzero(xx > 0.0))
    counter["neg"] += float(np.count_nonzero(xx < 0.0))


def _finalize_sign(row: Dict[str, float | str | int], prefix: str, counter: Dict[str, float]) -> None:
    row[f"{prefix}_nz"] = counter["nz"]
    row[f"{prefix}_pos_frac"] = _safe_ratio(counter["pos"], counter["nz"])
    row[f"{prefix}_neg_frac"] = _safe_ratio(counter["neg"], counter["nz"])


def _read_header_tree_compat(dat_path: Path) -> Tuple[Dict, np.ndarray, np.ndarray, np.ndarray]:
    try:
        hdr = get_header(str(dat_path))
        tree = get_tree_info(str(dat_path))
    except Exception:
        with dat_path.open("rb") as fh:
            hdr = get_header(fh)
            tree = get_tree_info(fh)

    if isinstance(tree, dict):
        offset_blocks = np.asarray(tree["offset_blocks"])
        node_info = np.asarray(tree["node_info"])
        levels = node_info[:, 0]
        indices = node_info[:, 1:4]
    elif isinstance(tree, tuple) and len(tree) >= 3:
        levels = np.asarray(tree[0])
        indices = np.asarray(tree[1])
        offset_blocks = np.asarray(tree[2])
    else:
        raise RuntimeError(f"Unsupported tree structure from datfile_io for {dat_path}")
    return hdr, levels, indices, offset_blocks


def analyze_frame(
    dat_path: Path,
    par_cfg: Dict[str, float],
    base_edges: np.ndarray,
    vfix: float,
    adaptive_threshold: bool,
    min_comp_cells: int,
    root_rmin: float,
    root_rmax: float,
    gravity_g0: float,
    sradius: float,
    sample_max: int,
    rng: np.random.Generator,
) -> Dict[str, float | str | int]:
    hdr, levels_raw, indices_raw, offset_blocks_raw = _read_header_tree_compat(dat_path)

    offset_blocks = np.asarray(offset_blocks_raw).astype(np.int64)
    levels = np.asarray(levels_raw).astype(np.int64)
    indices = np.asarray(indices_raw).astype(np.int64)
    block_shape = hdr["block_nx"].astype(np.int64)
    ndim = int(hdr["ndim"])
    nw = int(hdr["nw"])

    if ndim != 3:
        raise RuntimeError(f"{dat_path}: expected ndim=3, got ndim={ndim}")

    nr = int(par_cfg["domain_nx1"])
    nt = int(par_cfg["domain_nx2"])
    np_ = int(par_cfg["domain_nx3"])
    gamma = float(par_cfg["mhd_gamma"])

    shell_thr = np.full(nr, float(vfix), dtype=np.float64)
    if adaptive_threshold:
        shell_count = np.zeros(nr, dtype=np.float64)
        shell_sum = np.zeros(nr, dtype=np.float64)
        shell_sumsq = np.zeros(nr, dtype=np.float64)

        with dat_path.open("rb") as f:
            for ib in range(offset_blocks.shape[0]):
                lvl = int(levels[ib])
                ix1, ix2, ix3 = [int(v) for v in indices[ib]]
                r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
                r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
                base_r_idx = np.clip(np.floor((r_cent - base_edges[0]) / (base_edges[1] - base_edges[0])).astype(int), 0, nr - 1)

                fields = read_block_fields(f, int(offset_blocks[ib]), block_shape, ndim, nw)
                rho = fields[0]
                v_r, v_t, v_p, _pth = conservative_internal_to_primitive(
                    rho=rho, m1=fields[1], m2=fields[2], m3=fields[3], eint=fields[4], gamma=gamma
                )
                vabs = np.sqrt(v_r * v_r + v_t * v_t + v_p * v_p)

                for i_local, ridx in enumerate(base_r_idx):
                    vals = vabs[i_local].ravel()
                    vals = vals[np.isfinite(vals)]
                    if vals.size == 0:
                        continue
                    shell_count[ridx] += float(vals.size)
                    shell_sum[ridx] += float(vals.sum())
                    shell_sumsq[ridx] += float((vals * vals).sum())

        shell_mean = np.divide(shell_sum, np.maximum(shell_count, 1.0))
        shell_var = np.divide(shell_sumsq, np.maximum(shell_count, 1.0)) - shell_mean * shell_mean
        shell_std = np.sqrt(np.maximum(shell_var, 0.0))
        shell_thr = np.maximum(vfix, shell_mean + 3.0 * shell_std)

    base_mask = np.zeros((nr, nt, np_), dtype=np.bool_)
    r_centers = np.zeros(nr, dtype=np.float64)

    with dat_path.open("rb") as f:
        for ib in range(offset_blocks.shape[0]):
            lvl = int(levels[ib])
            ix1, ix2, ix3 = [int(v) for v in indices[ib]]
            r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
            t_edges = block_uniform_edges(
                par_cfg["xprobmin2"], par_cfg["xprobmax2"], nt, lvl, ix2, int(block_shape[1])
            )
            p_edges = block_uniform_edges(
                par_cfg["xprobmin3"], par_cfg["xprobmax3"], np_, lvl, ix3, int(block_shape[2])
            )
            r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
            t_cent = 0.5 * (t_edges[:-1] + t_edges[1:])
            p_cent = 0.5 * (p_edges[:-1] + p_edges[1:])
            base_r_idx = np.clip(np.floor((r_cent - base_edges[0]) / (base_edges[1] - base_edges[0])).astype(int), 0, nr - 1)
            base_t_idx = np.clip(
                np.floor((t_cent - par_cfg["xprobmin2"]) / (par_cfg["xprobmax2"] - par_cfg["xprobmin2"]) * nt).astype(int),
                0,
                nt - 1,
            )
            base_p_idx = np.mod(
                np.floor((p_cent - par_cfg["xprobmin3"]) / (par_cfg["xprobmax3"] - par_cfg["xprobmin3"]) * np_).astype(int),
                np_,
            )

            fields = read_block_fields(f, int(offset_blocks[ib]), block_shape, ndim, nw)
            rho = fields[0]
            v_r, v_t, v_p, _pth = conservative_internal_to_primitive(
                rho=rho, m1=fields[1], m2=fields[2], m3=fields[3], eint=fields[4], gamma=gamma
            )
            vabs = np.sqrt(v_r * v_r + v_t * v_t + v_p * v_p)
            thr_3d = shell_thr[base_r_idx][:, None, None]
            high = np.isfinite(vabs) & (vabs > thr_3d)
            base_mask[np.ix_(base_r_idx, base_t_idx, base_p_idx)] |= high
            r_centers[base_r_idx] = r_cent

    comp = largest_component(base_mask)
    comp_mask = np.zeros((nr, nt, np_), dtype=np.bool_)
    if comp.size > 0:
        comp_mask[comp[:, 0], comp[:, 1], comp[:, 2]] = True

    if comp.shape[0] >= min_comp_cells:
        rvals = r_centers[comp[:, 0]]
        comp_rmin = float(np.min(rvals))
        comp_r95 = float(np.quantile(rvals, 0.95))
        comp_size = int(comp.shape[0])
    else:
        comp_rmin = math.nan
        comp_r95 = math.nan
        comp_size = int(comp.shape[0])

    samplers: Dict[str, PrioritySampler] = {}
    sign_keys = [
        "vr",
        "f_pres_r",
        "jxb_r",
        "rho_g_r",
        "f_res_r",
        "a_pres_r",
        "a_jxb_r",
        "a_g_r",
        "a_res_r",
    ]
    sign_comp = {k: _init_sign_counter() for k in sign_keys}
    sign_bg = {k: _init_sign_counter() for k in sign_keys}

    root_cells = 0
    with dat_path.open("rb") as f:
        for ib in range(offset_blocks.shape[0]):
            lvl = int(levels[ib])
            ix1, ix2, ix3 = [int(v) for v in indices[ib]]
            r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
            t_edges = block_uniform_edges(
                par_cfg["xprobmin2"], par_cfg["xprobmax2"], nt, lvl, ix2, int(block_shape[1])
            )
            p_edges = block_uniform_edges(
                par_cfg["xprobmin3"], par_cfg["xprobmax3"], np_, lvl, ix3, int(block_shape[2])
            )
            r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
            t_cent = 0.5 * (t_edges[:-1] + t_edges[1:])
            p_cent = 0.5 * (p_edges[:-1] + p_edges[1:])
            dr = np.gradient(r_cent)
            dt = np.gradient(t_cent)
            dp = np.gradient(p_cent)

            base_r_idx = np.clip(np.floor((r_cent - base_edges[0]) / (base_edges[1] - base_edges[0])).astype(int), 0, nr - 1)
            base_t_idx = np.clip(
                np.floor((t_cent - par_cfg["xprobmin2"]) / (par_cfg["xprobmax2"] - par_cfg["xprobmin2"]) * nt).astype(int),
                0,
                nt - 1,
            )
            base_p_idx = np.mod(
                np.floor((p_cent - par_cfg["xprobmin3"]) / (par_cfg["xprobmax3"] - par_cfg["xprobmin3"]) * np_).astype(int),
                np_,
            )
            in_comp = comp_mask[np.ix_(base_r_idx, base_t_idx, base_p_idx)]

            fields = read_block_fields(f, int(offset_blocks[ib]), block_shape, ndim, nw)
            rho = fields[0]
            m1 = fields[1]
            m2 = fields[2]
            m3 = fields[3]
            eint = fields[4]
            b_r = fields[5]
            b_t = fields[6]
            b_p = fields[7]

            v_r, v_t, v_p, pth = conservative_internal_to_primitive(
                rho=rho, m1=m1, m2=m2, m3=m3, eint=eint, gamma=gamma
            )
            vabs = np.sqrt(v_r * v_r + v_t * v_t + v_p * v_p)
            b2 = b_r * b_r + b_t * b_t + b_p * b_p
            babs = np.sqrt(np.maximum(b2, 0.0))
            beta = 2.0 * pth / np.maximum(b2, EPS)
            va = babs / np.sqrt(np.maximum(rho, EPS))
            cs = np.sqrt(np.maximum(gamma * pth / np.maximum(rho, EPS), 0.0))
            align_vr = v_r / np.maximum(vabs, EPS)

            rr = r_cent[:, None, None]
            th = t_cent[None, :, None]
            sin_t = np.maximum(np.sin(th), 1.0e-8)
            r2 = rr * rr

            dpth_dr = _gradient_axis(pth, r_cent, axis=0)
            f_pres_r = -dpth_dr

            rbt = rr * b_t
            rbp = rr * b_p
            d_rbp_dr = _gradient_axis(rbp, r_cent, axis=0)
            d_rbt_dr = _gradient_axis(rbt, r_cent, axis=0)
            d_br_dth = _gradient_axis(b_r, t_cent, axis=1)
            d_br_dph = _gradient_axis(b_r, p_cent, axis=2)
            d_bt_dph = _gradient_axis(b_t, p_cent, axis=2)
            d_sinbp_dth = _gradient_axis(sin_t * b_p, t_cent, axis=1)

            curl_r = (d_sinbp_dth - d_bt_dph) / np.maximum(rr * sin_t, EPS)
            curl_t = (d_br_dph / np.maximum(sin_t, EPS) - d_rbp_dr) / np.maximum(rr, EPS)
            curl_p = (d_rbt_dr - d_br_dth) / np.maximum(rr, EPS)

            jxb_r = curl_t * b_p - curl_p * b_t
            g_r = gravity_g0 * (sradius * sradius) / np.maximum(rr * rr, EPS)
            rho_g_r = rho * g_r
            f_res_r = f_pres_r + jxb_r + rho_g_r

            rho_safe = np.maximum(rho, EPS)
            a_pres_r = f_pres_r / rho_safe
            a_jxb_r = jxb_r / rho_safe
            a_g_r = rho_g_r / rho_safe
            a_res_r = f_res_r / rho_safe

            d_r2vr_dr = _gradient_axis(rr * rr * v_r, r_cent, axis=0)
            d_sinvth_dth = _gradient_axis(sin_t * v_t, t_cent, axis=1)
            d_vp_dph = _gradient_axis(v_p, p_cent, axis=2)
            divv = d_r2vr_dr / np.maximum(rr * rr, EPS) + d_sinvth_dth / np.maximum(rr * sin_t, EPS) + d_vp_dph / np.maximum(
                rr * sin_t, EPS
            )

            if math.isfinite(comp_rmin) and math.isfinite(comp_r95):
                env = (rr >= comp_rmin) & (rr <= comp_r95)
            else:
                env = np.ones_like(in_comp, dtype=np.bool_)

            valid = np.isfinite(rho) & np.isfinite(pth) & np.isfinite(vabs) & np.isfinite(babs)
            mask_comp = valid & in_comp
            mask_bg = valid & (~in_comp) & env
            mask_root = valid & (rr >= root_rmin) & (rr <= root_rmax)
            mask_root_comp = mask_root & in_comp

            arrays = {
                "beta": beta,
                "va": va,
                "cs": cs,
                "rho": rho,
                "pth": pth,
                "e": eint,
                "babs": babs,
                "vr": v_r,
                "vabs": vabs,
                "align_vr": align_vr,
                "f_pres_r": f_pres_r,
                "jxb_r": jxb_r,
                "rho_g_r": rho_g_r,
                "f_res_r": f_res_r,
                "a_pres_r": a_pres_r,
                "a_jxb_r": a_jxb_r,
                "a_g_r": a_g_r,
                "a_res_r": a_res_r,
                "divv": divv,
            }

            for k, arr in arrays.items():
                _ensure_sampler(samplers, f"comp.{k}", sample_max, rng).add(arr[mask_comp])
                _ensure_sampler(samplers, f"bg.{k}", sample_max, rng).add(arr[mask_bg])

            root_cells += int(np.count_nonzero(mask_root_comp))

            for k in sign_keys:
                _update_sign_counter(sign_comp[k], arrays[k][mask_comp])
                _update_sign_counter(sign_bg[k], arrays[k][mask_bg])

    row: Dict[str, float | str | int] = {
        "frame": dat_path.stem,
        "idx": int(frame_id(dat_path.stem) or -1),
        "time": float(hdr["time"]),
        "it": int(hdr["it"]),
        "comp_size": comp_size,
        "comp_rmin": comp_rmin,
        "comp_r95": comp_r95,
        "root_comp_cells": root_cells,
    }

    keys = [
        "beta",
        "va",
        "cs",
        "rho",
        "pth",
        "e",
        "babs",
        "vr",
        "vabs",
        "align_vr",
        "f_pres_r",
        "jxb_r",
        "rho_g_r",
        "f_res_r",
        "a_pres_r",
        "a_jxb_r",
        "a_g_r",
        "a_res_r",
        "divv",
    ]
    for grp in ["comp", "bg"]:
        for k in keys:
            vals = samplers.get(f"{grp}.{k}", PrioritySampler(1, rng)).values() if f"{grp}.{k}" in samplers else np.empty(0)
            st = _qstats(vals)
            _add_group_stats(row, f"{grp}_{k}", st)

    for k in ["beta", "va", "vr", "vabs", "pth", "rho", "jxb_r", "f_res_r", "a_res_r", "divv"]:
        c = float(row.get(f"comp_{k}_p50", math.nan))
        b = float(row.get(f"bg_{k}_p50", math.nan))
        row[f"ratio_comp_bg_{k}_p50"] = _safe_ratio(c, b)

    j = abs(float(row.get("comp_jxb_r_p50", math.nan)))
    p = abs(float(row.get("comp_f_pres_r_p50", math.nan)))
    g = abs(float(row.get("comp_rho_g_r_p50", math.nan)))
    row["jxb_dom_ratio_p50"] = _safe_ratio(j, p + g)

    for k in sign_keys:
        _finalize_sign(row, f"comp_{k}", sign_comp[k])
        _finalize_sign(row, f"bg_{k}", sign_bg[k])

    return row


def write_csv(path: Path, rows: List[Dict[str, float | str | int]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows for {path}")
    keys: List[str] = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _arr(rows: Sequence[Dict[str, float | str | int]], key: str) -> np.ndarray:
    out = []
    for r in rows:
        v = r.get(key, math.nan)
        try:
            out.append(float(v))
        except Exception:
            out.append(math.nan)
    return np.array(out, dtype=np.float64)


def maybe_make_plots(rows: List[Dict[str, float | str | int]], out_prefix: Path) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = _arr(rows, "time")
    paths: List[Path] = []

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    axes[0].plot(t, _arr(rows, "comp_beta_p50"), label="comp beta P50")
    axes[0].plot(t, _arr(rows, "bg_beta_p50"), label="bg beta P50")
    axes[0].axhline(1.0e-3, color="k", linestyle="--", linewidth=1, label="beta=1e-3")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("beta")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(t, _arr(rows, "comp_va_p50"), label="comp vA P50")
    axes[1].plot(t, _arr(rows, "bg_va_p50"), label="bg vA P50")
    axes[1].plot(t, _arr(rows, "comp_vabs_p95"), label="comp |v| P95")
    axes[1].set_ylabel("speed")
    axes[1].legend(loc="best", fontsize=8)
    axes[1].grid(alpha=0.3)

    axes[2].plot(t, _arr(rows, "ratio_comp_bg_beta_p50"), label="beta comp/bg")
    axes[2].plot(t, _arr(rows, "ratio_comp_bg_va_p50"), label="vA comp/bg")
    axes[2].plot(t, _arr(rows, "ratio_comp_bg_vr_p50"), label="vr comp/bg")
    axes[2].set_xlabel("time")
    axes[2].set_ylabel("ratio")
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(alpha=0.3)
    fig.suptitle("Anomalous Component: plasma beta & Alfvenic evolution")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p1 = Path(str(out_prefix) + "_beta_va.png")
    fig.savefig(p1, dpi=170)
    plt.close(fig)
    paths.append(p1)

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(t, _arr(rows, "comp_f_pres_r_p50"), label="-(dp/dr)")
    axes[0].plot(t, _arr(rows, "comp_jxb_r_p50"), label="(JxB)_r")
    axes[0].plot(t, _arr(rows, "comp_rho_g_r_p50"), label="rho g_r")
    axes[0].plot(t, _arr(rows, "comp_f_res_r_p50"), label="F_res,r", linewidth=2)
    axes[0].set_ylabel("force")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(t, _arr(rows, "comp_a_pres_r_p50"), label="a_pres,r")
    axes[1].plot(t, _arr(rows, "comp_a_jxb_r_p50"), label="a_jxb,r")
    axes[1].plot(t, _arr(rows, "comp_a_g_r_p50"), label="a_g,r")
    axes[1].plot(t, _arr(rows, "comp_a_res_r_p50"), label="a_res,r", linewidth=2)
    axes[1].set_ylabel("acceleration")
    axes[1].legend(loc="best", fontsize=8)
    axes[1].grid(alpha=0.3)

    axes[2].plot(t, _arr(rows, "comp_vr_pos_frac"), label="frac(vr>0)")
    axes[2].plot(t, _arr(rows, "comp_f_res_r_pos_frac"), label="frac(F_res,r>0)")
    axes[2].plot(t, _arr(rows, "comp_a_res_r_pos_frac"), label="frac(a_res,r>0)")
    axes[2].plot(t, _arr(rows, "comp_jxb_r_pos_frac"), label="frac((JxB)_r>0)")
    axes[2].set_ylim(-0.02, 1.02)
    axes[2].set_xlabel("time")
    axes[2].set_ylabel("outward fraction (+r)")
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(alpha=0.3)
    fig.suptitle("Momentum budget, residual and direction (+r outward)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p2 = Path(str(out_prefix) + "_momentum_direction.png")
    fig.savefig(p2, dpi=170)
    plt.close(fig)
    paths.append(p2)
    return paths


def write_markdown_summary(rows: List[Dict[str, float | str | int]], out_prefix: Path) -> Path:
    if not rows:
        raise RuntimeError("No rows to summarize.")

    arr_v = _arr(rows, "comp_vabs_p95")
    arr_beta = _arr(rows, "comp_beta_p50")
    arr_va = _arr(rows, "comp_va_p50")
    arr_jdom = _arr(rows, "jxb_dom_ratio_p50")
    arr_ares = _arr(rows, "comp_a_res_r_p50")
    arr_ares_out = _arr(rows, "comp_a_res_r_pos_frac")
    arr_div = _arr(rows, "comp_divv_p50")

    def _imax_finite(x: np.ndarray) -> int:
        ok = np.where(np.isfinite(x))[0]
        if ok.size == 0:
            return -1
        return int(ok[np.nanargmax(x[ok])])

    def _imin_finite(x: np.ndarray) -> int:
        ok = np.where(np.isfinite(x))[0]
        if ok.size == 0:
            return -1
        return int(ok[np.nanargmin(x[ok])])

    i_v = _imax_finite(arr_v)
    i_b = _imin_finite(arr_beta)
    i_va = _imax_finite(arr_va)

    top_idx = np.argsort(np.where(np.isfinite(arr_v), arr_v, -np.inf))[-10:][::-1]
    top_idx = [int(i) for i in top_idx if np.isfinite(arr_v[i])]

    md_path = Path(str(out_prefix) + "_summary.md")
    with md_path.open("w") as f:
        f.write("# off_2270: anomalous high-speed component diagnosis\n\n")
        f.write("- Sign convention: `+r` is outward.\n")
        f.write("- Component definition: largest connected high-speed component in each frame.\n")
        f.write("- Focus: plasma beta, Alfvenic speed, momentum force budget, residual, acceleration and direction.\n\n")

        f.write("## Global stats\n\n")
        f.write(f"- Frames analyzed: `{len(rows)}`\n")
        f.write(
            f"- Mean `jxb_dom_ratio_p50`: `{np.nanmean(arr_jdom):.4g}` "
            f"(>1 means |(JxB)_r| dominates |-(dp/dr)|+|rho g_r| on median)\n"
        )
        f.write(f"- Mean `comp_a_res_r_pos_frac`: `{np.nanmean(arr_ares_out):.4g}`\n")
        f.write(f"- Mean `comp_divv_p50`: `{np.nanmean(arr_div):.4g}`\n\n")

        if i_v >= 0:
            r = rows[i_v]
            f.write(
                f"- Peak high-speed (`comp_vabs_p95`) at `{r['frame']}`: "
                f"`{float(r['comp_vabs_p95']):.6g}`, beta_p50=`{float(r['comp_beta_p50']):.6g}`, "
                f"vA_p50=`{float(r['comp_va_p50']):.6g}`\n"
            )
        if i_b >= 0:
            r = rows[i_b]
            f.write(
                f"- Lowest beta (`comp_beta_p50`) at `{r['frame']}`: "
                f"`{float(r['comp_beta_p50']):.6g}`, vabs_p95=`{float(r['comp_vabs_p95']):.6g}`\n"
            )
        if i_va >= 0:
            r = rows[i_va]
            f.write(
                f"- Highest Alfven speed (`comp_va_p50`) at `{r['frame']}`: "
                f"`{float(r['comp_va_p50']):.6g}`, beta_p50=`{float(r['comp_beta_p50']):.6g}`\n"
            )

        f.write("\n## Top 10 high-speed frames\n\n")
        f.write("| frame | time | vabs_p95 | beta_p50 | va_p50 | jxb_p50 | f_res_p50 | a_res_p50 | frac(a_res>0) |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for i in top_idx:
            r = rows[i]
            f.write(
                "| {frame} | {time:.6e} | {v:.6g} | {beta:.6g} | {va:.6g} | {jxb:.6g} | {fres:.6g} | {ares:.6g} | {aout:.4f} |\n".format(
                    frame=r["frame"],
                    time=float(r["time"]),
                    v=float(r["comp_vabs_p95"]),
                    beta=float(r["comp_beta_p50"]),
                    va=float(r["comp_va_p50"]),
                    jxb=float(r["comp_jxb_r_p50"]),
                    fres=float(r["comp_f_res_r_p50"]),
                    ares=float(r["comp_a_res_r_p50"]),
                    aout=float(r["comp_a_res_r_pos_frac"]),
                )
            )
    return md_path


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    par_path = case_dir / args.par
    if not par_path.exists():
        raise FileNotFoundError(f"Missing par file: {par_path}")

    cfg = read_par_config(par_path)
    base_edges = build_base_radial_edges(
        xmin=float(cfg["xprobmin1"]),
        xmax=float(cfg["xprobmax1"]),
        ncells=int(cfg["domain_nx1"]),
        q=float(cfg["qstretch_baselevel"]),
    )

    paths = select_paths(sorted(case_dir.glob(args.pattern)), args.start, args.end, args.stride)
    if not paths:
        raise RuntimeError(f"No dat matched pattern {args.pattern} in {case_dir}")

    rows: List[Dict[str, float | str | int]] = []
    rng = np.random.default_rng(args.seed)
    for i, p in enumerate(paths, start=1):
        row = analyze_frame(
            dat_path=p,
            par_cfg=cfg,
            base_edges=base_edges,
            vfix=args.vfix,
            adaptive_threshold=args.adaptive_threshold,
            min_comp_cells=args.min_comp_cells,
            root_rmin=args.root_rmin,
            root_rmax=args.root_rmax,
            gravity_g0=args.gravity_g0,
            sradius=args.sradius,
            sample_max=args.sample_max,
            rng=rng,
        )
        rows.append(row)
        print(
            f"[{i}/{len(paths)}] {row['frame']} "
            f"size={row['comp_size']} rmin={float(row['comp_rmin']):.3f} "
            f"beta50={float(row['comp_beta_p50']):.3e} vA50={float(row['comp_va_p50']):.3f} "
            f"|v|95={float(row['comp_vabs_p95']):.3f}"
        )

    rows.sort(key=lambda r: int(r["idx"]))

    out_prefix = case_dir / args.output_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(str(out_prefix) + "_timeseries.dat")
    write_csv(csv_path, rows)
    print(f"Wrote {csv_path}")

    plot_paths: List[Path] = []
    if args.make_plots:
        plot_paths = maybe_make_plots(rows, out_prefix)
        for pp in plot_paths:
            print(f"Wrote {pp}")

    md_path = write_markdown_summary(rows, out_prefix)
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
