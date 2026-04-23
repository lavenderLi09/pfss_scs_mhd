#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
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
    _gradient_axis,
    block_radial_edges,
    block_uniform_edges,
    build_base_radial_edges,
    conservative_internal_to_primitive,
    frame_id,
    read_block_fields,
    read_par_config,
    select_paths,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Probe a fixed onset angular location (theta,phi), then for each frame "
            "find max-speed point along that radial line within r<r_probe_max."
        )
    )
    p.add_argument("--case-dir", type=Path, default=CASE_DIR)
    p.add_argument("--par", default="amrvac.par")
    p.add_argument("--pattern", default="data/off*.dat")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=181)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--r-probe-max", type=float, default=5.0)
    p.add_argument("--r-onset-max", type=float, default=2.0)
    p.add_argument("--onset-start", type=int, default=40)
    p.add_argument("--onset-end", type=int, default=50)
    p.add_argument("--onset-vmin", type=float, default=6.0)
    p.add_argument("--n-lines", type=int, default=1, help="How many onset lines to track (top-N by onset-window speed).")
    p.add_argument(
        "--line-sep-cells",
        type=int,
        default=2,
        help="Minimum separation in base-grid cells between selected onset lines (theta/phi).",
    )
    p.add_argument("--gravity-g0", type=float, default=-14.072)
    p.add_argument("--sradius", type=float, default=1.0)

    p.add_argument("--unit-v-kms", type=float, default=116.448846777562)
    p.add_argument("--unit-t-sec", type=float, default=5972.5794)
    p.add_argument("--unit-l-rs", type=float, default=1.0)

    p.add_argument("--output-prefix", default="analysis/highstream_radial_line_probe")
    p.add_argument("--make-plot", action="store_true", default=True)
    p.add_argument("--no-make-plot", dest="make_plot", action="store_false")
    return p.parse_args()


def _read_header_tree_compat(dat_path: Path):
    try:
        hdr = get_header(str(dat_path))
        tree = get_tree_info(str(dat_path))
    except Exception:
        with dat_path.open("rb") as fh:
            hdr = get_header(fh)
            tree = get_tree_info(fh)

    if isinstance(tree, dict):
        offset_blocks = np.asarray(tree["offset_blocks"]).astype(np.int64)
        node_info = np.asarray(tree["node_info"])
        levels = node_info[:, 0].astype(np.int64)
        indices = node_info[:, 1:4].astype(np.int64)
    elif isinstance(tree, tuple) and len(tree) >= 3:
        levels = np.asarray(tree[0]).astype(np.int64)
        indices = np.asarray(tree[1]).astype(np.int64)
        offset_blocks = np.asarray(tree[2]).astype(np.int64)
    else:
        raise RuntimeError(f"Unsupported tree structure from datfile_io for {dat_path}")
    return hdr, levels, indices, offset_blocks


def _theta_phi_indices(theta: float, phi: float, cfg: Dict[str, float]) -> Tuple[int, int]:
    nt = int(cfg["domain_nx2"])
    np_ = int(cfg["domain_nx3"])
    tmin = float(cfg["xprobmin2"])
    tmax = float(cfg["xprobmax2"])
    pmin = float(cfg["xprobmin3"])
    pmax = float(cfg["xprobmax3"])
    jt = int(np.clip(np.floor((theta - tmin) / (tmax - tmin) * nt), 0, nt - 1))
    kp = int(np.mod(np.floor((phi - pmin) / (pmax - pmin) * np_), np_))
    return jt, kp


def _scan_frame_line(
    dat_path: Path,
    cfg: Dict[str, float],
    base_edges: np.ndarray,
    jt_fixed: int,
    kp_fixed: int,
    theta_fixed: float,
    phi_fixed: float,
    r_probe_max: float,
    gravity_g0: float,
    sradius: float,
) -> Dict[str, float]:
    hdr, levels, indices, offset_blocks = _read_header_tree_compat(dat_path)
    ndim = int(hdr["ndim"])
    nw = int(hdr["nw"])
    gamma = float(cfg["mhd_gamma"])
    nt = int(cfg["domain_nx2"])
    np_ = int(cfg["domain_nx3"])
    block_shape = hdr["block_nx"].astype(np.int64)

    rr_all: List[float] = []
    vv_all: List[float] = []
    vr_all: List[float] = []
    beta_all: List[float] = []
    ares_all: List[float] = []
    fres_all: List[float] = []
    fpres_all: List[float] = []
    fjxb_all: List[float] = []
    frhog_all: List[float] = []
    apres_all: List[float] = []
    ajxb_all: List[float] = []
    ag_all: List[float] = []

    with dat_path.open("rb") as f:
        for ib in range(offset_blocks.shape[0]):
            lvl = int(levels[ib])
            ix1, ix2, ix3 = [int(v) for v in indices[ib]]
            r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
            t_edges = block_uniform_edges(
                float(cfg["xprobmin2"]), float(cfg["xprobmax2"]), nt, lvl, ix2, int(block_shape[1])
            )
            p_edges = block_uniform_edges(
                float(cfg["xprobmin3"]), float(cfg["xprobmax3"]), np_, lvl, ix3, int(block_shape[2])
            )
            r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
            t_cent = 0.5 * (t_edges[:-1] + t_edges[1:])
            p_cent = 0.5 * (p_edges[:-1] + p_edges[1:])

            base_t_idx = np.clip(
                np.floor((t_cent - float(cfg["xprobmin2"])) / (float(cfg["xprobmax2"]) - float(cfg["xprobmin2"])) * nt).astype(int),
                0,
                nt - 1,
            )
            base_p_idx = np.mod(
                np.floor((p_cent - float(cfg["xprobmin3"])) / (float(cfg["xprobmax3"]) - float(cfg["xprobmin3"])) * np_).astype(int),
                np_,
            )

            it_loc = np.where(base_t_idx == jt_fixed)[0]
            ip_loc = np.where(base_p_idx == kp_fixed)[0]
            if it_loc.size == 0 or ip_loc.size == 0:
                continue

            # Robustly lock onto the target angular position inside this block,
            # instead of taking the first matching local index.
            it0 = int(it_loc[np.argmin(np.abs(t_cent[it_loc] - theta_fixed))])
            dphi = np.abs((p_cent[ip_loc] - phi_fixed + np.pi) % (2.0 * np.pi) - np.pi)
            ip0 = int(ip_loc[np.argmin(dphi)])

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
            beta = 2.0 * pth / np.maximum(b2, EPS)

            rr = r_cent[:, None, None]
            th = t_cent[None, :, None]
            sin_t = np.maximum(np.sin(th), 1.0e-8)

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

            curl_t = (d_br_dph / np.maximum(sin_t, EPS) - d_rbp_dr) / np.maximum(rr, EPS)
            curl_p = (d_rbt_dr - d_br_dth) / np.maximum(rr, EPS)
            jxb_r = curl_t * b_p - curl_p * b_t

            g_r = gravity_g0 * (sradius * sradius) / np.maximum(rr * rr, EPS)
            rho_g_r = rho * g_r
            f_res_r = f_pres_r + jxb_r + rho_g_r
            rho_safe = np.maximum(rho, EPS)
            a_res_r = f_res_r / rho_safe
            a_pres_r = f_pres_r / rho_safe
            a_jxb_r = jxb_r / rho_safe
            a_g_r = rho_g_r / rho_safe

            r_line = r_cent
            v_line = vabs[:, it0, ip0]
            vr_line = v_r[:, it0, ip0]
            beta_line = beta[:, it0, ip0]
            fres_line = f_res_r[:, it0, ip0]
            ares_line = a_res_r[:, it0, ip0]
            fpres_line = f_pres_r[:, it0, ip0]
            fjxb_line = jxb_r[:, it0, ip0]
            frhog_line = rho_g_r[:, it0, ip0]
            apres_line = a_pres_r[:, it0, ip0]
            ajxb_line = a_jxb_r[:, it0, ip0]
            ag_line = a_g_r[:, it0, ip0]

            mask = np.isfinite(r_line) & np.isfinite(v_line) & np.isfinite(ares_line) & (r_line <= r_probe_max)
            if not np.any(mask):
                continue

            rr_all.extend(r_line[mask].tolist())
            vv_all.extend(v_line[mask].tolist())
            vr_all.extend(vr_line[mask].tolist())
            beta_all.extend(beta_line[mask].tolist())
            ares_all.extend(ares_line[mask].tolist())
            fres_all.extend(fres_line[mask].tolist())
            fpres_all.extend(fpres_line[mask].tolist())
            fjxb_all.extend(fjxb_line[mask].tolist())
            frhog_all.extend(frhog_line[mask].tolist())
            apres_all.extend(apres_line[mask].tolist())
            ajxb_all.extend(ajxb_line[mask].tolist())
            ag_all.extend(ag_line[mask].tolist())

    row = {
        "frame": dat_path.stem,
        "idx": float(frame_id(dat_path.stem) or -1),
        "time_code": float(hdr["time"]),
        "it": float(hdr["it"]),
        "n_line_samples": float(len(rr_all)),
    }
    if len(rr_all) == 0:
        row.update(
            {
                "r_vmax_code": math.nan,
                "vmax_code": math.nan,
                "vr_at_vmax_code": math.nan,
                "beta_at_vmax": math.nan,
                "fpres_at_vmax_code": math.nan,
                "fjxb_at_vmax_code": math.nan,
                "frhog_at_vmax_code": math.nan,
                "fres_at_vmax_code": math.nan,
                "apres_at_vmax_code": math.nan,
                "ajxb_at_vmax_code": math.nan,
                "ag_at_vmax_code": math.nan,
                "ares_at_vmax_code": math.nan,
                "force_abs_frac_pres": math.nan,
                "force_abs_frac_jxb": math.nan,
                "force_abs_frac_grav": math.nan,
                "acc_abs_frac_pres": math.nan,
                "acc_abs_frac_jxb": math.nan,
                "acc_abs_frac_grav": math.nan,
                "force_jxb_dom_ratio": math.nan,
                "acc_jxb_dom_ratio": math.nan,
            }
        )
        return row

    rr_np = np.asarray(rr_all, dtype=np.float64)
    vv_np = np.asarray(vv_all, dtype=np.float64)
    vr_np = np.asarray(vr_all, dtype=np.float64)
    beta_np = np.asarray(beta_all, dtype=np.float64)
    fres_np = np.asarray(fres_all, dtype=np.float64)
    ares_np = np.asarray(ares_all, dtype=np.float64)
    fpres_np = np.asarray(fpres_all, dtype=np.float64)
    fjxb_np = np.asarray(fjxb_all, dtype=np.float64)
    frhog_np = np.asarray(frhog_all, dtype=np.float64)
    apres_np = np.asarray(apres_all, dtype=np.float64)
    ajxb_np = np.asarray(ajxb_all, dtype=np.float64)
    ag_np = np.asarray(ag_all, dtype=np.float64)

    imax = int(np.nanargmax(vv_np))
    force_abs_sum = abs(fpres_np[imax]) + abs(fjxb_np[imax]) + abs(frhog_np[imax]) + EPS
    acc_abs_sum = abs(apres_np[imax]) + abs(ajxb_np[imax]) + abs(ag_np[imax]) + EPS
    row.update(
        {
            "r_vmax_code": float(rr_np[imax]),
            "vmax_code": float(vv_np[imax]),
            "vr_at_vmax_code": float(vr_np[imax]),
            "beta_at_vmax": float(beta_np[imax]),
            "fpres_at_vmax_code": float(fpres_np[imax]),
            "fjxb_at_vmax_code": float(fjxb_np[imax]),
            "frhog_at_vmax_code": float(frhog_np[imax]),
            "fres_at_vmax_code": float(fres_np[imax]),
            "apres_at_vmax_code": float(apres_np[imax]),
            "ajxb_at_vmax_code": float(ajxb_np[imax]),
            "ag_at_vmax_code": float(ag_np[imax]),
            "ares_at_vmax_code": float(ares_np[imax]),
            "force_abs_frac_pres": float(abs(fpres_np[imax]) / force_abs_sum),
            "force_abs_frac_jxb": float(abs(fjxb_np[imax]) / force_abs_sum),
            "force_abs_frac_grav": float(abs(frhog_np[imax]) / force_abs_sum),
            "acc_abs_frac_pres": float(abs(apres_np[imax]) / acc_abs_sum),
            "acc_abs_frac_jxb": float(abs(ajxb_np[imax]) / acc_abs_sum),
            "acc_abs_frac_grav": float(abs(ag_np[imax]) / acc_abs_sum),
            "force_jxb_dom_ratio": float(abs(fjxb_np[imax]) / (abs(fpres_np[imax]) + abs(frhog_np[imax]) + EPS)),
            "acc_jxb_dom_ratio": float(abs(ajxb_np[imax]) / (abs(apres_np[imax]) + abs(ag_np[imax]) + EPS)),
        }
    )
    return row


def _phi_cell_dist(a: int, b: int, np_: int) -> int:
    d = abs(a - b)
    return min(d, np_ - d)


def _far_enough(jt: int, kp: int, chosen: List[Tuple[int, int]], sep: int, np_: int) -> bool:
    for j0, k0 in chosen:
        if abs(jt - j0) <= sep and _phi_cell_dist(kp, k0, np_) <= sep:
            return False
    return True


def _pick_onset_lines(
    ref_paths: Sequence[Path],
    cfg: Dict[str, float],
    base_edges: np.ndarray,
    r_onset_max: float,
    onset_vmin: float,
    n_lines: int,
    line_sep_cells: int,
) -> List[Dict[str, float]]:
    candidates: List[Tuple[float, float, float, float, float, int, int]] = []
    np_ = int(cfg["domain_nx3"])
    for p in ref_paths:
        hdr, levels, indices, offset_blocks = _read_header_tree_compat(p)
        ndim = int(hdr["ndim"])
        nw = int(hdr["nw"])
        gamma = float(cfg["mhd_gamma"])
        nt = int(cfg["domain_nx2"])
        np_ = int(cfg["domain_nx3"])
        block_shape = hdr["block_nx"].astype(np.int64)

        with p.open("rb") as f:
            for ib in range(offset_blocks.shape[0]):
                lvl = int(levels[ib])
                ix1, ix2, ix3 = [int(v) for v in indices[ib]]
                r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
                t_edges = block_uniform_edges(
                    float(cfg["xprobmin2"]), float(cfg["xprobmax2"]), nt, lvl, ix2, int(block_shape[1])
                )
                p_edges = block_uniform_edges(
                    float(cfg["xprobmin3"]), float(cfg["xprobmax3"]), np_, lvl, ix3, int(block_shape[2])
                )
                r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
                t_cent = 0.5 * (t_edges[:-1] + t_edges[1:])
                p_cent = 0.5 * (p_edges[:-1] + p_edges[1:])

                fields = read_block_fields(f, int(offset_blocks[ib]), block_shape, ndim, nw)
                rho = fields[0]
                v_r, v_t, v_p, _pth = conservative_internal_to_primitive(
                    rho=rho, m1=fields[1], m2=fields[2], m3=fields[3], eint=fields[4], gamma=gamma
                )
                vabs = np.sqrt(v_r * v_r + v_t * v_t + v_p * v_p)

                rr = r_cent[:, None, None]
                mask = np.isfinite(vabs) & (rr <= r_onset_max)
                if not np.any(mask):
                    continue

                cand = np.where(mask & (vabs >= onset_vmin))
                if cand[0].size == 0:
                    continue

                vals = vabs[cand]
                kpick = min(16, vals.size)
                top = np.argpartition(vals, -kpick)[-kpick:]
                for ii in top:
                    ir = int(cand[0][ii])
                    it = int(cand[1][ii])
                    ip = int(cand[2][ii])
                    vcode = float(vabs[ir, it, ip])
                    theta = float(t_cent[it])
                    phi = float(p_cent[ip])
                    rr0 = float(r_cent[ir])
                    fid = float(frame_id(p.stem) or -1)
                    jt, kp = _theta_phi_indices(theta, phi, cfg)
                    candidates.append((vcode, theta, phi, rr0, fid, jt, kp))

    if not candidates:
        raise RuntimeError("Failed to find onset angular location in reference window.")

    candidates.sort(key=lambda x: x[0], reverse=True)
    chosen_idx: List[Tuple[int, int]] = []
    out: List[Dict[str, float]] = []
    for vcode, theta, phi, rr0, fid, jt, kp in candidates:
        if not _far_enough(jt, kp, chosen_idx, max(0, int(line_sep_cells)), np_):
            continue
        chosen_idx.append((jt, kp))
        out.append(
            {
                "theta_rad": theta,
                "phi_rad": phi,
                "onset_ref_r_Rs": rr0,
                "onset_ref_frame_idx": fid,
                "jt": float(jt),
                "kp": float(kp),
                "onset_vmax_code": vcode,
            }
        )
        if len(out) >= max(1, int(n_lines)):
            break

    if not out:
        raise RuntimeError("No onset lines selected after separation filtering.")
    return out


def write_csv(path: Path, rows: List[Dict[str, float]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    keys = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def maybe_plot(rows: List[Dict[str, float]], out_prefix: Path) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_line: Dict[int, List[Dict[str, float]]] = {}
    for r in rows:
        lid = int(round(float(r.get("line_id", 0.0))))
        by_line.setdefault(lid, []).append(r)

    out_paths: List[Path] = []
    for lid in sorted(by_line.keys()):
        rr = sorted(by_line[lid], key=lambda x: float(x["idx"]))
        tmin = np.array([r["time_min"] for r in rr], dtype=np.float64)
        vmax = np.array([r["vmax_km_s"] for r in rr], dtype=np.float64)
        vrmax = np.array([r["vr_at_vmax_km_s"] for r in rr], dtype=np.float64)
        ares = np.array([r["ares_at_vmax_km_s2"] for r in rr], dtype=np.float64)
        apres = np.array([r["apres_at_vmax_km_s2"] for r in rr], dtype=np.float64)
        ajxb = np.array([r["ajxb_at_vmax_km_s2"] for r in rr], dtype=np.float64)
        agrav = np.array([r["ag_at_vmax_km_s2"] for r in rr], dtype=np.float64)
        fpres_frac = np.array([r["force_abs_frac_pres"] for r in rr], dtype=np.float64)
        fjxb_frac = np.array([r["force_abs_frac_jxb"] for r in rr], dtype=np.float64)
        fgrav_frac = np.array([r["force_abs_frac_grav"] for r in rr], dtype=np.float64)
        beta = np.array([r["beta_at_vmax"] for r in rr], dtype=np.float64)
        rv = np.array([r["r_vmax_Rs"] for r in rr], dtype=np.float64)
        th = float(rr[0]["theta_fixed_rad"])
        ph = float(rr[0]["phi_fixed_rad"])

        fig, axes = plt.subplots(5, 1, figsize=(10, 14), sharex=True)
        axes[0].plot(tmin, vmax, "o-", label="|v| at radial-line max")
        axes[0].plot(tmin, vrmax, "s--", label="v_r at same point")
        axes[0].set_ylabel("km/s")
        axes[0].legend(loc="best")
        axes[0].grid(alpha=0.3)

        axes[1].plot(tmin, ares, "o-", color="tab:red", label="a_res,r")
        axes[1].plot(tmin, ajxb, "-", color="tab:blue", label="a_jxb,r")
        axes[1].plot(tmin, apres, "-", color="tab:orange", label="a_pres,r")
        axes[1].plot(tmin, agrav, "-", color="tab:gray", label="a_grav,r")
        axes[1].axhline(0.0, color="k", lw=1)
        axes[1].set_ylabel("km/s^2")
        axes[1].set_yscale("symlog", linthresh=1e-3)
        axes[1].legend(loc="best")
        axes[1].grid(alpha=0.3)

        axes[2].plot(tmin, fjxb_frac, "o-", color="tab:blue", label="|JxB| / sum|F_i|")
        axes[2].plot(tmin, fpres_frac, "s-", color="tab:orange", label="|-(dp/dr)| / sum|F_i|")
        axes[2].plot(tmin, fgrav_frac, "^-", color="tab:gray", label="|rho g| / sum|F_i|")
        axes[2].set_ylim(-0.02, 1.02)
        axes[2].set_ylabel("force abs fraction")
        axes[2].legend(loc="best")
        axes[2].grid(alpha=0.3)

        axes[3].plot(tmin, beta, "o-", color="tab:green")
        axes[3].set_yscale("log")
        axes[3].set_ylabel("beta")
        axes[3].grid(alpha=0.3)

        axes[4].plot(tmin, rv, "o-", color="tab:purple")
        axes[4].set_ylabel("r at line-max (Rs)")
        axes[4].set_xlabel("time (minutes)")
        axes[4].grid(alpha=0.3)

        fig.suptitle(f"Line {lid}: theta={th:.4f} rad, phi={ph:.4f} rad")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        p = Path(str(out_prefix) + f"_line{lid:02d}_timeseries.png")
        fig.savefig(p, dpi=180)
        plt.close(fig)
        out_paths.append(p)
    return out_paths


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    par_path = case_dir / args.par
    cfg = read_par_config(par_path)
    base_edges = build_base_radial_edges(
        xmin=float(cfg["xprobmin1"]),
        xmax=float(cfg["xprobmax1"]),
        ncells=int(cfg["domain_nx1"]),
        q=float(cfg["qstretch_baselevel"]),
    )

    paths = select_paths(sorted(case_dir.glob(args.pattern)), args.start, args.end, args.stride)
    if not paths:
        raise RuntimeError("No dat files selected.")

    ref_paths = [p for p in paths if args.onset_start <= int(frame_id(p.stem) or -1) <= args.onset_end]
    if not ref_paths:
        raise RuntimeError("No onset-window frames selected.")

    onset_lines = _pick_onset_lines(
        ref_paths=ref_paths,
        cfg=cfg,
        base_edges=base_edges,
        r_onset_max=args.r_onset_max,
        onset_vmin=args.onset_vmin,
        n_lines=args.n_lines,
        line_sep_cells=args.line_sep_cells,
    )

    rows: List[Dict[str, float]] = []
    unit_a = args.unit_v_kms / args.unit_t_sec
    npaths = len(paths)
    nlines = len(onset_lines)
    for lid, line in enumerate(onset_lines, start=1):
        theta0 = float(line["theta_rad"])
        phi0 = float(line["phi_rad"])
        jt0 = int(round(float(line["jt"])))
        kp0 = int(round(float(line["kp"])))
        r0 = float(line["onset_ref_r_Rs"])
        ref_idx = float(line["onset_ref_frame_idx"])
        v0 = float(line["onset_vmax_code"])
        print(
            f"LINE {lid}/{nlines}: theta={theta0:.6f}, phi={phi0:.6f}, "
            f"base(jt={jt0},kp={kp0}), onset_frame={ref_idx:.0f}, onset_r={r0:.3f}Rs, onset_v={v0*args.unit_v_kms:.2f} km/s"
        )

        for i, p in enumerate(paths, start=1):
            r = _scan_frame_line(
                dat_path=p,
                cfg=cfg,
                base_edges=base_edges,
                jt_fixed=jt0,
                kp_fixed=kp0,
                theta_fixed=theta0,
                phi_fixed=phi0,
                r_probe_max=args.r_probe_max,
                gravity_g0=args.gravity_g0,
                sradius=args.sradius,
            )
            r["line_id"] = float(lid)
            r["time_min"] = r["time_code"] * args.unit_t_sec / 60.0
            r["r_vmax_Rs"] = r["r_vmax_code"] * args.unit_l_rs
            r["vmax_km_s"] = r["vmax_code"] * args.unit_v_kms
            r["vr_at_vmax_km_s"] = r["vr_at_vmax_code"] * args.unit_v_kms
            r["fpres_at_vmax_code"] = r["fpres_at_vmax_code"]
            r["fjxb_at_vmax_code"] = r["fjxb_at_vmax_code"]
            r["frhog_at_vmax_code"] = r["frhog_at_vmax_code"]
            r["fres_at_vmax_code"] = r["fres_at_vmax_code"]
            r["apres_at_vmax_km_s2"] = r["apres_at_vmax_code"] * unit_a
            r["ajxb_at_vmax_km_s2"] = r["ajxb_at_vmax_code"] * unit_a
            r["ag_at_vmax_km_s2"] = r["ag_at_vmax_code"] * unit_a
            r["ares_at_vmax_km_s2"] = r["ares_at_vmax_code"] * unit_a
            r["theta_fixed_rad"] = theta0
            r["phi_fixed_rad"] = phi0
            r["theta_fixed_base_index"] = float(jt0)
            r["phi_fixed_base_index"] = float(kp0)
            r["onset_ref_r_Rs"] = r0
            r["onset_ref_frame_idx"] = ref_idx
            r["onset_ref_v_km_s"] = v0 * args.unit_v_kms
            rows.append(r)
            print(
                f"[line {lid}/{nlines} | {i}/{npaths}] {p.stem} "
                f"r={r['r_vmax_Rs']:.3f}Rs v={r['vmax_km_s']:.2f}km/s "
                f"a_res={r['ares_at_vmax_km_s2']:.3e} "
                f"(a_jxb={r['ajxb_at_vmax_km_s2']:.3e}, a_pres={r['apres_at_vmax_km_s2']:.3e}, a_g={r['ag_at_vmax_km_s2']:.3e}) "
                f"beta={r['beta_at_vmax']:.3e}"
            )

    rows = sorted(rows, key=lambda x: (x["line_id"], x["idx"]))
    out_prefix = case_dir / args.output_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_csv = Path(str(out_prefix) + "_timeseries.csv")
    write_csv(out_csv, rows)
    print(f"Wrote {out_csv}")
    if args.make_plot:
        out_pngs = maybe_plot(rows, out_prefix)
        for pp in out_pngs:
            print(f"Wrote {pp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
