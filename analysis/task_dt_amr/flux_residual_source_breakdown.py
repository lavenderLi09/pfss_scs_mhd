# Copied to project analysis task toolbox
# Source case path: off_2270_initwind/off_2270_initwind_hpc/analysis/flux_residual_source_breakdown.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_initwind/off_2270_initwind_hpc/analysis/flux_residual_source_breakdown.py

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List

import numpy as np

from dt_cfl_common import (
    EPS,
    block_radial_edges,
    block_uniform_edges,
    build_base_radial_edges,
    get_header,
    get_tree_info,
    parse_log,
    read_block_fields_selected,
    read_par_config,
)


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _phi_delta(a: np.ndarray, b: float) -> np.ndarray:
    d = np.abs(a - b)
    return np.minimum(d, 2.0 * math.pi - d)


def _compute_profile_metrics(
    dat_path: Path,
    qstretch: float,
    th0: float,
    ph0: float,
    nbin: int,
    win_deg: float,
) -> Dict[str, float]:
    with dat_path.open("rb") as handle:
        header = get_header(handle)
        block_lvls, block_ixs, block_offsets = get_tree_info(handle)

        names = {n.strip(): i for i, n in enumerate(header["w_names"])}
        fmap = {"rho": names["rho"], "m1": names["m1"]}

        block_shape = header["block_nx"].astype(int)
        domain_nx = header["domain_nx"].astype(int)
        xmins = [float(v) for v in np.asarray(header["xmin"], dtype=np.float64)]
        xmaxs = [float(v) for v in np.asarray(header["xmax"], dtype=np.float64)]
        ndim = int(header["ndim"])

        rmin = xmins[0]
        rmax = xmaxs[0]
        dr_bin = (rmax - rmin) / float(nbin)
        r_cent = rmin + (np.arange(nbin, dtype=np.float64) + 0.5) * dr_bin

        sum_g = np.zeros(nbin, dtype=np.float64)
        w_g = np.zeros(nbin, dtype=np.float64)
        sum_l = np.zeros(nbin, dtype=np.float64)
        w_l = np.zeros(nbin, dtype=np.float64)

        base_edges = build_base_radial_edges(rmin, rmax, int(domain_nx[0]), qstretch)
        win = math.radians(win_deg)

        for level, bix, offset in zip(block_lvls.astype(int), block_ixs.astype(int), block_offsets.astype(int)):
            fld = read_block_fields_selected(handle, int(offset), block_shape, ndim, fmap)
            rho = fld["rho"]
            m1 = fld["m1"]
            vr = m1 / np.maximum(rho, EPS)

            r_edges = block_radial_edges(base_edges, int(level), int(bix[0]), int(block_shape[0]))
            th_edges = block_uniform_edges(
                xmins[1], xmaxs[1], int(domain_nx[1]), int(level), int(bix[1]), int(block_shape[1])
            )
            ph_edges = block_uniform_edges(
                xmins[2], xmaxs[2], int(domain_nx[2]), int(level), int(bix[2]), int(block_shape[2])
            )
            rc = 0.5 * (r_edges[:-1] + r_edges[1:])
            tc = 0.5 * (th_edges[:-1] + th_edges[1:])
            pc = 0.5 * (ph_edges[:-1] + ph_edges[1:])
            dth = np.diff(th_edges)
            dph = np.diff(ph_edges)

            rr, tt, pp = np.meshgrid(rc, tc, pc, indexing="ij")
            flux = rho * vr * rr * rr
            wt = np.maximum(np.sin(tt), EPS) * dth[None, :, None] * dph[None, None, :]

            bidx = np.floor((rr - rmin) / max(dr_bin, EPS)).astype(np.int64)
            bidx = np.clip(bidx, 0, nbin - 1)

            num = np.bincount(bidx.ravel(), weights=(flux * wt).ravel(), minlength=nbin)
            den = np.bincount(bidx.ravel(), weights=wt.ravel(), minlength=nbin)
            sum_g += num
            w_g += den

            mask_local = (np.abs(tt - th0) <= win) & (_phi_delta(pp, ph0) <= win)
            if np.any(mask_local):
                num_l = np.bincount(
                    bidx[mask_local].ravel(),
                    weights=(flux[mask_local] * wt[mask_local]).ravel(),
                    minlength=nbin,
                )
                den_l = np.bincount(
                    bidx[mask_local].ravel(),
                    weights=wt[mask_local].ravel(),
                    minlength=nbin,
                )
                sum_l += num_l
                w_l += den_l

        Fg = np.full(nbin, np.nan, dtype=np.float64)
        Fl = np.full(nbin, np.nan, dtype=np.float64)
        mg = w_g > 0.0
        ml = w_l > 0.0
        Fg[mg] = sum_g[mg] / w_g[mg]
        Fl[ml] = sum_l[ml] / w_l[ml]

        Jg = np.full(nbin - 1, np.nan, dtype=np.float64)
        Jl = np.full(nbin - 1, np.nan, dtype=np.float64)
        vg = np.isfinite(Fg[:-1]) & np.isfinite(Fg[1:])
        vl = np.isfinite(Fl[:-1]) & np.isfinite(Fl[1:])
        Jg[vg] = np.abs(Fg[1:][vg] - Fg[:-1][vg]) / np.maximum(np.abs(Fg[:-1][vg]), EPS)
        Jl[vl] = np.abs(Fl[1:][vl] - Fl[:-1][vl]) / np.maximum(np.abs(Fl[:-1][vl]), EPS)

        if np.any(np.isfinite(Jg)):
            ig = int(np.nanargmax(Jg))
            jg_max = float(Jg[ig])
            rg_max = float(r_cent[ig])
        else:
            ig = -1
            jg_max = math.nan
            rg_max = math.nan

        if np.any(np.isfinite(Jl)):
            il = int(np.nanargmax(Jl))
            jl_max = float(Jl[il])
            rl_max = float(r_cent[il])
        else:
            il = -1
            jl_max = math.nan
            rl_max = math.nan

        frac_g = float(np.mean(Jg[np.isfinite(Jg)] > 0.30)) if np.any(np.isfinite(Jg)) else math.nan
        frac_l = float(np.mean(Jl[np.isfinite(Jl)] > 0.30)) if np.any(np.isfinite(Jl)) else math.nan

        return {
            "jg_max": jg_max,
            "rg_max": rg_max,
            "jl_max": jl_max,
            "rl_max": rl_max,
            "frac_g_gt_0p30": frac_g,
            "frac_l_gt_0p30": frac_l,
        }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--nbin", type=int, default=256)
    ap.add_argument("--win-deg", type=float, default=5.0)
    args = ap.parse_args()

    case = args.case_dir.resolve()
    out_dir = (args.out_dir or (case / "analysis/dt_cfl_hpc")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    per = _read_csv(case / "analysis/dt_cfl_hpc/rho_bottleneck_diagnostic_per_frame.csv")
    gate = {r["frame"]: r for r in _read_csv(case / "analysis/dt_cfl_hpc/rho_bottleneck_gate_report.csv")}
    top = _read_csv(case / "analysis/dt_cfl_hpc/dt_cfl_topcells.csv")
    top1 = {}
    for r in top:
        if int(float(r.get("rank", "1"))) == 1 and r["frame"] not in top1:
            top1[r["frame"]] = r

    par = read_par_config(case / "amrvac.par")
    qstretch = float(par["qstretch_baselevel"])
    log_by_it = parse_log(case / "data/amr_probe.log")

    rows = []
    for r in per:
        fr = r["frame"]
        gr = gate.get(fr, {})
        if gr.get("flux_residual_flag", "False").lower() != "true":
            continue
        t = top1.get(fr)
        if t is None:
            continue
        dat = case / "data" / f"{fr}.dat"
        if not dat.exists():
            continue
        th0 = float(t["theta"])
        ph0 = float(t["phi"])
        m = _compute_profile_metrics(dat, qstretch, th0, ph0, args.nbin, args.win_deg)
        it = int(float(r["it"]))
        dtlog = log_by_it.get(it, {}).get("dt", math.nan)
        row = {
            "frame": fr,
            "it": it,
            "time": float(r["time"]),
            "dt_log": dtlog,
            "flux_residual_flag": gr.get("flux_residual_flag", ""),
            "flux_residual_used": float(gr.get("flux_residual_flag", "False").lower() == "true"),
            **m,
        }
        rows.append(row)

    csv_out = out_dir / "flux_residual_source_breakdown.csv"
    with csv_out.open("w", newline="") as f:
        if rows:
            keys = list(rows[0].keys())
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        else:
            f.write("")

    md_out = out_dir / "flux_residual_source_breakdown.md"
    lines = ["# Flux Residual Source Breakdown", ""]
    lines.append(f"- analyzed flagged frames: {len(rows)}")
    if rows:
        rg = [r["rg_max"] for r in rows if math.isfinite(r["rg_max"])]
        rl = [r["rl_max"] for r in rows if math.isfinite(r["rl_max"])]
        fg = [r["frac_g_gt_0p30"] for r in rows if math.isfinite(r["frac_g_gt_0p30"])]
        fl = [r["frac_l_gt_0p30"] for r in rows if math.isfinite(r["frac_l_gt_0p30"])]
        if rg:
            lines.append(f"- global max-jump radius mean/std: {np.mean(rg):.6f} / {np.std(rg):.6f}")
        if rl:
            lines.append(f"- local max-jump radius mean/std: {np.mean(rl):.6f} / {np.std(rl):.6f}")
        if fg:
            lines.append(f"- global jump-width mean(frac J>0.30): {np.mean(fg):.6f}")
        if fl:
            lines.append(f"- local jump-width mean(frac J>0.30): {np.mean(fl):.6f}")
        lines.append("")
        lines.append("## Interpretation Rule")
        lines.append("- fixed-shell jump: low std of rg_max and low global jump-width")
        lines.append("- global discontinuity: high global jump-width across many bins")
        lines.append("- local driver: local jump-width >> global jump-width")
    md_out.write_text("\n".join(lines) + "\n")

    print(csv_out)
    print(md_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
