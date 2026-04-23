#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
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
            "Fixed-point (Eulerian) diagnostics on selected lines: choose per-line fixed radius "
            "as the most frequent r-location of line-wise vmax, then output v/beta/force split/dvdt."
        )
    )
    p.add_argument("--case-dir", type=Path, default=CASE_DIR)
    p.add_argument("--par", default="amrvac.par")
    p.add_argument("--pattern", default="data/off*.dat")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=181)
    p.add_argument("--stride", type=int, default=2)

    p.add_argument(
        "--line-source-csv",
        default="analysis/highstream_radial_line_off0000_0181_s2_n3_timeseries.csv",
        help="CSV from multi-line radial-line run to extract line geometry and r_mode.",
    )
    p.add_argument("--line-ids", default="", help="Comma-separated line ids, empty means all.")
    p.add_argument("--r-mode-round", type=int, default=3, help="Round digits for r_vmax mode counting.")

    p.add_argument("--gravity-g0", type=float, default=-14.072)
    p.add_argument("--sradius", type=float, default=1.0)
    p.add_argument("--unit-v-kms", type=float, default=116.448846777562)
    p.add_argument("--unit-t-sec", type=float, default=5972.5794)
    p.add_argument("--output-prefix", default="analysis/highstream_fixedpoint_lines")
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


def _parse_line_ids(s: str) -> set[int]:
    s = s.strip()
    if not s:
        return set()
    out = set()
    for t in s.split(","):
        tt = t.strip()
        if not tt:
            continue
        out.add(int(tt))
    return out


def _load_line_defs(csv_path: Path, select_ids: set[int], r_mode_round: int) -> List[Dict[str, float]]:
    by_line: Dict[int, List[Dict[str, float]]] = defaultdict(list)
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            lid = int(round(float(row["line_id"])))
            if select_ids and lid not in select_ids:
                continue
            by_line[lid].append(row)

    if not by_line:
        raise RuntimeError(f"No lines found in {csv_path} for selection={sorted(select_ids)}")

    defs: List[Dict[str, float]] = []
    for lid in sorted(by_line.keys()):
        rows = by_line[lid]
        theta = float(rows[0]["theta_fixed_rad"])
        phi = float(rows[0]["phi_fixed_rad"])
        jt = int(round(float(rows[0]["theta_fixed_base_index"])))
        kp = int(round(float(rows[0]["phi_fixed_base_index"])))

        rv = []
        v_for_r = defaultdict(list)
        for r in rows:
            try:
                rr = float(r["r_vmax_Rs"])
                vv = float(r["vmax_km_s"])
            except Exception:
                continue
            if math.isfinite(rr):
                key = round(rr, r_mode_round)
                rv.append(key)
                if math.isfinite(vv):
                    v_for_r[key].append(vv)
        if not rv:
            raise RuntimeError(f"line {lid}: no finite r_vmax in {csv_path}")

        c = Counter(rv)
        top_count = max(c.values())
        top_keys = [k for k, v in c.items() if v == top_count]
        if len(top_keys) == 1:
            r_mode = float(top_keys[0])
        else:
            # tie-break by larger mean vmax on that radius bucket
            r_mode = float(sorted(top_keys, key=lambda k: (np.mean(v_for_r[k]) if v_for_r[k] else -np.inf), reverse=True)[0])

        defs.append(
            {
                "line_id": float(lid),
                "theta_rad": theta,
                "phi_rad": phi,
                "jt": float(jt),
                "kp": float(kp),
                "r_fixed_Rs": r_mode,
                "r_mode_count": float(top_count),
            }
        )
    return defs


def _probe_fixed_point(
    dat_path: Path,
    cfg: Dict[str, float],
    base_edges: np.ndarray,
    line_def: Dict[str, float],
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

    theta = float(line_def["theta_rad"])
    phi = float(line_def["phi_rad"])
    jt = int(round(float(line_def["jt"])))
    kp = int(round(float(line_def["kp"])))
    r_target = float(line_def["r_fixed_Rs"])

    best = None
    best_dr = float("inf")
    best_level = -1

    with dat_path.open("rb") as f:
        for ib in range(offset_blocks.shape[0]):
            lvl = int(levels[ib])
            ix1, ix2, ix3 = [int(v) for v in indices[ib]]
            r_edges = block_radial_edges(base_edges, lvl, ix1, int(block_shape[0]))
            t_edges = block_uniform_edges(float(cfg["xprobmin2"]), float(cfg["xprobmax2"]), nt, lvl, ix2, int(block_shape[1]))
            p_edges = block_uniform_edges(float(cfg["xprobmin3"]), float(cfg["xprobmax3"]), np_, lvl, ix3, int(block_shape[2]))
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

            it_loc = np.where(base_t_idx == jt)[0]
            ip_loc = np.where(base_p_idx == kp)[0]
            if it_loc.size == 0 or ip_loc.size == 0:
                continue

            it0 = int(it_loc[np.argmin(np.abs(t_cent[it_loc] - theta))])
            dphi = np.abs((p_cent[ip_loc] - phi + np.pi) % (2.0 * np.pi) - np.pi)
            ip0 = int(ip_loc[np.argmin(dphi)])

            ir0 = int(np.argmin(np.abs(r_cent - r_target)))
            dr = abs(float(r_cent[ir0]) - r_target)

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
            a_pres_r = f_pres_r / rho_safe
            a_jxb_r = jxb_r / rho_safe
            a_g_r = rho_g_r / rho_safe
            a_res_r = f_res_r / rho_safe

            # Eulerian radial momentum closure term:
            # d_t v_r = a_res_r - (v · grad v_r)
            dvr_dr = _gradient_axis(v_r, r_cent, axis=0)
            dvr_dth = _gradient_axis(v_r, t_cent, axis=1)
            dvr_dph = _gradient_axis(v_r, p_cent, axis=2)
            adv_r = (
                v_r * dvr_dr
                + v_t / np.maximum(rr, EPS) * dvr_dth
                + v_p / np.maximum(rr * sin_t, EPS) * dvr_dph
                - (v_t * v_t + v_p * v_p) / np.maximum(rr, EPS)
            )
            rhs_dvrdt = a_res_r - adv_r

            cand = {
                "r_sample_Rs": float(r_cent[ir0]),
                "dr_abs_Rs": dr,
                "v_km_s": float(vabs[ir0, it0, ip0]),
                "vr_km_s": float(v_r[ir0, it0, ip0]),
                "beta": float(beta[ir0, it0, ip0]),
                "f_pres_code": float(f_pres_r[ir0, it0, ip0]),
                "f_jxb_code": float(jxb_r[ir0, it0, ip0]),
                "f_grav_code": float(rho_g_r[ir0, it0, ip0]),
                "f_res_code": float(f_res_r[ir0, it0, ip0]),
                "a_pres_code": float(a_pres_r[ir0, it0, ip0]),
                "a_jxb_code": float(a_jxb_r[ir0, it0, ip0]),
                "a_grav_code": float(a_g_r[ir0, it0, ip0]),
                "a_res_code": float(a_res_r[ir0, it0, ip0]),
                "adv_r_code": float(adv_r[ir0, it0, ip0]),
                "rhs_dvrdt_code": float(rhs_dvrdt[ir0, it0, ip0]),
                "level": float(lvl),
            }

            if (dr < best_dr - 1.0e-12) or (abs(dr - best_dr) <= 1.0e-12 and lvl > best_level):
                best = cand
                best_dr = dr
                best_level = lvl

    out = {
        "frame": dat_path.stem,
        "idx": float(frame_id(dat_path.stem) or -1),
        "time_code": float(hdr["time"]),
        "it": float(hdr["it"]),
        "line_id": float(line_def["line_id"]),
        "theta_rad": float(line_def["theta_rad"]),
        "phi_rad": float(line_def["phi_rad"]),
        "jt": float(line_def["jt"]),
        "kp": float(line_def["kp"]),
        "r_fixed_Rs": float(line_def["r_fixed_Rs"]),
        "r_mode_count": float(line_def["r_mode_count"]),
    }
    if best is None:
        for k in [
            "r_sample_Rs",
            "dr_abs_Rs",
            "v_km_s",
            "vr_km_s",
            "beta",
            "f_pres_code",
            "f_jxb_code",
            "f_grav_code",
            "f_res_code",
            "a_pres_code",
            "a_jxb_code",
            "a_grav_code",
            "a_res_code",
            "adv_r_code",
            "rhs_dvrdt_code",
            "level",
        ]:
            out[k] = math.nan
        return out

    out.update(best)
    return out


def write_csv(path: Path, rows: List[Dict[str, float]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    keys = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _append_time_derivatives(rows: List[Dict[str, float]], unit_t_sec: float) -> None:
    by_line: Dict[int, List[Dict[str, float]]] = defaultdict(list)
    for r in rows:
        lid = int(round(float(r["line_id"])))
        by_line[lid].append(r)

    for lid, rr in by_line.items():
        rr.sort(key=lambda x: float(x["idx"]))
        t = np.array([float(x["time_code"]) * unit_t_sec for x in rr], dtype=np.float64)
        v = np.array([float(x["v_km_s"]) for x in rr], dtype=np.float64)
        vr = np.array([float(x["vr_km_s"]) for x in rr], dtype=np.float64)

        dvdt = np.full_like(v, np.nan)
        dvrdt = np.full_like(vr, np.nan)
        if rr and np.count_nonzero(np.isfinite(v)) >= 2:
            dvdt = np.gradient(v, t, edge_order=1)
        if rr and np.count_nonzero(np.isfinite(vr)) >= 2:
            dvrdt = np.gradient(vr, t, edge_order=1)

        for i in range(len(rr)):
            rr[i]["dv_dt_km_s2"] = float(dvdt[i]) if np.isfinite(dvdt[i]) else math.nan
            rr[i]["dvr_dt_km_s2"] = float(dvrdt[i]) if np.isfinite(dvrdt[i]) else math.nan


def maybe_plot(rows: List[Dict[str, float]], out_prefix: Path, unit_t_sec: float) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_paths: List[Path] = []
    by_line: Dict[int, List[Dict[str, float]]] = defaultdict(list)
    for r in rows:
        lid = int(round(float(r["line_id"])))
        by_line[lid].append(r)

    for lid in sorted(by_line.keys()):
        rr = sorted(by_line[lid], key=lambda x: float(x["idx"]))
        tmin = np.array([float(r["time_code"]) * unit_t_sec / 60.0 for r in rr], dtype=np.float64)
        v = np.array([float(r["v_km_s"]) for r in rr], dtype=np.float64)
        vr = np.array([float(r["vr_km_s"]) for r in rr], dtype=np.float64)
        ares = np.array([float(r["a_res_km_s2"]) for r in rr], dtype=np.float64)
        apres = np.array([float(r["a_pres_km_s2"]) for r in rr], dtype=np.float64)
        ajxb = np.array([float(r["a_jxb_km_s2"]) for r in rr], dtype=np.float64)
        agrav = np.array([float(r["a_grav_km_s2"]) for r in rr], dtype=np.float64)
        adv = np.array([float(r.get("adv_r_km_s2", math.nan)) for r in rr], dtype=np.float64)
        rhs = np.array([float(r.get("rhs_dvrdt_km_s2", math.nan)) for r in rr], dtype=np.float64)
        dvdt = np.array([float(r.get("dv_dt_km_s2", math.nan)) for r in rr], dtype=np.float64)
        dvrdt = np.array([float(r.get("dvr_dt_km_s2", math.nan)) for r in rr], dtype=np.float64)
        beta = np.array([float(r["beta"]) for r in rr], dtype=np.float64)
        dr = np.array([float(r["dr_abs_Rs"]) for r in rr], dtype=np.float64)
        r_fix = float(rr[0]["r_fixed_Rs"])
        th = float(rr[0]["theta_rad"])
        ph = float(rr[0]["phi_rad"])

        fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
        axes[0].plot(tmin, v, "o-", label="|v| fixed-point")
        axes[0].plot(tmin, vr, "s--", label="v_r fixed-point")
        axes[0].set_ylabel("km/s")
        axes[0].legend(loc="best")
        axes[0].grid(alpha=0.3)

        axes[1].plot(tmin, ares, "o-", label="a_res")
        axes[1].plot(tmin, adv, "-", label="adv_r=(v·∇)v_r")
        axes[1].plot(tmin, rhs, "-", label="rhs=a_res-adv")
        axes[1].plot(tmin, ajxb, "-", label="a_jxb")
        axes[1].plot(tmin, apres, "-", label="a_pres")
        axes[1].plot(tmin, agrav, "-", label="a_grav")
        axes[1].plot(tmin, dvrdt, "k--", label="dvr/dt")
        axes[1].plot(tmin, dvdt, "k:", label="d|v|/dt")
        axes[1].axhline(0.0, color="k", lw=1)
        axes[1].set_ylabel("km/s^2")
        axes[1].set_yscale("symlog", linthresh=1e-3)
        axes[1].legend(loc="best", fontsize=8)
        axes[1].grid(alpha=0.3)

        axes[2].plot(tmin, beta, "o-", color="tab:green")
        axes[2].set_ylabel("beta")
        axes[2].set_yscale("log")
        axes[2].grid(alpha=0.3)

        axes[3].plot(tmin, dr, "o-", color="tab:purple")
        axes[3].set_ylabel("|r_sample-r_fixed| (Rs)")
        axes[3].set_xlabel("time (minutes)")
        axes[3].grid(alpha=0.3)

        fig.suptitle(f"Line {lid}: fixed point r={r_fix:.3f}Rs, theta={th:.4f}, phi={ph:.4f}")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        pp = Path(str(out_prefix) + f"_line{lid:02d}_fixedpoint.png")
        fig.savefig(pp, dpi=180)
        plt.close(fig)
        out_paths.append(pp)

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

    select_ids = _parse_line_ids(args.line_ids)
    line_csv = case_dir / args.line_source_csv
    line_defs = _load_line_defs(line_csv, select_ids=select_ids, r_mode_round=args.r_mode_round)

    rows: List[Dict[str, float]] = []
    unit_a = args.unit_v_kms / args.unit_t_sec
    nlines = len(line_defs)
    npaths = len(paths)
    for li, ld in enumerate(line_defs, start=1):
        print(
            f"LINE {li}/{nlines}: id={int(ld['line_id'])}, theta={ld['theta_rad']:.6f}, phi={ld['phi_rad']:.6f}, "
            f"r_fixed(mode)={ld['r_fixed_Rs']:.3f}Rs, mode_count={int(ld['r_mode_count'])}"
        )
        for i, p in enumerate(paths, start=1):
            r = _probe_fixed_point(
                dat_path=p,
                cfg=cfg,
                base_edges=base_edges,
                line_def=ld,
                gravity_g0=args.gravity_g0,
                sradius=args.sradius,
            )
            r["v_km_s"] = r["v_km_s"] * args.unit_v_kms
            r["vr_km_s"] = r["vr_km_s"] * args.unit_v_kms
            r["a_pres_km_s2"] = r["a_pres_code"] * unit_a
            r["a_jxb_km_s2"] = r["a_jxb_code"] * unit_a
            r["a_grav_km_s2"] = r["a_grav_code"] * unit_a
            r["a_res_km_s2"] = r["a_res_code"] * unit_a
            r["adv_r_km_s2"] = r["adv_r_code"] * unit_a
            r["rhs_dvrdt_km_s2"] = r["rhs_dvrdt_code"] * unit_a
            rows.append(r)
            print(
                f"[line {int(ld['line_id'])} | {i}/{npaths}] {p.stem} "
                f"r_sample={r['r_sample_Rs']:.3f} (target={r['r_fixed_Rs']:.3f}) "
                f"v={r['v_km_s']:.1f} a_res={r['a_res_km_s2']:.3e} adv={r['adv_r_km_s2']:.3e} "
                f"rhs={r['rhs_dvrdt_km_s2']:.3e} beta={r['beta']:.3e}"
            )

    rows.sort(key=lambda x: (int(round(float(x["line_id"]))), int(round(float(x["idx"])))))
    _append_time_derivatives(rows, unit_t_sec=args.unit_t_sec)

    out_prefix = case_dir / args.output_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_csv = Path(str(out_prefix) + "_timeseries.csv")
    write_csv(out_csv, rows)
    print(f"Wrote {out_csv}")

    if args.make_plot:
        outs = maybe_plot(rows, out_prefix=out_prefix, unit_t_sec=args.unit_t_sec)
        for pp in outs:
            print(f"Wrote {pp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
