#!/usr/bin/env python3
"""Per-snapshot DT/CFL bottleneck analysis for off_2270_initwind_hpc-like cases."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from dt_cfl_common import (
    EPS,
    block_radial_edges,
    block_uniform_edges,
    build_base_radial_edges,
    classify_bottleneck,
    conservative_to_primitive,
    get_angle_scales,
    get_header,
    get_tree_info,
    parse_log,
    read_block_fields_selected,
    read_par_config,
    write_csv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze minimum dt location and CFL decomposition per DAT frame.")
    parser.add_argument("--case-dir", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--pattern", default="data/amr_probe*.dat")
    parser.add_argument("--log", default="data/amr_probe.log")
    parser.add_argument("--par", default="amrvac.par")
    parser.add_argument("--output-dir", default="analysis/dt_cfl_hpc")
    parser.add_argument("--frames", default="", help="Comma-separated stems, e.g. amr_probe0000,amr_probe0001")
    parser.add_argument("--topn", type=int, default=10)
    parser.add_argument("--angle-mode", choices=["auto", "normalized_2pi", "radian"], default="auto")
    parser.add_argument("--no-plots", action="store_true", help="Kept for interface compatibility; no plots are produced.")
    return parser.parse_args()


def select_frames(paths: Sequence[Path], frames_arg: str) -> List[Path]:
    if not frames_arg.strip():
        return list(paths)
    wanted = {token.strip() for token in frames_arg.split(",") if token.strip()}
    return [p for p in paths if p.stem in wanted]


def extract_gamma(header: Dict[str, object], fallback_gamma: float) -> float:
    names = list(header.get("param_names", []))
    params = np.array(header.get("params", []), dtype=np.float64)
    if "gamma" in names:
        return float(params[names.index("gamma")])
    return float(fallback_gamma)


def analyze_snapshot(
    dat_path: Path,
    par_cfg: Dict[str, float],
    log_by_it: Dict[int, Dict[str, float]],
    topn: int,
    angle_mode: str,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    with dat_path.open("rb") as handle:
        header = get_header(handle)
        block_lvls, block_ixs, block_offsets = get_tree_info(handle)

        name_to_idx = {name.strip(): idx for idx, name in enumerate(header["w_names"]) }
        required = ["rho", "m1", "m2", "m3", "e", "b1", "b2", "b3"]
        missing = [name for name in required if name not in name_to_idx]
        if missing:
            raise KeyError(f"{dat_path.name}: missing fields {missing}")
        field_map = {name: name_to_idx[name] for name in required}

        # Use coordinates from the DAT header, i.e. the *actual* simulation coordinates.
        # In spherical geometry, AMRVAC may transform par-domain angles before runtime.
        coord_cfg = {
            "xprobmax2": float(header["xmax"][1]),
            "xprobmax3": float(header["xmax"][2]),
        }
        mode_used, theta_scale, phi_scale = get_angle_scales(coord_cfg, angle_mode)
        gamma = extract_gamma(header, par_cfg["mhd_gamma"])
        it = int(header["it"])
        log_row = log_by_it.get(it, {})
        dt_log = float(log_row.get("dt", math.nan))
        time_log = float(log_row.get("time", math.nan))

        block_shape = header["block_nx"].astype(int)
        domain_nx = header["domain_nx"].astype(int)
        ndim = int(header["ndim"])
        xmins = [float(v) for v in np.asarray(header["xmin"], dtype=np.float64)]
        xmaxs = [float(v) for v in np.asarray(header["xmax"], dtype=np.float64)]
        courant = float(par_cfg["courantpar"])
        base_edges = build_base_radial_edges(
            xmins[0],
            xmaxs[0],
            int(domain_nx[0]),
            float(par_cfg["qstretch_baselevel"]),
        )

        top_rows: List[Dict[str, object]] = []

        for block_id, (level, bix, offset) in enumerate(zip(block_lvls.astype(int), block_ixs.astype(int), block_offsets.astype(int)), start=1):
            fields = read_block_fields_selected(handle, int(offset), block_shape, ndim, field_map)

            r_edges = block_radial_edges(base_edges, int(level), int(bix[0]), int(block_shape[0]))
            th_edges = block_uniform_edges(xmins[1], xmaxs[1], int(domain_nx[1]), int(level), int(bix[1]), int(block_shape[1]))
            ph_edges = block_uniform_edges(xmins[2], xmaxs[2], int(domain_nx[2]), int(level), int(bix[2]), int(block_shape[2]))

            r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
            th_centers = 0.5 * (th_edges[:-1] + th_edges[1:]) * theta_scale
            ph_centers = 0.5 * (ph_edges[:-1] + ph_edges[1:]) * phi_scale

            dr = np.diff(r_edges)
            dth_phys = np.diff(th_edges) * theta_scale
            dph_phys = np.diff(ph_edges) * phi_scale

            rr, tt, pp = np.meshgrid(r_centers, th_centers, ph_centers, indexing="ij")
            dr3 = dr[:, None, None]
            ds_th = rr * dth_phys[None, :, None]
            ds_ph = rr * np.maximum(np.sin(tt), EPS) * dph_phys[None, None, :]

            rho = fields["rho"]
            m1 = fields["m1"]
            m2 = fields["m2"]
            m3 = fields["m3"]
            eint = fields["e"]
            b_r = fields["b1"]
            b_th = fields["b2"]
            b_ph = fields["b3"]

            v_r, v_th, v_ph, pth = conservative_to_primitive(rho, m1, m2, m3, eint, gamma)
            valid = (rho > 0.0) & (pth > 0.0)
            if not np.any(valid):
                continue

            rho_safe = np.maximum(rho, EPS)
            b_sq = b_r * b_r + b_th * b_th + b_ph * b_ph
            c_s_sq = np.where(valid, gamma * pth / rho_safe, np.nan)
            v_a_sq = np.where(valid, b_sq / rho_safe, np.nan)
            c_s = np.sqrt(np.maximum(c_s_sq, 0.0))
            v_a = np.sqrt(np.maximum(v_a_sq, 0.0))

            def fast_speed(normal_b: np.ndarray) -> np.ndarray:
                b_n_sq = normal_b * normal_b
                disc = (c_s_sq + v_a_sq) ** 2 - 4.0 * c_s_sq * b_n_sq / rho_safe
                disc = np.maximum(disc, 0.0)
                cf_sq = 0.5 * ((c_s_sq + v_a_sq) + np.sqrt(disc))
                return np.sqrt(np.maximum(cf_sq, 0.0))

            cf_r = fast_speed(b_r)
            cf_th = fast_speed(b_th)
            cf_ph = fast_speed(b_ph)
            cf_max = np.maximum(np.maximum(cf_r, cf_th), cf_ph)

            c_r = np.where(valid, (np.abs(v_r) + cf_r) / np.maximum(dr3, EPS), np.nan)
            c_th = np.where(valid, (np.abs(v_th) + cf_th) / np.maximum(ds_th, EPS), np.nan)
            c_ph = np.where(valid, (np.abs(v_ph) + cf_ph) / np.maximum(ds_ph, EPS), np.nan)
            c_sum = c_r + c_th + c_ph

            finite = np.isfinite(c_sum) & (c_sum > 0.0)
            if not np.any(finite):
                continue

            dt_pred = np.where(finite, courant / c_sum, np.nan)
            flat_candidates = np.flatnonzero(finite)
            if flat_candidates.size == 0:
                continue
            keep = min(topn, int(flat_candidates.size))
            cand_vals = dt_pred.ravel()[flat_candidates]
            local_pick = np.argpartition(cand_vals, keep - 1)[:keep]
            chosen = flat_candidates[local_pick]
            chosen = chosen[np.argsort(dt_pred.ravel()[chosen])]

            refine_ratio = 2 ** (int(level) - 1)
            start_r = (int(bix[0]) - 1) * int(block_shape[0])
            start_t = (int(bix[1]) - 1) * int(block_shape[1])
            start_p = (int(bix[2]) - 1) * int(block_shape[2])
            base_r_idx = (np.arange(start_r, start_r + int(block_shape[0])) // refine_ratio).astype(np.int64)
            base_t_idx = (np.arange(start_t, start_t + int(block_shape[1])) // refine_ratio).astype(np.int64)
            base_p_idx = (np.arange(start_p, start_p + int(block_shape[2])) // refine_ratio).astype(np.int64)

            for flat_idx in chosen:
                i_cell, j_cell, k_cell = np.unravel_index(int(flat_idx), c_sum.shape)
                c_rv = float(c_r[i_cell, j_cell, k_cell])
                c_thv = float(c_th[i_cell, j_cell, k_cell])
                c_phv = float(c_ph[i_cell, j_cell, k_cell])
                csumv = float(c_sum[i_cell, j_cell, k_cell])
                dtv = float(dt_pred[i_cell, j_cell, k_cell])

                vrv = float(v_r[i_cell, j_cell, k_cell])
                vtv = float(v_th[i_cell, j_cell, k_cell])
                vpv = float(v_ph[i_cell, j_cell, k_cell])
                cfrv = float(cf_r[i_cell, j_cell, k_cell])
                cfthv = float(cf_th[i_cell, j_cell, k_cell])
                cfphv = float(cf_ph[i_cell, j_cell, k_cell])
                csv = float(c_s[i_cell, j_cell, k_cell])
                vav = float(v_a[i_cell, j_cell, k_cell])
                dr_v = float(dr3[i_cell, 0, 0])
                ds_th_v = float(ds_th[i_cell, j_cell, k_cell])
                ds_ph_v = float(ds_ph[i_cell, j_cell, k_cell])

                dom_idx = int(np.argmax([c_rv, c_thv, c_phv]))
                v_n = abs([vrv, vtv, vpv][dom_idx])
                cf_n = [cfrv, cfthv, cfphv][dom_idx]
                dom_dir, dom_phys = classify_bottleneck(
                    c_rv,
                    c_thv,
                    c_phv,
                    dr_v,
                    ds_th_v,
                    ds_ph_v,
                    v_n,
                    cf_n,
                    csv,
                    vav,
                )

                rho_v = float(rho[i_cell, j_cell, k_cell])
                p_v = float(pth[i_cell, j_cell, k_cell])
                brv = float(b_r[i_cell, j_cell, k_cell])
                btv = float(b_th[i_cell, j_cell, k_cell])
                bpv = float(b_ph[i_cell, j_cell, k_cell])
                vabs = math.sqrt(vrv * vrv + vtv * vtv + vpv * vpv)
                babs = math.sqrt(brv * brv + btv * btv + bpv * bpv)

                theta_rad = float(tt[i_cell, j_cell, k_cell])
                phi_rad = float(pp[i_cell, j_cell, k_cell])
                theta_deg = theta_rad * 180.0 / math.pi
                phi_deg = phi_rad * 180.0 / math.pi

                row = {
                    "frame": dat_path.stem,
                    "it": it,
                    "time": float(header["time"]),
                    "time_log": time_log,
                    "dt_log": dt_log,
                    "dt_pred": dtv,
                    "dt_ratio": (dtv / dt_log) if (math.isfinite(dt_log) and dt_log > 0.0) else math.nan,
                    "angle_mode_used": mode_used,
                    "theta_scale": theta_scale,
                    "phi_scale": phi_scale,
                    "gamma": gamma,
                    "level": int(level),
                    "block_id": int(block_id),
                    "i": int(i_cell),
                    "j": int(j_cell),
                    "k": int(k_cell),
                    "i_global": int(start_r + i_cell),
                    "j_global": int(start_t + j_cell),
                    "k_global": int(start_p + k_cell),
                    "r": float(rr[i_cell, j_cell, k_cell]),
                    "theta": theta_rad,
                    "phi": phi_rad,
                    "theta_deg": theta_deg,
                    "phi_deg": phi_deg,
                    "rho": rho_v,
                    "p": p_v,
                    "v_abs": vabs,
                    "b_abs": babs,
                    "v_r": vrv,
                    "v_theta": vtv,
                    "v_phi": vpv,
                    "b_r": brv,
                    "b_theta": btv,
                    "b_phi": bpv,
                    "c_s": csv,
                    "v_A": vav,
                    "cf_r": cfrv,
                    "cf_th": cfthv,
                    "cf_ph": cfphv,
                    "cf_max": float(cf_max[i_cell, j_cell, k_cell]),
                    "dr": dr_v,
                    "r_dtheta": ds_th_v,
                    "r_sintheta_dphi": ds_ph_v,
                    "shell_index": int(base_r_idx[i_cell]),
                    "theta_bin": int(base_t_idx[j_cell]),
                    "phi_bin": int(base_p_idx[k_cell]),
                    "C_r": c_rv,
                    "C_th": c_thv,
                    "C_ph": c_phv,
                    "C_sum": csumv,
                    "dominant_direction": dom_dir,
                    "dominant_physics": dom_phys,
                    "nleafs": int(header["nleafs"]),
                    "n1_log": int(log_row.get("n1", math.nan)) if "n1" in log_row else "",
                    "n2_log": int(log_row.get("n2", math.nan)) if "n2" in log_row else "",
                }
                top_rows.append(row)

    if not top_rows:
        raise RuntimeError(f"No valid CFL cells found in {dat_path}")

    top_rows.sort(key=lambda r: float(r["dt_pred"]))
    top_rows = top_rows[:topn]
    for rank, row in enumerate(top_rows, start=1):
        row["rank"] = rank

    best = dict(top_rows[0])
    best["dt_pred_min"] = best["dt_pred"]
    best["dt_ratio"] = (float(best["dt_pred_min"]) / dt_log) if (math.isfinite(dt_log) and dt_log > 0.0) else math.nan
    best["in_pole_guard_1deg"] = bool(best["theta_deg"] <= 1.0 or best["theta_deg"] >= 179.0)
    summary = {
        "frame": best["frame"],
        "it": best["it"],
        "time": best["time"],
        "time_log": best["time_log"],
        "dt_log": best["dt_log"],
        "dt_pred_min": best["dt_pred_min"],
        "dt_ratio": best["dt_ratio"],
        "angle_mode_used": best["angle_mode_used"],
        "theta_scale": best["theta_scale"],
        "phi_scale": best["phi_scale"],
        "level": best["level"],
        "block_id": best["block_id"],
        "i": best["i"],
        "j": best["j"],
        "k": best["k"],
        "i_global": best["i_global"],
        "j_global": best["j_global"],
        "k_global": best["k_global"],
        "r": best["r"],
        "theta": best["theta"],
        "phi": best["phi"],
        "theta_deg": best["theta_deg"],
        "phi_deg": best["phi_deg"],
        "shell_index": best["shell_index"],
        "theta_bin": best["theta_bin"],
        "phi_bin": best["phi_bin"],
        "C_r": best["C_r"],
        "C_th": best["C_th"],
        "C_ph": best["C_ph"],
        "C_sum": best["C_sum"],
        "dominant_direction": best["dominant_direction"],
        "dominant_physics": best["dominant_physics"],
        "v_abs": best["v_abs"],
        "c_s": best["c_s"],
        "v_A": best["v_A"],
        "cf_r": best["cf_r"],
        "cf_th": best["cf_th"],
        "cf_ph": best["cf_ph"],
        "cf_max": best["cf_max"],
        "dr": best["dr"],
        "r_dtheta": best["r_dtheta"],
        "r_sintheta_dphi": best["r_sintheta_dphi"],
        "rho": best["rho"],
        "p": best["p"],
        "nleafs": best["nleafs"],
        "n1_log": best["n1_log"],
        "n2_log": best["n2_log"],
        "in_pole_guard_1deg": best["in_pole_guard_1deg"],
        "gamma": best["gamma"],
    }
    return summary, top_rows


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    output_dir = (case_dir / args.output_dir).resolve()

    par_cfg = read_par_config(case_dir / args.par)
    log_by_it = parse_log(case_dir / args.log)
    dat_paths = sorted(case_dir.glob(args.pattern))
    dat_paths = select_frames(dat_paths, args.frames)
    if not dat_paths:
        raise FileNotFoundError(f"No DAT files matched pattern '{args.pattern}' under {case_dir}")

    per_frame: List[Dict[str, object]] = []
    topcells: List[Dict[str, object]] = []

    for dat_path in dat_paths:
        summary, top_rows = analyze_snapshot(
            dat_path=dat_path,
            par_cfg=par_cfg,
            log_by_it=log_by_it,
            topn=int(args.topn),
            angle_mode=args.angle_mode,
        )
        per_frame.append(summary)
        topcells.extend(top_rows)

    per_frame.sort(key=lambda r: int(r["it"]))
    topcells.sort(key=lambda r: (str(r["frame"]), float(r["dt_pred"])))

    write_csv(output_dir / "dt_cfl_per_frame.csv", per_frame)
    write_csv(output_dir / "dt_cfl_topcells.csv", topcells)

    print(f"Wrote {output_dir / 'dt_cfl_per_frame.csv'}")
    print(f"Wrote {output_dir / 'dt_cfl_topcells.csv'}")
    print("\nframe        it      dt_log      dt_pred_min   ratio      dir     physics    r      theta_deg  level")
    for row in per_frame:
        print(
            f"{str(row['frame']):<12} {int(row['it']):>6d} "
            f"{float(row['dt_log']):>11.3e} {float(row['dt_pred_min']):>12.3e} "
            f"{float(row['dt_ratio']):>8.3f} {str(row['dominant_direction']):>8} "
            f"{str(row['dominant_physics']):>10} {float(row['r']):>7.3f} "
            f"{float(row['theta_deg']):>10.4f} {int(row['level']):>6d}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
