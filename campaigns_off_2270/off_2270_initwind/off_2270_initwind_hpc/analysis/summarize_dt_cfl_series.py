#!/usr/bin/env python3
"""Summarize DT/CFL per-frame diagnostics into evolution tables and markdown report."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np

from dt_cfl_common import parse_log, safe_float, write_csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build evolution summary and markdown report from per-frame DT/CFL results.")
    parser.add_argument("--case-dir", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--input", default="analysis/dt_cfl_hpc/dt_cfl_per_frame.csv")
    parser.add_argument("--log", default="data/amr_probe.log")
    parser.add_argument("--output", default="analysis/dt_cfl_hpc")
    return parser.parse_args()


def load_rows(path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(dict(row))
    if not rows:
        raise RuntimeError(f"No rows found in {path}")
    return rows


def to_bool(raw: object) -> bool:
    s = str(raw).strip().lower()
    return s in {"1", "true", "t", "yes", "y"}


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    input_path = (case_dir / args.input).resolve()
    out_dir = (case_dir / args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(input_path)
    log_by_it = parse_log(case_dir / args.log)

    normalized: List[Dict[str, object]] = []
    for row in rows:
        it = int(float(row.get("it", 0)))
        log_row = log_by_it.get(it, {})
        dt_log = safe_float(row.get("dt_log"))
        if (not math.isfinite(dt_log)) and ("dt" in log_row):
            dt_log = float(log_row["dt"])

        n1 = safe_float(row.get("n1_log"))
        n2 = safe_float(row.get("n2_log"))
        if not math.isfinite(n1) and ("n1" in log_row):
            n1 = float(log_row["n1"])
        if not math.isfinite(n2) and ("n2" in log_row):
            n2 = float(log_row["n2"])

        dt_pred = safe_float(row.get("dt_pred_min"))
        if not math.isfinite(dt_pred):
            dt_pred = safe_float(row.get("dt_pred"))

        c_r = safe_float(row.get("C_r"))
        c_th = safe_float(row.get("C_th"))
        c_ph = safe_float(row.get("C_ph"))
        c_sum = safe_float(row.get("C_sum"))

        frac_r = c_r / c_sum if (math.isfinite(c_r) and math.isfinite(c_sum) and c_sum > 0.0) else math.nan
        frac_th = c_th / c_sum if (math.isfinite(c_th) and math.isfinite(c_sum) and c_sum > 0.0) else math.nan
        frac_ph = c_ph / c_sum if (math.isfinite(c_ph) and math.isfinite(c_sum) and c_sum > 0.0) else math.nan

        v_abs = safe_float(row.get("v_abs"))
        cf_max = safe_float(row.get("cf_max"))
        cf_over_v = cf_max / v_abs if (math.isfinite(cf_max) and math.isfinite(v_abs) and v_abs > 0.0) else math.nan

        theta_deg = safe_float(row.get("theta_deg"))
        in_pole_guard = bool(theta_deg <= 1.0 or theta_deg >= 179.0) if math.isfinite(theta_deg) else to_bool(row.get("in_pole_guard_1deg", ""))

        normalized.append(
            {
                "frame": row.get("frame", ""),
                "it": it,
                "time": safe_float(row.get("time")),
                "dt_log": dt_log,
                "dt_pred_min": dt_pred,
                "dt_ratio": (dt_pred / dt_log) if (math.isfinite(dt_pred) and math.isfinite(dt_log) and dt_log > 0.0) else math.nan,
                "r": safe_float(row.get("r")),
                "theta_deg": theta_deg,
                "phi_deg": safe_float(row.get("phi_deg")),
                "level": int(float(row.get("level", 0))),
                "block_id": int(float(row.get("block_id", 0))),
                "dominant_direction": row.get("dominant_direction", ""),
                "dominant_physics": row.get("dominant_physics", ""),
                "C_r": c_r,
                "C_th": c_th,
                "C_ph": c_ph,
                "C_sum": c_sum,
                "frac_C_r": frac_r,
                "frac_C_th": frac_th,
                "frac_C_ph": frac_ph,
                "v_abs": v_abs,
                "c_s": safe_float(row.get("c_s")),
                "v_A": safe_float(row.get("v_A")),
                "cf_max": cf_max,
                "cf_over_v": cf_over_v,
                "n1": n1,
                "n2": n2,
                "nleafs": safe_float(row.get("nleafs")),
                "in_pole_guard_1deg": in_pole_guard,
                "angle_mode_used": row.get("angle_mode_used", ""),
            }
        )

    normalized.sort(key=lambda r: int(r["it"]))

    write_csv(out_dir / "dt_cfl_evolution_summary.csv", normalized)

    finite_dt_log = np.array([float(r["dt_log"]) for r in normalized if math.isfinite(float(r["dt_log"])) and float(r["dt_log"]) > 0.0])
    finite_dt_pred = np.array([float(r["dt_pred_min"]) for r in normalized if math.isfinite(float(r["dt_pred_min"])) and float(r["dt_pred_min"]) > 0.0])
    finite_ratio = np.array([float(r["dt_ratio"]) for r in normalized if math.isfinite(float(r["dt_ratio"]))])
    finite_cfov = np.array([float(r["cf_over_v"]) for r in normalized if math.isfinite(float(r["cf_over_v"]))])

    dirs = Counter(str(r["dominant_direction"]) for r in normalized)
    phys = Counter(str(r["dominant_physics"]) for r in normalized)

    frac_ph_vals = np.array([float(r["frac_C_ph"]) for r in normalized if math.isfinite(float(r["frac_C_ph"]))])
    frac_r_vals = np.array([float(r["frac_C_r"]) for r in normalized if math.isfinite(float(r["frac_C_r"]))])
    frac_th_vals = np.array([float(r["frac_C_th"]) for r in normalized if math.isfinite(float(r["frac_C_th"]))])

    theta_vals = np.array([float(r["theta_deg"]) for r in normalized if math.isfinite(float(r["theta_deg"]))])
    r_vals = np.array([float(r["r"]) for r in normalized if math.isfinite(float(r["r"]))])

    n2_vals = np.array([float(r["n2"]) for r in normalized if math.isfinite(float(r["n2"]))])
    dt_vals_for_corr = np.array([float(r["dt_log"]) for r in normalized if math.isfinite(float(r["n2"])) and math.isfinite(float(r["dt_log"]))])
    n2_vals_for_corr = np.array([float(r["n2"]) for r in normalized if math.isfinite(float(r["n2"])) and math.isfinite(float(r["dt_log"]))])
    corr_n2_dt = math.nan
    if n2_vals_for_corr.size >= 3:
        corr_n2_dt = float(np.corrcoef(n2_vals_for_corr, dt_vals_for_corr)[0, 1])

    pole_hits = sum(1 for r in normalized if bool(r["in_pole_guard_1deg"]))

    top_small_dt = sorted(
        [r for r in normalized if math.isfinite(float(r["dt_pred_min"]))],
        key=lambda r: float(r["dt_pred_min"]),
    )[:5]

    today = datetime.now().strftime("%Y-%m-%d")
    md_path = out_dir / f"dt_cfl_audit_summary_{today}.md"

    lines: List[str] = []
    lines.append(f"# DT/CFL 审查总结 ({today})")
    lines.append("")
    lines.append(f"输入: `{input_path}`")
    lines.append(f"样本帧数: `{len(normalized)}`")
    lines.append("")

    lines.append("## 最小 dt 空间轨迹")
    if r_vals.size > 0 and theta_vals.size > 0:
        lines.append(f"- r 范围: `{np.min(r_vals):.6f} -> {np.max(r_vals):.6f}`")
        lines.append(f"- theta 范围(deg): `{np.min(theta_vals):.6f} -> {np.max(theta_vals):.6f}`")
    if normalized:
        first = normalized[0]
        last = normalized[-1]
        lines.append(
            f"- 起点: `{first['frame']}` it={first['it']} r={float(first['r']):.6f} theta={float(first['theta_deg']):.6f}deg"
        )
        lines.append(
            f"- 终点: `{last['frame']}` it={last['it']} r={float(last['r']):.6f} theta={float(last['theta_deg']):.6f}deg"
        )
    lines.append("")

    lines.append("## CFL 主导项演化")
    if frac_ph_vals.size > 0:
        lines.append(f"- 平均 `C_ph/C_sum`: `{np.mean(frac_ph_vals):.4f}`")
    if frac_r_vals.size > 0:
        lines.append(f"- 平均 `C_r/C_sum`: `{np.mean(frac_r_vals):.4f}`")
    if frac_th_vals.size > 0:
        lines.append(f"- 平均 `C_th/C_sum`: `{np.mean(frac_th_vals):.4f}`")
    lines.append(f"- `dominant_direction` 计数: `{dict(dirs)}`")
    lines.append(f"- `dominant_physics` 计数: `{dict(phys)}`")
    lines.append("")

    lines.append("## cf_max 与 |v| 关系")
    if finite_cfov.size > 0:
        lines.append(
            f"- `cf_max/|v|`: p50={np.percentile(finite_cfov,50):.4f}, p90={np.percentile(finite_cfov,90):.4f}, max={np.max(finite_cfov):.4f}"
        )
    else:
        lines.append("- 无有效 `cf_max/|v|` 样本")
    lines.append("")

    lines.append("## 与 AMR(n1/n2) 的对应关系")
    if n2_vals.size > 0:
        lines.append(f"- n2 范围: `{np.min(n2_vals):.0f} -> {np.max(n2_vals):.0f}`")
    if math.isfinite(corr_n2_dt):
        lines.append(f"- corr(n2, dt_log): `{corr_n2_dt:.4f}`")
    lines.append(f"- 最小dt位置落入 1deg 极区保护带的帧数: `{pole_hits}/{len(normalized)}`")
    lines.append("")

    lines.append("## 时间步与一致性")
    if finite_dt_log.size > 0:
        lines.append(
            f"- dt_log: min={np.min(finite_dt_log):.6e}, p50={np.percentile(finite_dt_log,50):.6e}, max={np.max(finite_dt_log):.6e}"
        )
    if finite_dt_pred.size > 0:
        lines.append(
            f"- dt_pred_min: min={np.min(finite_dt_pred):.6e}, p50={np.percentile(finite_dt_pred,50):.6e}, max={np.max(finite_dt_pred):.6e}"
        )
    if finite_ratio.size > 0:
        lines.append(
            f"- dt_ratio(dt_pred/dt_log): p50={np.percentile(finite_ratio,50):.4f}, p90={np.percentile(finite_ratio,90):.4f}"
        )
    lines.append("")

    lines.append("## 最小 dt 前5帧")
    for row in top_small_dt:
        lines.append(
            f"- `{row['frame']}` it={row['it']} dt_pred_min={float(row['dt_pred_min']):.6e} "
            f"(r={float(row['r']):.6f}, theta={float(row['theta_deg']):.6f}deg, dir={row['dominant_direction']}, phys={row['dominant_physics']})"
        )
    lines.append("")

    md_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {out_dir / 'dt_cfl_evolution_summary.csv'}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
