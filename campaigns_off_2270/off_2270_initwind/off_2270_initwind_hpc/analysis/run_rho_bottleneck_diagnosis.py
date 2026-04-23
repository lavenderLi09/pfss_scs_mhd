#!/usr/bin/env python3
"""Run a minimal rho bottleneck diagnosis pipeline from DT/CFL CSV inputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from dt_cfl_common import (
    EPS,
    block_radial_edges,
    block_uniform_edges,
    build_base_radial_edges,
    get_header,
    get_tree_info,
    read_block_fields_selected,
    read_par_config,
)
from rho_bottleneck_diagnostics import (
    classify_final_label,
    evaluate_numerical_gate,
    physical_score,
)

PER_FRAME_COLUMNS = [
    "frame",
    "it",
    "time",
    "dt_pred_min",
    "rho",
    "p",
    "v_abs",
    "b_abs",
    "v_A",
    "rho_ratio_to_initial",
    "shell_percentile",
    "flux_jump",
    "persistence_frames",
    "score",
    "S1",
    "S2",
    "S3",
    "S4",
    "gate_failed",
    "final_label",
]

GATE_REPORT_COLUMNS = [
    "frame",
    "it",
    "time",
    "neg_rho_or_p",
    "nan_inf",
    "divb_flag",
    "flux_residual_flag",
    "gate_failed",
    "divb_median",
    "divb_p999",
]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run rho bottleneck diagnosis on per-frame DT/CFL outputs.")
    parser.add_argument("--case-dir", required=True, help="Case directory containing the analysis inputs.")
    parser.add_argument("--per-frame", required=True, help="Path to the per-frame CSV file.")
    parser.add_argument("--topcells", required=True, help="Path to the top-cells CSV file.")
    parser.add_argument("--output-dir", required=True, help="Directory for diagnosis outputs.")
    parser.add_argument("--frames", default="", help="Optional comma-separated frame filter.")
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation.")
    return parser


def output_paths(out_dir: Path) -> dict[str, Path]:
    return {
        "per_frame": out_dir / "rho_bottleneck_diagnostic_per_frame.csv",
        "gate": out_dir / "rho_bottleneck_gate_report.csv",
        "summary": out_dir / "rho_bottleneck_diagnosis_summary.md",
        "plot": out_dir / "rho_bottleneck_trend.png",
    }


def _resolve_existing_path(case_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path.resolve()
    case_candidate = (case_dir / path).resolve()
    if case_candidate.exists():
        return case_candidate
    return path.resolve()


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _parse_frame_filter(raw_frames: str) -> list[str] | None:
    frames = [frame.strip() for frame in raw_frames.split(",") if frame.strip()]
    return frames or None


def _to_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    raw_value = row.get(key, "")
    if raw_value in ("", None):
        return default
    return float(raw_value)


def _to_int(row: dict[str, str], key: str, default: int = 0) -> int:
    raw_value = row.get(key, "")
    if raw_value in ("", None):
        return default
    return int(float(raw_value))


def _optional_float(row: dict[str, str], *keys: str) -> float:
    for key in keys:
        raw_value = row.get(key, "")
        if raw_value in ("", None):
            continue
        return float(raw_value)
    return float("nan")


def _relative_jump(current_value: float, previous_value: float | None) -> float:
    if previous_value is None:
        return 0.0
    if previous_value == 0.0:
        return 0.0 if current_value == 0.0 else 1.0
    return abs(current_value - previous_value) / abs(previous_value)


def _shell_percentile(shell_index: int, max_shell_index: int) -> float:
    if max_shell_index <= 0:
        return 0.0
    return 100.0 * float(shell_index) / float(max_shell_index)


def _select_rank1_topcells(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        if _to_int(row, "rank", default=-1) != 1:
            continue
        frame = row["frame"]
        if frame not in selected:
            selected[frame] = row
    return selected


def _merged_timeline_rows(
    per_frame_rows: list[dict[str, str]],
    rank1_topcells: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    merged_rows: list[dict[str, str]] = []
    missing_frames: list[str] = []
    for row in per_frame_rows:
        frame = row["frame"]
        topcell = rank1_topcells.get(frame)
        if topcell is None:
            missing_frames.append(frame)
            continue
        merged = dict(row)
        merged.update(topcell)
        merged_rows.append(merged)
    if missing_frames:
        missing = ", ".join(missing_frames)
        raise ValueError(f"Missing rank==1 topcell rows for frames: {missing}")
    return merged_rows


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_summary(path: Path, per_frame_rows: list[dict[str, Any]]) -> None:
    label_counts = Counter(row["final_label"] for row in per_frame_rows)
    latest = per_frame_rows[-1]
    lines = [
        "# Rho Bottleneck Diagnosis Summary",
        "",
        f"- Frames analyzed: {len(per_frame_rows)}",
        "- Assumptions:",
        "  - flux_jump is the frame-to-frame relative change in rank-1 b_abs over the full merged timeline.",
        "  - shell_percentile is 100 * shell_index / max(shell_index) over the full merged timeline.",
        "  - divB and flux-residual are computed from DAT when available.",
        "  - if DAT diagnostics are unavailable for a frame, non-finite values are passed to the numerical gate as a conservative safe default.",
        "",
        "## Final Label Distribution",
    ]
    for label in sorted(label_counts):
        lines.append(f"- {label}: {label_counts[label]}")
    lines.extend(
        [
            "",
            "## Latest Frame",
            f"- frame={latest['frame']} it={latest['it']} time={latest['time']}",
            f"- Latest frame label: {latest['final_label']}",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def _safe_dat_diagnostics(
    dat_path: Path,
    qstretch: float,
    bottleneck_row: dict[str, str],
) -> tuple[float, float]:
    if not dat_path.exists():
        return float("nan"), float("nan")

    try:
        block_id = _to_int(bottleneck_row, "block_id", default=-1)
        i_cell = _to_int(bottleneck_row, "i", default=-1)
        j_cell = _to_int(bottleneck_row, "j", default=-1)
        k_cell = _to_int(bottleneck_row, "k", default=-1)
        r0 = _to_float(bottleneck_row, "r", default=math.nan)
    except Exception:
        return float("nan"), float("nan")

    if not math.isfinite(r0):
        return float("nan"), float("nan")

    try:
        with dat_path.open("rb") as handle:
            header = get_header(handle)
            block_lvls, block_ixs, block_offsets = get_tree_info(handle)

            name_to_idx = {name.strip(): idx for idx, name in enumerate(header["w_names"])}
            required = ["rho", "m1", "b1", "b2", "b3"]
            if any(name not in name_to_idx for name in required):
                return float("nan"), float("nan")
            field_map = {name: name_to_idx[name] for name in required}

            block_shape = header["block_nx"].astype(int)
            domain_nx = header["domain_nx"].astype(int)
            ndim = int(header["ndim"])
            xmins = [float(v) for v in np.asarray(header["xmin"], dtype=np.float64)]
            xmaxs = [float(v) for v in np.asarray(header["xmax"], dtype=np.float64)]
            base_edges = build_base_radial_edges(
                xmins[0],
                xmaxs[0],
                int(domain_nx[0]),
                float(qstretch),
            )

            targets = [0.98 * r0, r0, 1.02 * r0]
            sum_flux = np.zeros(3, dtype=np.float64)
            sum_w = np.zeros(3, dtype=np.float64)
            divb_value = float("nan")

            for idx, (level, bix, offset) in enumerate(
                zip(block_lvls.astype(int), block_ixs.astype(int), block_offsets.astype(int)),
                start=1,
            ):
                fields = read_block_fields_selected(handle, int(offset), block_shape, ndim, field_map)
                rho = fields["rho"]
                m1 = fields["m1"]
                b_r = fields["b1"]
                b_th = fields["b2"]
                b_ph = fields["b3"]

                r_edges = block_radial_edges(base_edges, int(level), int(bix[0]), int(block_shape[0]))
                th_edges = block_uniform_edges(
                    xmins[1],
                    xmaxs[1],
                    int(domain_nx[1]),
                    int(level),
                    int(bix[1]),
                    int(block_shape[1]),
                )
                ph_edges = block_uniform_edges(
                    xmins[2],
                    xmaxs[2],
                    int(domain_nx[2]),
                    int(level),
                    int(bix[2]),
                    int(block_shape[2]),
                )

                r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
                th_centers = 0.5 * (th_edges[:-1] + th_edges[1:])
                ph_centers = 0.5 * (ph_edges[:-1] + ph_edges[1:])
                dr = np.diff(r_edges)
                dth = np.diff(th_edges)
                dph = np.diff(ph_edges)

                rr, tt, _ = np.meshgrid(r_centers, th_centers, ph_centers, indexing="ij")
                rho_safe = np.maximum(rho, EPS)
                vr = m1 / rho_safe
                flux = rho * vr * rr * rr
                weight = np.maximum(np.sin(tt), EPS) * dth[None, :, None] * dph[None, None, :]
                dr3 = dr[:, None, None]

                for t_idx, target in enumerate(targets):
                    mask = np.abs(rr - target) <= 0.5 * dr3
                    if np.any(mask):
                        sum_flux[t_idx] += float(np.sum(flux[mask] * weight[mask]))
                        sum_w[t_idx] += float(np.sum(weight[mask]))

                if idx == block_id and i_cell >= 0 and j_cell >= 0 and k_cell >= 0:
                    ni, nj, nk = b_r.shape
                    if i_cell < ni and j_cell < nj and k_cell < nk and ni >= 3 and nj >= 3 and nk >= 3:
                        rr_b, tt_b, _ = np.meshgrid(r_centers, th_centers, ph_centers, indexing="ij")
                        sin_t = np.maximum(np.sin(tt_b), EPS)
                        d_r_term = np.gradient(rr_b * rr_b * b_r, r_centers, axis=0, edge_order=2)
                        d_th_term = np.gradient(sin_t * b_th, th_centers, axis=1, edge_order=2)
                        d_ph_term = np.gradient(b_ph, ph_centers, axis=2, edge_order=2)
                        div_b = (
                            d_r_term / np.maximum(rr_b * rr_b, EPS)
                            + d_th_term / np.maximum(rr_b * sin_t, EPS)
                            + d_ph_term / np.maximum(rr_b * sin_t, EPS)
                        )
                        val = float(div_b[i_cell, j_cell, k_cell])
                        divb_value = abs(val) if math.isfinite(val) else float("nan")

            F = np.full(3, np.nan, dtype=np.float64)
            valid = sum_w > 0.0
            F[valid] = sum_flux[valid] / sum_w[valid]
            if np.all(np.isfinite(F)):
                denom = max(abs(float(F[1])), EPS)
                flux_residual = max(abs(float(F[2] - F[1])), abs(float(F[1] - F[0]))) / denom
            else:
                flux_residual = float("nan")
            return float(divb_value), float(flux_residual)
    except Exception:
        return float("nan"), float("nan")


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    case_dir = Path(args.case_dir).resolve()
    per_frame_path = _resolve_existing_path(case_dir, args.per_frame)
    topcells_path = _resolve_existing_path(case_dir, args.topcells)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        qstretch = float(read_par_config(case_dir / "amrvac.par")["qstretch_baselevel"])
    except Exception:
        qstretch = float("nan")

    frame_filter = _parse_frame_filter(args.frames)
    selected_frames = set(frame_filter) if frame_filter else None

    per_frame_rows_all = _load_csv_rows(per_frame_path)
    topcell_rows_all = _load_csv_rows(topcells_path)
    rank1_topcells = _select_rank1_topcells(topcell_rows_all)

    try:
        merged_rows = _merged_timeline_rows(per_frame_rows_all, rank1_topcells)
    except ValueError as exc:
        parser.error(str(exc))

    if not merged_rows:
        parser.error("No merged per-frame rows were available for diagnosis.")

    initial_rho = _to_float(merged_rows[0], "rho")
    max_shell_index = max(_to_int(row, "shell_index") for row in merged_rows)
    previous_b_abs: float | None = None
    low_rho_streak = 0

    all_per_frame_out: list[dict[str, Any]] = []
    all_gate_report_out: list[dict[str, Any]] = []

    for row in merged_rows:
        rho = _to_float(row, "rho")
        p = _to_float(row, "p")
        v_abs = _to_float(row, "v_abs")
        b_abs = _to_float(row, "b_abs")
        v_a = _to_float(row, "v_A")
        shell_index = _to_int(row, "shell_index")
        divb_value = _optional_float(row, "divb", "divB")
        flux_residual = _optional_float(row, "flux_residual")
        if not math.isfinite(divb_value) or not math.isfinite(flux_residual):
            dat_path = case_dir / "data" / f"{row['frame']}.dat"
            if math.isfinite(qstretch):
                dat_divb, dat_flux_residual = _safe_dat_diagnostics(dat_path, qstretch, row)
                if math.isfinite(dat_divb):
                    divb_value = dat_divb
                if math.isfinite(dat_flux_residual):
                    flux_residual = dat_flux_residual

        rho_ratio_to_initial = 0.0 if initial_rho == 0.0 else rho / initial_rho
        shell_percentile = _shell_percentile(shell_index, max_shell_index)
        flux_jump = _relative_jump(b_abs, previous_b_abs)

        if rho_ratio_to_initial < 0.4:
            low_rho_streak += 1
        else:
            low_rho_streak = 0

        gate_result = evaluate_numerical_gate(
            rho=[rho],
            p=[p],
            divb=[divb_value],
            flux_residual=flux_residual,
        )
        score_result = physical_score(
            rho_ratio_to_initial=rho_ratio_to_initial,
            shell_percentile=shell_percentile,
            flux_jump=flux_jump,
            persistence_frames=low_rho_streak,
        )
        final_label = classify_final_label(gate_result["gate_failed"], score_result["score"])

        all_per_frame_out.append(
            {
                "frame": row["frame"],
                "it": row["it"],
                "time": row["time"],
                "dt_pred_min": row["dt_pred_min"],
                "rho": rho,
                "p": p,
                "v_abs": v_abs,
                "b_abs": b_abs,
                "v_A": v_a,
                "rho_ratio_to_initial": rho_ratio_to_initial,
                "shell_percentile": shell_percentile,
                "flux_jump": flux_jump,
                "persistence_frames": low_rho_streak,
                "score": score_result["score"],
                "S1": score_result["S1"],
                "S2": score_result["S2"],
                "S3": score_result["S3"],
                "S4": score_result["S4"],
                "gate_failed": gate_result["gate_failed"],
                "final_label": final_label,
            }
        )
        all_gate_report_out.append(
            {
                "frame": row["frame"],
                "it": row["it"],
                "time": row["time"],
                "neg_rho_or_p": gate_result["neg_rho_or_p"],
                "nan_inf": gate_result["nan_inf"],
                "divb_flag": gate_result["divb_flag"],
                "flux_residual_flag": gate_result["flux_residual_flag"],
                "gate_failed": gate_result["gate_failed"],
                "divb_median": gate_result["divb_median"],
                "divb_p999": gate_result["divb_p999"],
            }
        )

        previous_b_abs = b_abs if math.isfinite(b_abs) else previous_b_abs

    if selected_frames is None:
        per_frame_out = all_per_frame_out
        gate_report_out = all_gate_report_out
    else:
        per_frame_out = [row for row in all_per_frame_out if row["frame"] in selected_frames]
        gate_report_out = [row for row in all_gate_report_out if row["frame"] in selected_frames]

    if not per_frame_out:
        parser.error("No per-frame rows matched the requested frame selection.")

    outputs = output_paths(output_dir)
    _write_csv(outputs["per_frame"], PER_FRAME_COLUMNS, per_frame_out)
    _write_csv(outputs["gate"], GATE_REPORT_COLUMNS, gate_report_out)
    _write_summary(outputs["summary"], per_frame_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
