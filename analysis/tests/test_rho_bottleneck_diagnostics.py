# Copied to project analysis task toolbox
# Source case path: off_2270_initwind/off_2270_initwind_hpc/analysis/tests/test_rho_bottleneck_diagnostics.py
# Original file: /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_initwind/off_2270_initwind_hpc/analysis/tests/test_rho_bottleneck_diagnostics.py

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from run_rho_bottleneck_diagnosis import build_arg_parser, main, output_paths
from rho_bottleneck_diagnostics import (
    classify_final_label,
    evaluate_numerical_gate,
    physical_score,
)


@pytest.mark.parametrize(
    ("gate_failed", "score", "expected"),
    [
        (True, 60, "NUMERICAL_DOMINANT"),
        (False, 60, "PHYSICAL_HIGH_RISK"),
        (False, 30, "MIXED_UNCERTAIN"),
        (False, 10, "PHYSICAL_ACCEPTABLE"),
    ],
)
def test_classify_final_label(gate_failed: bool, score: float, expected: str) -> None:
    assert classify_final_label(gate_failed, score) == expected


def test_evaluate_numerical_gate_flags_negative_density() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([-1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([0.0, 0.0]),
        flux_residual=0.0,
    )

    assert result["gate_failed"] is True
    assert result["neg_rho_or_p"] is True


def test_evaluate_numerical_gate_flags_nonpositive_pressure() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 0.0]),
        divb=np.array([0.0, 0.0]),
        flux_residual=0.0,
    )

    assert result["neg_rho_or_p"] is True
    assert result["gate_failed"] is True


def test_evaluate_numerical_gate_flags_divb_outlier() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([1.0] * 998 + [1000.0] * 2),
        flux_residual=0.0,
    )

    assert result["divb_flag"] is True
    assert result["divb_median"] == 1.0
    assert result["divb_p999"] == 1000.0


def test_evaluate_numerical_gate_flags_zero_median_divb_outlier() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([0.0] * 998 + [1000.0] * 2),
        flux_residual=0.0,
    )

    assert result["divb_flag"] is True
    assert result["divb_median"] == 0.0
    assert result["divb_p999"] == 1000.0


def test_evaluate_numerical_gate_no_finite_divb_returns_inf_metrics() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([np.nan, np.inf]),
        flux_residual=0.0,
    )

    assert result["nan_inf"] is True
    assert result["gate_failed"] is True
    assert result["divb_median"] == np.inf
    assert result["divb_p999"] == np.inf


def test_evaluate_numerical_gate_flags_nonfinite_rho_or_p() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([np.nan, 1.0]),
        p=np.array([1.0, np.inf]),
        divb=np.array([0.0, 0.0]),
        flux_residual=0.0,
    )

    assert result["nan_inf"] is True
    assert result["gate_failed"] is True


def test_evaluate_numerical_gate_flags_flux_residual_boundary() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([0.0, 0.0]),
        flux_residual=0.30,
    )

    assert result["flux_residual_flag"] is False


def test_evaluate_numerical_gate_flags_flux_residual() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([0.0, 0.0]),
        flux_residual=0.31,
    )

    assert result["flux_residual_flag"] is True


def test_evaluate_numerical_gate_flags_nonfinite_flux_residual() -> None:
    result = evaluate_numerical_gate(
        rho=np.array([1.0, 1.0]),
        p=np.array([1.0, 1.0]),
        divb=np.array([0.0, 0.0]),
        flux_residual=np.inf,
    )

    assert result["flux_residual_flag"] is True
    assert result["gate_failed"] is True


def test_physical_score_thresholds() -> None:
    result = physical_score(
        rho_ratio_to_initial=0.1,
        shell_percentile=0.05,
        flux_jump=0.35,
        persistence_frames=6,
    )

    assert result == {
        "S1": 40,
        "S2": 30,
        "S3": 20,
        "S4": 10,
        "score": 100,
    }


@pytest.mark.parametrize(
    ("rho_ratio_to_initial", "expected_s1"),
    [
        (0.2, 20),
        (0.4, 0),
    ],
)
def test_physical_score_s1_boundaries(
    rho_ratio_to_initial: float,
    expected_s1: int,
) -> None:
    result = physical_score(
        rho_ratio_to_initial=rho_ratio_to_initial,
        shell_percentile=0.05,
        flux_jump=0.35,
        persistence_frames=6,
    )

    assert result["S1"] == expected_s1


@pytest.mark.parametrize(
    ("shell_percentile", "expected_s2"),
    [
        (0.1, 15),
        (1.0, 0),
    ],
)
def test_physical_score_s2_boundaries(
    shell_percentile: float,
    expected_s2: int,
) -> None:
    result = physical_score(
        rho_ratio_to_initial=0.1,
        shell_percentile=shell_percentile,
        flux_jump=0.35,
        persistence_frames=6,
    )

    assert result["S2"] == expected_s2


@pytest.mark.parametrize(
    ("flux_jump", "expected_s3"),
    [
        (0.30, 10),
        (0.15, 10),
        (0.14999999999999997, 0),
    ],
)
def test_physical_score_s3_boundaries(flux_jump: float, expected_s3: int) -> None:
    result = physical_score(
        rho_ratio_to_initial=0.1,
        shell_percentile=0.05,
        flux_jump=flux_jump,
        persistence_frames=6,
    )

    assert result["S3"] == expected_s3


def test_run_rho_bottleneck_diagnosis_cli_parser_smoke() -> None:
    parser = build_arg_parser()
    ns = parser.parse_args(
        [
            "--case-dir",
            ".",
            "--per-frame",
            "analysis/dt_cfl_hpc/dt_cfl_per_frame.csv",
            "--topcells",
            "analysis/dt_cfl_hpc/dt_cfl_topcells.csv",
            "--output-dir",
            "analysis/dt_cfl_hpc",
        ]
    )

    assert ns.case_dir == "."


def test_output_paths_names(tmp_path: Path) -> None:
    outputs = output_paths(tmp_path)

    assert outputs == {
        "per_frame": tmp_path / "rho_bottleneck_diagnostic_per_frame.csv",
        "gate": tmp_path / "rho_bottleneck_gate_report.csv",
        "summary": tmp_path / "rho_bottleneck_diagnosis_summary.md",
        "plot": tmp_path / "rho_bottleneck_trend.png",
    }


def test_run_rho_bottleneck_diagnosis_writes_outputs(tmp_path: Path) -> None:
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    per_frame_path = case_dir / "per_frame.csv"
    topcells_path = case_dir / "topcells.csv"
    out_dir = tmp_path / "out"

    per_frame_rows = [
        {"frame": "amr_probe0000", "it": "0", "time": "0.0", "dt_pred_min": "1.0"},
        {"frame": "amr_probe0001", "it": "10", "time": "0.1", "dt_pred_min": "0.8"},
        {"frame": "amr_probe0002", "it": "20", "time": "0.2", "dt_pred_min": "0.7"},
    ]
    with per_frame_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_frame_rows[0]))
        writer.writeheader()
        writer.writerows(per_frame_rows)

    topcells_rows = [
        {
            "frame": "amr_probe0000",
            "it": "0",
            "time": "0.0",
            "rho": "10.0",
            "p": "5.0",
            "v_abs": "1.0",
            "b_abs": "10.0",
            "v_A": "2.0",
            "shell_index": "100",
            "rank": "1",
        },
        {
            "frame": "amr_probe0000",
            "it": "0",
            "time": "0.0",
            "rho": "999.0",
            "p": "999.0",
            "v_abs": "999.0",
            "b_abs": "999.0",
            "v_A": "999.0",
            "shell_index": "999",
            "rank": "2",
        },
        {
            "frame": "amr_probe0001",
            "it": "10",
            "time": "0.1",
            "rho": "3.0",
            "p": "5.0",
            "v_abs": "2.0",
            "b_abs": "14.0",
            "v_A": "4.0",
            "shell_index": "100",
            "rank": "1",
        },
        {
            "frame": "amr_probe0002",
            "it": "20",
            "time": "0.2",
            "rho": "-1.0",
            "p": "5.0",
            "v_abs": "3.0",
            "b_abs": "14.0",
            "v_A": "5.0",
            "shell_index": "200",
            "rank": "1",
        },
    ]
    with topcells_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(topcells_rows[0]))
        writer.writeheader()
        writer.writerows(topcells_rows)

    rc = main(
        [
            "--case-dir",
            str(case_dir),
            "--per-frame",
            "per_frame.csv",
            "--topcells",
            "topcells.csv",
            "--output-dir",
            str(out_dir),
            "--no-plots",
        ]
    )

    assert rc == 0

    outputs = output_paths(out_dir)
    assert outputs["per_frame"].exists()
    assert outputs["gate"].exists()
    assert outputs["summary"].exists()
    assert outputs["plot"] == out_dir / "rho_bottleneck_trend.png"

    with outputs["per_frame"].open() as handle:
        rows = list(csv.DictReader(handle))

    assert [row["frame"] for row in rows] == ["amr_probe0000", "amr_probe0001", "amr_probe0002"]
    assert rows[0]["rho"] == "10.0"
    assert rows[0]["gate_failed"] == "True"
    assert rows[0]["final_label"] == "NUMERICAL_DOMINANT"
    assert rows[1]["final_label"] == "NUMERICAL_DOMINANT"
    assert rows[2]["gate_failed"] == "True"
    assert rows[2]["final_label"] == "NUMERICAL_DOMINANT"
    assert set(rows[0]) >= {
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
    }

    with outputs["gate"].open() as handle:
        gate_rows = list(csv.DictReader(handle))

    assert len(gate_rows) == 3
    assert set(gate_rows[0]) == {
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
    }
    assert gate_rows[2]["neg_rho_or_p"] == "True"

    summary_text = outputs["summary"].read_text()
    assert "NUMERICAL_DOMINANT" in summary_text
    assert "Latest frame label: NUMERICAL_DOMINANT" in summary_text


def test_frames_filter_does_not_change_rho_ratio_baseline(tmp_path: Path) -> None:
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    per_frame_path = case_dir / "per_frame.csv"
    topcells_path = case_dir / "topcells.csv"
    unfiltered_out = tmp_path / "out_all"
    filtered_out = tmp_path / "out_filtered"

    per_frame_rows = [
        {"frame": "amr_probe0000", "it": "0", "time": "0.0", "dt_pred_min": "1.0"},
        {"frame": "amr_probe0001", "it": "10", "time": "0.1", "dt_pred_min": "0.8"},
        {"frame": "amr_probe0002", "it": "20", "time": "0.2", "dt_pred_min": "0.7"},
    ]
    with per_frame_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_frame_rows[0]))
        writer.writeheader()
        writer.writerows(per_frame_rows)

    topcells_rows = [
        {"frame": "amr_probe0000", "it": "0", "time": "0.0", "rho": "10.0", "p": "5.0", "v_abs": "1.0", "b_abs": "10.0", "v_A": "2.0", "shell_index": "10", "rank": "1", "divb": "1.0", "flux_residual": "0.1"},
        {"frame": "amr_probe0001", "it": "10", "time": "0.1", "rho": "5.0", "p": "5.0", "v_abs": "2.0", "b_abs": "12.0", "v_A": "4.0", "shell_index": "20", "rank": "1", "divb": "1.0", "flux_residual": "0.1"},
        {"frame": "amr_probe0002", "it": "20", "time": "0.2", "rho": "2.5", "p": "5.0", "v_abs": "3.0", "b_abs": "14.0", "v_A": "5.0", "shell_index": "30", "rank": "1", "divb": "1.0", "flux_residual": "0.1"},
    ]
    with topcells_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(topcells_rows[0]))
        writer.writeheader()
        writer.writerows(topcells_rows)

    assert main(
        [
            "--case-dir",
            str(case_dir),
            "--per-frame",
            "per_frame.csv",
            "--topcells",
            "topcells.csv",
            "--output-dir",
            str(unfiltered_out),
            "--no-plots",
        ]
    ) == 0
    assert main(
        [
            "--case-dir",
            str(case_dir),
            "--per-frame",
            "per_frame.csv",
            "--topcells",
            "topcells.csv",
            "--output-dir",
            str(filtered_out),
            "--frames",
            "amr_probe0002",
            "--no-plots",
        ]
    ) == 0

    with output_paths(unfiltered_out)["per_frame"].open() as handle:
        unfiltered_rows = {row["frame"]: row for row in csv.DictReader(handle)}
    with output_paths(filtered_out)["per_frame"].open() as handle:
        filtered_rows = list(csv.DictReader(handle))

    assert [row["frame"] for row in filtered_rows] == ["amr_probe0002"]
    assert filtered_rows[0]["rho_ratio_to_initial"] == unfiltered_rows["amr_probe0002"]["rho_ratio_to_initial"]
    assert filtered_rows[0]["rho_ratio_to_initial"] == "0.25"


def test_missing_diagnostics_force_safe_gate_failure(tmp_path: Path) -> None:
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    per_frame_path = case_dir / "per_frame.csv"
    topcells_path = case_dir / "topcells.csv"
    out_dir = tmp_path / "out"

    per_frame_rows = [
        {"frame": "amr_probe0000", "it": "0", "time": "0.0", "dt_pred_min": "1.0"},
    ]
    with per_frame_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_frame_rows[0]))
        writer.writeheader()
        writer.writerows(per_frame_rows)

    topcells_rows = [
        {"frame": "amr_probe0000", "it": "0", "time": "0.0", "rho": "10.0", "p": "5.0", "v_abs": "1.0", "b_abs": "10.0", "v_A": "2.0", "shell_index": "10", "rank": "1"},
    ]
    with topcells_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(topcells_rows[0]))
        writer.writeheader()
        writer.writerows(topcells_rows)

    assert main(
        [
            "--case-dir",
            str(case_dir),
            "--per-frame",
            "per_frame.csv",
            "--topcells",
            "topcells.csv",
            "--output-dir",
            str(out_dir),
            "--no-plots",
        ]
    ) == 0

    with output_paths(out_dir)["per_frame"].open() as handle:
        row = next(csv.DictReader(handle))
    with output_paths(out_dir)["gate"].open() as handle:
        gate_row = next(csv.DictReader(handle))

    assert row["gate_failed"] == "True"
    assert row["final_label"] == "NUMERICAL_DOMINANT"
    assert gate_row["nan_inf"] == "True"
    assert gate_row["flux_residual_flag"] == "True"


def test_case_dir_relative_paths_take_precedence_over_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cwd = tmp_path / "cwd"
    case_dir = tmp_path / "case"
    cwd.mkdir()
    case_dir.mkdir()
    out_dir = tmp_path / "out"

    for base, marker in [(cwd, "99.0"), (case_dir, "10.0")]:
        per_frame_path = base / "shared_per_frame.csv"
        topcells_path = base / "shared_topcells.csv"
        per_frame_rows = [{"frame": "amr_probe0000", "it": "0", "time": "0.0", "dt_pred_min": "1.0"}]
        topcells_rows = [
            {
                "frame": "amr_probe0000",
                "it": "0",
                "time": "0.0",
                "rho": marker,
                "p": "5.0",
                "v_abs": "1.0",
                "b_abs": "10.0",
                "v_A": "2.0",
                "shell_index": "10",
                "rank": "1",
                "divb": "1.0",
                "flux_residual": "0.1",
            }
        ]
        with per_frame_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(per_frame_rows[0]))
            writer.writeheader()
            writer.writerows(per_frame_rows)
        with topcells_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(topcells_rows[0]))
            writer.writeheader()
            writer.writerows(topcells_rows)

    monkeypatch.chdir(cwd)

    assert main(
        [
            "--case-dir",
            str(case_dir),
            "--per-frame",
            "shared_per_frame.csv",
            "--topcells",
            "shared_topcells.csv",
            "--output-dir",
            str(out_dir),
            "--no-plots",
        ]
    ) == 0

    with output_paths(out_dir)["per_frame"].open() as handle:
        row = next(csv.DictReader(handle))

    assert row["rho"] == "10.0"


@pytest.mark.parametrize(
    ("persistence_frames", "expected_s4"),
    [
        (5, 10),
        (4, 0),
    ],
)
def test_physical_score_s4_boundaries(
    persistence_frames: int,
    expected_s4: int,
) -> None:
    result = physical_score(
        rho_ratio_to_initial=0.1,
        shell_percentile=0.05,
        flux_jump=0.35,
        persistence_frames=persistence_frames,
    )

    assert result["S4"] == expected_s4
