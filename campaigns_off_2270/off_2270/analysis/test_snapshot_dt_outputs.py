#!/usr/bin/env python3
"""Smoke tests for snapshot dt analysis outputs."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


REPO = Path("/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd")
CASE_DIR = REPO / "amrvac_polytropic" / "off_2270"
ANALYSIS_DIR = CASE_DIR / "analysis"
SCRIPT = ANALYSIS_DIR / "analyze_snapshot_dt_bottleneck.py"
OUTPUT_DIR = Path("/tmp/off_2270_analysis_test")

if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))


def test_internal_energy_conversion() -> None:
    from analyze_snapshot_dt_bottleneck import conservative_internal_to_primitive

    rho = np.array([2.0, 4.0], dtype=np.float64)
    m1 = np.array([6.0, 8.0], dtype=np.float64)
    m2 = np.array([0.0, -4.0], dtype=np.float64)
    m3 = np.array([2.0, 0.0], dtype=np.float64)
    eint = np.array([10.0, 20.0], dtype=np.float64)
    gamma = 1.05

    v_r, v_th, v_ph, pth = conservative_internal_to_primitive(rho, m1, m2, m3, eint, gamma)

    np.testing.assert_allclose(v_r, np.array([3.0, 2.0]))
    np.testing.assert_allclose(v_th, np.array([0.0, -1.0]))
    np.testing.assert_allclose(v_ph, np.array([1.0, 0.0]))
    np.testing.assert_allclose(pth, np.array([0.5, 1.0]))


def main() -> int:
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    test_internal_energy_conversion()

    cmd = [
        sys.executable,
        str(SCRIPT),
        "--case-dir",
        str(CASE_DIR),
        "--frames",
        "off0000",
        "--topn",
        "3",
        "--output-dir",
        str(OUTPUT_DIR),
    ]
    subprocess.run(cmd, check=True)

    expected = [
        OUTPUT_DIR / "dt_bottleneck_summary.csv",
        OUTPUT_DIR / "bottleneck_hotspots.csv",
        OUTPUT_DIR / "bottleneck_hotspots.md",
        OUTPUT_DIR / "physical_speed_summary.csv",
        OUTPUT_DIR / "surface_shell_summary.csv",
        OUTPUT_DIR / "top_vabs_cells.csv",
        OUTPUT_DIR / "top_vabs_inner_cells.csv",
        OUTPUT_DIR / "off0000_rtheta_maps.png",
        OUTPUT_DIR / "off0000_bottleneck_scatter.png",
        OUTPUT_DIR / "off0000_top_vabs_br_overlay.png",
        OUTPUT_DIR / "off0000_top_vabs_inner_br_overlay.png",
    ]
    missing = [path for path in expected if not path.exists()]
    if missing:
        raise AssertionError(f"Missing expected outputs: {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
