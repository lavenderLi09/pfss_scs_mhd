#!/usr/bin/env python3
"""
Export stretched spherical grid Brtp to VTK for ParaView / VisIt.

Uses pfss.Brtp2vts (structured .vts, fast) or Brtp2vtu (.vtu, very slow for large grids).

Typical workflow:
  1. Run initial_magnetic.ipynb (saves Brtp_combined.npy + rtp_stretched.npz).
  2. python export_brtp_vtk.py --rtp rtp_stretched.npz --brtp Brtp_combined.npy --out Brtp_1_15 --format vts

If rtp_stretched.npz is missing, pass the same R_min, R_max, q, nr, nt, nphi as in the notebook.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "codes"))

from pfss.pfss_module import Brtp2vts, Brtp2vtu


def build_stretched_mesh(R_min: float, R_max: float, q: float, nr: int, nt: int, nphi: int):
    dr0 = (R_max - R_min) * (q - 1) / (q**nr - 1)
    r_edges = R_min + dr0 * (q ** np.arange(nr + 1) - 1) / (q - 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    dth = np.pi / nt
    dph = 2 * np.pi / nphi
    t_centers = np.linspace(dth / 2, np.pi - dth / 2, nt)
    p_centers = np.linspace(dph / 2, 2 * np.pi - dph / 2, nphi)
    rr, tt, pp = np.meshgrid(r_centers, t_centers, p_centers, indexing="ij")
    return rr, tt, pp


def main() -> None:
    ap = argparse.ArgumentParser(description="Export Brtp on stretched grid to VTK")
    ap.add_argument("--brtp", default="Brtp_combined.npy", help="numpy array shape (3, nr, nt, nphi)")
    ap.add_argument(
        "--rtp",
        default="rtp_stretched.npz",
        help="npz with rr, tt, pp (from notebook). If missing, use --rebuild-mesh",
    )
    ap.add_argument("--rebuild-mesh", action="store_true", help="ignore --rtp; rebuild mesh from parameters")
    ap.add_argument("--out", default="Brtp_stretched", help="output base name (no extension)")
    ap.add_argument("--format", choices=["vts", "vtu"], default="vts")
    ap.add_argument("--R-min", type=float, default=1.0)
    ap.add_argument("--R-max", type=float, default=15.0)
    ap.add_argument("--q", type=float, default=1.005)
    ap.add_argument("--nr", type=int, default=500)
    ap.add_argument("--nt", type=int, default=96)
    ap.add_argument("--nphi", type=int, default=192)
    args = ap.parse_args()

    brtp_path = ROOT / args.brtp
    if not brtp_path.is_file():
        raise SystemExit(f"Not found: {brtp_path}")

    Brtp = np.load(brtp_path)
    if Brtp.ndim != 4 or Brtp.shape[0] != 3:
        raise SystemExit(f"Expected Brtp shape (3, nr, nt, nphi), got {Brtp.shape}")

    rtp_file = ROOT / args.rtp
    if args.rebuild_mesh or not rtp_file.is_file():
        if not args.rebuild_mesh and not rtp_file.is_file():
            print(f"Warning: {rtp_file} not found, rebuilding mesh from CLI parameters.")
        rr, tt, pp = build_stretched_mesh(
            args.R_min, args.R_max, args.q, args.nr, args.nt, args.nphi
        )
    else:
        d = np.load(rtp_file)
        rr, tt, pp = d["rr"], d["tt"], d["pp"]

    if tuple(Brtp.shape[1:]) != rr.shape:
        raise SystemExit(
            f"Shape mismatch: Brtp[1:]={Brtp.shape[1:]} vs grid {rr.shape}. "
            "Use matching --nr/--nt/--nphi/--R-max/--q or save rtp_stretched.npz from the same run."
        )

    out_base = ROOT / args.out
    rtp = [rr, tt, pp]

    if args.format == "vts":
        Brtp2vts(Brtp, rtp, vts_name=str(out_base))
        print(f"Wrote {out_base}.vts")
    else:
        npts = rr.size
        if npts > 500_000:
            print(
                f"Warning: {npts} points — Brtp2vtu uses Python loops and may take a long time. Prefer --format vts."
            )
        Brtp2vtu(Brtp, rtp, vtu_name=str(out_base))
        print(f"Wrote {out_base}.vtu")


if __name__ == "__main__":
    main()
