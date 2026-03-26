#!/usr/bin/env python3
"""Plot an HMI synoptic Br map from a FITS file."""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

_cache_root = Path(tempfile.gettempdir()) / "mpl-cache-pfss-scs"
_cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "config"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root / "xdg"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def load_map(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path) as hdul:
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None and hdu.data.ndim == 2:
                return np.asarray(hdu.data, dtype=float), hdu.header
    raise ValueError(f"No 2D image HDU found in {path}")


def latitude_edges_deg(header: fits.Header, ny: int) -> np.ndarray:
    sinlat_edges = ((np.arange(ny + 1) + 0.5) - header["CRPIX2"]) * header["CDELT2"] + header["CRVAL2"]
    sinlat_edges = np.clip(sinlat_edges, -1.0, 1.0)
    return np.degrees(np.arcsin(sinlat_edges))


def longitude_edges_deg(nx: int) -> np.ndarray:
    return np.linspace(0.0, 360.0, nx + 1)


def plot_br(
    data: np.ndarray,
    header: fits.Header,
    output: Path,
    percentile: float = 99.5,
) -> None:
    ny, nx = data.shape
    lon_edges = longitude_edges_deg(nx)
    lat_edges = latitude_edges_deg(header, ny)
    clip = np.nanpercentile(np.abs(data), percentile)

    fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)
    mesh = ax.pcolormesh(
        lon_edges,
        lat_edges,
        data,
        cmap="RdBu_r",
        shading="auto",
        vmin=-clip,
        vmax=clip,
    )

    rotation = header.get("CAR_ROT", "unknown")
    content = header.get("CONTENT", "HMI Synoptic Br Map")
    fig.colorbar(mesh, ax=ax, pad=0.02, label=r"$B_r$ [Mx cm$^{-2}$]")
    ax.set_title(f"{content} (CR {rotation})")
    ax.set_xlabel("Carrington longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.grid(color="0.7", alpha=0.25, linewidth=0.5)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fits_path", type=Path, help="Path to the HMI synoptic FITS file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="PNG output path (default: alongside FITS file)",
    )
    parser.add_argument(
        "--percentile",
        type=float,
        default=99.5,
        help="Symmetric clip percentile for the color scale",
    )
    args = parser.parse_args()

    data, header = load_map(args.fits_path)
    output = args.output or args.fits_path.with_suffix(".png")
    plot_br(data, header, output, percentile=args.percentile)
    print(output)


if __name__ == "__main__":
    main()
