#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import vtk
from PIL import Image, ImageDraw
from vtk.util.numpy_support import vtk_to_numpy


CASE_DIR = Path(
    "/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/"
    "amrvac_polytropic/analytic_bipolar_1to20_stretched"
)
DATA_DIR = CASE_DIR / "data"
OUT_DIR = CASE_DIR / "analysis" / "pressure_anomaly"
INITIAL_VTU = DATA_DIR / "pc_smoke_offpole0000.vtu"
FINAL_VTU = DATA_DIR / "pc_smoke_offpole0001.vtu"


def read_vtu(path: Path) -> tuple[float, dict[str, np.ndarray], np.ndarray]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()

    centers = vtk.vtkCellCenters()
    centers.SetInputData(grid)
    centers.Update()
    points = vtk_to_numpy(centers.GetOutput().GetPoints().GetData()).astype(float)

    cell_data = grid.GetCellData()
    arrays = {}
    for name in ["p", "rho", "Vr", "Vt", "Vp"]:
        arrays[name] = vtk_to_numpy(cell_data.GetArray(name)).astype(float)

    time_array = grid.GetFieldData().GetArray("TIME")
    time_value = float(vtk_to_numpy(time_array)[0])
    return time_value, arrays, points


def spherical_from_cart(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    radius = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(np.clip(z / radius, -1.0, 1.0))
    phi = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    return radius, theta, phi


def unique_index(values: np.ndarray, decimals: int = 12) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(values, decimals=decimals)
    uniq = np.unique(rounded)
    index = np.searchsorted(uniq, rounded)
    return uniq, index


def build_max_grid(
    radial_index: np.ndarray,
    theta_index: np.ndarray,
    values: np.ndarray,
    nr: int,
    nt: int,
) -> np.ndarray:
    grid = np.full((nr, nt), -np.inf)
    np.maximum.at(grid, (radial_index, theta_index), values)
    return grid


def build_shell_grid(
    shell_mask: np.ndarray,
    theta_index: np.ndarray,
    phi_index: np.ndarray,
    values: np.ndarray,
    nt: int,
    np_: int,
) -> np.ndarray:
    grid = np.full((nt, np_), np.nan)
    tsel = theta_index[shell_mask]
    psel = phi_index[shell_mask]
    vsel = values[shell_mask]
    for t_idx, p_idx, val in zip(tsel, psel, vsel):
        grid[t_idx, p_idx] = val
    return grid


def sequential_rgb(data: np.ndarray, log_scale: bool = False) -> np.ndarray:
    finite = np.isfinite(data)
    vals = data.copy()
    vals[~finite] = np.nan
    if log_scale:
        vals[finite] = np.log10(np.maximum(vals[finite], 1.0e-30))
    vmin = np.nanmin(vals)
    vmax = np.nanmax(vals)
    norm = np.zeros_like(vals)
    if vmax > vmin:
        norm[finite] = (vals[finite] - vmin) / (vmax - vmin)
    rgb = np.zeros(data.shape + (3,), dtype=np.uint8)
    rgb[..., 0] = (255.0 * np.clip(norm ** 0.7, 0.0, 1.0)).astype(np.uint8)
    rgb[..., 1] = (255.0 * np.clip(0.2 + 0.8 * norm, 0.0, 1.0)).astype(np.uint8)
    rgb[..., 2] = (255.0 * np.clip(0.05 + 0.35 * (1.0 - norm), 0.0, 1.0)).astype(np.uint8)
    rgb[~finite] = np.array([255, 255, 255], dtype=np.uint8)
    return rgb


def diverging_rgb(data: np.ndarray) -> np.ndarray:
    finite = np.isfinite(data)
    vmax = np.nanmax(np.abs(data[finite]))
    if vmax <= 0.0:
        vmax = 1.0
    norm = np.zeros_like(data)
    norm[finite] = np.clip(data[finite] / vmax, -1.0, 1.0)
    rgb = np.full(data.shape + (3,), 255, dtype=np.uint8)
    pos = finite & (norm >= 0.0)
    neg = finite & (norm < 0.0)
    rgb[pos, 0] = 255
    rgb[pos, 1] = (255.0 * (1.0 - norm[pos])).astype(np.uint8)
    rgb[pos, 2] = (255.0 * (1.0 - norm[pos])).astype(np.uint8)
    rgb[neg, 0] = (255.0 * (1.0 + norm[neg])).astype(np.uint8)
    rgb[neg, 1] = (255.0 * (1.0 + norm[neg])).astype(np.uint8)
    rgb[neg, 2] = 255
    rgb[~finite] = np.array([255, 255, 255], dtype=np.uint8)
    return rgb


def to_image(rgb: np.ndarray, scale_x: int = 8, scale_y: int = 8) -> Image.Image:
    img = Image.fromarray(np.flipud(rgb), mode="RGB")
    return img.resize((img.width * scale_x, img.height * scale_y), Image.Resampling.NEAREST)


def add_caption(image: Image.Image, caption: str) -> Image.Image:
    pad = 30
    canvas = Image.new("RGB", (image.width, image.height + pad), color=(255, 255, 255))
    canvas.paste(image, (0, pad))
    draw = ImageDraw.Draw(canvas)
    draw.text((10, 8), caption, fill=(0, 0, 0))
    return canvas


def format_sci(value: float) -> str:
    return f"{value:.6e}"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t0, init_data, points = read_vtu(INITIAL_VTU)
    t1, final_data, _ = read_vtu(FINAL_VTU)
    radius, theta, phi = spherical_from_cart(points)
    delta_p = final_data["p"] - init_data["p"]

    radial_values, radial_index = unique_index(radius, decimals=10)
    theta_values, theta_index = unique_index(theta, decimals=10)
    phi_values, phi_index = unique_index(phi, decimals=10)

    p_max_grid = build_max_grid(
        radial_index, theta_index, final_data["p"], radial_values.size, theta_values.size
    )
    dp_max_grid = build_max_grid(
        radial_index, theta_index, delta_p, radial_values.size, theta_values.size
    )

    inner_shell_index = 0
    inner_shell_mask = radial_index == inner_shell_index
    inner_shell_p = build_shell_grid(
        inner_shell_mask,
        theta_index,
        phi_index,
        final_data["p"],
        theta_values.size,
        phi_values.size,
    )
    inner_shell_dp = build_shell_grid(
        inner_shell_mask,
        theta_index,
        phi_index,
        delta_p,
        theta_values.size,
        phi_values.size,
    )

    rt_p_image = add_caption(
        to_image(sequential_rgb(p_max_grid, log_scale=True)),
        "Final pressure, phi-max over each (r, theta); bottom is r_min, left is theta_min",
    )
    rt_dp_image = add_caption(
        to_image(diverging_rgb(dp_max_grid)),
        "Pressure change, phi-max over each (r, theta); red = increase, blue = decrease",
    )
    inner_p_image = add_caption(
        to_image(sequential_rgb(inner_shell_p, log_scale=True)),
        f"Final pressure on first radial shell r={radial_values[inner_shell_index]:.6f}",
    )
    inner_dp_image = add_caption(
        to_image(diverging_rgb(inner_shell_dp)),
        "Pressure change on first radial shell; left/right edges are phi = 0/2pi",
    )

    two_panel_1 = Image.new(
        "RGB",
        (rt_p_image.width + rt_dp_image.width, max(rt_p_image.height, rt_dp_image.height)),
        color=(255, 255, 255),
    )
    two_panel_1.paste(rt_p_image, (0, 0))
    two_panel_1.paste(rt_dp_image, (rt_p_image.width, 0))
    two_panel_1.save(OUT_DIR / "pressure_r_theta_summary.png")

    two_panel_2 = Image.new(
        "RGB",
        (inner_p_image.width + inner_dp_image.width, max(inner_p_image.height, inner_dp_image.height)),
        color=(255, 255, 255),
    )
    two_panel_2.paste(inner_p_image, (0, 0))
    two_panel_2.paste(inner_dp_image, (inner_p_image.width, 0))
    two_panel_2.save(OUT_DIR / "pressure_inner_shell_summary.png")

    pmax_idx = int(np.argmax(final_data["p"]))
    theta_min = float(theta_values.min())
    theta_max = float(theta_values.max())
    dtheta_to_cut = np.minimum(theta - theta_min, theta_max - theta)
    near_theta_cut = dtheta_to_cut < 0.05
    near_inner_r = radius < 1.3
    p_top1_mask = final_data["p"] >= np.percentile(final_data["p"], 99.0)
    p_gt_1e2_mask = final_data["p"] > 1.0e-2
    p_gt_1e3_mask = final_data["p"] > 1.0e-3

    top20_idx = np.argsort(final_data["p"])[-20:][::-1]
    table_lines = []
    for idx in top20_idx:
        table_lines.append(
            "| "
            + " | ".join(
                [
                    str(int(idx)),
                    format_sci(final_data["p"][idx]),
                    f"{radius[idx]:.6f}",
                    f"{theta[idx]:.6f}",
                    f"{phi[idx]:.6f}",
                    f"{final_data['Vr'][idx]:.3f}",
                    f"{final_data['Vt'][idx]:.3f}",
                    f"{final_data['Vp'][idx]:.3f}",
                    format_sci(final_data["rho"][idx]),
                ]
            )
            + " |"
        )

    summary = f"""Pressure anomaly summary for pc_smoke_offpole

inputs
initial_vtu={INITIAL_VTU}
final_vtu={FINAL_VTU}
time_start={t0:.9e}
time_end={t1:.9e}

peak_pressure
p_max={format_sci(final_data["p"][pmax_idx])}
r={radius[pmax_idx]:.6f}
theta={theta[pmax_idx]:.6f}
phi={phi[pmax_idx]:.6f}
Vr={final_data["Vr"][pmax_idx]:.3f}
Vt={final_data["Vt"][pmax_idx]:.3f}
Vp={final_data["Vp"][pmax_idx]:.3f}

boundary_localization
p_gt_1e2_count={int(np.sum(p_gt_1e2_mask))}
p_gt_1e2_near_theta_cut={int(np.sum(p_gt_1e2_mask & near_theta_cut))}
p_gt_1e2_near_inner_r={int(np.sum(p_gt_1e2_mask & near_inner_r))}
p_gt_1e2_near_both={int(np.sum(p_gt_1e2_mask & near_theta_cut & near_inner_r))}
p_gt_1e3_count={int(np.sum(p_gt_1e3_mask))}
p_gt_1e3_near_theta_cut={int(np.sum(p_gt_1e3_mask & near_theta_cut))}
p_gt_1e3_near_inner_r={int(np.sum(p_gt_1e3_mask & near_inner_r))}
p_gt_1e3_near_both={int(np.sum(p_gt_1e3_mask & near_theta_cut & near_inner_r))}
p_top1_count={int(np.sum(p_top1_mask))}
p_top1_near_theta_cut={int(np.sum(p_top1_mask & near_theta_cut))}
p_top1_near_inner_r={int(np.sum(p_top1_mask & near_inner_r))}
p_top1_near_both={int(np.sum(p_top1_mask & near_theta_cut & near_inner_r))}
first_radial_shell_r={radial_values[inner_shell_index]:.6f}

top20_pressure_cells
{chr(10).join(table_lines)}
"""

    (OUT_DIR / "pressure_anomaly_summary.txt").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
