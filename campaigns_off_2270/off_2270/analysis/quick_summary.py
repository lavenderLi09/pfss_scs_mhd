#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quick high-speed screening (inner-origin and outward propagation)."
    )
    parser.add_argument("--case-dir", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--pattern", default="data/off*.dat")
    parser.add_argument("--par", default="amrvac.par")
    parser.add_argument("--log", default="data/off.log")
    parser.add_argument("--topn", type=int, default=10)
    parser.add_argument("--start", type=int, default=None)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--stride", type=int, default=1, help="Sample every N frames, e.g. 10.")
    parser.add_argument("--frames", default="", help="Comma list (off0030,off0040,...) overrides start/end/stride.")
    parser.add_argument("--outer-rmin", type=float, default=3.0, help="Outer propagation threshold radius.")
    parser.add_argument("--output", default="analysis/highstream_quick_summary.csv")
    return parser.parse_args()


def load_analyzer(case_dir: Path):
    if str(case_dir) not in sys.path:
        sys.path.insert(0, str(case_dir))
    candidates = [
        case_dir / "analysis" / "analyze_snapshot_dt_bottleneck.py",
        case_dir / "analyze_snapshot_dt_bottleneck.py",
    ]
    script = next((p for p in candidates if p.exists()), None)
    if script is None:
        raise FileNotFoundError("Cannot find analyze_snapshot_dt_bottleneck.py")
    spec = importlib.util.spec_from_file_location("dtmod", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to import analyzer from {script}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def frame_id(stem: str) -> int | None:
    m = re.search(r"(\d+)$", stem)
    return int(m.group(1)) if m else None


def select_paths(paths: Sequence[Path], frames_arg: str, start: int | None, end: int | None, stride: int) -> List[Path]:
    if frames_arg.strip():
        wanted = {x.strip() for x in frames_arg.split(",") if x.strip()}
        return [p for p in paths if p.stem in wanted]

    selected: List[Path] = []
    for p in paths:
        n = frame_id(p.stem)
        if n is None:
            continue
        if start is not None and n < start:
            continue
        if end is not None and n > end:
            continue
        if stride > 1:
            base = start if start is not None else 0
            if (n - base) % stride != 0:
                continue
        selected.append(p)
    return selected


def fval(row: Dict[str, object], key: str) -> float:
    try:
        return float(row.get(key, math.nan))
    except (TypeError, ValueError):
        return math.nan


def best_by_speed(rows: Iterable[Dict[str, object]], rmin: float | None = None) -> Dict[str, object] | None:
    best = None
    best_v = -math.inf
    for row in rows:
        r = fval(row, "r")
        v = fval(row, "v_abs")
        if not math.isfinite(v) or not math.isfinite(r):
            continue
        if rmin is not None and r < rmin:
            continue
        if v > best_v:
            best = row
            best_v = v
    return best


def p90_radius(rows: Sequence[Dict[str, object]]) -> float:
    rs = sorted([fval(r, "r") for r in rows if math.isfinite(fval(r, "r"))])
    if not rs:
        return math.nan
    idx = int(round(0.9 * (len(rs) - 1)))
    return rs[idx]


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.resolve()
    mod = load_analyzer(case_dir)

    par_cfg = mod.read_par_config(case_dir / args.par)
    log_by_it = mod.parse_log(case_dir / args.log)
    base_edges = mod.build_base_radial_edges(
        par_cfg["xprobmin1"],
        par_cfg["xprobmax1"],
        int(par_cfg["domain_nx1"]),
        par_cfg["qstretch_baselevel"],
    )

    paths = sorted(case_dir.glob(args.pattern))
    paths = select_paths(paths, args.frames, args.start, args.end, args.stride)
    if not paths:
        raise FileNotFoundError("No frames selected.")

    rows: List[Dict[str, object]] = []
    for i, dat_path in enumerate(paths, start=1):
        summary, _top_rows, top_vabs_rows, top_vabs_inner_rows, physical_summary, *_ = mod.analyze_snapshot(
            dat_path, par_cfg, log_by_it, args.topn, base_edges
        )
        inner_top = best_by_speed(top_vabs_inner_rows)
        global_top = best_by_speed(top_vabs_rows)
        outer_top = best_by_speed(top_vabs_rows, rmin=args.outer_rmin)

        row = {
            "frame": summary["frame"],
            "it": summary["it"],
            "time": summary["time"],
            "dt_log": summary["dt_log"],
            "dominant_physics": summary["dominant_physics"],
            "vabs_p99": physical_summary["vabs_p99"],
            "vabs_max": physical_summary["vabs_max"],
            "inner_top_vabs": fval(inner_top or {}, "v_abs"),
            "inner_top_r": fval(inner_top or {}, "r"),
            "inner_top_rho": fval(inner_top or {}, "rho"),
            "global_top_vabs": fval(global_top or {}, "v_abs"),
            "global_top_r": fval(global_top or {}, "r"),
            "global_top_rho": fval(global_top or {}, "rho"),
            "outer_top_vabs": fval(outer_top or {}, "v_abs"),
            "outer_top_r": fval(outer_top or {}, "r"),
            "outer_top_rho": fval(outer_top or {}, "rho"),
            "top_vabs_r_p90": p90_radius(top_vabs_rows),
            "outer_present": int(outer_top is not None),
        }
        rows.append(row)
        print(
            f"[{i}/{len(paths)}] {row['frame']} "
            f"inner_v={row['inner_top_vabs']:.3e} "
            f"outer_v={row['outer_top_vabs']:.3e} "
            f"r_top={row['global_top_r']:.3f}",
            flush=True,
        )

    out_path = case_dir / args.output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
