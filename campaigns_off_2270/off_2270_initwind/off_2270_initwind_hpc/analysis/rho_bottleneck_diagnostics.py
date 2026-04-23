from __future__ import annotations

from typing import Any

import numpy as np

NUMERICAL_DOMINANT = "NUMERICAL_DOMINANT"
PHYSICAL_HIGH_RISK = "PHYSICAL_HIGH_RISK"
MIXED_UNCERTAIN = "MIXED_UNCERTAIN"
PHYSICAL_ACCEPTABLE = "PHYSICAL_ACCEPTABLE"


def evaluate_numerical_gate(rho, p, divb, flux_residual: float) -> dict[str, Any]:
    rho_arr = np.asarray(rho)
    p_arr = np.asarray(p)
    divb_arr = np.abs(np.asarray(divb))

    neg_rho_or_p = bool(np.any(rho_arr <= 0) or np.any(p_arr <= 0))

    nan_inf = bool(
        (not np.all(np.isfinite(rho_arr)))
        or (not np.all(np.isfinite(p_arr)))
        or (not np.all(np.isfinite(divb_arr)))
    )

    finite_mask = np.isfinite(divb_arr)
    finite_divb = divb_arr[finite_mask]
    if finite_divb.size:
        divb_median = float(np.median(finite_divb))
        divb_p999 = float(np.percentile(finite_divb, 99.9))
    else:
        divb_median = float("inf")
        divb_p999 = float("inf")

    divb_flag = bool(
        np.isfinite(divb_median)
        and divb_p999 > 20 * divb_median
    )

    flux_residual_flag = bool((not np.isfinite(flux_residual)) or flux_residual > 0.30)
    gate_failed = bool(neg_rho_or_p or nan_inf or divb_flag or flux_residual_flag)

    return {
        "neg_rho_or_p": neg_rho_or_p,
        "nan_inf": nan_inf,
        "divb_flag": divb_flag,
        "divb_median": divb_median,
        "divb_p999": divb_p999,
        "flux_residual_flag": flux_residual_flag,
        "gate_failed": gate_failed,
    }


def physical_score(
    *,
    rho_ratio_to_initial: float,
    shell_percentile: float,
    flux_jump: float,
    persistence_frames: int,
) -> dict[str, int]:
    s1 = 40 if rho_ratio_to_initial < 0.2 else 20 if rho_ratio_to_initial < 0.4 else 0
    s2 = 30 if shell_percentile < 0.1 else 15 if shell_percentile < 1.0 else 0
    s3 = 20 if flux_jump > 0.30 else 10 if flux_jump >= 0.15 else 0
    s4 = 10 if persistence_frames >= 5 else 0
    score = s1 + s2 + s3 + s4
    return {"S1": s1, "S2": s2, "S3": s3, "S4": s4, "score": score}


def classify_final_label(gate_failed: bool, score: float) -> str:
    if gate_failed:
        return NUMERICAL_DOMINANT
    if score >= 60:
        return PHYSICAL_HIGH_RISK
    if score >= 30:
        return MIXED_UNCERTAIN
    return PHYSICAL_ACCEPTABLE
