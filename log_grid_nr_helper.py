"""Helpers for a pure logarithmic/geometric radial grid."""

from __future__ import annotations

import math


def compute_log_grid_nr(r_min: float, r_out: float, q: float) -> int:
    """Return the number of radial grid points for a log grid.

    The grid is defined by:
        r_i = r_min * q**i

    The returned ``nr`` counts grid points, including both endpoints. The
    growth factor ``q`` is kept exact, and ``nr`` is chosen so that the last
    grid point is the smallest one that reaches or exceeds ``r_out``.
    """

    if r_min <= 0.0:
        raise ValueError("r_min must be positive.")
    if r_out <= r_min:
        raise ValueError("r_out must be greater than r_min.")
    if q <= 1.0:
        raise ValueError("q must be greater than 1.")

    n_cells = math.ceil(math.log(r_out / r_min) / math.log(q))
    return n_cells + 1


def radial_point(r_min: float, q: float, i: int) -> float:
    """Return the i-th radial grid point for the log grid."""

    if i < 0:
        raise ValueError("i must be non-negative.")
    return r_min * (q**i)


if __name__ == "__main__":
    R_MIN = 1.0
    R_OUT = 12.0
    Q = 1.005

    nr = compute_log_grid_nr(R_MIN, R_OUT, Q)
    n_cells = nr - 1
    r_prev = radial_point(R_MIN, Q, nr - 2)
    r_last = radial_point(R_MIN, Q, nr - 1)

    print(f"r_min={R_MIN}")
    print(f"r_out={R_OUT}")
    print(f"q={Q}")
    print(f"nr={nr}")
    print(f"n_cells={n_cells}")
    print(f"r[nr-2]={r_prev:.12f}")
    print(f"r[nr-1]={r_last:.12f}")
