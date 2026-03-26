"""Basic verification for the log-grid nr helper."""

from __future__ import annotations

import math

from log_grid_nr_helper import compute_log_grid_nr, radial_point


def main() -> None:
    r_min = 1.0
    r_out = 12.0
    q = 1.005

    nr = compute_log_grid_nr(r_min, r_out, q)
    n_cells = nr - 1
    raw_ratio = math.log(r_out / r_min) / math.log(q)
    r_prev = radial_point(r_min, q, nr - 2)
    r_last = radial_point(r_min, q, nr - 1)

    assert nr == 500, f"expected nr=500, got {nr}"
    assert n_cells == 499, f"expected n_cells=499, got {n_cells}"
    assert r_last >= r_out, "last grid point must reach or exceed r_out"
    assert r_prev < r_out, "previous grid point must stay below r_out"

    print(f"log(12.0)/log(1.005) = {raw_ratio:.12f}")
    print(f"nr = {nr}")
    print(f"n_cells = {n_cells}")
    print(f"r[nr-2] = {r_prev:.12f}")
    print(f"r[nr-1] = {r_last:.12f}")
    print("verification passed")


if __name__ == "__main__":
    main()
