# Analytic Bipolar 1-20 Rsun: PC Stability Ladder

This folder is the local stability-test variant of the analytical bipolar AMRVAC case. It keeps the same physics and grid as the deeper-rooted setup, but adds staged one-core runs that are meant to expose genuine instabilities on a PC instead of masking them with AMRVAC's small-value replacement.

Key setup choices:

- Geometry: `spherical_3D`
- Domain: `r = 1.0..20.0 Rsun`, full theta range, periodic phi
- Base mesh: `72 x 72 x 72`
- Radial stretching: geometric, `stretch_dim(1)='uni'`, `qstretch_baselevel = 1.005`
- MHD model: polytropic with `mhd_gamma = 1.05`
- Divergence control: constrained transport (`typedivbfix='ct'`)
- Small-value handling in staged runs: `fix_small_values = .false.`

The magnetic field is the analytical buried-charge bipolar field implemented in `bipolar_field()`. This case does not include the Titov flux rope, electric driving, or PFSS/SCS input files.

Parker background constants in `mod_usr.t`:

- `rc = 3.45d0`
- `Vs = 117.54d0`
- `rhob = 5.0d9/1000.0d0/unit_numberdensity`
- `Tiso = 2.0d0`

Buried-bipole defaults used by all staged runs:

- `f_q = 1.0d4`
- `f_d = 2.5`
- `f_L = 0.5`

## Staged PC runs

- `pc_init.par`: initialization only, `time_max = 0`
- `pc_smoke.par`: short evolution, `it_max = 50`
- `pc_smoke_offpole.par`: same as smoke, but trims away the exact poles for pole-sensitivity diagnosis
- `pc_medium.par`: moderate evolution, `it_max = 200`
- `pc_long.par`: extended local run, `it_max = 1000`

Each stage writes VTU output under `data/` and keeps log/stdout files with the same stage prefix.

## How to run

Build with your AMRVAC 3.1 tree:

```bash
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 make
```

Run one stage:

```bash
OMP_NUM_THREADS=1 ./run_pc_stability.sh pc_smoke
```

Run the full ladder:

```bash
OMP_NUM_THREADS=1 ./run_pc_stability.sh all
```

The runner stops on the first hard-failure marker it finds in stdout or the AMRVAC log, including `NaN`, `Inf`, `negative`, `dtmin`, `abort`, `mpistop`, `fatal`, or a segmentation fault.
