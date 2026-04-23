# Bipolar HPC Test

This directory is an isolated HPC-only sibling of `analytic_bipolar_1to20_stretched`.

It exists so that:

- the parent PC stability ladder stays unchanged
- HPC mesh and AMR tuning can use a different `mod_usr.t`
- HPC runs can have separate `.par`, `.sh`, and job files

The baseline source for this case is the off-pole branch of the parent directory at the time this folder was created. After creation, changes inside this folder should not be mirrored back into the parent case unless a later review explicitly decides to do that.

## Pilot configuration

This isolated case uses:

- off-pole theta range: `0.005 .. 0.495`
- base mesh: `96 x 96 x 96`
- block size: `12 x 12 x 12`
- AMR ceiling: `refine_max_level=2`
- hybrid AMR rule in this folder's `mod_usr.t`

## Local preflight

```bash
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 /Users/zhaoyan/Documents/codes/amrvac-amrvac3.1/setup.pl -d=3
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 make
OMP_NUM_THREADS=1 ./amrvac -i hpc_pilot_init_offpole.par
./run_hpc_pilot.sh hpc_pilot_init_offpole 1
```

## Cluster submission

```bash
bsub < main_hpc_pilot.job
```
