# Off-Pole HPC Pilot Mesh And AMR Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a first cheap-but-meaningful HPC pilot for the off-pole analytical bipolar case with a higher base mesh, `refine_max_level=2`, and a controlled hybrid AMR policy that refines the inner corona and equatorial current-sheet-like region without feeding the known theta-cut artifact.

**Architecture:** Keep the current off-pole theta-boundary branch and local PC ladder intact. Add a separate HPC parameter file pair and job runner, then replace the blanket `qt==0` refinement rule in `mod_usr.t` with user-configurable geometric-plus-field AMR controls read from `&usr_list`. Use `specialthreshold` only to suppress refinement near the theta cuts and outer radial shell; do not penalize periodic phi boundaries.

**Tech Stack:** MPI-AMRVAC 3.1, preprocessed Fortran in `mod_usr.t`, AMRVAC namelist `.par` files, Bash runner scripts, LSF `bsub` job script.

---

## File Structure

- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t`
  Responsibility: add user-tunable AMR controls, implement the hybrid forced-refinement logic in `special_refine_grid()`, and narrow `specialthreshold()` so it suppresses refinement only near the theta cuts and outer radial boundary.
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_init_offpole.par`
  Responsibility: one-shot initialization verification for the new HPC mesh and AMR namelist.
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_offpole.par`
  Responsibility: the actual first HPC pilot run configuration.
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/run_hpc_pilot.sh`
  Responsibility: run `init`, `pilot`, or `all` with MPI and fail fast on AMRVAC hard-failure markers.
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/main_hpc_pilot.job`
  Responsibility: submit the HPC pilot through the same LSF family already used by `main.job`.
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/README.md`
  Responsibility: document the new off-pole HPC pilot, its mesh/AMR choices, and how to run it locally versus on the cluster.

## Assumptions Locked In

- The first HPC campaign stays **off-pole only** and does not attempt to rescue the known full-pole crash.
- The cluster scheduler is still **LSF** because the existing job file uses `#BSUB`.
- The first pilot uses a **higher base mesh first**: `96 x 96 x 96`, `block_nx*=12`, `refine_max_level=2`.
- The first pilot keeps the current **analytical-continuation theta boundary** unchanged; this plan only changes mesh and AMR behavior.
- The AMR “current sheet area” is approximated by a hybrid rule: broad equatorial window plus a low-`|Br|/|B|` trigger, not a full current-density calculation inside `special_refine_grid()`.

### Task 1: Add HPC Pilot Parameter Files First

**Files:**
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_init_offpole.par`
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_offpole.par`
- Test: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_init_offpole.par`

- [ ] **Step 1: Create the initialization-only HPC pilot file with the new AMR knobs before `mod_usr.t` knows about them**

```fortran
$AMRVAC_DIR/setup.pl -d=3

 &filelist
        base_filename='data/hpc_pilot_init_offpole'
        saveprim=.true.
        autoconvert=.true.
        convert_type='vtuBCCmpi'
        nwauxio=11
 /

 &savelist
        ditsave_log=5
        itsave(1,1)=0
        itsave(1,2)=0
/

 &stoplist
        dtmin=1.D-8
        time_max=0.0d0
        wall_time_max=0.15d0
 /

 &methodlist
        time_stepper='threestep'
        flux_scheme=20*'hll'
        limiter=20*'vanleer'
        tvdlfeps=0.d0
        small_density = 2497.9540699435920d0
        small_values_method = 'error'
        fix_small_values=.false.
 /

 &boundlist
        typeboundary_min1 = 8*'special'
        typeboundary_max1 = 8*'special'
        typeboundary_min2 = 8*'special'
        typeboundary_max2 = 8*'special'
        typeboundary_min3 = 8*'periodic'
        typeboundary_max3 = 8*'periodic'
        nghostcells=4
 /

 &meshlist
        refine_criterion=3
        refine_max_level=2
        refine_threshold=0.18d0,0.24d0,18*0.18d0
        derefine_ratio=20*0.08d0
        w_refine_weight(4)=0.5d0
        w_refine_weight(5)=0.25d0
        w_refine_weight(6)=0.25d0
        block_nx1=12
        block_nx2=12
        block_nx3=12
        domain_nx1=96
        domain_nx2=96
        domain_nx3=96
        xprobmin1=1.0d0
        xprobmax1=20.0d0
        xprobmin2=0.005d0
        xprobmax2=0.495d0
        xprobmin3=0.0d0
        xprobmax3=1.0d0
        stretch_dim(1)='uni'
        qstretch_baselevel=1.005d0
        ditregrid=2
        max_blocks=40000
 /

 &paramlist
        typecourant='maxsum'
        courantpar=0.8d0
 /

 &mhd_list
        mhd_energy=.true.
        mhd_internal_e=.true.
        mhd_viscosity=.false.
        mhd_gravity=.true.
        mhd_gamma=1.05d0
        typedivbfix='ct'
/

 &usr_list
        f_q=1.0d4
        f_d=2.5d0
        f_L=0.5d0
        amr_r_core_max=4.0d0
        amr_r_sheet_max=8.0d0
        amr_theta_guard=0.02d0
        amr_sheet_halfwidth=0.35d0
        amr_br_ratio_max=0.20d0
        amr_force_max_level=2
        amr_outer_relax_frac=0.25d0
        amr_theta_relax_frac=0.10d0
/
```

- [ ] **Step 2: Create the actual HPC pilot run file with the same mesh/AMR setup and a short evolution budget**

```fortran
$AMRVAC_DIR/setup.pl -d=3

 &filelist
        base_filename='data/hpc_pilot_offpole'
        saveprim=.true.
        autoconvert=.true.
        convert_type='vtuBCCmpi'
        nwauxio=11
 /

 &savelist
        ditsave_log=5
        dtsave_custom=0.01d0
        dtsave_dat=0.02d0
        itsave(1,1)=0
        itsave(1,2)=0
/

 &stoplist
        dtmin=1.D-8
        it_max=120
        wall_time_max=0.40d0
 /

 &methodlist
        time_stepper='threestep'
        flux_scheme=20*'hll'
        limiter=20*'vanleer'
        tvdlfeps=0.d0
        small_density = 2497.9540699435920d0
        small_values_method = 'error'
        fix_small_values=.false.
 /

 &boundlist
        typeboundary_min1 = 8*'special'
        typeboundary_max1 = 8*'special'
        typeboundary_min2 = 8*'special'
        typeboundary_max2 = 8*'special'
        typeboundary_min3 = 8*'periodic'
        typeboundary_max3 = 8*'periodic'
        nghostcells=4
 /

 &meshlist
        refine_criterion=3
        refine_max_level=2
        refine_threshold=0.18d0,0.24d0,18*0.18d0
        derefine_ratio=20*0.08d0
        w_refine_weight(4)=0.5d0
        w_refine_weight(5)=0.25d0
        w_refine_weight(6)=0.25d0
        block_nx1=12
        block_nx2=12
        block_nx3=12
        domain_nx1=96
        domain_nx2=96
        domain_nx3=96
        xprobmin1=1.0d0
        xprobmax1=20.0d0
        xprobmin2=0.005d0
        xprobmax2=0.495d0
        xprobmin3=0.0d0
        xprobmax3=1.0d0
        stretch_dim(1)='uni'
        qstretch_baselevel=1.005d0
        ditregrid=2
        max_blocks=40000
 /

 &paramlist
        typecourant='maxsum'
        courantpar=0.8d0
 /

 &mhd_list
        mhd_energy=.true.
        mhd_internal_e=.true.
        mhd_viscosity=.false.
        mhd_gravity=.true.
        mhd_gamma=1.05d0
        typedivbfix='ct'
/

 &usr_list
        f_q=1.0d4
        f_d=2.5d0
        f_L=0.5d0
        amr_r_core_max=4.0d0
        amr_r_sheet_max=8.0d0
        amr_theta_guard=0.02d0
        amr_sheet_halfwidth=0.35d0
        amr_br_ratio_max=0.20d0
        amr_force_max_level=2
        amr_outer_relax_frac=0.25d0
        amr_theta_relax_frac=0.10d0
/
```

- [ ] **Step 3: Run the new init file before changing `mod_usr.t` to verify the namelist fails for the new AMR keys**

Run:

```bash
cd /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched
OMP_NUM_THREADS=1 ./amrvac -i hpc_pilot_init_offpole.par
```

Expected: FAIL during parameter read because `amr_r_core_max` and the other new `&usr_list` entries are not defined yet.

- [ ] **Step 4: Commit the parameter-file skeletons**

```bash
git add \
  /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_init_offpole.par \
  /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_offpole.par
git commit -m "test: add hpc pilot parameter skeletons"
```

### Task 2: Implement The Hybrid AMR Controls In `mod_usr.t`

**Files:**
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t:5-12`
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t:53-67`
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t:791-806`
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t:901-931`
- Test: `amrvac_polytropic/analytic_bipolar_1to20_stretched/hpc_pilot_init_offpole.par`

- [ ] **Step 1: Extend the AMR control declarations near the top of `mod_usr.t`**

Replace the declaration block near the existing `f_q`, `f_d`, `f_L` fields with:

```fortran
  double precision :: f_q, f_d, f_L
  double precision :: amr_r_core_max, amr_r_sheet_max, amr_theta_guard
  double precision :: amr_sheet_halfwidth, amr_br_ratio_max
  double precision :: amr_outer_relax_frac, amr_theta_relax_frac
  integer :: amr_force_max_level
```

- [ ] **Step 2: Make `usr_params_read()` accept the new `&usr_list` keys and give them safe defaults**

Replace the current namelist and default block in `usr_params_read()` with:

```fortran
    namelist /usr_list/ f_q, f_d, f_L, &
         amr_r_core_max, amr_r_sheet_max, amr_theta_guard, &
         amr_sheet_halfwidth, amr_br_ratio_max, amr_force_max_level, &
         amr_outer_relax_frac, amr_theta_relax_frac

    f_q=1.d0
    f_d=1.d0
    f_L=1.d0
    amr_r_core_max=4.d0
    amr_r_sheet_max=8.d0
    amr_theta_guard=0.02d0
    amr_sheet_halfwidth=0.35d0
    amr_br_ratio_max=0.20d0
    amr_force_max_level=1
    amr_outer_relax_frac=0.25d0
    amr_theta_relax_frac=0.10d0
```

- [ ] **Step 3: Replace the blanket `qt==0` refine-all rule with the hybrid inner-corona and sheet-window rule**

Replace `special_refine_grid()` with:

```fortran
  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: rmin, thetamin, thetamax, thetac, br_ratio_min
    double precision :: bmag(ixI^S), br_ratio(ixI^S)
    logical :: safe_theta, force_inner, force_sheet

    refine=0
    coarsen=0

    rmin=minval(x(ixO^S,1))
    thetamin=minval(x(ixO^S,2))
    thetamax=maxval(x(ixO^S,2))
    thetac=0.5d0*(thetamin+thetamax)

    safe_theta = (thetamin >= xprobmin2 + amr_theta_guard) .and. &
                 (thetamax <= xprobmax2 - amr_theta_guard)

    bmag(ixI^S)=sqrt(w(ixI^S,mag(1))**2 + w(ixI^S,mag(2))**2 + w(ixI^S,mag(3))**2)
    br_ratio(ixI^S)=abs(w(ixI^S,mag(1))) / (bmag(ixI^S) + smalldouble)
    br_ratio_min=minval(br_ratio(ixO^S))

    force_inner = safe_theta .and. (rmin <= amr_r_core_max)
    force_sheet = safe_theta .and. (rmin <= amr_r_sheet_max) .and. &
                  (abs(thetac-0.5d0*dpi) <= amr_sheet_halfwidth) .and. &
                  (br_ratio_min <= amr_br_ratio_max)

    if(level < amr_force_max_level .and. (force_inner .or. force_sheet)) then
      refine=1
      coarsen=-1
    end if

  end subroutine special_refine_grid
```

- [ ] **Step 4: Narrow `specialthreshold()` so it suppresses AMR near the theta cuts and the outer shell only**

Replace the body of `specialthreshold()` with:

```fortran
  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level

    double precision :: theta_zone, outer_zone, tol_add, dist_theta_cut

    tol_add=0.d0
    theta_zone=max(amr_theta_guard, amr_theta_relax_frac*(xprobmax2-xprobmin2))
    outer_zone=amr_outer_relax_frac*(xprobmax1-xprobmin1)
    dist_theta_cut=min(xlocal(2)-xprobmin2, xprobmax2-xlocal(2))

    if(dist_theta_cut < theta_zone) then
      tol_add=tol_add + (1.d0-dist_theta_cut/theta_zone)*0.5d0
    end if

    if(xprobmax1-xlocal(1) < outer_zone) then
      tol_add=tol_add + (1.d0-(xprobmax1-xlocal(1))/outer_zone)*0.5d0
    end if

    tolerance=tolerance+tol_add

  end subroutine specialthreshold
```

- [ ] **Step 5: Rebuild AMRVAC and rerun the init file to verify the new namelist and refinement logic compile cleanly**

Run:

```bash
cd /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 make
OMP_NUM_THREADS=1 ./amrvac -i hpc_pilot_init_offpole.par
```

Expected:
- `make` exits `0`
- the init run exits `0`
- `data/hpc_pilot_init_offpole.log` exists
- no `negative`, `NaN`, `Inf`, or `mpistop` appears in stdout/log

- [ ] **Step 6: Commit the AMR implementation**

```bash
git add /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/mod_usr.t
git commit -m "feat: add configurable off-pole hpc amr controls"
```

### Task 3: Add The HPC Runner, LSF Job, And README Handoff

**Files:**
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/run_hpc_pilot.sh`
- Create: `amrvac_polytropic/analytic_bipolar_1to20_stretched/main_hpc_pilot.job`
- Modify: `amrvac_polytropic/analytic_bipolar_1to20_stretched/README.md`
- Test: `amrvac_polytropic/analytic_bipolar_1to20_stretched/run_hpc_pilot.sh`

- [ ] **Step 1: Add a dedicated runner for the new HPC pilot stages**

Create `run_hpc_pilot.sh` with:

```bash
#!/bin/bash
set -euo pipefail

case_dir="/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched"
cd "$case_dir"

stages=(hpc_pilot_init_offpole hpc_pilot_offpole)

if [[ $# -gt 2 ]]; then
  echo "Usage: $0 [all|hpc_pilot_init_offpole|hpc_pilot_offpole] [np]"
  exit 2
fi

target="${1:-all}"
np="${2:-${LSB_DJOB_NUMPROC:-32}}"

if [[ "$target" != "all" ]]; then
  stages=("$target")
fi

run_stage() {
  local stage="$1"
  local par_file="${stage}.par"
  local stdout_file="data/${stage}.stdout"
  local log_file="data/${stage}.log"

  echo "== running ${stage} with np=${np} =="
  rm -f "data/${stage}"*.dat "data/${stage}"*.vtu "data/${stage}.log" "data/${stage}.stdout"

  OMP_NUM_THREADS=1 mpirun -np "$np" ./amrvac -i "$par_file" > "$stdout_file" 2>&1

  if [[ ! -f "$log_file" ]]; then
    echo "Missing log file for ${stage}: ${log_file}"
    return 1
  fi

  if grep -Eqi '(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|negative|dtmin|abort|mpistop|fatal|segmentation fault' "$stdout_file" "$log_file"; then
    echo "Detected a hard-failure marker in ${stage}; inspect ${stdout_file} and ${log_file}."
    return 1
  fi

  echo "Last log line for ${stage}:"
  tail -n 1 "$log_file"
  echo
}

for stage in "${stages[@]}"; do
  run_stage "$stage"
done
```

- [ ] **Step 2: Add an LSF submission script for the cheap first HPC attempt**

Create `main_hpc_pilot.job` with:

```bash
#!/bin/bash
#BSUB -q mpi
#BSUB -J bipolar_hpc_pilot
#BSUB -n 32
#BSUB -o ./data/hpc_pilot.out
#BSUB -e ./data/hpc_pilot.err

set -euo pipefail

cd /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched
echo -n "==job started:"
date
echo
./run_hpc_pilot.sh all "${LSB_DJOB_NUMPROC:-32}"
echo
echo -n "==job finished:"
date
echo
```

- [ ] **Step 3: Update `README.md` so the next operator can see the PC ladder versus the new HPC pilot**

Append this section to `README.md`:

```markdown
## Off-pole HPC pilot

The first HPC campaign should use the off-pole branch only. It keeps the current theta-boundary implementation, raises the base mesh from `72 x 72 x 72` to `96 x 96 x 96`, and increases `refine_max_level` from `1` to `2`.

The AMR policy is intentionally not a full-domain refine-all rule. `special_refine_grid()` should force extra resolution only when both of these conditions are met:

- the block is safely away from the truncated theta caps by `amr_theta_guard`
- the block is either inside the inner-corona window (`r <= amr_r_core_max`) or inside the equatorial sheet window (`r <= amr_r_sheet_max`, `|theta-pi/2| <= amr_sheet_halfwidth`) and contains low `|Br|/|B|`

The new files are:

- `hpc_pilot_init_offpole.par`: initialization-only verification of the new mesh and AMR controls
- `hpc_pilot_offpole.par`: the actual cheap first HPC pilot
- `run_hpc_pilot.sh`: runs `init`, `pilot`, or `all`
- `main_hpc_pilot.job`: LSF submission script with `32` MPI ranks by default

Local verification:

```bash
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 make
OMP_NUM_THREADS=1 ./amrvac -i hpc_pilot_init_offpole.par
```

Cluster submission:

```bash
bsub < main_hpc_pilot.job
```
```

- [ ] **Step 4: Verify the new runner and job file are syntactically clean**

Run:

```bash
cd /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched
bash -n run_hpc_pilot.sh
bash -n main_hpc_pilot.job
./run_hpc_pilot.sh hpc_pilot_init_offpole 1
```

Expected:
- both `bash -n` commands exit `0`
- the one-rank init stage exits `0`
- `data/hpc_pilot_init_offpole.stdout` and `data/hpc_pilot_init_offpole.log` are created

- [ ] **Step 5: Commit the HPC handoff files**

```bash
git add \
  /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/run_hpc_pilot.sh \
  /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/main_hpc_pilot.job \
  /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched/README.md
git commit -m "docs: add off-pole hpc pilot runbook"
```

## Final Verification

- Rebuild once more:

```bash
cd /Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/analytic_bipolar_1to20_stretched
AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac-amrvac3.1 make
```

- Run the complete local preflight:

```bash
./run_hpc_pilot.sh all 1
```

- Inspect the outputs:

```bash
grep -Eqi '(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|negative|dtmin|abort|mpistop|fatal|segmentation fault' \
  data/hpc_pilot_init_offpole.stdout data/hpc_pilot_init_offpole.log \
  data/hpc_pilot_offpole.stdout data/hpc_pilot_offpole.log
```

Expected:
- build exits `0`
- `init` exits `0`
- `pilot` exits `0`; if it exits non-zero or writes any hard-failure marker, stop here and do **not** submit `main_hpc_pilot.job`
- the `grep` command returns non-zero, meaning no hard-failure markers were found

## Self-Review Notes

- Spec coverage: this plan covers the three requested changes directly: higher mesh, modest AMR increase, and a specific `special_refine_grid()` policy.
- Placeholder scan: no `TODO`, `TBD`, or deferred policy decisions remain.
- Type consistency: the AMR parameter names used in the `.par` files match the `usr_params_read()` namelist and the logic in both `special_refine_grid()` and `specialthreshold()`.
