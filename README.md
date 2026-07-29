# LBM Simple Solver

This repository contains a coarray Fortran D2Q9 solver for five cases:
`shear`, `density`, `couette`, `poiseuille`, and `lid`.

## Requirements

For local runs, install:

- Bash, `awk`, and a Fortran compiler;
- GNU Fortran for one-image runs, Intel `ifx` for shared-memory multi-image
  runs, or OpenCoarrays (`caf` and `cafrun`) for portable multi-image runs;
- Fortran Package Manager (`fpm`) only if you prefer to use the project
  manifest directly;
- Python 3 with NumPy and Matplotlib only if plots or output validation are
  required.

If Intel oneAPI is installed but `ifx` is not yet in `PATH`, initialize it
before running a script:

```bash
source /opt/intel/oneapi/setvars.sh
```

## Test the code

Run the Fortran unit tests:

```bash
bash test.sh
```

Run the unit tests followed by a 20-step smoke test of all five cases:

```bash
SMOKE=1 bash test.sh
```

GNU Fortran is selected automatically when Intel Fortran is unavailable. A
compiler can also be selected explicitly:

```bash
FC=gfortran bash test.sh
FC=ifx bash test.sh
```

Validate result files from a completed physical run:

```bash
python3 test/validate_outputs.py --data-dir runs/my_run/shear
```

## Run locally

The common interface is:

```text
bash run_local.sh CASE [OMEGA] [WALL_VELOCITY] [MAX_STEPS]
```

For example:

```bash
bash run_local.sh shear
bash run_local.sh density 1.2
bash run_local.sh couette 1.0 0.05
bash run_local.sh poiseuille 1.0
bash run_local.sh lid 1.0 0.0555555555555556
```

The scripts create a time-stamped directory below `runs/`. Set environment
variables to change the build or output:

```bash
NX=300 NY=300 IMAGES=4 FC=ifx RUN_DIR=runs/lid_300 \
  bash run_local.sh lid 1.0 0.0555555555555556
```

Supported settings are:

- `NX`, `NY`: lattice dimensions, default `128 x 128`;
- `IMAGES`: coarray images, default `1`;
- `FC`: `ifx`, `caf`, or `gfortran`, detected automatically;
- `BUILD_TYPE`: `release` or `debug`;
- `RUN_DIR`: destination for the executable, log, and result files.

Plain GNU Fortran uses `-fcoarray=single` and therefore supports one image.
Use Intel Fortran or OpenCoarrays for `IMAGES` greater than one.

To run every case with the same build, use `all`. Supplying `MAX_STEPS` makes
each case run exactly that many steps:

```bash
NX=64 NY=64 bash run_local.sh all 1.0 0.02 100
```

Omit `MAX_STEPS` for a physical run with normal diagnostics and convergence
checks. To cap such a run without disabling those checks, set
`LBM_STEP_LIMIT`:

```bash
LBM_STEP_LIMIT=100000 bash run_local.sh lid 1.0 0.05
```

Convenience drivers are also provided:

```bash
bash run_local_lid.sh
bash run_local_poiseuille.sh
bash run_local_shear_sweep.sh
bash run_lid_reynolds_sweep.sh
```

## Measure local performance

The fourth argument enables fixed-step timing. `PERFORMANCE_ONLY=1` also
removes field and diagnostic file output from the compiled executable:

```bash
NX=512 NY=512 IMAGES=4 FC=ifx PERFORMANCE_ONLY=1 \
  bash run_local.sh lid 1.0 0.05 5000
```

Each run writes `performance.csv` with the lattice size, image decomposition,
steps, wall time, and MLUPS.

## Run on a Slurm cluster

`bwunicluster.slurm` builds inside the allocation and supports the same five
case names. Request one Slurm task per coarray image:

```bash
sbatch --nodes=1 --ntasks=16 \
  --export=ALL,CASE=shear,NX=512,NY=512,NIMAGES=16 \
  bwunicluster.slurm
```

Run all cases in one allocation:

```bash
sbatch --nodes=1 --ntasks=16 \
  --export=ALL,CASE=all,NX=256,NY=256,NIMAGES=16,STEPS=1000 \
  bwunicluster.slurm
```

The job defaults to Intel distributed coarrays. It also accepts `FC=caf` for
OpenCoarrays and `FC=gfortran` for a one-task job. Cluster software names
differ, so export the modules or setup script required at the target site
before submission. For example:

```bash
export MODULES="compiler/intel mpi/impi"
export FC=ifx
sbatch --export=ALL --nodes=2 --ntasks=64 bwunicluster.slurm
```

Alternatively, set `ENV_SETUP` to a shell setup file. Other useful overrides
are `COARRAY_MODE`, `OMEGA`, `WALL_VELOCITY`, `STEPS`, `RUN_DIR`, and
`PERFORMANCE_ONLY`.

## Measure cluster performance

Submit one output-free fixed-step measurement:

```bash
sbatch --nodes=4 --ntasks=128 --ntasks-per-node=32 \
  --export=ALL,CASE=lid,NIMAGES=128,NODES=4,TASKS_PER_NODE=32,NX=1200,NY=1200,STEPS=5000 \
  lid_performance2.slurm
```

Submit repeated placement studies with the helper script:

```bash
DRY_RUN=1 CASE=lid IMAGES=512 NODE_COUNTS="16 32 64" \
  SIZES="5000x5000 10000x10000" bash submit_lid_scaling2.sh

CASE=lid IMAGES=512 NODE_COUNTS="16 32 64" \
  SIZES="5000x5000 10000x10000" bash submit_lid_scaling2.sh
```

The performance jobs write unique directories below `cluster_runs/` containing
`performance.csv`, `allocation.csv`, `placement.csv`, and `simulation.log`.

Create performance plots after the jobs finish:

```bash
python3 performance_plots.py --data-dir cluster_runs
python3 placement_plots.py --data-dir cluster_runs
```

Create figures for any available physical-case outputs:

```bash
python3 plots.py --data-dir runs/my_run --output-dir runs/my_run/figures
```
