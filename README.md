# LBM Simple Solver: two-dimensional lattice Boltzmann simulations

The executable supports five D2Q9/BGK cases: shear-wave decay, density-wave
propagation, Couette flow, force-driven Poiseuille flow, and a moving-lid
cavity. The command-line syntax is

```text
milestone4 CASE [OMEGA] [WALL_VELOCITY] [MAX_STEPS]
```

The implementation is separated into `lbm_params` (lattice constants),
`domain_decomposition` (process grid and halo exchange), `lbm_solver`
(initialization, collision, streaming, and boundaries), and `lbm_output`
(diagnostics and field/profile output). The main program contains the run
configuration and timestep orchestration.

`OMEGA` defaults to `1/0.7`, and the program rejects values outside
`0 < omega < 2`. The lattice viscosity is

```text
nu = (1/3) * (1/omega - 1/2).
```

For the lid cavity, the Reynolds number is printed at startup and is computed
as `Re = WALL_VELOCITY * NX / nu`. Compile that case on a 300 x 300 grid to
match the assignment.

The optional `MAX_STEPS` argument is intended for reproducible performance
runs: it forces exactly that many steps and disables early convergence and
in-loop diagnostic/output work, so the measurement covers the LBM kernels and
halo communication rather than file I/O. Every
run reports synchronized wall-clock time and MLUPS and appends one row to
`performance.csv`. MLUPS uses exactly `NX * NY * completed_steps`; ghost cells
are not counted.

Set `LBM_STEP_LIMIT` to cap a physical run while retaining residual checks,
diagnostic output, and early convergence. This is separate from positional
`MAX_STEPS`, which remains the fixed-step performance mode.

The two-dimensional decomposition supports arbitrary lattice dimensions. If a
dimension is not divisible by its process-grid dimension, the first blocks get
one extra cell, and global output coordinates are reconstructed from explicit
block offsets.

## Local runs

For a resource-conscious 300 x 300 moving-lid run, use the supplied local
driver. It uses four coarray images, lowers the process priority, defaults to
`Re=100`, creates a fresh output directory, and generates the streamplot after
the solver finishes:

```bash
bash run_local_lid.sh
```

The defaults can be changed without editing the script, for example
`IMAGES=2 OMEGA=1.0 WALL_VELOCITY=0.0555556 bash run_local_lid.sh`.

Run the standard 300 x 300 lid-cavity Reynolds-number sweep with:

```bash
bash run_lid_reynolds_sweep.sh
```

It uses `Re = 100 500 1000 2500 5000`, a lid speed of 0.1, and a
convergence-aware limit of 1,600,000 steps by default. Override the cap with,
for example, `STEP_LIMIT=500000 bash run_lid_reynolds_sweep.sh`.

The force-driven Poiseuille case has a separate resource-conscious local
driver. It defaults to a 64 x 64 grid and runs the solver, plot generation,
and output validation in one command:

```bash
bash run_local_poiseuille.sh
```

Its grid and runtime parameters can also be overridden, for example
`NX=128 NY=128 IMAGES=4 OMEGA=1.0 bash run_local_poiseuille.sh`.

With Intel oneAPI initialized, a four-image run is for example:

```bash
FOR_COARRAY_NUM_IMAGES=4 fpm run --compiler ifx \
  --flag '-coarray=shared -O3' -- shear 1.4
```

Run the shear case at several omega values to generate the measured-viscosity
points. Each run creates a separate `viscosity_omega_*.txt` file.

Generate every available report figure with the existing Conda environment:

```bash
mkdir -p /tmp/matplotlib
MPLCONFIGDIR=/tmp/matplotlib conda run -n meep \
  python plots.py --output-dir figures
```

The plotting script skips cases that have not yet been run. It generates:

- density and velocity evolution for both wave experiments;
- shear-wave amplitude decay and measured viscosity versus omega;
- transient Couette and Poiseuille profiles, analytical steady profiles, and
  the half-way wall positions; and
- streamplots at selected times for the moving-lid cavity.

For pressure/body-force-driven plane Poiseuille flow, the steady analytical
solution used in the output and plots is

```text
u_x(y) = g_x y (L-y) / (2 nu),   0 <= y <= L,
```

where `g_x` is the imposed acceleration, `L = NY`, and the lattice nodes are at
`y = j - 1/2`. Thus the physical walls are at `y=0` and `y=L`.

## Cluster runs

`bwunicluster.slurm` fixes the lid case at 300 x 300 and accepts `OMEGA` and
`WALL_VELOCITY` through `sbatch --export`, so the Reynolds number changes
without changing the box. Other cases also accept `NX` and `NY`. The script has
a viscosity sweep mode:

```bash
sbatch --export=ALL,CASE=viscosity-sweep,OMEGAS="0.6 0.8 1.0 1.2 1.4 1.6 1.7" \
  bwunicluster.slurm
```

It can run any individual case (`shear`, `density`, `couette`, `poiseuille`, or
`lid`). Set `STEPS` to request a fixed-step timed run.

For a fixed-image placement study, the supplied sweep keeps 512 coarray images
and a 10000x10000 lattice fixed while comparing 16, 32, and 64 nodes:

```bash
DRY_RUN=1 bash submit_lid_scaling2.sh
bash submit_lid_scaling2.sh
```

This corresponds to 32, 16, and 8 images per node. Three repetitions are
submitted for every placement. To use the smaller communication-sensitive
lattice instead, set `SIZES=5000x5000`. Each job writes `performance.csv`,
`allocation.csv`, and a combined `placement.csv` to a unique directory under
`cluster_runs/`. Full-field text output is disabled in these performance jobs
because it is prohibitively large at these lattice sizes.

After the jobs finish, compare node placements using the median of the repeats:

```bash
MPLCONFIGDIR=/tmp/matplotlib conda run -n meep \
  python placement_plots.py --data-dir cluster_runs
```

For the conventional strong-scaling results already present under
`cluster_runs/`, combine measurements into the image-count performance figures
with:

```bash
mkdir -p /tmp/matplotlib
MPLCONFIGDIR=/tmp/matplotlib conda run -n meep \
  python performance_plots.py --data-dir cluster_runs
```

This produces runtime, MLUPS, speedup, and parallel-efficiency plots, a combined
four-panel dashboard, and `performance_summary.csv` under
`cluster_runs/performance_figures/`. Repeated measurements at the same lattice
size and image count are combined using their median by default.

## Tests

`fpm test` checks the D2Q9 quadrature weights, isotropic moments, opposite
directions, and the default omega range. The simulation output itself is
validated by comparison with shear-wave decay and the analytical Couette and
Poiseuille profiles, as well as mass conservation and lid initial conditions:

```bash
conda run -n meep python test/validate_outputs.py --data-dir .
```
