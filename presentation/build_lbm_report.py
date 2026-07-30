"""Build an executable, results-first Jupyter report for LBM Simple Solver."""

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "presentation" / "lbm_results_report.ipynb"


def markdown(source: str):
    """Create a report markdown cell."""
    return nbf.v4.new_markdown_cell(source.strip())


def code(source: str):
    """Create a report code cell."""
    return nbf.v4.new_code_cell(source.strip())


nb = nbf.v4.new_notebook()
nb.metadata.update(
    {
        "kernelspec": {
            "display_name": "Python 3 (meep)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.13"},
    }
)

nb.cells = [
    markdown(
        r"""
# Results Report: Parallel D2Q9 Lattice Boltzmann Solver

**LBM Simple Solver · coarray Fortran · BGK collision · D2Q9 lattice**

This report evaluates two questions:

1. Does the solver reproduce the expected physical behaviour?
2. How well does the current coarray implementation perform and scale?

The evidence is presented before the numerical theory. All tables and plots below
are generated from archived result files in this repository, so the notebook can be
rerun when new measurements are added.
"""
    ),
    markdown(
        r"""
## Executive summary

The solver passes its main validation checks and reaches a measured peak of
**3.16 GLUPS** in the available cluster data. The most important findings are:

- The shear-wave sweep recovers the imposed kinematic viscosity to within
  **0.0264% in the worst case**, well inside the 1% validation threshold.
- Couette and Poiseuille flow reproduce the expected linear and parabolic velocity
  profiles. The converged Poiseuille run has a relative L2 error of
  **$1.11\times10^{-4}$**.
- Four of five 300×300 lid-driven-cavity runs converge to a residual below
  $10^{-8}$. The Re=5000 case reaches its 1.6-million-step limit and must be
  treated as a bounded, non-converged result.
- The largest recorded median throughput is **3157 MLUPS** for 1200×1200 on
  256 coarray images. A 2000×2000 grid reaches **2596 MLUPS** at the same image
  count.
- The 300×300 problem has already saturated: throughput falls from
  **790.7 MLUPS at 128 images to 768.3 MLUPS at 256**.
- Several repeated measurements have wide ranges, and some very short parallel
  runs finish in only a few seconds. Longer timed regions and a placement study
  are the highest-value next experiments.
"""
    ),
    code(
        r"""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

root = Path.cwd()
if root.name == "presentation":
    root = root.parent

plt.style.use("seaborn-v0_8-whitegrid")
COLORS = ["#176B87", "#2A9D8F", "#E9C46A", "#F4A261", "#D1495B"]

def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

perf = read_csv(
    root / "cluster_runs" / "performance_figures_latest" / "performance_summary.csv"
)
viscosity = read_csv(
    root / "shear_viscosity_sweep_300x300" / "figures" / "viscosity_summary.csv"
)
lid = read_csv(root / "lid_reynolds_results_300x300" / "summary.csv")

len(perf), len(viscosity), len(lid)
"""
    ),
    markdown(
        r"""
## 1. Cluster performance

The benchmark path runs a fixed 5,000 steps with convergence checks, diagnostics,
and field output disabled. The reported wall time is the maximum across images, and
throughput is

$$
\mathrm{MLUPS} = \frac{N_xN_yN_{\mathrm{steps}}}{10^6t_{\mathrm{wall}}}.
$$

Each plotted point is the median of three runs; the shaded band spans the measured
minimum and maximum.
"""
    ),
    code(
        r"""
sizes = sorted({(int(row["nx"]), int(row["ny"])) for row in perf})
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

for index, (nx, ny) in enumerate(sizes):
    rows = sorted(
        (row for row in perf if int(row["nx"]) == nx and int(row["ny"]) == ny),
        key=lambda row: int(row["images"]),
    )
    images = np.array([int(row["images"]) for row in rows])
    mlups = np.array([float(row["mlups"]) for row in rows])
    low = np.array([float(row["mlups_min"]) for row in rows])
    high = np.array([float(row["mlups_max"]) for row in rows])
    efficiency = np.array([float(row["efficiency_percent"]) for row in rows])
    color = COLORS[index % len(COLORS)]
    label = f"{nx}×{ny}"
    axes[0].plot(images, mlups, "o-", color=color, label=label)
    axes[0].fill_between(images, low, high, color=color, alpha=0.12)
    axes[1].plot(images, efficiency, "o-", color=color, label=label)

for ax in axes:
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted({int(row["images"]) for row in perf}))
    ax.set_xticklabels(sorted({int(row["images"]) for row in perf}))
    ax.set_xlabel("Coarray images")
axes[0].set_yscale("log")
axes[0].set_ylabel("Median throughput [MLUPS]")
axes[0].set_title("Throughput and repeat range")
axes[1].axhline(100, color="0.35", ls="--", lw=1.2, label="ideal efficiency")
axes[1].set_ylabel("Parallel efficiency [%]")
axes[1].set_title("Efficiency relative to each one-image baseline")
axes[0].legend(title="Grid")
axes[1].legend(title="Grid", ncol=2)
fig.tight_layout()
plt.show()
"""
    ),
    code(
        r"""
peak = max(perf, key=lambda row: float(row["mlups"]))
rows_256 = sorted(
    (row for row in perf if int(row["images"]) == 256),
    key=lambda row: int(row["nx"]),
)

print(
    f"Peak: {float(peak['mlups']):.1f} MLUPS on "
    f"{peak['nx']}×{peak['ny']} with {peak['images']} images."
)
print("\n256-image comparison")
print("grid       runtime [s]   MLUPS    efficiency")
for row in rows_256:
    print(
        f"{row['nx']}×{row['ny']:<5} "
        f"{float(row['wall_seconds']):11.3f} "
        f"{float(row['mlups']):8.1f} "
        f"{float(row['efficiency_percent']):9.1f}%"
    )
"""
    ),
    markdown(
        r"""
### Performance interpretation

The data show a clear problem-size effect. A 300×300 grid gives each image too
little work at high image counts, so halo synchronization and communication dominate.
The larger grids continue to benefit from additional images through 256.

Efficiencies above 100% should not be read as physically superlinear scaling without
qualification. The one-image baselines and several intermediate points contain large
run-to-run variation, while cache effects can also improve per-core performance after
decomposition. In particular, some high-image-count runs last only 1–8 seconds, making
startup effects and system noise a significant fraction of the measurement.
"""
    ),
    code(
        r"""
variability = []
for row in perf:
    median = float(row["mlups"])
    spread = (float(row["mlups_max"]) - float(row["mlups_min"])) / median
    variability.append((spread, row))

print("Five widest repeat ranges, normalized by median MLUPS")
print("grid       images   runtime [s]   relative range")
for spread, row in sorted(variability, reverse=True)[:5]:
    print(
        f"{row['nx']}×{row['ny']:<5} {int(row['images']):6d} "
        f"{float(row['wall_seconds']):11.3f} {100*spread:13.1f}%"
    )
"""
    ),
    markdown(
        r"""
## 2. Physical validation

### 2.1 Shear-wave viscosity

A periodic sinusoidal velocity mode should decay according to
$u_x(t)=u_0\exp(-\nu k^2t)$. Fitting the simulated decay therefore measures the
viscosity independently of the input relation
$\nu=\frac{1}{3}(1/\omega-1/2)$.
"""
    ),
    code(
        r"""
omega = np.array([float(row["omega"]) for row in viscosity])
applied = np.array([float(row["applied_viscosity"]) for row in viscosity])
measured = np.array([float(row["measured_viscosity"]) for row in viscosity])
error_pct = 100 * np.array([float(row["relative_error"]) for row in viscosity])

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5))
axes[0].plot(applied, applied, "--", color="0.35", label="exact agreement")
axes[0].plot(applied, measured, "o", ms=7, color=COLORS[0], label="measured")
axes[0].set(
    xlabel="Applied viscosity",
    ylabel="Measured viscosity",
    title="Recovered viscosity",
)
axes[0].legend()
axes[1].bar([f"{value:.1f}" for value in omega], error_pct, color=COLORS[1])
axes[1].axhline(1.0, color=COLORS[4], ls="--", label="1% threshold")
axes[1].set_yscale("log")
axes[1].set(
    xlabel=r"Relaxation rate $\omega$",
    ylabel="Relative error [%]",
    title="All sweep points pass",
)
axes[1].legend()
fig.tight_layout()
plt.show()

print(f"Worst relative viscosity error: {error_pct.max():.5f}%")
"""
    ),
    markdown(
        r"""
### 2.2 Channel-flow profiles

Couette flow tests the moving-wall boundary and should produce a linear velocity
profile. Poiseuille flow tests stationary bounce-back plus Guo body forcing and
should produce a parabola. Agreement in both cases checks not only collision and
streaming, but also the physical half-way wall location.
"""
    ),
    code(
        r"""
couette = np.loadtxt(root / "couette_profile.txt", comments="#")
poiseuille = np.loadtxt(
    root / "local_poiseuille_output_20260716_101555" / "poiseuille_profile.txt",
    comments="#",
)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.7))
for ax, data, title in [
    (axes[0], couette, "Couette flow"),
    (axes[1], poiseuille, "Poiseuille flow"),
]:
    ax.plot(data[:, 1], data[:, 0], lw=3, color=COLORS[0], label="numerical")
    ax.plot(data[:, 2], data[:, 0], "--", lw=2, color=COLORS[4], label="analytical")
    ax.set(xlabel=r"$u_x$", ylabel="y", title=title)
    ax.legend()
fig.tight_layout()
plt.show()
"""
    ),
    markdown(
        r"""
### 2.3 Lid-driven cavity convergence

Increasing Reynolds number reduces viscosity and greatly increases the time required
to reach a steady residual. Re=5000 did not meet the convergence criterion within the
step limit, so its flow field is useful as a finite-time result but not as a converged
steady benchmark.
"""
    ),
    code(
        r"""
reynolds = np.array([int(row["reynolds"]) for row in lid])
steps = np.array([int(row["completed_steps"]) for row in lid])
residual = np.array([float(row["final_residual"]) for row in lid])
converged = np.array([row["converged"].lower() == "true" for row in lid])

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.5))
bar_colors = [COLORS[1] if value else COLORS[4] for value in converged]
axes[0].bar([str(value) for value in reynolds], steps / 1e6, color=bar_colors)
axes[0].set(
    xlabel="Reynolds number",
    ylabel="Completed steps [millions]",
    title="Cost of reaching steady state",
)
axes[1].semilogy(reynolds, residual, "o-", color=COLORS[0])
axes[1].axhline(1e-8, color=COLORS[4], ls="--", label="convergence criterion")
axes[1].set(
    xlabel="Reynolds number",
    ylabel="Final velocity residual",
    title="Final residual",
)
axes[1].legend()
fig.tight_layout()
plt.show()
"""
    ),
    markdown(
        r"""
## 3. What was tested

| Case | Principal test | Result |
|---|---|---|
| Shear wave | Viscous decay and $\omega\leftrightarrow\nu$ relation | Worst viscosity error 0.0264% |
| Density wave | Propagation, damping, and mass conservation | Mean-density error below $10^{-10}$ at saved slices |
| Couette | Moving wall and half-way bounce-back | Numerical and linear profiles agree |
| Poiseuille | Body force and stationary walls | Relative L2 error $1.11\times10^{-4}$ |
| Lid cavity | Four-wall boundary interaction and steady convergence | Re≤2500 converged; Re=5000 bounded |

Together these cases exercise periodic and physical boundaries, moving and stationary
walls, body forcing, distributed halo exchange, convergence checks, and global output
assembly.
"""
    ),
    markdown(
        r"""
## 4. Numerical method in brief

The solver stores nine D2Q9 populations at each lattice cell. A time step consists of:

1. computing density and velocity;
2. local BGK collision, with Guo forcing for Poiseuille flow;
3. exchanging one-cell post-collision halos between coarray images;
4. pull streaming from local or halo source cells;
5. reconstructing missing wall populations with half-way bounce-back.

The global grid is decomposed into a near-square image grid. Pull streaming gives each
destination cell a single local writer, while remote source values are read from halo
cells. This keeps the time-step kernel uniform across the interior and image boundaries.

The kinematic viscosity is controlled by the relaxation rate:

$$
\nu = \frac{1}{3}\left(\frac{1}{\omega}-\frac{1}{2}\right), \qquad 0<\omega<2.
$$

These details are sufficient to interpret the validation results; a fuller derivation
and annotated Fortran excerpts remain in `lbm_code_presentation.ipynb`.
"""
    ),
    markdown(
        r"""
## 5. Limitations

- Strong-scaling curves currently stop at 256 images and use only three repeats.
- The same 5,000-step duration produces very different wall times across grid sizes
  and image counts; the shortest measurements are vulnerable to startup and scheduler
  noise.
- No completed fixed-image placement study is present, so the cost of spreading the
  same number of images across more nodes is not yet quantified.
- The performance result measures the present implementation, including separate
  collision/streaming arrays and a full `f = f_new` copy. It is not a hardware roofline.
- The Re=5000 cavity case is not a converged steady solution.
- There is no external reference-profile comparison for the lid cavity (for example,
  centerline velocity data), so the lid result is presently a convergence and qualitative
  flow-structure check.
"""
    ),
    markdown(
        r"""
## 6. Recommended next cluster measurements

If cluster time is limited, run these in order:

1. **Fixed-image placement study (highest value).** Use 384 total images on a
   5000×5000 or 10000×10000 grid, spread across 16, 24, and 48 nodes, with at least
   three exclusive repeats. These placements use 24, 16, and 8 images per node and
   directly measure the trade-off between memory
   bandwidth/contention within a node and network communication between nodes.
2. **Longer, cleaner timing at the scaling frontier.** Extend the timed region until
   each sample lasts roughly 20–60 seconds. Repeat 5–7 times, especially at 128, 256,
   and 384 images. This will make efficiency claims much more defensible.
3. **Large-grid strong scaling.** Measure one larger fixed problem at 64, 128, 256,
   and 384 images. The existing 300×300 case is already saturated and should not be
   extended further.
4. **Cluster A/B test of halo synchronization.** Compare the old `sync all` version
   with the current neighbor-scoped `sync images` version using the same allocation
   and randomized run order. The local result suggests a benefit, but it does not
   establish the multi-node effect.

The repository already contains `submit_lid_scaling2.sh` and
`lid_performance2.slurm` for the first experiment. A dry run is shown in the appendix.
"""
    ),
    markdown(
        r"""
## Appendix: reproducibility

Rebuild and execute this report from the repository root:

```bash
conda run -n meep python presentation/build_lbm_report.py
conda run -n meep jupyter nbconvert \
  --to notebook --execute --inplace \
  presentation/lbm_results_report.ipynb
```

Preview the recommended placement submissions without sending jobs:

```bash
DRY_RUN=1 CASE=lid IMAGES=384 NODE_COUNTS="16 24 48" \
  SIZES="5000x5000 10000x10000" STEPS=5000 REPEATS=3 \
  bash submit_lid_scaling2.sh
```

After new jobs finish, regenerate the performance and placement summaries:

```bash
python3 performance_plots.py --data-dir cluster_runs \
  --output-dir cluster_runs/performance_figures_latest
python3 placement_plots.py --data-dir cluster_runs
```
"""
    ),
]

OUT.write_text(nbf.writes(nb), encoding="utf-8")
print(f"Wrote {OUT}")
