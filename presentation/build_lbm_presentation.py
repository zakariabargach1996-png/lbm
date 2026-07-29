"""Build the executable Jupyter/Reveal.js presentation for LBM Simple Solver."""

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "presentation" / "lbm_code_presentation.ipynb"


def slide(source: str, kind: str = "markdown", slide_type: str = "slide"):
    """Create one notebook cell with Reveal.js slide metadata."""
    cell = (
        nbf.v4.new_markdown_cell(source)
        if kind == "markdown"
        else nbf.v4.new_code_cell(source)
    )
    cell.metadata["slideshow"] = {"slide_type": slide_type}
    if kind == "code":
        cell.metadata["jupyter"] = {"source_hidden": True}
    return cell


nb = nbf.v4.new_notebook()
nb.metadata.update(
    {
        "kernelspec": {
            "display_name": "Python 3 (meep)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.13"},
        "celltoolbar": "Slideshow",
        "rise": {
            "autolaunch": False,
            "scroll": True,
            "theme": "black",
            "transition": "fade",
        },
    }
)

cells = [
    slide(
        r"""
# Parallel D2Q9 Lattice Boltzmann Solver

### Algorithms, boundary conditions, validation, and performance

**LBM Simple Solver · Fortran coarrays · BGK + Guo forcing**

<small>Use **Space** / arrow keys in slideshow view. The notebook is also an executable report.</small>
"""
    ),
    slide(
        r"""
## What matters most

- A D2Q9, single-relaxation-time (BGK) solver supports five physical cases.
- **Pull streaming** makes every owned-cell write local and deterministic.
- One-cell **coarray halos** turn remote sources into ordinary local reads.
- **Half-way bounce-back** reconstructs only populations missing at physical walls.
- Validation compares against analytical decay/profile solutions and conservation laws.
- Fixed-step benchmark mode removes diagnostics and I/O from kernel timing.

> The central design idea: communicate post-collision edge values, then let every image
> execute the same local pull-streaming kernel.
"""
    ),
    slide(
        r"""
## D2Q9 state and lattice

Every lattice node stores nine populations $f_q$:

```
      7  3  6
       \ | /
    4 -- 1 -- 2
       / | \
      8  5  9
```

| Quantity | Definition |
|---|---|
| Density | $\rho = \sum_q f_q$ |
| Momentum | $\rho\mathbf{u} = \sum_q f_q\mathbf{c}_q$ |
| Weights | $4/9$ (rest), $1/9$ (axes), $1/36$ (diagonals) |
| Sound speed | $c_s^2=1/3$ |
| Viscosity | $\nu=c_s^2(1/\omega-1/2)$, with $0<\omega<2$ |

The arrays use the layout `f(q, x, y)` and a one-cell ghost layer on every side.
"""
    ),
    slide(
        r"""
## One complete time step

```text
f(t)
 │
 ├─ 1. moments: rho, ux, uy
 ├─ 2. local BGK collision + optional Guo force
 ├─ 3. exchange post-collision edge populations into halos
 ├─ 4. pull-stream from local or halo source cells
 ├─ 5. reconstruct missing wall populations (bounce-back)
 └─ 6. swap/promote f_new → f(t+1)
```

The ordering is essential: halo exchange must publish **new post-collision**
populations, not values from the previous time level.

Source: `app/main.f90`, main time-step loop.
"""
    ),
    slide(
        r"""
## Collision: equilibrium + forcing

The second-order low-Mach equilibrium is

$$
f_q^{eq}=w_q\rho\left[
1+\frac{\mathbf{c}_q\cdot\mathbf{u}}{c_s^2}
+\frac{(\mathbf{c}_q\cdot\mathbf{u})^2}{2c_s^4}
-\frac{\mathbf{u}^2}{2c_s^2}\right].
$$

BGK relaxes each population locally:

$$f_q^*=f_q-\omega(f_q-f_q^{eq})+(1-\omega/2)S_q.$$

- `S_q` is the Guo forcing term.
- The macroscopic velocity includes the matching half-force correction.
- Only Poiseuille uses a body force; all other cases use zero forcing.
- Collision has no communication and is therefore the natural optimization target.
"""
    ),
    slide(
        r"""
## Fortran excerpt: BGK collision

This is the central update inside `collide` (names shortened only by the slide layout):

```fortran
do jj = 1, ny_loc
  do ii = 1, nx_loc
    u2 = ux(ii,jj)**2 + uy(ii,jj)**2
    do qq = 1, q
      cu = real(cx(qq),dp)*ux(ii,jj) + real(cy(qq),dp)*uy(ii,jj)
      feq = w(qq)*rho(ii,jj) * (1.0_dp + cu/cs2 + &
        cu**2/(2.0_dp*cs2**2) - u2/(2.0_dp*cs2))

      f(qq,ii,jj) = f(qq,ii,jj) - &
        relaxation_omega*(f(qq,ii,jj)-feq) + &
        (1.0_dp-0.5_dp*relaxation_omega)*source
    end do
  end do
end do
```

The three nested loops expose the performance-critical structure: lattice cells outside,
nine populations inside, with no remote access.
""",
        slide_type="skip",
    ),
    slide(
        r"""
## Why pull streaming instead of push?

### Pull (implemented)

```fortran
f_next(q,i,j) = f_post(q,i-cx(q),j-cy(q))
```

- one writer owns every destination cell;
- contiguous destination traversal;
- remote values are already present in ghost layers;
- wall directions are visibly “missing” and repaired afterward.

### Push (alternative)

```fortran
f_next(q,i+cx(q),j+cy(q)) = f_post(q,i,j)
```

- scatters writes to neighboring destinations;
- awkward at image boundaries and easier to race in threaded variants;
- needs special handling or communication for remote destinations.

Pull does not reduce the mathematical work, but it simplifies ownership,
parallel correctness, and later OpenMP/SIMD work.
"""
    ),
    slide(
        r"""
## Fortran excerpt: pull streaming

The complete streaming kernel is deliberately small:

```fortran
f_next = 0.0_dp
do jj = 1, ny_loc
  do ii = 1, nx_loc
    do qq = 1, q
      f_next(qq,ii,jj) = &
        f_post(qq,ii-cx(qq),jj-cy(qq))
    end do
  end do
end do
```

Because `cx` and `cy` are only −1, 0, or +1:

- a one-cell halo is sufficient;
- the destination is always owned by the current image;
- a source may be local, a halo value, or absent at a physical wall.
""",
        slide_type="skip",
    ),
    slide(
        r"""
## Coarray domain decomposition

The global lattice is factored into a near-square `px × py` image grid.

- quotient/remainder partitioning supports non-divisible dimensions;
- `x_offset`, `y_offset` preserve global coordinates;
- `0` and `n_local+1` are ghost indices;
- only `f` is a coarray; macroscopic fields and `f_new` stay local;
- periodic global edges point to the image on the opposite side;
- physical walls use neighbor ID `0`.

Halo exchange proceeds in two stages:

1. left/right columns;
2. top/bottom rows **including the new column ghosts**, which propagates corners.

The exchange uses three dependency-matched `sync images` phases: horizontal
neighbors before column reads, all cardinal neighbors before row reads, and
vertical neighbors after row reads.
"""
    ),
    slide(
        r"""
## Fortran excerpt: coarray halo reads

Square brackets select storage exposed by another coarray image:

```fortran
sync images(horizontal_images)

if (left_img /= 0) then
  f(:,0,1:ny_loc) = &
    f(:,remote_nx,1:ny_loc)[left_img]
end if

if (right_img /= 0) then
  f(:,nx_loc+1,1:ny_loc) = &
    f(:,1,1:ny_loc)[right_img]
end if

sync images(cardinal_images)

if (up_img /= 0) then
  f(:,0:nx_loc+1,ny_loc+1) = &
    f(:,0:nx_loc+1,1)[up_img]
end if
```

The second-stage row includes column ghosts, so diagonal corner data travels without
a separate corner-message phase.
""",
        slide_type="skip",
    ),
    slide(
        r"""
## Stationary half-way bounce-back

After pull streaming, populations that would have come from outside a wall are unset.
They are reconstructed from the outgoing post-collision population at the same node:

$$f_q(\mathbf{x},t+\Delta t)=f_{\bar q}^{*}(\mathbf{x},t).$$

Example at the left wall:

```fortran
do jj = 1, ny_loc
  f_next(2,1,jj) = f_post(4,1,jj)  ! E  <- W
  f_next(6,1,jj) = f_post(8,1,jj)  ! NE <- SW
  f_next(9,1,jj) = f_post(7,1,jj)  ! SE <- NW
end do
```

- The physical wall lies halfway between the boundary fluid node and exterior node.
- This gives no-slip without storing solid nodes.
- Only images touching a global physical edge apply the boundary routine.
- Couette, Poiseuille, and lid side/bottom walls use this rule.
"""
    ),
    slide(
        r"""
## Moving-lid bounce-back

At the top lid, normal momentum is reflected and diagonal populations receive
tangential momentum:

```fortran
f_next(S ) = f_post(N)
f_next(SW) = f_post(NE) - rho_wall*u_wall/6
f_next(SE) = f_post(NW) + rho_wall*u_wall/6
```

For D2Q9, the moving-wall correction reduces to $\rho u_w/6$.

Implementation choices:

- $\rho_{wall}$ is estimated from local post-collision populations;
- top-wall handling runs after side walls, so it defines the corner convention;
- the simple corner treatment is stable here, but higher-Re work may benefit from
  a more carefully validated corner scheme.
"""
    ),
    slide(
        r"""
## Five cases, one solver

| Case | Initial/boundary condition | Physical check |
|---|---|---|
| Shear wave | sinusoidal $u_x(y)$, periodic | exponential viscous decay |
| Density wave | sinusoidal $\rho(x)$, periodic | propagation + mass conservation |
| Couette | stationary bottom, moving top | linear velocity profile |
| Poiseuille | stationary walls + body force | parabolic velocity profile |
| Lid cavity | moving top, three fixed walls | vortical flow + convergence |

The case selection changes initialization, boundaries, force, run length, and output;
the collision/communication/stream kernels remain common.
"""
    ),
    slide(
        r"""
## Shear-wave result: viscosity is recovered

The mode should decay as

$$u_x(t)=u_0\exp(-\nu k^2t), \qquad k=2\pi/N_y.$$

Across the 300×300 sweep ($\omega=0.6$ to $1.6$):

- measured relative viscosity error ranges from about $4\times10^{-10}$ to
  $2.63\times10^{-4}$;
- every point is far below the 1% validation threshold;
- the result directly verifies the $\omega\leftrightarrow\nu$ relation.

![Measured and theoretical shear-wave decay](../shear_viscosity_sweep_300x300/figures/measured_vs_theory_decay.png)
"""
    ),
    slide(
        r"""
## Density-wave result

The periodic 128×128 density-wave case was regenerated for this report:

- 1,109 time steps;
- $\rho(x,0)=1+0.01\sin(2\pi x/N_x)$;
- the disturbance propagates and damps while its spatial mean remains one.

The validation rule checks every saved time slice:

$$|\langle\rho\rangle-1| < 10^{-10}.$$
""",
    ),
    slide(
        r"""
## Viscosity controls density-wave damping

<div style="display:flex; gap:1.5rem; justify-content:center; align-items:flex-start;">
  <figure style="width:43%; margin:0;">
    <img src="../density_wavegifs/density.gif"
         alt="Density wave with faster viscous damping"
         style="width:100%; max-height:430px; object-fit:contain; margin:0;">
    <figcaption><strong>Higher viscosity</strong><br>
      The density contrast damps rapidly and the field becomes uniform.
    </figcaption>
  </figure>
  <figure style="width:43%; margin:0;">
    <img src="../density_wavegifs/density7.gif"
         alt="Density wave with slower viscous damping"
         style="width:100%; max-height:430px; object-fit:contain; margin:0;">
    <figcaption><strong>Lower viscosity</strong><br>
      The propagating wave remains visible for many more time steps.
    </figcaption>
  </figure>
</div>

<small>Both animations use the same density color scale (0.29–0.31); compare how
long spatial contrast survives.</small>
"""
    ),
    slide(
        r"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path.cwd()
if root.name == "presentation":
    root = root.parent
density = np.loadtxt(root / "presentation" / "density_profiles.txt", comments="#")
steps = np.unique(density[:, 0]).astype(int)
chosen = [steps[0], steps[len(steps)//3], steps[2*len(steps)//3], steps[-1]]
fig, ax = plt.subplots(figsize=(10, 4.6))
for t in chosen:
    block = density[density[:, 0] == t]
    ax.plot(block[:, 1], block[:, 2] - 1.0, label=f"t={t}")
ax.set(xlabel="x", ylabel=r"$\rho-1$", title="Density-wave evolution")
ax.grid(alpha=.25); ax.legend(ncol=2)
plt.show()
""",
        "code",
        "fragment",
    ),
    slide(
        r"""
## Channel-flow results

### Couette flow

- 128×128 root result;
- numerical profile follows $u_x(y)=u_w y/L$;
- half-way wall positions are $y=0$ and $y=L$.

### Poiseuille flow

- 64×64 local result;
- Guo body forcing recovers
  $u_x(y)=g_x y(L-y)/(2\nu)$;
- numerical/theoretical difference is visually negligible.

These two cases jointly validate moving-wall bounce-back, stationary bounce-back,
forcing, and the chosen physical wall location.
"""
    ),
    slide(
        r"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path.cwd()
if root.name == "presentation":
    root = root.parent
couette = np.loadtxt(root / "couette_profile.txt", comments="#")
poiseuille = np.loadtxt(
    root / "local_poiseuille_output_20260716_101555" / "poiseuille_profile.txt",
    comments="#",
)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, data, title in [
    (axes[0], couette, "Couette: linear profile"),
    (axes[1], poiseuille, "Poiseuille: parabolic profile"),
]:
    ax.plot(data[:, 1], data[:, 0], lw=3, label="numerical")
    ax.plot(data[:, 2], data[:, 0], "--", lw=2, label="analytical")
    ax.set(xlabel=r"$u_x$", ylabel="y", title=title)
    ax.grid(alpha=.25); ax.legend()
plt.tight_layout()
plt.show()
""",
        "code",
        "fragment",
    ),
    slide(
        r"""
## Poiseuille transient: zero flow → stable parabola

The 64×64, four-image run starts from rest and saves a full profile every 622 steps.
It converges at step **42,500**:

- final velocity residual: $8.59\times10^{-9}$;
- relative L2 profile error: $1.11\times10^{-4}$;
- maximum absolute error: $4.07\times10^{-6}$.

The curves below show representative saved profiles; the dashed curve is the steady
analytical solution. Early evolution is plug-like in the interior, while viscous
momentum diffusion enforces the wall condition and builds the parabola.
"""
    ),
    slide(
        r"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path.cwd()
if root.name == "presentation":
    root = root.parent
history = np.loadtxt(
    root / "local_poiseuille_output_20260716_101555" / "poiseuille_evolution.txt",
    comments="#",
)
saved_steps = np.unique(history[:, 0]).astype(int)
indices = np.unique(np.linspace(0, len(saved_steps)-1, 9).astype(int))
colors = plt.cm.viridis(np.linspace(0, 0.9, len(indices)))

fig, ax = plt.subplots(figsize=(10.5, 5.2))
for color, index in zip(colors, indices):
    t = saved_steps[index]
    block = history[history[:, 0] == t]
    ax.plot(block[:, 2], block[:, 1], color=color, lw=1.8, label=f"t={t:,}")
final = history[history[:, 0] == saved_steps[-1]]
ax.plot(final[:, 3], final[:, 1], "k--", lw=2.5, label="steady analytical")
ax.set(xlabel=r"$u_x$", ylabel="y", title="Poiseuille profiles approaching steady state")
ax.grid(alpha=.25)
ax.legend(ncol=2, fontsize=9)
plt.tight_layout()
plt.show()
""",
        "code",
        "fragment",
    ),
    slide(
        r"""
## Every saved Poiseuille profile

This animation steps through all **70 recorded profiles**, including $t=0$ and the
converged result at $t=42{,}500$. Use the playback controls below.
"""
    ),
    slide(
        r"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

root = Path.cwd()
if root.name == "presentation":
    root = root.parent
history = np.loadtxt(
    root / "local_poiseuille_output_20260716_101555" / "poiseuille_evolution.txt",
    comments="#",
)
saved_steps = np.unique(history[:, 0]).astype(int)
reference = history[history[:, 0] == saved_steps[-1]]

fig, ax = plt.subplots(figsize=(8.8, 5.0))
numerical_line, = ax.plot([], [], lw=3, color="#2a9d8f", label="numerical")
ax.plot(reference[:, 3], reference[:, 1], "k--", lw=2, label="steady analytical")
ax.set(xlim=(0, 0.052), ylim=(0, 64), xlabel=r"$u_x$", ylabel="y")
ax.grid(alpha=.25); ax.legend(loc="upper right")
title = ax.set_title("")

def update(frame):
    step = saved_steps[frame]
    block = history[history[:, 0] == step]
    numerical_line.set_data(block[:, 2], block[:, 1])
    title.set_text(f"Poiseuille transient · step {step:,} / {saved_steps[-1]:,}")
    return numerical_line, title

animation = FuncAnimation(fig, update, frames=len(saved_steps), interval=130, blit=False)
plt.close(fig)
HTML(animation.to_jshtml(fps=8, default_mode="once"))
""",
        "code",
        "fragment",
    ),
    slide(
        r"""
## Moving-lid Reynolds sweep

All runs use a 300×300 cavity and $u_{lid}=0.1$.

| Re | $\omega$ | Steps | Residual | Status |
|---:|---:|---:|---:|:---|
| 100 | 0.7143 | 114,500 | $9.73\times10^{-9}$ | converged |
| 500 | 1.4706 | 292,000 | $9.66\times10^{-9}$ | converged |
| 1,000 | 1.6949 | 488,500 | $9.96\times10^{-9}$ | converged |
| 2,500 | 1.8657 | 1,291,500 | $9.99\times10^{-9}$ | converged |
| 5,000 | 1.9305 | 1,600,000 | $6.74\times10^{-7}$ | hit limit |

Higher Reynolds number means lower viscosity, slower steady convergence, and
increasingly fine vortical structure. Re=5000 is a **bounded run**, not a converged result.
"""
    ),
    slide(
        r"""
## Lid flow field: Re = 1,000

![Moving-lid cavity evolution](../lid_reynolds_results_300x300/Re_1000/figures/moving_lid_evolution.png)

The streamplot shows the developing primary recirculation and corner structures.
The full field is assembled only for output; the time-stepping state stays distributed.
"""
    ),
    slide(
        r"""
## Performance methodology

The optional `MAX_STEPS` argument activates a reproducible benchmark path:

- exactly the requested number of steps;
- early convergence disabled;
- in-loop diagnostics and output disabled;
- global work count is exactly `NX × NY × steps`;
- wall time is the slowest-image time (`co_max`);
- throughput is reported as

$$\mathrm{MLUPS}=\frac{N_xN_yN_{steps}}{10^6t_{wall}}.$$

Cluster results use 5,000 steps and three repeats; the plotted central value is
the median. This isolates kernels + halo communication, not file I/O.
"""
    ),
    slide(
        r"""
## Cluster strong scaling

![Runtime, throughput, speedup, and efficiency](../cluster_runs/performance_figures_latest/performance_dashboard.png)

Highlights from the latest aggregate:

- peak median throughput is about **3,157 MLUPS** for 1200×1200 on 256 images;
- 2000×2000 reaches about **2,596 MLUPS** on 256 images;
- the 300×300 problem saturates: 790.7 MLUPS at 128 images, then 768.3 at 256;
- larger grids amortize communication better;
- superlinear efficiencies at some points signal baseline/run variability and cache effects,
  so they should not be interpreted as algorithmic scaling beyond 100%.
"""
    ),
    slide(
        r"""
## Local A/B test: halo synchronization

The same 300×300 moving-lid benchmark was run locally for 5,000 steps with
three repeats per image count. Values below are median MLUPS.

| Images | Three `sync all` barriers | Neighbor-scoped `sync images` | Median change |
|---:|---:|---:|---:|
| 1 | 18.62 | 17.11 | −8.1% |
| 2 | 32.79 | **38.49** | **+17.4%** |
| 3 | 56.11 | **57.58** | +2.6% |
| 4 | 72.27 | **73.64** | +1.9% |

- The neighbor-scoped version synchronizes left/right, then all cardinal
  neighbors, then up/down—matching the two-stage halo dependencies.
- One image has no communication partner, so its difference is timing noise,
  not a synchronization effect.
- Repeat ranges overlap at 2 and 4 images; the result supports a modest local
  improvement, but more randomized repeats are needed for a precise estimate.

> The main value is removing unrelated images from each halo synchronization;
> the benefit should be reassessed in the larger cluster scaling sweep.
"""
    ),
    slide(
        r"""
## Performance context: an optimized LBM framework

| System | Peak shown here | Execution and workload |
|---|---:|---|
| This solver | **3.157 GLUPS** | 256 CPU coarray images, double-precision D2Q9 |
| OpenLB 1.5 | **≈1.33 TLUPS** | 512 NVIDIA A100 GPUs, single-precision D3Q19 |

This is a **scale reference, not a head-to-head benchmark**:

- the lattice, precision, hardware, boundaries, and problem sizes differ;
- OpenLB uses mature CUDA, AVX2/AVX-512, MPI/OpenMP, and optimized streaming paths;
- the comparison shows the performance ceiling enabled by accelerator-aware data layout,
  vectorization, communication overlap, and years of framework optimization.

<small>Source: [OpenLB performance report, HoreKa 2022](https://www.openlb.net/performance/).
OpenLB reports ≈1.33 TLUPS on 512 A100 GPUs for single-precision D3Q19 BGK with
Periodic Shift streaming.</small>
"""
    ),
    slide(
        r"""
## Where time is likely going

Each lattice update moves two full population fields through memory:

- read 9 populations for moments/collision;
- write 9 collided populations;
- read 9 again during streaming;
- write 9 into `f_new`;
- then `f = f_new` copies the full buffer once more.

At scale, communication adds:

- edge/corner halo traffic proportional to subdomain perimeter;
- three neighbor-scoped synchronization phases per exchange;
- remote coarray reads and synchronization latency.

The implementation is therefore both **memory-bandwidth bound** locally and
**latency/synchronization sensitive** at high image counts.
"""
    ),
    slide(
        r"""
## Optimization roadmap

### Highest value, lowest risk

1. Replace `f = f_new` with pointer/allocatable buffer swapping.
2. Fuse moments + collision where practical; keep scalars in registers.
3. Add `do concurrent` / `!$omp simd` after checking compiler vectorization reports.
4. Benchmark loop order and population layout on the target CPU.

### Parallel scaling

5. Evaluate coarray events or asynchronous exchange beyond the current
   neighbor-scoped synchronization.
6. Pack/send only populations that actually cross each face, not all nine.
7. Overlap interior collision/streaming with halo communication.
8. Use shared-memory threading inside a node and fewer coarray images across nodes.

### Numerical robustness

9. Consider MRT/TRT collision for improved stability near $\omega\to2$.
10. Validate a refined lid-corner treatment at high Reynolds number.
11. Add a solid/fluid geometry mask and link-wise bounce-back; use interpolated
    bounce-back or an immersed-boundary method for curved surfaces.

Every change should be checked against the five-case validation suite and MLUPS benchmark.
"""
    ),
    slide(
        r"""
## Beyond BGK: what are TRT and MRT?

### BGK — one relaxation time

All non-conserved behavior relaxes with the same rate $\omega$. It is simple and fast,
but couples viscosity to the damping of less physical modes.

### TRT — two relaxation times

Populations are split into **symmetric** and **antisymmetric** pairs. One rate sets
viscosity; the second can tune wall accuracy and damp non-hydrodynamic error.
TRT adds modest cost and is often the practical next step.

### MRT — multiple relaxation times

Populations are transformed into physical moments such as density, momentum, stress,
and higher-order modes. Each non-conserved moment receives its own relaxation rate.
This gives more stability and control, especially near $\omega\to2$, at the cost of
extra transforms, parameters, and tuning.

**Both retain the same streaming and halo structure; only the local collision changes.**
"""
    ),
    slide(
        r"""
## Correctness evidence and caveats

### Evidence

- unit tests verify D2Q9 weights, moments, opposite directions, and omega range;
- shear viscosity agrees with theory to better than 0.03% in the 300×300 sweep;
- Couette and Poiseuille profiles match their analytical solutions;
- density-wave and lid checks cover mass conservation and finite fields;
- global coordinates and reductions are decomposition independent.

### Caveats

- BGK becomes fragile as $\omega$ approaches 2;
- Re=5000 did not satisfy the steady residual threshold within 1.6M steps;
- boundaries are limited to axis-aligned outer walls: there are no curved or internal
  solid geometries yet;
- scaling points show substantial run variability and apparent superlinearity;
- neighbor-scoped synchronization still incurs latency at high image counts.
"""
    ),
    slide(
        r"""
# Takeaways

1. **Pull streaming + one-cell halos** gives a clear ownership model and simple local kernel.
2. **Post-collision bounce-back** correctly fills only the populations missing at walls.
3. A single solver reproduces decay, wave, channel, and cavity physics.
4. Validation is strong for the analytical cases; the highest-Re lid result needs more runtime
   or a more robust collision model.
5. The clearest next wins are removing the full-array copy, improving vectorization, and overlapping
   useful work with the neighbor-aware communication.

### Run the presentation

```bash
conda run -n meep jupyter notebook presentation/lbm_code_presentation.ipynb
```

Or open the exported `presentation/lbm_code_presentation.slides.html`.
"""
    ),
]

nb["cells"] = cells
OUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, OUT)
print(f"Wrote {OUT} with {len(cells)} cells")
