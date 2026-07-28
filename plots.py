"""Create the report figures produced by the milestone-4 LBM simulations."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def selected_steps(data: np.ndarray, maximum: int = 7) -> np.ndarray:
    """Return evenly spaced available time steps, including first and last."""
    steps = np.unique(data[:, 0].astype(int))
    if len(steps) <= maximum:
        return steps
    indices = np.linspace(0, len(steps) - 1, maximum).round().astype(int)
    return steps[np.unique(indices)]


def plot_wave_profiles(
    data_file: Path, output_file: Path, title: str, coordinate: str
) -> None:
    """Plot density and longitudinal-velocity profiles at selected times."""
    data = np.loadtxt(data_file, comments="#")
    fig, (ax_rho, ax_u) = plt.subplots(
        2, 1, figsize=(7.2, 7.2), sharex=True, constrained_layout=True
    )
    for step in selected_steps(data):
        snapshot = data[data[:, 0] == step]
        label = rf"$t={step}$"
        ax_rho.plot(snapshot[:, 1], snapshot[:, 2] - 1.0, label=label)
        ax_u.plot(snapshot[:, 1], snapshot[:, 3], label=label)

    ax_rho.set_ylabel(r"Density perturbation $\rho-1$")
    ax_rho.set_title(title)
    ax_rho.grid(True, alpha=0.3)
    ax_rho.legend(ncol=2, fontsize="small")
    ax_u.set_xlabel(coordinate)
    ax_u.set_ylabel(r"Velocity $u_x$")
    ax_u.grid(True, alpha=0.3)
    fig.savefig(output_file, dpi=200)


def plot_shear_decay(data_file: Path, output_file: Path) -> None:
    """Plot numerical and theoretical shear-wave amplitude decay."""
    data = np.loadtxt(data_file, comments="#")
    step, numerical, theoretical, relative_error = data.T
    fig, (ax_decay, ax_error) = plt.subplots(
        2,
        1,
        figsize=(7.0, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        constrained_layout=True,
    )
    ax_decay.semilogy(step, theoretical, color="black", label="Analytical decay")
    ax_decay.semilogy(step, numerical, "o-", ms=3, label="LBM")
    ax_decay.set_ylabel("Shear-wave amplitude")
    ax_decay.set_title("Shear-wave decay")
    ax_decay.grid(True, which="both", alpha=0.3)
    ax_decay.legend()
    ax_error.semilogy(step, np.maximum(relative_error, 1e-16), color="tab:orange")
    ax_error.set_xlabel("Time step")
    ax_error.set_ylabel("Relative error")
    ax_error.grid(True, which="both", alpha=0.3)
    fig.savefig(output_file, dpi=200)


def plot_viscosity(files: list[Path], output_file: Path) -> None:
    """Plot measured viscosity versus omega and the BGK prediction."""
    measurements = np.vstack([np.loadtxt(path, comments="#", ndmin=2) for path in files])
    measurements = measurements[np.argsort(measurements[:, 0])]
    omega_min = max(0.05, 0.95 * measurements[:, 0].min())
    omega_max = min(1.7, 1.05 * measurements[:, 0].max())
    if np.isclose(omega_min, omega_max):
        omega_min = max(0.05, omega_min - 0.2)
        omega_max = min(1.7, omega_max + 0.2)
    omega_curve = np.linspace(omega_min, omega_max, 400)
    analytical_curve = (1.0 / omega_curve - 0.5) / 3.0

    fig, ax = plt.subplots(figsize=(7.0, 5.2), constrained_layout=True)
    ax.plot(
        omega_curve,
        analytical_curve,
        color="black",
        label=r"Analytical $\nu=\frac{1}{3}(\omega^{-1}-\frac{1}{2})$",
    )
    ax.scatter(
        measurements[:, 0], measurements[:, 1], color="tab:blue", zorder=3,
        label="Measured from shear-wave decay"
    )
    ax.set_xlabel(r"Relaxation parameter $\omega$")
    ax.set_ylabel(r"Kinematic viscosity $\nu$")
    ax.set_title("Measured viscosity")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(output_file, dpi=200)


def plot_channel_evolution(
    data_file: Path, output_file: Path, title: str
) -> None:
    """Plot transient channel profiles, steady solution, and wall locations."""
    data = np.loadtxt(data_file, comments="#")
    fig, ax = plt.subplots(figsize=(6.8, 6.2), constrained_layout=True)
    for step in selected_steps(data):
        snapshot = data[data[:, 0] == step]
        ax.plot(snapshot[:, 2], snapshot[:, 1], label=rf"$t={step}$")

    final_snapshot = data[data[:, 0] == np.max(data[:, 0])]
    ax.plot(
        final_snapshot[:, 3], final_snapshot[:, 1], "k--", lw=2,
        label="Steady analytical solution"
    )
    bottom_wall, top_wall = final_snapshot[0, 1], final_snapshot[-1, 1]
    ax.axhline(bottom_wall, color="tab:red", ls=":", lw=2, label="Walls")
    ax.axhline(top_wall, color="tab:red", ls=":", lw=2)
    ax.set_ylim(bottom_wall, top_wall)
    ax.set_xlabel(r"Streamwise velocity $u_x$")
    ax.set_ylabel(r"Channel coordinate $y$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")
    fig.savefig(output_file, dpi=200)


def load_flow_field(path: Path) -> tuple[np.ndarray, ...]:
    """Load one rectangular flow field and return mesh-shaped arrays."""
    data = np.loadtxt(path, comments="#")
    x_values, y_values = np.unique(data[:, 0]), np.unique(data[:, 1])
    shape = (len(y_values), len(x_values))
    x = data[:, 0].reshape(shape)
    y = data[:, 1].reshape(shape)
    rho = data[:, 2].reshape(shape)
    ux = data[:, 3].reshape(shape)
    uy = data[:, 4].reshape(shape)
    return x, y, rho, ux, uy


def plot_lid_evolution(files: list[Path], output_file: Path) -> None:
    """Plot selected lid-cavity snapshots as velocity streamplots."""
    indices = np.linspace(0, len(files) - 1, min(6, len(files))).round().astype(int)
    chosen = [files[index] for index in np.unique(indices)]
    columns = min(3, len(chosen))
    rows = int(np.ceil(len(chosen) / columns))
    fig, axes = plt.subplots(
        rows, columns, figsize=(5.0 * columns, 4.5 * rows),
        squeeze=False, constrained_layout=True
    )
    for ax, path in zip(axes.flat, chosen):
        x, y, _rho, ux, uy = load_flow_field(path)
        stride = max(1, min(x.shape) // 75)
        speed = np.hypot(ux, uy)
        if np.max(speed) > np.finfo(float).eps:
            stream = ax.streamplot(
                x[0, ::stride], y[::stride, 0],
                ux[::stride, ::stride], uy[::stride, ::stride],
                color=speed[::stride, ::stride], cmap="viridis", density=1.2
            )
            fig.colorbar(stream.lines, ax=ax, label="Speed")
        else:
            ax.text(
                0.5, 0.5, "Fluid initially at rest", ha="center", va="center",
                transform=ax.transAxes
            )
        step = int(path.stem.rsplit("_", 1)[-1])
        ax.set_title(rf"$t={step}$")
        ax.set_xlim(0.0, x.max() + 0.5)
        ax.set_ylim(0.0, y.max() + 0.5)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    for ax in axes.flat[len(chosen):]:
        ax.set_visible(False)
    fig.suptitle("Moving-lid cavity: temporal evolution")
    fig.savefig(output_file, dpi=180)


def plot_performance(files: list[Path], output_file: Path) -> None:
    """Plot lid-cavity MLUPS against the number of coarray images."""
    records = []
    for path in files:
        data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
        for row in np.atleast_1d(data):
            if row["case"] == "lid":
                records.append((int(row["nx"]), int(row["ny"]),
                                int(row["images"]), float(row["mlups"])))
    if not records:
        return

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for nx_value, ny_value in sorted({(row[0], row[1]) for row in records}):
        subset = [row for row in records if row[:2] == (nx_value, ny_value)]
        images = sorted({row[2] for row in subset})
        # Average repeated measurements at the same process count.
        mlups = [np.mean([row[3] for row in subset if row[2] == count])
                 for count in images]
        ax.plot(images, mlups, "o-", label=f"{nx_value}×{ny_value}")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    process_ticks = sorted({row[2] for row in records})
    ax.set_xticks(process_ticks, labels=[str(value) for value in process_ticks])
    ax.set_xlabel("MPI processes / coarray images")
    ax.set_ylabel("Performance [MLUPS]")
    ax.set_title("Moving-lid LBM strong scaling")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(title="Lattice")
    fig.savefig(output_file, dpi=200)


def generate_available_figures(data_dir: Path, output_dir: Path) -> list[Path]:
    """Generate every figure for which the corresponding data exist."""
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    jobs = [
        ("shear_profiles.txt", "shear_profiles.png", plot_wave_profiles,
         ("Shear-wave profiles", r"Lattice coordinate $y$")),
        ("density_profiles.txt", "density_profiles.png", plot_wave_profiles,
         ("Density-wave profiles", r"Lattice coordinate $x$")),
        ("shear_wave_decay.txt", "shear_wave_decay.png", plot_shear_decay, ()),
        ("couette_evolution.txt", "couette_evolution.png", plot_channel_evolution,
         ("Couette-flow evolution",)),
        ("poiseuille_evolution.txt", "poiseuille_evolution.png", plot_channel_evolution,
         ("Poiseuille-flow evolution",)),
    ]
    for input_name, output_name, function, extra_args in jobs:
        input_path = data_dir / input_name
        if input_path.exists():
            output_path = output_dir / output_name
            function(input_path, output_path, *extra_args)
            written.append(output_path)

    viscosity_files = sorted(data_dir.glob("viscosity_omega_*.txt"))
    if viscosity_files:
        output_path = output_dir / "viscosity_vs_omega.png"
        plot_viscosity(viscosity_files, output_path)
        written.append(output_path)

    lid_files = sorted(data_dir.glob("moving_lid_frame_*.txt"))
    if lid_files:
        output_path = output_dir / "moving_lid_evolution.png"
        plot_lid_evolution(lid_files, output_path)
        written.append(output_path)

    performance_files = sorted(data_dir.rglob("performance.csv"))
    if performance_files:
        output_path = output_dir / "lid_performance.png"
        plot_performance(performance_files, output_path)
        if output_path.exists():
            written.append(output_path)
    return written


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    written = generate_available_figures(args.data_dir, args.output_dir)
    if not written:
        raise SystemExit(f"No recognized simulation outputs found in {args.data_dir}")
    for path in written:
        print(f"Wrote {path}")
    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
