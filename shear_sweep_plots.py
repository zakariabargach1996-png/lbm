"""Plot a multi-omega shear-wave sweep using linear (normal) axes."""

from __future__ import annotations

from argparse import ArgumentParser
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


INITIAL_AMPLITUDE = 0.05


def load_viscosities(data_dir: Path) -> np.ndarray:
    """Load and sort omega, measured nu, applied nu, and relative error."""
    files = sorted(data_dir.glob("omega_*/viscosity_omega_*.txt"))
    if not files:
        raise FileNotFoundError(f"No viscosity measurements found below {data_dir}")
    values = np.vstack(
        [np.loadtxt(path, comments="#", ndmin=2) for path in files]
    )
    return values[np.argsort(values[:, 2])]


def plot_viscosity_comparison(values: np.ndarray, output_path: Path) -> None:
    """Compare measured viscosity directly with the applied BGK viscosity."""
    applied = values[:, 2]
    measured = values[:, 1]
    lower = min(applied.min(), measured.min())
    upper = max(applied.max(), measured.max())
    margin = 0.06 * (upper - lower)

    fig, ax = plt.subplots(figsize=(6.5, 5.8), constrained_layout=True)
    ax.plot(
        [lower - margin, upper + margin],
        [lower - margin, upper + margin],
        "k--",
        label="Ideal: measured = applied",
    )
    scatter = ax.scatter(
        applied,
        measured,
        c=values[:, 0],
        cmap="viridis",
        s=65,
        edgecolor="black",
        linewidth=0.5,
        zorder=3,
    )
    for omega, measured_nu, applied_nu, _error in values:
        ax.annotate(
            rf"$\omega={omega:.1f}$",
            (applied_nu, measured_nu),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlim(lower - margin, upper + margin)
    ax.set_ylim(lower - margin, upper + margin)
    ax.set_xlabel(r"Applied kinematic viscosity $\nu$")
    ax.set_ylabel(r"Measured kinematic viscosity $\nu_{\mathrm{measured}}$")
    ax.set_title("Applied versus measured viscosity")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.colorbar(scatter, ax=ax, label=r"Relaxation parameter $\omega$")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def load_decay_files(data_dir: Path) -> list[tuple[float, Path, np.ndarray]]:
    """Load each omega's measured and theoretical amplitude history."""
    histories: list[tuple[float, Path, np.ndarray]] = []
    for path in sorted(data_dir.glob("omega_*/shear_wave_decay.txt")):
        omega = float(path.parent.name.removeprefix("omega_"))
        histories.append((omega, path, np.loadtxt(path, comments="#", ndmin=2)))
    if not histories:
        raise FileNotFoundError(f"No shear-wave decay files found below {data_dir}")
    return histories


def measured_viscosity_at_steps(
    histories: list[tuple[float, Path, np.ndarray]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return applied and measured viscosities at shared sampled time steps."""
    common_steps: set[int] | None = None
    for _omega, _path, data in histories:
        positive_steps = set(data[data[:, 0] > 0, 0].astype(int))
        common_steps = (
            positive_steps
            if common_steps is None
            else common_steps.intersection(positive_steps)
        )
    if not common_steps:
        raise ValueError("The decay histories have no common positive time steps")

    available = np.asarray(sorted(common_steps), dtype=int)
    fractions = np.asarray([0.0, 0.1, 0.25, 0.5, 0.75, 1.0])
    indices = np.unique(
        np.rint(fractions * (len(available) - 1)).astype(int)
    )
    selected_steps = available[indices]

    ordered = sorted(histories, key=lambda item: (1.0 / item[0] - 0.5) / 3.0)
    applied = np.asarray([(1.0 / omega - 0.5) / 3.0 for omega, _, _ in ordered])
    measured = np.empty((len(selected_steps), len(ordered)))
    for column, (_omega, _path, data) in enumerate(ordered):
        amplitude_by_step = {
            int(row[0]): float(row[1])
            for row in data
        }
        wave_number = 2.0 * np.pi / 300.0
        positive = data[data[:, 0] > 0]
        applied_nu = applied[column]
        estimates = -np.log(
            np.abs(positive[:, 2] / INITIAL_AMPLITUDE)
        ) / (applied_nu * positive[:, 0])
        inferred_k2 = float(np.median(estimates))
        if not np.isfinite(inferred_k2) or inferred_k2 <= 0.0:
            inferred_k2 = wave_number**2
        for row, step in enumerate(selected_steps):
            amplitude = amplitude_by_step[int(step)]
            measured[row, column] = -np.log(
                abs(amplitude / INITIAL_AMPLITUDE)
            ) / (inferred_k2 * step)
    return selected_steps, applied, measured


def plot_viscosity_by_timestep(
    histories: list[tuple[float, Path, np.ndarray]], output_path: Path
) -> None:
    """Compare applied and measured viscosity at several simulation times."""
    steps, applied, measured = measured_viscosity_at_steps(histories)
    lower = min(applied.min(), measured.min())
    upper = max(applied.max(), measured.max())
    margin = 0.06 * (upper - lower)
    colors = plt.colormaps["viridis"](
        np.linspace(0.08, 0.92, len(steps))
    )

    fig, (ax, ax_error) = plt.subplots(
        2,
        1,
        figsize=(7.4, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        constrained_layout=True,
    )
    ax.plot(
        [lower - margin, upper + margin],
        [lower - margin, upper + margin],
        "k--",
        lw=1.8,
        label="Ideal: measured = applied",
    )
    for step, values, color in zip(steps, measured, colors):
        ax.plot(
            applied,
            values,
            "o-",
            color=color,
            ms=5,
            lw=1.4,
            alpha=0.9,
            label=f"t = {step}",
        )
        ax_error.plot(
            applied,
            100.0 * (values - applied) / applied,
            "o-",
            color=color,
            ms=4,
            lw=1.4,
            alpha=0.9,
        )

    ax.set_xlim(lower - margin, upper + margin)
    ax.set_ylim(lower - margin, upper + margin)
    ax.set_ylabel(r"Measured kinematic viscosity $\nu_{\mathrm{measured}}(t)$")
    ax.set_title("Applied versus measured viscosity during decay")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", ncol=2)
    ax_error.axhline(0.0, color="black", ls="--", lw=1.2)
    ax_error.set_xlabel(r"Applied kinematic viscosity $\nu$")
    ax_error.set_ylabel("Signed relative error [%]")
    ax_error.grid(True, alpha=0.3)
    scalar_map = plt.cm.ScalarMappable(
        norm=Normalize(vmin=float(steps.min()), vmax=float(steps.max())),
        cmap="viridis",
    )
    fig.colorbar(scalar_map, ax=[ax, ax_error], label="Time step")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_viscosity_time_history(
    histories: list[tuple[float, Path, np.ndarray]], output_path: Path
) -> None:
    """Plot inferred viscosity over time with the BGK value for each omega."""
    common_final_step = min(float(data[-1, 0]) for _, _, data in histories)
    colors = plt.colormaps["tab10"](np.linspace(0.0, 0.5, len(histories)))

    fig, ax = plt.subplots(figsize=(13.5, 7.0), constrained_layout=True)
    for color, (omega, _path, data) in zip(colors, histories):
        applied_nu = (1.0 / omega - 0.5) / 3.0
        positive = data[
            (data[:, 0] > 0.0)
            & (data[:, 0] <= common_final_step)
            & (np.abs(data[:, 1]) > 0.0)
        ]
        step = positive[:, 0]

        k2_samples = -np.log(
            np.abs(positive[:, 2] / INITIAL_AMPLITUDE)
        ) / (applied_nu * step)
        k2 = float(np.median(k2_samples[np.isfinite(k2_samples)]))
        measured_nu = -np.log(
            np.abs(positive[:, 1] / INITIAL_AMPLITUDE)
        ) / (k2 * step)

        ax.plot(
            step,
            measured_nu,
            color=color,
            lw=1.8,
            label=rf"Measured, $\omega={omega:.1f}$",
        )
        ax.axhline(
            applied_nu,
            color=color,
            ls="--",
            lw=1.4,
            alpha=0.9,
            label=rf"Theory, $\omega={omega:.1f}$ ($\nu={applied_nu:.5f}$)",
        )

    ax.set_xlim(left=0.0, right=common_final_step)
    ax.set_xlabel("Time step")
    ax.set_ylabel(r"Measured kinematic viscosity $\nu$")
    ax.set_title(
        r"Measured viscosity for different $\omega$ values (300$\times$300 grid)"
    )
    ax.grid(True, alpha=0.35)
    ax.legend(loc="center right", fontsize=8.5, ncol=2)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_decay_comparison(
    histories: list[tuple[float, Path, np.ndarray]], output_path: Path
) -> None:
    """Plot measured and theoretical exponential decay on normal axes."""
    columns = 3
    rows = math.ceil(len(histories) / columns)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(5.0 * columns, 3.8 * rows),
        squeeze=False,
        constrained_layout=True,
    )

    for ax, (omega, _path, data) in zip(axes.flat, histories):
        step, measured, theoretical, relative_error = data.T
        applied_nu = (1.0 / omega - 0.5) / 3.0
        ax.plot(step, theoretical, color="black", lw=2, label="Theory")
        ax.plot(
            step,
            measured,
            "o",
            color="tab:blue",
            ms=3,
            markevery=max(1, len(step) // 24),
            label="Measured",
        )
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel("Time step")
        ax.set_ylabel("Shear-wave amplitude")
        ax.set_title(
            rf"$\omega={omega:.1f}$, $\nu={applied_nu:.5f}$"
            f"\nmax relative error = {100.0 * relative_error.max():.2f}%"
        )
        ax.grid(True, alpha=0.3)
        ax.legend()

    for ax in axes.flat[len(histories):]:
        ax.set_visible(False)
    fig.suptitle("Measured versus theoretical exponential shear-wave decay")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_summary(values: np.ndarray, output_path: Path) -> None:
    """Write a convenient combined table for reports and further analysis."""
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["omega", "applied_viscosity", "measured_viscosity", "relative_error"]
        )
        for omega, measured, applied, relative_error in values:
            writer.writerow([omega, applied, measured, relative_error])


def main() -> None:
    """Parse arguments and create the shear-sweep reports."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    values = load_viscosities(args.data_dir)
    histories = load_decay_files(args.data_dir)
    viscosity_plot = args.output_dir / "applied_vs_measured_viscosity.png"
    timestep_plot = (
        args.output_dir / "applied_vs_measured_viscosity_by_timestep.png"
    )
    time_history_plot = (
        args.output_dir / "measured_viscosity_vs_theory_by_omega.png"
    )
    decay_plot = args.output_dir / "measured_vs_theory_decay.png"
    summary = args.output_dir / "viscosity_summary.csv"

    plot_viscosity_comparison(values, viscosity_plot)
    plot_viscosity_by_timestep(histories, timestep_plot)
    plot_viscosity_time_history(histories, time_history_plot)
    plot_decay_comparison(histories, decay_plot)
    write_summary(values, summary)
    print(f"Wrote {viscosity_plot}")
    print(f"Wrote {timestep_plot}")
    print(f"Wrote {time_history_plot}")
    print(f"Wrote {decay_plot}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
