"""Generate strong-scaling plots from LBM ``performance.csv`` files.

The solver writes one performance row per completed run.  This script searches
recursively below ``--data-dir``, combines repeated measurements, and creates:

* wall-clock time versus coarray images;
* throughput in MLUPS versus coarray images;
* strong-scaling speedup versus coarray images;
* parallel efficiency versus coarray images;
* a four-panel dashboard containing all of the above; and
* a combined CSV table containing the aggregated measurements.

Speedup and efficiency are calculated separately for every lattice size.  The
smallest available image count is used as that lattice's baseline.  For the
usual strong-scaling study this should be the one-image measurement.
"""

from __future__ import annotations

from argparse import ArgumentParser
from collections import defaultdict
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
REQUIRED_COLUMNS = {
    "case",
    "nx",
    "ny",
    "images",
    "px",
    "py",
    "steps",
    "wall_seconds",
    "mlups",
}


@dataclass(frozen=True)
class Measurement:
    """One row produced by ``report_performance`` in ``app/main.f90``."""

    case: str
    nx: int
    ny: int
    images: int
    px: int
    py: int
    steps: int
    wall_seconds: float
    mlups: float
    source: Path


@dataclass(frozen=True)
class Summary:
    """Repeated measurements aggregated at one strong-scaling point."""

    case: str
    nx: int
    ny: int
    steps: int
    images: int
    px: int
    py: int
    repeats: int
    wall_seconds: float
    wall_min: float
    wall_max: float
    mlups: float
    mlups_min: float
    mlups_max: float

    @property
    def series_key(self) -> tuple[str, int, int, int]:
        """Identify measurements that belong to one strong-scaling curve."""
        return self.case, self.nx, self.ny, self.steps


def parse_integer(row: dict[str, str], column: str, path: Path) -> int:
    """Parse an integer column and include the source path in any error."""
    try:
        return int(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid {column!r} value") from error


def parse_float(row: dict[str, str], column: str, path: Path) -> float:
    """Parse a finite floating-point column."""
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid {column!r} value") from error
    if not np.isfinite(value):
        raise ValueError(f"{path}: {column!r} is not finite")
    return value


def read_measurements(data_dir: Path, selected_case: str) -> list[Measurement]:
    """Read all matching performance rows below ``data_dir``."""
    paths = sorted(data_dir.rglob("performance.csv"))
    if not paths:
        raise FileNotFoundError(
            f"No performance.csv files found below {data_dir.resolve()}"
        )

    measurements: list[Measurement] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            columns = set(reader.fieldnames or ())
            missing = REQUIRED_COLUMNS - columns
            if missing:
                names = ", ".join(sorted(missing))
                print(f"Warning: skipping {path}; missing columns: {names}",
                      file=sys.stderr)
                continue

            for row in reader:
                case = row["case"].strip().lower()
                if selected_case != "all" and case != selected_case:
                    continue
                measurement = Measurement(
                    case=case,
                    nx=parse_integer(row, "nx", path),
                    ny=parse_integer(row, "ny", path),
                    images=parse_integer(row, "images", path),
                    px=parse_integer(row, "px", path),
                    py=parse_integer(row, "py", path),
                    steps=parse_integer(row, "steps", path),
                    wall_seconds=parse_float(row, "wall_seconds", path),
                    mlups=parse_float(row, "mlups", path),
                    source=path,
                )
                if measurement.images < 1 or measurement.steps < 1:
                    raise ValueError(
                        f"{path}: images and steps must both be positive"
                    )
                if measurement.wall_seconds <= 0.0 or measurement.mlups <= 0.0:
                    raise ValueError(
                        f"{path}: wall_seconds and mlups must both be positive"
                    )
                measurements.append(measurement)

    if not measurements:
        raise ValueError(
            f"No case={selected_case!r} performance rows found below "
            f"{data_dir.resolve()}"
        )
    return measurements


def aggregate_measurements(
    measurements: list[Measurement], statistic: str
) -> list[Summary]:
    """Combine repeated runs with the requested central statistic."""
    groups: dict[
        tuple[str, int, int, int, int], list[Measurement]
    ] = defaultdict(list)
    for item in measurements:
        key = item.case, item.nx, item.ny, item.steps, item.images
        groups[key].append(item)

    if statistic == "median":
        reduce_values = np.median
    elif statistic == "mean":
        reduce_values = np.mean
    else:
        # "best" means the shortest time and, equivalently, greatest MLUPS.
        reduce_values = None

    summaries: list[Summary] = []
    for (case, nx, ny, steps, images), items in groups.items():
        times = np.asarray([item.wall_seconds for item in items])
        rates = np.asarray([item.mlups for item in items])
        if statistic == "best":
            selected = int(np.argmin(times))
            central_time = float(times[selected])
            central_rate = float(rates[selected])
        else:
            central_time = float(reduce_values(times))
            central_rate = float(reduce_values(rates))

        # The decomposition should be deterministic.  Retaining the most recent
        # px/py is sufficient for the table and highlights inconsistent input.
        decompositions = {(item.px, item.py) for item in items}
        if len(decompositions) > 1:
            print(
                f"Warning: {nx}x{ny}, p={images} has multiple decompositions: "
                f"{sorted(decompositions)}",
                file=sys.stderr,
            )
        px, py = items[-1].px, items[-1].py
        summaries.append(
            Summary(
                case=case,
                nx=nx,
                ny=ny,
                steps=steps,
                images=images,
                px=px,
                py=py,
                repeats=len(items),
                wall_seconds=central_time,
                wall_min=float(np.min(times)),
                wall_max=float(np.max(times)),
                mlups=central_rate,
                mlups_min=float(np.min(rates)),
                mlups_max=float(np.max(rates)),
            )
        )
    return sorted(
        summaries,
        key=lambda item: (item.case, item.nx, item.ny, item.steps, item.images),
    )


def organize_series(
    summaries: list[Summary],
) -> dict[tuple[str, int, int, int], list[Summary]]:
    """Group aggregated points into strong-scaling curves."""
    series: dict[tuple[str, int, int, int], list[Summary]] = defaultdict(list)
    for item in summaries:
        series[item.series_key].append(item)
    for points in series.values():
        points.sort(key=lambda item: item.images)
    return dict(series)


def series_label(key: tuple[str, int, int, int], show_steps: bool) -> str:
    """Create a concise legend label."""
    _case, nx, ny, steps = key
    label = f"{nx}×{ny}"
    if show_steps:
        label += f", {steps} steps"
    return label


def configure_image_axis(axis: plt.Axes, all_images: list[int]) -> None:
    """Use a base-two log axis when possible and label every measured count."""
    axis.set_xscale("log", base=2)
    axis.set_xticks(all_images, labels=[str(value) for value in all_images])
    axis.set_xlabel("Coarray images")
    axis.grid(True, which="both", alpha=0.3)


def draw_runtime(
    axis: plt.Axes,
    series: dict[tuple[str, int, int, int], list[Summary]],
    show_steps: bool,
) -> None:
    """Draw wall-clock time for each fixed-size problem."""
    for key, points in series.items():
        images = [point.images for point in points]
        times = [point.wall_seconds for point in points]
        axis.plot(images, times, "o-", label=series_label(key, show_steps))
    axis.set_yscale("log")
    axis.set_ylabel("Wall-clock time [s]")
    axis.set_title("Runtime")


def draw_mlups(
    axis: plt.Axes,
    series: dict[tuple[str, int, int, int], list[Summary]],
    show_steps: bool,
) -> None:
    """Draw lattice-update throughput."""
    for key, points in series.items():
        images = [point.images for point in points]
        rates = [point.mlups for point in points]
        axis.plot(images, rates, "o-", label=series_label(key, show_steps))
    axis.set_yscale("log")
    axis.set_ylabel("Performance [MLUPS]")
    axis.set_title("Throughput")


def scaling_values(points: list[Summary]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return image counts, speedup, and efficiency for one curve."""
    images = np.asarray([point.images for point in points], dtype=float)
    times = np.asarray([point.wall_seconds for point in points], dtype=float)
    baseline_images = images[0]
    speedup = times[0] / times
    ideal_speedup = images / baseline_images
    efficiency = 100.0 * speedup / ideal_speedup
    return images, speedup, efficiency


def draw_speedup(
    axis: plt.Axes,
    series: dict[tuple[str, int, int, int], list[Summary]],
    show_steps: bool,
) -> None:
    """Draw measured speedup and the ideal linear reference."""
    maximum_ratio = 1.0
    baselines: set[int] = set()
    for key, points in series.items():
        images, speedup, _efficiency = scaling_values(points)
        ratio = images / images[0]
        baselines.add(int(images[0]))
        maximum_ratio = max(maximum_ratio, float(np.max(ratio)))
        axis.plot(images, speedup, "o-", label=series_label(key, show_steps))

    # Normally every curve has a one-image baseline.  If incomplete data leave
    # different baselines, draw a correct ideal reference for each one.
    all_images = sorted(
        {point.images for points in series.values() for point in points}
    )
    for baseline in sorted(baselines):
        reference_images = np.asarray(
            [value for value in all_images if value >= baseline], dtype=float
        )
        label = "Ideal linear speedup"
        if len(baselines) > 1:
            label += f" (baseline p={baseline})"
        axis.plot(
            reference_images,
            reference_images / baseline,
            "k--",
            linewidth=1.5,
            alpha=0.8,
            label=label,
        )
    axis.set_ylabel("Speedup")
    axis.set_title("Strong-scaling speedup")
    axis.set_ylim(bottom=0.0, top=max(1.1, 1.08 * maximum_ratio))


def draw_efficiency(
    axis: plt.Axes,
    series: dict[tuple[str, int, int, int], list[Summary]],
    show_steps: bool,
) -> None:
    """Draw strong-scaling parallel efficiency."""
    for key, points in series.items():
        images, _speedup, efficiency = scaling_values(points)
        axis.plot(images, efficiency, "o-", label=series_label(key, show_steps))
    axis.axhline(100.0, color="black", linestyle="--", linewidth=1.5,
                 label="Ideal efficiency")
    axis.set_ylabel("Parallel efficiency [%]")
    axis.set_title("Parallel efficiency")
    axis.set_ylim(bottom=0.0)


def write_single_plot(
    output_path: Path,
    draw_function,
    series: dict[tuple[str, int, int, int], list[Summary]],
    all_images: list[int],
    show_steps: bool,
) -> None:
    """Write one presentation-ready performance plot."""
    figure, axis = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    draw_function(axis, series, show_steps)
    configure_image_axis(axis, all_images)
    axis.legend(title="Lattice", fontsize="small")
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def write_dashboard(
    output_path: Path,
    series: dict[tuple[str, int, int, int], list[Summary]],
    all_images: list[int],
    show_steps: bool,
) -> None:
    """Write a compact four-panel summary figure."""
    figure, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.5), constrained_layout=True
    )
    drawers = (draw_runtime, draw_mlups, draw_speedup, draw_efficiency)
    for axis, draw_function in zip(axes.flat, drawers):
        draw_function(axis, series, show_steps)
        configure_image_axis(axis, all_images)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        title="Lattice",
        loc="outside right center",
    )
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def write_summary(path: Path, summaries: list[Summary]) -> None:
    """Write the aggregated data and derived scaling metrics."""
    series = organize_series(summaries)
    scaling: dict[tuple[str, int, int, int, int], tuple[float, float]] = {}
    for key, points in series.items():
        images, speedup, efficiency = scaling_values(points)
        for image_count, speedup_value, efficiency_value in zip(
            images, speedup, efficiency
        ):
            scaling[(*key, int(image_count))] = (
                float(speedup_value),
                float(efficiency_value),
            )

    fieldnames = [
        "case", "nx", "ny", "steps", "images", "px", "py", "repeats",
        "wall_seconds", "wall_min", "wall_max",
        "mlups", "mlups_min", "mlups_max", "speedup", "efficiency_percent",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for item in summaries:
            speedup, efficiency = scaling[(*item.series_key, item.images)]
            writer.writerow(
                {
                    "case": item.case,
                    "nx": item.nx,
                    "ny": item.ny,
                    "steps": item.steps,
                    "images": item.images,
                    "px": item.px,
                    "py": item.py,
                    "repeats": item.repeats,
                    "wall_seconds": f"{item.wall_seconds:.12g}",
                    "wall_min": f"{item.wall_min:.12g}",
                    "wall_max": f"{item.wall_max:.12g}",
                    "mlups": f"{item.mlups:.12g}",
                    "mlups_min": f"{item.mlups_min:.12g}",
                    "mlups_max": f"{item.mlups_max:.12g}",
                    "speedup": f"{speedup:.12g}",
                    "efficiency_percent": f"{efficiency:.12g}",
                }
            )


def main() -> None:
    """Parse arguments and generate all performance-analysis artifacts."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "cluster_runs",
        help="directory searched recursively for performance.csv files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="figure directory (default: DATA_DIR/performance_figures)",
    )
    parser.add_argument(
        "--case",
        default="lid",
        help="case to plot, or 'all' (default: lid)",
    )
    parser.add_argument(
        "--statistic",
        choices=("median", "mean", "best"),
        default="median",
        help="how repeated runs are combined (default: median)",
    )
    args = parser.parse_args()

    data_dir = args.data_dir.expanduser()
    output_dir = (
        args.output_dir.expanduser()
        if args.output_dir is not None
        else data_dir / "performance_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    measurements = read_measurements(data_dir, args.case.lower())
    summaries = aggregate_measurements(measurements, args.statistic)
    series = organize_series(summaries)
    all_images = sorted({item.images for item in summaries})
    # Include steps in the legend only if a lattice was benchmarked for more
    # than one fixed step count.
    lattice_steps: dict[tuple[str, int, int], set[int]] = defaultdict(set)
    for item in summaries:
        lattice_steps[(item.case, item.nx, item.ny)].add(item.steps)
    show_steps = any(len(values) > 1 for values in lattice_steps.values())

    outputs = {
        "runtime_vs_images.png": draw_runtime,
        "mlups_vs_images.png": draw_mlups,
        "speedup_vs_images.png": draw_speedup,
        "parallel_efficiency.png": draw_efficiency,
    }
    for filename, draw_function in outputs.items():
        write_single_plot(
            output_dir / filename,
            draw_function,
            series,
            all_images,
            show_steps,
        )

    write_dashboard(
        output_dir / "performance_dashboard.png",
        series,
        all_images,
        show_steps,
    )
    write_summary(output_dir / "performance_summary.csv", summaries)

    print(f"Read {len(measurements)} measurements from {data_dir.resolve()}")
    print(f"Created {len(summaries)} aggregated scaling points")
    for filename in [*outputs, "performance_dashboard.png",
                     "performance_summary.csv"]:
        print(f"Wrote {output_dir / filename}")


if __name__ == "__main__":
    main()
