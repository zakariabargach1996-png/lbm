"""Summarize fixed-image LBM runs made with different node counts."""

from __future__ import annotations

from argparse import ArgumentParser
from collections import defaultdict
import csv
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt


REQUIRED = {
    "case", "nx", "ny", "images", "steps", "wall_seconds", "mlups",
    "nodes", "tasks_per_node", "repeat",
}


def read_rows(data_dir: Path) -> list[dict[str, str]]:
    """Read placement measurements recursively from a directory."""
    rows: list[dict[str, str]] = []
    for path in sorted(data_dir.rglob("placement.csv")):
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            missing = REQUIRED - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            rows.extend(reader)
    if not rows:
        raise FileNotFoundError(f"No placement.csv files found below {data_dir}")
    return rows


def summarize(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    """Combine repeated placement measurements by their median."""
    groups: dict[tuple[object, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            row["case"], int(row["nx"]), int(row["ny"]), int(row["images"]),
            int(row["steps"]), int(row["nodes"]), int(row["tasks_per_node"]),
        )
        groups[key].append(row)

    result: list[dict[str, object]] = []
    for key, samples in groups.items():
        case, nx, ny, images, steps, nodes, tasks_per_node = key
        times = [float(sample["wall_seconds"]) for sample in samples]
        rates = [float(sample["mlups"]) for sample in samples]
        result.append({
            "case": case,
            "nx": nx,
            "ny": ny,
            "images": images,
            "steps": steps,
            "nodes": nodes,
            "tasks_per_node": tasks_per_node,
            "repeats": len(samples),
            "median_wall_seconds": median(times),
            "median_mlups": median(rates),
            "min_wall_seconds": min(times),
            "max_wall_seconds": max(times),
        })
    return sorted(
        result,
        key=lambda row: (
            row["case"], row["nx"], row["ny"], row["images"], row["steps"],
            row["nodes"],
        ),
    )


def write_summary(rows: list[dict[str, object]], path: Path) -> None:
    """Write aggregated placement measurements as CSV."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def draw(rows: list[dict[str, object]], path: Path) -> None:
    """Plot throughput against nodes and tasks per node."""
    series: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = row["nx"], row["ny"], row["images"], row["steps"]
        series[key].append(row)

    figure, (runtime_axis, rate_axis) = plt.subplots(1, 2, figsize=(11, 4.5))
    for (nx, ny, images, steps), points in series.items():
        points.sort(key=lambda row: int(row["nodes"]))
        nodes = [int(row["nodes"]) for row in points]
        label = f"{nx}×{ny}, {images} images, {steps} steps"
        runtime_axis.plot(
            nodes, [float(row["median_wall_seconds"]) for row in points],
            marker="o", label=label,
        )
        rate_axis.plot(
            nodes, [float(row["median_mlups"]) for row in points],
            marker="o", label=label,
        )

    for axis in (runtime_axis, rate_axis):
        axis.set_xscale("log", base=2)
        axis.set_xlabel("Nodes (fixed total image count)")
        axis.grid(True, which="both", alpha=0.3)
        axis.legend()
    runtime_axis.set_ylabel("Median wall time [s] (lower is better)")
    rate_axis.set_ylabel("Median throughput [MLUPS] (higher is better)")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    """Parse arguments and generate placement summaries."""
    parser = ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("cluster_runs"))
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("cluster_runs/placement_figures"),
    )
    args = parser.parse_args()

    rows = summarize(read_rows(args.data_dir))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "placement_summary.csv"
    figure_path = args.output_dir / "placement_comparison.png"
    write_summary(rows, summary_path)
    draw(rows, figure_path)

    best = min(rows, key=lambda row: float(row["median_wall_seconds"]))
    print(f"Wrote {summary_path}")
    print(f"Wrote {figure_path}")
    print(
        "Fastest measured placement: "
        f"{best['nodes']} nodes × {best['tasks_per_node']} images/node "
        f"({best['median_wall_seconds']:.3f} s, {best['median_mlups']:.3f} MLUPS)"
    )


if __name__ == "__main__":
    main()
