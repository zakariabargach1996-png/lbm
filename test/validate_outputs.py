"""Validate available LBM result files against physical expectations."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np


def relative_l2(numerical: np.ndarray, analytical: np.ndarray) -> float:
    """Return the relative L2 difference between two arrays."""
    return float(np.linalg.norm(numerical - analytical) / np.linalg.norm(analytical))


def validate(data_dir: Path) -> list[str]:
    """Validate every recognized solver output in a directory."""
    checks: list[str] = []

    shear_decay = data_dir / "shear_wave_decay.txt"
    if shear_decay.exists():
        data = np.loadtxt(shear_decay, comments="#")
        assert np.isfinite(data).all(), "shear decay contains non-finite values"
        assert data[-1, 3] < 1.0e-2, "shear-wave decay error exceeds 1%"
        checks.append(f"shear decay error {data[-1, 3]:.3e}")

    density_profiles = data_dir / "density_profiles.txt"
    if density_profiles.exists():
        data = np.loadtxt(density_profiles, comments="#")
        assert np.isfinite(data).all(), "density profiles contain non-finite values"
        for step in np.unique(data[:, 0]):
            mean_density = data[data[:, 0] == step, 2].mean()
            assert abs(mean_density - 1.0) < 1.0e-10, "density-wave mass drift"
        checks.append("density-wave mass conservation")

    for name, label, tolerance in (
        ("couette_profile.txt", "Couette", 1.0e-2),
        ("poiseuille_profile.txt", "Poiseuille", 2.0e-2),
    ):
        path = data_dir / name
        if path.exists():
            data = np.loadtxt(path, comments="#")
            error = relative_l2(data[1:-1, 1], data[1:-1, 2])
            assert error < tolerance, f"{label} profile error exceeds {tolerance:.1%}"
            checks.append(f"{label} relative L2 error {error:.3e}")

    viscosity_files = sorted(data_dir.glob("viscosity_omega_*.txt"))
    if viscosity_files:
        data = np.vstack(
            [np.loadtxt(path, comments="#", ndmin=2) for path in viscosity_files]
        )
        assert np.all(data[:, 0] <= 1.7), "viscosity sweep contains omega > 1.7"
        assert np.all(data[:, 3] < 1.0e-2), "measured viscosity error exceeds 1%"
        checks.append(f"{len(data)} viscosity measurement(s) within 1%")

    lid_files = sorted(data_dir.glob("moving_lid_frame_*.txt"))
    if lid_files:
        initial = np.loadtxt(lid_files[0], comments="#")
        final = np.loadtxt(lid_files[-1], comments="#")
        assert np.isfinite(final).all(), "lid field contains non-finite values"
        assert np.max(np.abs(initial[:, 3:5])) < 1.0e-14, "lid initial velocity is nonzero"
        assert abs(final[:, 2].mean() - 1.0) < 1.0e-8, "lid-cavity mass drift"
        checks.append("lid initial conditions, finiteness, and mass conservation")

    if not checks:
        raise AssertionError(f"no recognized result files in {data_dir}")
    return checks


def main() -> None:
    """Parse arguments and print completed validation checks."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    for check in validate(args.data_dir):
        print(f"PASS: {check}")


if __name__ == "__main__":
    main()
