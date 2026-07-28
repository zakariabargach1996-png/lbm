"""Reduce full density-wave frames to compact x-profiles for the presentation."""

from pathlib import Path
import re

import numpy as np


SOURCE = Path("/tmp/lbm_density_presentation")
DESTINATION = Path(__file__).with_name("density_profiles.txt")

available = sorted(SOURCE.glob("density_frame_*.txt"))
if not available:
    raise FileNotFoundError(f"No density frames found below {SOURCE}")

selected_indices = [0, len(available) // 3, 2 * len(available) // 3, len(available) - 1]
rows = []
for path in [available[index] for index in selected_indices]:
    match = re.search(r"(\d+)$", path.stem)
    step = int(match.group(1))
    field = np.loadtxt(path, comments="#")
    for x in np.unique(field[:, 0]).astype(int):
        mean_density = field[field[:, 0] == x, 2].mean()
        rows.append((step, x, mean_density))

data = np.asarray(rows)
np.savetxt(
    DESTINATION,
    data,
    header="step x mean_density",
    fmt=["%d", "%d", "%.16e"],
)
print(f"Wrote {DESTINATION} from {len(selected_indices)} full frames")
