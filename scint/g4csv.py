"""Reading the CSV files written by Geant4's analysis manager.

Geant4's CSV writer does not emit a conventional header row. Column names and
types live in ``#column <type> <name>`` comment lines, and the data follows with
no header at all, so a naive ``read_csv`` silently treats the first event as
column names. This module parses the real header.

CSV is preferred over ROOT for pipeline output because it is diffable, needs no
ROOT installation to read, and can be checksummed meaningfully for the run record.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

COLUMN_PREFIX = "#column "


class G4CsvError(ValueError):
    """Raised when a file does not look like Geant4 analysis-manager CSV."""


def read_ntuple(path: str | Path) -> dict[str, np.ndarray]:
    """Read a Geant4 CSV ntuple into {column name: array}.

    Empty files (a run that produced no rows) return empty arrays rather than
    raising, so that a sweep point with no hits is still a valid result.
    """
    path = Path(path)
    names: list[str] = []
    rows: list[str] = []

    with path.open() as handle:
        for line in handle:
            if line.startswith(COLUMN_PREFIX):
                parts = line[len(COLUMN_PREFIX):].split()
                if len(parts) < 2:
                    raise G4CsvError(f"{path}: malformed column header: {line.strip()}")
                names.append(parts[-1])
            elif line.startswith("#"):
                continue
            elif line.strip():
                rows.append(line)

    if not names:
        raise G4CsvError(
            f"{path}: no '#column' headers found -- is this a Geant4 analysis-manager CSV?"
        )

    if not rows:
        return {name: np.empty(0) for name in names}

    data = np.genfromtxt(rows, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(names):
        raise G4CsvError(
            f"{path}: header declares {len(names)} columns but rows carry {data.shape[1]}"
        )
    return {name: data[:, index] for index, name in enumerate(names)}


def resolve_output(stem: str | Path, ntuple: str = "events") -> Path:
    """Return the file the analysis manager actually wrote for an output stem.

    Geant4 appends ``_nt_<ntuple>`` to the stem, which is easy to trip over.
    """
    stem = Path(stem)
    return stem.with_name(f"{stem.name}_nt_{ntuple}.csv")
