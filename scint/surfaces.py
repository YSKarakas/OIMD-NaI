"""The surface-treatment / wrapping / coupling grid, as Geant4 supports it.

Geant4 flattens three physically independent choices -- mechanical treatment,
reflector material and optical coupling -- into a single enumerator such as
``polishedteflonair``. A light-collection study needs to vary them separately,
so ``sim/src/SurfaceCatalogue.cc`` re-expresses them as three axes and refuses
combinations Geant4 has no look-up table for. This module reads that catalogue
rather than duplicating it, so Python and C++ cannot drift apart.
"""

from __future__ import annotations

import csv
import subprocess
from dataclasses import dataclass
from functools import lru_cache

from scint.geant4 import require_g4data

TREATMENTS = ("polished", "etched", "ground")
WRAPPINGS = ("none", "lumirror", "teflon", "tio", "tyvek", "esr")
COUPLINGS = ("air", "glue")


@dataclass(frozen=True)
class SurfaceOption:
    """One point of the grid, and whether Geant4 can model it."""

    treatment: str
    wrapping: str
    coupling: str
    supported: bool
    g4_finish: str
    reason: str

    @property
    def label(self) -> str:
        return f"{self.treatment}/{self.wrapping}/{self.coupling}"


@lru_cache(maxsize=1)
def catalogue() -> tuple[SurfaceOption, ...]:
    """Every combination of the three axes, annotated with Geant4 support."""
    binary = require_g4data()
    result = subprocess.run([str(binary), "surfaces"], capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise RuntimeError(f"g4data surfaces failed: {result.stderr[-1000:]}")
    rows = csv.DictReader(
        line for line in result.stdout.splitlines() if not line.startswith("#")
    )
    return tuple(
        SurfaceOption(
            treatment=row["treatment"],
            wrapping=row["wrapping"],
            coupling=row["coupling"],
            supported=row["supported"] == "true",
            g4_finish=row["g4_finish"],
            reason=row["reason"],
        )
        for row in rows
    )


def supported_options() -> tuple[SurfaceOption, ...]:
    """The combinations that can actually be simulated."""
    return tuple(option for option in catalogue() if option.supported)


def resolve(treatment: str, wrapping: str, coupling: str) -> SurfaceOption:
    """One point of the grid by name, so a caller can read its refusal reason."""
    for option in catalogue():
        if (option.treatment, option.wrapping, option.coupling) == (
            treatment, wrapping, coupling
        ):
            return option
    raise KeyError(f"{treatment}/{wrapping}/{coupling} is not on the grid")


def sweep_axis() -> list[str]:
    """Supported combinations as ``treatment/wrapping/coupling`` sweep values."""
    return [option.label for option in supported_options()]
