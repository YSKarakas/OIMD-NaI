"""Photon attenuation via Geant4, with a version-keyed on-disk cache.

Attenuation coefficients are not computed in Python. They come from the
``g4data`` helper (``sim/tools/g4data.cc``), which queries Geant4's own EPICS2017
cross sections through ``G4EmCalculator``. Two reasons:

1. Licence. The EPICS2017 data files shipped with Geant4 are marked "NOT FOR
   COMMERCIAL USE AND MUST BE USED WITHIN GEANT4", so they cannot be vendored
   into a redistributable Python package.
2. Consistency. The screening layer and the transport simulation then agree by
   construction; a discrepancy between them can never be a data-version artefact.

The cache key includes the Geant4 version tag, so results from one build are
never served to another -- the exact mistake that a ``geant4-config`` version
string, which cannot distinguish a beta from a release, would let through.
"""

from __future__ import annotations

import csv
import hashlib
import subprocess
from collections.abc import Sequence
from pathlib import Path

from scint.geant4 import Geant4NotBuiltError, geant4_env, require_g4data, version_key
from scint.materials import Material

__all__ = ["Geant4NotBuiltError", "attenuation_lengths_cm", "geant4_env"]

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CACHE_DIR = _REPO_ROOT / "data" / "derived" / "attenuation"


def _cache_key(material: Material, energies: Sequence[float]) -> str:
    payload = "|".join(
        [
            version_key(),
            material.formula,
            f"{material.density_g_cm3!r}",
            material.massfrac_arg(),
            ",".join(f"{e!r}" for e in energies),
        ]
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def attenuation_lengths_cm(
    material: Material,
    energies_mev: Sequence[float],
    *,
    refresh: bool = False,
) -> dict[float, float]:
    """Return {energy in MeV: total attenuation length in cm} for ``material``.

    The attenuation length is the total photon mean free path, 1/mu, including
    coherent (Rayleigh) scattering -- the same convention as standard tabulations.
    """
    binary = require_g4data()
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    safe_name = material.name.replace(" ", "_").replace("/", "_")
    cache_file = _CACHE_DIR / f"{safe_name}_{_cache_key(material, energies_mev)}.csv"

    if refresh or not cache_file.exists():
        cmd = [
            str(binary), "attenuation",
            "--name", safe_name,
            "--density", repr(material.density_g_cm3),
            "--massfrac", material.massfrac_arg(),
            "--energies", ",".join(repr(e) for e in energies_mev),
            "--out", str(cache_file),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, env=geant4_env())
        if result.returncode != 0:
            raise RuntimeError(
                f"g4data failed for {material.name} (exit {result.returncode}):\n"
                f"{result.stderr[-2000:]}"
            )

    with cache_file.open() as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        return {float(r["energy_MeV"]): float(r["attenuation_length_cm"]) for r in rows}
