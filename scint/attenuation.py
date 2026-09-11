"""Photon attenuation via Geant4, with an on-disk cache.

Attenuation coefficients are not computed in Python. They are obtained from the
``g4data`` helper (``sim/tools/g4data.cc``), which queries Geant4's own EPICS2017
cross sections through ``G4EmCalculator``. Two reasons:

1. Licence. The EPICS2017 data files shipped with Geant4 are marked "NOT FOR
   COMMERCIAL USE AND MUST BE USED WITHIN GEANT4", so they cannot be vendored
   into a redistributable Python package.
2. Consistency. The screening layer and the transport simulation then agree by
   construction; a discrepancy between them can never be a data-version artefact.

Results are cached under ``data/derived/attenuation/`` keyed by a hash of the
exact query, so repeated screening runs do not re-invoke Geant4.
"""

from __future__ import annotations

import csv
import hashlib
import os
import subprocess
from collections.abc import Sequence
from pathlib import Path

from scint.materials import Material

_REPO_ROOT = Path(__file__).resolve().parent.parent
_G4DATA = _REPO_ROOT / "build" / "sim" / "g4data"
_CACHE_DIR = _REPO_ROOT / "data" / "derived" / "attenuation"
_GEANT4_DATA_ROOT = Path("/usr/local/share/Geant4/data")

# Geant4 locates its datasets through these variables. geant4.sh does not survive
# a non-interactive shell reliably, so they are set explicitly for the subprocess.
_DATASET_ENV = {
    "G4LEDATA": "G4EMLOW8.7",
    "G4LEVELGAMMADATA": "PhotonEvaporation6.1",
    "G4RADIOACTIVEDATA": "RadioactiveDecay6.1.2",
    "G4PARTICLEXSDATA": "G4PARTICLEXS4.1",
    "G4PIIDATA": "G4PII1.3",
    "G4REALSURFACEDATA": "RealSurface2.2",
    "G4SAIDXSDATA": "G4SAIDDATA2.0",
    "G4ABLADATA": "G4ABLA3.3",
    "G4INCLDATA": "G4INCL1.2",
    "G4ENSDFSTATEDATA": "G4ENSDFSTATE3.0",
    "G4NEUTRONHPDATA": "G4NDL4.7.1",
    "G4CHANNELINGDATA": "G4CHANNELING1.0",
}


class Geant4NotBuiltError(RuntimeError):
    """Raised when the g4data helper has not been compiled."""


def geant4_env(base: dict[str, str] | None = None, data_root: Path | None = None) -> dict[str, str]:
    """Return an environment with Geant4 dataset variables set explicitly."""
    root = data_root or _GEANT4_DATA_ROOT
    env = dict(base if base is not None else os.environ)
    for var, directory in _DATASET_ENV.items():
        env[var] = str(root / directory)
    return env


def _cache_key(material: Material, energies: Sequence[float]) -> str:
    payload = "|".join(
        [
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
    if not _G4DATA.exists():
        raise Geant4NotBuiltError(
            f"{_G4DATA} not found. Build it with:\n"
            "    cmake -S sim -B build/sim -DCMAKE_BUILD_TYPE=Release && cmake --build build/sim -j8"
        )

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = _CACHE_DIR / f"{material.name}_{_cache_key(material, energies_mev)}.csv"

    if refresh or not cache_file.exists():
        cmd = [
            str(_G4DATA), "attenuation",
            "--name", material.name.replace(" ", "_"),
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
