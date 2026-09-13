"""Locating, describing and invoking the Geant4 installation.

The version handling here exists because of a specific trap. ``geant4-config
--version`` reports ``11.4.0`` for both the 11.4 beta and the 11.4 release: the
beta already carries the target ``G4VERSION_NUMBER`` (1140). The only reliable
discriminator is ``G4VERSION_TAG``, which the ``g4data version`` command exports.

The development machine used for this project turned out to be running
``geant4-11-04-beta-01`` (G4Date 26-June-2025) while reporting 11.4.0. Since the
11.4 cycle refactored ``G4OpBoundaryProcess`` and moved scintillation onto a
model-based framework -- the very code the light-yield geometry study rests on --
that distinction decides whether a result is publishable.
"""

from __future__ import annotations

import json
import os
import subprocess
from functools import lru_cache
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
# The version tool must be the one built against the toolkit actually in use.
# Inside the container that is build/container/g4data, not the host build, and
# the launcher says so through the environment; the default is the host build.
G4DATA_BIN = Path(os.environ.get("SCINT_G4DATA_BIN", str(_REPO_ROOT / "build" / "sim" / "g4data")))
GEANT4_DATA_ROOT = Path(os.environ.get("SCINT_G4_DATA_ROOT", "/usr/local/share/Geant4/data"))

PRERELEASE_MARKERS = ("beta", "alpha", "rc", "cand", "ref")

# Geant4 locates its datasets through these variables. geant4.sh does not survive
# a non-interactive shell reliably, so they are set explicitly for subprocesses.
# Inside the container they are already present in the environment and are left
# untouched.
_DATASET_DIRS = {
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


def require_g4data() -> Path:
    """Return the path to the g4data binary, or explain how to build it."""
    if not G4DATA_BIN.exists():
        raise Geant4NotBuiltError(
            f"{G4DATA_BIN} not found. Build it with:\n"
            "    cmake -S sim -B build/sim -DCMAKE_BUILD_TYPE=Release\n"
            "    cmake --build build/sim -j8"
        )
    return G4DATA_BIN


def geant4_env(base: dict[str, str] | None = None) -> dict[str, str]:
    """Return an environment with Geant4 dataset variables set.

    Variables already present in the environment win, so a container that has
    them baked in is not overridden with host paths.
    """
    env = dict(base if base is not None else os.environ)
    for var, directory in _DATASET_DIRS.items():
        env.setdefault(var, str(GEANT4_DATA_ROOT / directory))
    return env


@lru_cache(maxsize=1)
def version_info() -> dict[str, object]:
    """Return the Geant4 version as reported by the library itself.

    Falls back to ``geant4-config`` when g4data has not been built, but that
    fallback cannot distinguish a beta from a release, so it is marked as such.
    """
    try:
        binary = require_g4data()
    except Geant4NotBuiltError:
        reported = _config_version()
        return {
            "version_number": None,
            "tag": None,
            "date": None,
            "config_version": reported,
            "source": "geant4-config",
            "is_prerelease": None,  # unknown: geant4-config cannot tell
        }

    result = subprocess.run([str(binary), "version"], capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        raise RuntimeError(f"g4data version failed: {result.stderr[-1000:]}")
    info = json.loads(result.stdout)
    info["source"] = "g4data"
    info["is_prerelease"] = is_prerelease(info.get("tag"))
    # Recorded alongside the tag precisely to show that the two disagree in
    # substance: this string is identical for a series' beta and its release.
    info["config_version"] = _config_version()
    return info


def _config_version() -> str | None:
    """What ``geant4-config --version`` reports. Informative only."""
    try:
        out = subprocess.run(
            ["geant4-config", "--version"], capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def is_prerelease(tag: str | None) -> bool | None:
    """Whether a Geant4 tag denotes a pre-release build."""
    if not tag:
        return None
    lowered = tag.lower()
    return any(marker in lowered for marker in PRERELEASE_MARKERS)


def version_key() -> str:
    """A short, stable string identifying the Geant4 build, for cache keys.

    Cross-version cache reuse is the failure mode this prevents: physics results
    computed by one Geant4 build must never be silently served to another.
    """
    info = version_info()
    return str(info.get("tag") or info.get("config_version") or "unknown-geant4")
