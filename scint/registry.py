"""Run provenance: deterministic identity and complete environment capture.

The contract is simple and strict: a run directory contains everything needed to
explain and repeat the run. ``config.yaml`` is the *entire* input -- if a result
depends on something not in that file, that is a bug. The run identifier is a
hash of the configuration's canonical form, so the same inputs always map to the
same directory and the input-to-output correspondence cannot be lost or
retroactively edited.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_RUNS_DIR = _REPO_ROOT / "runs"
RUN_ID_LENGTH = 16


class RunExistsError(RuntimeError):
    """Raised when a run directory for this configuration already exists."""


def canonical_yaml(config: dict[str, Any]) -> str:
    """Serialise a configuration to a stable, order-independent string.

    Mapping keys are sorted, so two configurations that differ only in key order
    produce identical text and therefore identical run identifiers.
    """
    return yaml.safe_dump(config, sort_keys=True, default_flow_style=False, allow_unicode=True)


def run_id(config: dict[str, Any]) -> str:
    """Return the deterministic identifier for a configuration."""
    digest = hashlib.sha256(canonical_yaml(config).encode("utf-8")).hexdigest()
    return digest[:RUN_ID_LENGTH]


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args], cwd=_REPO_ROOT, capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _command_version(*cmd: str) -> str | None:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def _geant4_version() -> dict[str, Any]:
    """Version details from the Geant4 library itself, never fatal."""
    try:
        from scint.geant4 import version_info

        return dict(version_info())
    except Exception as exc:  # a provenance record must never block a run
        return {"error": f"{type(exc).__name__}: {exc}"}


def _installed_packages() -> dict[str, str]:
    try:
        from importlib.metadata import distributions
    except ImportError:  # pragma: no cover - Python < 3.8 only
        return {}
    return {
        dist.metadata["Name"]: dist.version
        for dist in distributions()
        if dist.metadata.get("Name")
    }


def capture_environment(seed: int | None = None) -> dict[str, Any]:
    """Record everything outside ``config.yaml`` that could influence a result.

    ``git.dirty`` is the value to watch: a run made from an uncommitted tree is
    not reproducible from the recorded commit alone, and reports flag it.
    """
    status = _git("status", "--porcelain")
    git = {
        "commit": _git("rev-parse", "HEAD"),
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(status) if status is not None else None,
    }
    # The container has no git, so a run made inside it would record nothing
    # about the tree it came from -- the one fact this record exists to keep.
    # The host-side launcher therefore captures the state before entering the
    # container and passes it through the environment; the record says which
    # route it came by, so the two are never confused.
    if git["commit"] is None and os.environ.get("SCINT_GIT_COMMIT"):
        git = {
            "commit": os.environ["SCINT_GIT_COMMIT"],
            "branch": os.environ.get("SCINT_GIT_BRANCH"),
            "dirty": os.environ.get("SCINT_GIT_DIRTY") == "1",
            "source": "captured on the host before entering the container",
        }
    return {
        "captured_at": datetime.now(timezone.utc).isoformat(),
        "git": git,
        # geant4-config alone is not enough: it reports the same version string
        # for a series' beta and its release. scint.geant4.version_info reads the
        # tag out of the library itself and flags pre-release builds.
        "geant4": {
            **_geant4_version(),
            "config_version": _command_version("geant4-config", "--version"),
            "prefix": _command_version("geant4-config", "--prefix"),
        },
        # The dataset directories Geant4 will read, from the environment the
        # run inherits (geant4.sh sets them); a reproduction needs the versions
        # in these names, not only the toolkit tag.
        "geant4_datasets": {k: v for k, v in sorted(os.environ.items())
                            if k.startswith("G4") and k.endswith("DATA")},
        "root": {"version": _command_version("root-config", "--version")},
        "python": {
            "version": sys.version.split()[0],
            "executable": sys.executable,
            "packages": _installed_packages(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "in_container": Path("/.dockerenv").exists(),
        },
        "seed": seed,
    }


@dataclass(frozen=True)
class Run:
    """A materialised run directory."""

    id: str
    path: Path
    config: dict[str, Any]

    @property
    def config_path(self) -> Path:
        return self.path / "config.yaml"

    @property
    def env_path(self) -> Path:
        return self.path / "env.json"

    @property
    def status_path(self) -> Path:
        return self.path / "status.json"

    def status(self) -> dict[str, Any] | None:
        if not self.status_path.exists():
            return None
        return json.loads(self.status_path.read_text())


def create_run(
    config: dict[str, Any],
    *,
    runs_dir: Path | None = None,
    seed: int | None = None,
    exist_ok: bool = False,
) -> Run:
    """Materialise a run directory for ``config`` and record the environment.

    Refuses to overwrite an existing run unless ``exist_ok`` is set: because the
    identifier is a hash of the configuration, an existing directory means this
    exact run has been done before, and silently redoing it would destroy the
    earlier result.
    """
    base = runs_dir or DEFAULT_RUNS_DIR
    rid = run_id(config)
    path = base / rid

    if path.exists() and not exist_ok:
        raise RunExistsError(
            f"run {rid} already exists at {path}. Pass exist_ok=True to reuse it, "
            "or change the configuration."
        )

    path.mkdir(parents=True, exist_ok=True)
    path.joinpath("config.yaml").write_text(canonical_yaml(config))
    path.joinpath("env.json").write_text(json.dumps(capture_environment(seed), indent=2) + "\n")
    return Run(id=rid, path=path, config=config)


def record_status(run: Run, **fields: Any) -> None:
    """Write or update ``status.json`` for a run."""
    current = run.status() or {}
    current.update(fields)
    run.status_path.write_text(json.dumps(current, indent=2) + "\n")


def checksum(path: Path, chunk_size: int = 1 << 20) -> str:
    """Return the SHA-256 of a file, for output integrity records."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def iter_runs(runs_dir: Path | None = None):
    """Yield every complete run under ``runs_dir``, oldest identifier first."""
    base = runs_dir or DEFAULT_RUNS_DIR
    if not base.exists():
        return
    for path in sorted(base.iterdir()):
        config_file = path / "config.yaml"
        if path.is_dir() and config_file.exists():
            yield Run(id=path.name, path=path, config=yaml.safe_load(config_file.read_text()))


def load_run(rid: str, runs_dir: Path | None = None) -> Run:
    """Load a single run by identifier."""
    base = runs_dir or DEFAULT_RUNS_DIR
    path = base / rid
    config_file = path / "config.yaml"
    if not config_file.exists():
        raise FileNotFoundError(f"no run {rid} under {base}")
    return Run(id=rid, path=path, config=yaml.safe_load(config_file.read_text()))
